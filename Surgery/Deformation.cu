//#include "math\CVector3d.h"//包含会引起错误

#include "Deformation.cuh"

#include <thrust\device_ptr.h>
#include <thrust\for_each.h>
#include <thrust\iterator\zip_iterator.h>
#include <thrust\sort.h>

#include <iostream>


// error makro
#define CUERR {                                                              \
	cudaError_t err;                                                         \
	if ((err = cudaGetLastError()) != cudaSuccess) {                         \
		std::cout << "CUDA error: " << cudaGetErrorString(err) << " : "      \
					<< __FILE__ << ", line " << __LINE__ << std::endl;       \
		/*exit(1);                                                            */ \
	}                                                                        \
}

//---------------------------------------------------------------------------
//global variables
//---------------------------------------------------------------------------
const int MySIZE = 43250;

// simulation parameters in constant memory
__constant__ PhysParam params;


//---------------------------------------------------------------------------
//function defination
//---------------------------------------------------------------------------

extern "C" void RegisterCudaBuffer(struct cudaGraphicsResource** dst, GLuint srcvbo)
{
	cudaGraphicsGLRegisterBuffer(dst, srcvbo, cudaGraphicsMapFlagsNone); CUERR
}

extern "C" void UnRegisterGLBufferObject(struct cudaGraphicsResource *cuda_vbo_resource)
{
	cudaGraphicsUnregisterResource(cuda_vbo_resource);
}

extern "C" void AllocateArray(void **devPtr, size_t size)
{
	cudaMalloc(devPtr, size); CUERR
}

extern "C" void FreeArray(void *devPtr)
{
	cudaFree(devPtr); CUERR
}

extern "C" void CopyArrayToDevice(void *device, const void *host, int offset, int size)
{
	cudaMemcpy((char*)device + offset, host, size, cudaMemcpyHostToDevice); CUERR
}

extern "C" void CopyArrayFromDevice(void *host, const void *device, struct cudaGraphicsResource **cuda_vbo_resource, int size)
{
	if (cuda_vbo_resource)
	{
		device = MapGLBufferObject(cuda_vbo_resource);
	}

	cudaMemcpy(host, device, size, cudaMemcpyDeviceToHost); CUERR

	if (cuda_vbo_resource)
	{
		UnMapGLBufferObject(*cuda_vbo_resource);
	}
}

extern "C" void CopyArrayFromDeviceToDevice(void *dTarg, void *dSrc, struct cudaGraphicsResource **cuda_vbo_resource, int size)
{
	if (cuda_vbo_resource)
	{
		dSrc = MapGLBufferObject(cuda_vbo_resource);
	}

	cudaMemcpy(dTarg, dSrc, size, cudaMemcpyDeviceToDevice); CUERR

	if (cuda_vbo_resource)
	{
		UnMapGLBufferObject(*cuda_vbo_resource);
	}
}

extern "C" void SetValueToDevice(void *device, int value, int size)
{
	cudaMemset(device, value, size); CUERR
}

extern "C" void* MapGLBufferObject(struct cudaGraphicsResource **cuda_vbo_resource)
{
	void *ptr;
	cudaGraphicsMapResources(1, cuda_vbo_resource, 0); CUERR
	
	size_t num_bytes;
	cudaGraphicsResourceGetMappedPointer((void**)&ptr, &num_bytes, *cuda_vbo_resource); CUERR

	return ptr;
}

extern "C" void UnMapGLBufferObject(struct cudaGraphicsResource *cuda_vbo_resource)
{
	cudaGraphicsUnmapResources(1, &cuda_vbo_resource, 0); CUERR
}

extern "C" void SetParameters(PhysParam *hostParams)
{
	cudaMemcpyToSymbol(params, hostParams, sizeof(PhysParam)); CUERR
}


//--------------------------------------------------------------------
//Find neighbors
//--------------------------------------------------------------------

// calculate position in uniform grid
__device__ int3 CalcGridPos(double3 p)
{
	int3 gridPos;
	gridPos.x = floor((p.x - params.worldOrigin.x) / params.cellSize.x);
	gridPos.y = floor((p.y - params.worldOrigin.y) / params.cellSize.y);
	gridPos.z = floor((p.z - params.worldOrigin.z) / params.cellSize.z);
	return gridPos;
}

// calculate address in grid from position (clamping to edges)
__device__ uint CalcGridHash(int3 gridPos)
{
	gridPos.x = gridPos.x & (params.gridSize.x - 1);  // wrap grid, assumes size is power of 2
	gridPos.y = gridPos.y & (params.gridSize.y - 1);
	gridPos.z = gridPos.z & (params.gridSize.z - 1);

	return gridPos.z * params.gridSize.y * params.gridSize.x + gridPos.y * params.gridSize.x + gridPos.x;
}

__global__ void CalcHashD(uint *gridParticleHash, uint *gridParticleIndex, double3 *pos, uint numParticles)
{
	uint index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= numParticles)
		return;

	volatile double3 p = pos[index];

	//get address in grid
	int3 gridPos = CalcGridPos(make_double3(p.x, p.y, p.z));
	uint hash = CalcGridHash(gridPos);

	gridParticleHash[index] = hash;
	gridParticleIndex[index] = index;
}


//Round a / b to nearest higher integer value
uint iDivUp(uint a, uint b)
{
	return (a % b != 0) ? (a / b + 1) : (a / b);
}

extern "C" void ComputeGridSize(uint numParticles, uint defaultBlockSize, uint &numBlocks, uint &numThreads)
{
	numThreads = min(defaultBlockSize, numParticles);
	numThreads = nextPow2(numThreads);
	numBlocks = iDivUp(numParticles, numThreads);
}

extern "C" void CalcHash(uint *gridParticleHash, uint *gridParticleIndex, double3 *pos, int numParticles)
{
	uint numThreads, numBlocks;
	ComputeGridSize(numParticles, 256, numBlocks, numThreads);

	CalcHashD << <numBlocks, numThreads >> > (gridParticleHash, gridParticleIndex, pos, numParticles);
	cudaDeviceSynchronize(); CUERR
}

// sort particles based on hash
extern "C" void SortParticles(uint *dGridParticleHash, uint *dGridParticleIndex, uint numParticles)
{
	thrust::sort_by_key(thrust::device_ptr<uint>(dGridParticleHash), thrust::device_ptr<uint>(dGridParticleHash + numParticles), thrust::device_ptr<uint>(dGridParticleIndex));
}

//------------------------------------------------------------------------------------------------------
// find the start & end particle index of each cell in the sorted hash array
//------------------------------------------------------------------------------------------------------
__global__ void FindCellStartEndD(uint   *cellStart,        // output: cell start index
								uint   *cellEnd,          // output: cell end index
								uint   *gridParticleHash, // input: sorted grid hashes
								uint    numParticles)
{
	extern __shared__ uint sharedHash[];    // blockSize + 1 elements

	uint index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= numParticles)
		return;

	uint hash;
	hash = gridParticleHash[index]; //当前particle的cell索引

	// Load hash data into shared memory so that we can look
	// at neighboring particle's hash value without loading
	// two hash values per thread
	sharedHash[threadIdx.x + 1] = hash;

	if (index > 0 && threadIdx.x == 0) //block's 1st thread except first block
	{
		// first thread in block must load neighbor particle hash
		sharedHash[0] = gridParticleHash[index - 1];
	}
	__syncthreads();

	// If this particle has a different cell index to the previous
	// particle then it must be the first particle in the cell,
	// As it isn't the first particle, it must also be the cell end of
	// the previous particle's cell (all particle's index in one cell < cell end)
	if (index == 0 || hash != sharedHash[threadIdx.x])
	{
		cellStart[hash] = index;

		if (index > 0)
			cellEnd[sharedHash[threadIdx.x]] = index;
	}
	if (index == numParticles - 1)
	{
		cellEnd[hash] = index + 1;
	}
}

extern "C" void FindCellStartEnd(uint *cellStart, uint *cellEnd,        //output
							uint *gridParticleHash, 
							uint numParticles, uint numCells)
{
	uint numThreads, numBlocks;
	ComputeGridSize(numParticles, 256, numBlocks, numThreads);

	// set all cells to empty
	cudaMemset(cellStart, 0xffffffff, numCells * sizeof(uint)); CUERR

	uint smemSize = sizeof(uint) * (numThreads + 1);
	FindCellStartEndD << <numBlocks, numThreads, smemSize >> > (cellStart, cellEnd, gridParticleHash, numParticles);
	cudaDeviceSynchronize(); CUERR
}

//----------------------------------------------------------------------------------------
//reorder pos & vel base on sorted particle indices
//----------------------------------------------------------------------------------------
__global__ void ReorderDataD(double3 *sortedPos,        // output: sorted positions
							double3 *sortedVel,        // output: sorted velocities
							uint   *gridParticleIndex, // input: sorted particle indices
							double3 *oldPos,           // input: unsorted position array
							double3 *oldVel,           // input: unsorted velocity array
							uint    numParticles)
{
	uint index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= numParticles)
		return;

	uint sortedIndex = gridParticleIndex[index];

	sortedPos[index] = oldPos[sortedIndex];
	sortedVel[index] = oldVel[sortedIndex];
}

extern "C" void ReorderData(double3 *sortedPos, double3 *sortedVel, //output
							uint *gridParticleIndex,
							double3 *oldPos, double3 *oldVel,
							uint numParticles)
{
	uint numThreads, numBlocks;
	ComputeGridSize(numParticles, 256, numBlocks, numThreads);

	ReorderDataD << <numBlocks, numThreads >> > (sortedPos, sortedVel, gridParticleIndex, oldPos, oldVel, numParticles);
	cudaDeviceSynchronize(); CUERR
}


//-------------------------------------------------------------------------------
//find fixed distance neighbors from particles in neighboring cells
//-------------------------------------------------------------------------------
__device__ bool isTwoParticleWithinDst(double3 a, double3 b, double dst)
{
	double distance = (b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y) + (b.z - a.z) * (b.z - a.z);

	return distance <= (dst * dst) ? true : false;
}

__global__ void FindNeighborsWithinDstD(int* neighborIndex, //output
										uint *neighborCount, //output
										uint *gridParticleIndex,
										double3 *sortedPos, double dst,
										uint *cellStart, uint *cellEnd,
										uint numParticles)
{
	uint index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= numParticles)
		return;

	uint originalIndex = gridParticleIndex[index];
	double3 pos = sortedPos[index];
	int3 gridPos = CalcGridPos(pos);

	uint neighborId = 0;

	for (int z = -1; z <= 1; z++)
	{
		for (int y = -1; y <= 1; y++)
		{
			for (int x = -1; x <= 1; x++)
			{
				int3 neighbourGridPos = gridPos + make_int3(x, y, z);
				uint gridHash = CalcGridHash(neighbourGridPos);

				uint startIndex = cellStart[gridHash];
				if (startIndex != 0xffffffff) //cell is not empty
				{
					uint endIndex = cellEnd[gridHash];
					for (int j = startIndex; j < endIndex; j++)
					{
						if (j != index)
						{
							double3 targPos = sortedPos[j];
							if (isTwoParticleWithinDst(pos, targPos, dst))
							{
								if (neighborId >= 26)
									break;

								uint originalTargIndex = gridParticleIndex[j];

								neighborIndex[originalIndex * 26 + neighborId] = originalTargIndex; //默认每个particle周围最多有26个邻居
								neighborId++;
							}
						}
					}
					neighborCount[originalIndex] = neighborId;
				}
			}
		}
	}
}

extern "C" void FindNeighborsWithinDst(int* neighborIndx, //output: original unsorted particles neighbor indeices
										uint *neighborCount, //output
										uint *gridParticleIndex,
										double3 *sortedPos, double dst,
										uint *cellStart, uint *cellEnd,
										uint numParticles)
{
	uint numThreads, numBlocks;
	ComputeGridSize(numParticles, 256, numBlocks, numThreads);

	FindNeighborsWithinDstD << <numBlocks, numThreads >> > (neighborIndx, neighborCount, gridParticleIndex, sortedPos, dst, cellStart, cellEnd, numParticles);
	cudaDeviceSynchronize(); CUERR
}

//---------------------------------------------------------------------------
// build distance constraints based on neighborIdx
//---------------------------------------------------------------------------

__global__ void UpdateDistanceConstraintsD(int2 *dDistanceConstraints, //output
											double *dRestLength, //output
											uint *dNumDistanceConstraints, //output
											double3 *dPos,
											int *dNeighborIndex, uint *dNeighborCount,
											uint numParticles)
{
	uint index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= numParticles)
		return;

	uint constraintId = 0;
	for (uint j = 0; j < dNeighborCount[index]; j++)
	{
		int neighborIdx = dNeighborIndex[index * 26 + j];
		if (neighborIdx > index)
		{
			int2 pair = make_int2(index, neighborIdx);
			double3 dis = dPos[index] - dPos[neighborIdx];

			dDistanceConstraints[constraintId] = pair;
			dRestLength[constraintId] = YH::Length(dis);
			constraintId++;
		}
	}

	*dNumDistanceConstraints = constraintId;
}

extern "C" void UpdateDistanceConstraints(int2 *dDistanceConstraints, //output
											double *dRestLength, //output
											uint *dNumDistanceConstraints, //output
											double3 *dPos,
											int* dNeighborIndex, uint *dNeighborCount,
											uint numParticles)
{
	uint numThreads, numBlocks;
	ComputeGridSize(numParticles, 256, numBlocks, numThreads);

	UpdateDistanceConstraintsD << <1, 1 >> > (dDistanceConstraints, dRestLength, dNumDistanceConstraints, dPos, dNeighborIndex, dNeighborCount, numParticles);
	cudaDeviceSynchronize(); CUERR
}

//---------------------------------------------------------------------------
// find nearby particle indices from target position
//---------------------------------------------------------------------------

__global__ void FindParticlesWithinDstFromTargetD(int* nearbyIndex, //output
													uint *nearbyCount, //output
													uint *gridParticleIndex,
													double3 *sortedPos,
													double3 targ, double dst,
													uint *cellStart, uint *cellEnd)
{
	int3 gridPos = CalcGridPos(targ);

	uint nearbyId = 0;

	//当dst取得比较大时，这里可以放宽搜索范围
	for (int z = -2; z <= 2; z++)
	{
		for (int y = -2; y <= 2; y++)
		{
			for (int x = -2; x <= 2; x++)
			{
				int3 neighbourGridPos = gridPos + make_int3(x, y, z);
				uint gridHash = CalcGridHash(neighbourGridPos);

				uint startIndex = cellStart[gridHash];
				if (startIndex != 0xffffffff) //cell is not empty
				{
					uint endIndex = cellEnd[gridHash];
					for (int j = startIndex; j < endIndex; j++)
					{
						double3 testPos = sortedPos[j];
						if (isTwoParticleWithinDst(testPos, targ, dst))
						{
							uint originalTestPosIndex = gridParticleIndex[j];
							nearbyIndex[nearbyId] = originalTestPosIndex;
							nearbyId++;
						}
					}
					nearbyCount[0] = nearbyId;
				}
			}
		}
	}
}

extern "C" void FindParticlesWithinDstFromTarget(int* nearbyIndex, //output
												uint *nearbyCount, //output
												uint *gridParticleIndex, 
												double3 *sortedPos, 
												double3 targ, double dst, 
												uint *cellStart, uint *cellEnd)
{
	//！注意，这里不用每个粒子循环都求一遍，要更改
	FindParticlesWithinDstFromTargetD << <1, 1 >> > (nearbyIndex, nearbyCount, gridParticleIndex, sortedPos, targ, dst, cellStart, cellEnd);
	cudaDeviceSynchronize(); CUERR
}

//-----------------------------------------------------------------------------
// update velocity: v += f/m * dt
//-----------------------------------------------------------------------------

__global__ void UpdataVelocitiesD(double3 *dVel, double dt, 
								int* nearbyIndex, uint *nearbyCount,
								double3 hapticForce, double3 environmentForce, 
								double *dInvMasses, 
								uint numParticles)
{
	uint index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= numParticles)
		return;

	for (uint i = 0; i < nearbyCount[0]; i++)
	{
		if (index == nearbyIndex[i])
		{
			dVel[index] += hapticForce * dInvMasses[index] * dt;
		}
	}

	dVel[index] += environmentForce * dInvMasses[index] * dt;
}

extern "C" void UpdateVelocities(double3 *dVel, double dt, 
								int* nearbyIndex, uint* nearbyCount, 
								double3 hapticForce, double3 environmentForce, 
								double *dInvMasses, 
								uint numParticles)
{
	uint numThreads, numBlocks;
	ComputeGridSize(numParticles, 256, numBlocks, numThreads);

	UpdataVelocitiesD << <numBlocks, numThreads >> > (dVel, dt, nearbyIndex, nearbyCount, hapticForce, environmentForce, dInvMasses, numParticles);
	cudaDeviceSynchronize(); CUERR
}

//------------------------------------------------------------------------------
//semi implicit update vel & pos
//------------------------------------------------------------------------------

__global__ void SemiImplicitEulerD(double3 *dPredictPos, 
	double3 *dPos, double3 *dVel,
	double dt,
	int* nearbyIndex, uint *nearbyCount,
	double3 hapticForce, double3 environmentForce,
	double *dInvMasses,
	uint numParticles)
{
	uint index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= numParticles)
		return;

	dVel[index].x = dVel[index].y = dVel[index].z = 0;

	for (uint i = 0; i < nearbyCount[0]; i++)
	{
		if (index == nearbyIndex[i])
		{
			dVel[index] += hapticForce * dInvMasses[index] * dt;
		}
	}

	dVel[index] += environmentForce * dInvMasses[index] * dt;

	dPredictPos[index] = dPos[index] + dVel[index] * dt;
}

void SemiImplicitEuler(double3 *dPredictPos, 
	double3 *dPos, double3 *dVel,
	double dt, 
	int* nearbyIndex, uint *nearbyCount, 
	double3 hapticForce, double3 environmentForce, 
	double *dInvMasses, 
	uint numParticles)
{
	uint numThreads, numBlocks;
	ComputeGridSize(numParticles, 256, numBlocks, numThreads);

	SemiImplicitEulerD << <numBlocks, numThreads >> > (dPredictPos, dPos, dVel, dt, nearbyIndex, nearbyCount, hapticForce, environmentForce, dInvMasses, numParticles);
	cudaDeviceSynchronize(); CUERR
}

//------------------------------------------------------------------------------
// predict position
//------------------------------------------------------------------------------

// Utility class used to avoid linker errors with extern
// unsized shared memory arrays with templated type
template<class T>
struct SharedMemory
{
	__device__ inline operator T *()
	{
		extern __shared__ int __smem[];
		return (T *)__smem;
	}

	__device__ inline operator const T *() const
	{
		extern __shared__ int __smem[];
		return (T *)__smem;
	}
};

template <class T>
__global__ void Reduce_ArraySumD(T* dOut, T* dIn, uint numData)
{
	uint globalId = blockIdx.x * blockDim.x + threadIdx.x;
	uint localId = threadIdx.x;

	if (globalId >= numData)
		return;

	T *sdata = SharedMemory<T>();

	sdata[localId] = dIn[globalId];
	__syncthreads();

	for (unsigned i = blockDim.x / 2; i > 0; i >>= 1)
	{
		if (localId < i)
		{
			sdata[localId] += sdata[localId + i];
		}
		__syncthreads();
	}

	// write result for this block to global mem
	if (localId == 0)
	{
		dOut[blockIdx.x] = sdata[0];
	}
}


__global__ void AverageD(double3* dAver, double3* dSum, uint count)
{
	uint index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index != 0)
		return;

	(*dAver).x = (*dSum).x / count;
	(*dAver).y = (*dSum).y / count;
	(*dAver).z = (*dSum).z / count;
}

extern "C" unsigned int nextPow2(unsigned int x)
{
	--x;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return ++x;
}

template<class T>
void SumFromArray(T *dOut, T *dSrc, uint numParticles)
{
	uint numThreads, numBlocks;
	ComputeGridSize(numParticles, 256, numBlocks, numThreads);

	T *dTmpSum, *dSum;
	AllocateArray((void**)&dTmpSum, sizeof(T) * numBlocks);
	AllocateArray((void**)&dSum, sizeof(T) * numBlocks);

	uint sMemSize = sizeof(T) * numThreads;
	Reduce_ArraySumD<T> << <numBlocks, numThreads, sMemSize >> > (dTmpSum, dSrc, numParticles); //归约操作必须保证线程数是2的指数倍
	cudaDeviceSynchronize(); CUERR

	while (numBlocks != 1)
	{
		ComputeGridSize(numBlocks, 256, numBlocks, numThreads);

		sMemSize = sizeof(T) * numThreads;
		Reduce_ArraySumD<T> << <numBlocks, numThreads, sMemSize >> > (dSum, dTmpSum, numParticles);
		cudaDeviceSynchronize(); CUERR

		CopyArrayFromDeviceToDevice(dTmpSum, dSum, 0, sizeof(T) * numBlocks);
	}

	CopyArrayFromDeviceToDevice(dOut, dSum, 0, sizeof(T));
}

template 
void SumFromArray<double3>(double3 *dOut, double3 *dSrc, uint numParticles);

template
void SumFromArray<Mat3d>(Mat3d *dOut, Mat3d *dSrc, uint numParticles);


__global__ void Prepare_L_I_R_ArrayD(double3 *dLArray, Mat3d *dIArray, double3 *dR, double3 *dPos, double3 *dVel, double3 *dAverPos, uint numParticles)
{
	uint index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= numParticles)
		return;

	double3 ri = dPos[index] - dAverPos[0];
	dLArray[index] = YH::Cross(ri, dVel[index]);

	Mat3d tmp(0, ri.z, -ri.y,
			-ri.z, 0, ri.x,
			ri.y, -ri.x, 0);
	dIArray[index] = tmp * YH::Transpose(tmp);

	dR[index] = ri;
}

__global__ void PredictPosD(double3 *dPredictPos, //output
							double3 *dPos, 
							double3 *dVel, double3 *dAverVel, 
							Mat3d *dI, double3 *dL, double3 *dR, 
							double damping, double dt, 
							uint numParticles)
{
	uint index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= numParticles)
		return;

	double3 W = YH::Invert(dI[0]) * dL[0];
	double3 deltVi = dAverVel[0] + YH::Cross(W, dR[index]) - dVel[index];

	dVel[index] += damping * deltVi;
	dPredictPos[index] = dPos[index] + dVel[index] * dt;
}

extern "C" void IntegrateExplicitWithDamping(double3 *dPredictPos, double3 *dPos, double3 *dVel, double damping, double dt, uint numParticles)
{
	double3 *dSumPos, *dSumVel, *dAverPos, *dAverVel;
	AllocateArray((void**)&dSumPos, sizeof(double3));
	AllocateArray((void**)&dSumVel, sizeof(double3));
	AllocateArray((void**)&dAverPos, sizeof(double3));
	AllocateArray((void**)&dAverVel, sizeof(double3));

	//---------------------位置平均----------------------
	SumFromArray<double3>(dSumPos, dPos, numParticles);
	AverageD << <1, 1 >> > (dAverPos, dSumPos, numParticles);
	//-------------------------速度平均------------------
	SumFromArray<double3>(dSumVel, dVel, numParticles);
	AverageD << <1, 1 >> > (dAverVel, dSumVel, numParticles);

	double3 hSumPos;
	CopyArrayFromDevice(&hSumPos, dSumPos, 0, sizeof(double3));
	//printf("dSumPos: %f %f %f\n", hSumPos.x, hSumPos.y, hSumPos.z);
	
	//----------------------预测位置_准备---------------------
	Mat3d hI, *dI; 
	hI.Identity();

	double3 *dL, *dLArray;
	double3 *dR;
	Mat3d *dIArray;

	AllocateArray((void**)&dI, sizeof(Mat3d));
	AllocateArray((void**)&dL, sizeof(double3));
	AllocateArray((void**)&dLArray, sizeof(double3) * numParticles);
	AllocateArray((void**)&dR, sizeof(double3) * numParticles);
	AllocateArray((void**)&dIArray, sizeof(Mat3d) * numParticles);

	CopyArrayToDevice(dI, &hI, 0, sizeof(Mat3d));
	SetValueToDevice(dL, 0, sizeof(double3));
	SetValueToDevice(dLArray, 0, sizeof(double3) * numParticles);
	SetValueToDevice(dR, 0, sizeof(double3) * numParticles);
	SetValueToDevice(dIArray, 0, sizeof(Mat3d) * numParticles);

	uint numThreads, numBlocks;
	ComputeGridSize(numParticles, 256, numBlocks, numThreads);
	Prepare_L_I_R_ArrayD << <numBlocks, numThreads >> > (dLArray, dIArray, dR, dPos, dVel, dAverPos, numParticles);

	//做一次归约
	SumFromArray<double3>(dL, dLArray, numParticles);
	SumFromArray<Mat3d>(dI, dIArray, numParticles);

	//test output
	/*CopyArrayFromDevice(&hI, dI, 0, sizeof(Mat3d));

	printf("dI: \n%f %f %f \n %f %f %f \n %f %f %f\n", hI(0, 0), hI(0, 1), hI(0, 2),
													hI(1, 0), hI(1, 1), hI(1, 2),
													hI(2, 0), hI(2, 1), hI(2, 2));
	*/
	//-------------------------预测位置-----------------------------
	PredictPosD << <numBlocks, numThreads >> > (dPredictPos, dPos, dVel, dAverVel, dI, dL, dR, damping, dt, numParticles);
	cudaDeviceSynchronize(); CUERR
	

	FreeArray(dSumPos);
	FreeArray(dSumVel);
	FreeArray(dAverPos);
	FreeArray(dAverVel);
	FreeArray(dI);
	FreeArray(dL);
	FreeArray(dR);
	FreeArray(dLArray);
	FreeArray(dIArray);
}

//-----------------------------------------------------------------------------
// update all constraints
//-----------------------------------------------------------------------------

__global__ void UpdateDistanceConstraintsD(double3 *dPredictPos, 
	double* dInvMasses,
	int2 *dDistanceConstraint, double compressionStiffness, double stretchStiffness, double *dRestLength, int numDistanceConstraints)
{
	uint index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= numDistanceConstraints)
		return;

	int p0 = dDistanceConstraint[index].x;
	int p1 = dDistanceConstraint[index].y;
	double restLen = dRestLength[index];

	double invMass0 = dInvMasses[p0];
	double invMass1 = dInvMasses[p1];
	double invMassSum = invMass0 + invMass1;

	if (invMassSum <= 0.000001)
		return;

	double3 n = dPredictPos[p1] - dPredictPos[p0];
	double len = YH::Length(n);
	n = YH::Normalize(n);

	double3 corr;
	if (len < restLen)
		corr = compressionStiffness * n * (len - restLen) / invMassSum;
	else
		corr = stretchStiffness * n * (len - restLen) / invMassSum;

	dPredictPos[p0] -= invMass0 * corr;
	dPredictPos[p1] += invMass1 * corr;
}

void UpdateAllConstraints(double3 *dPredictPos, //output
	double *dInvMasses,
	int2 *dDistanceConstraint, double compressionStiffness, double stretchStiffness, double *dRestLength, int numDistanceConstraints,
	uint numIteration)
{
	uint numThreads0, numBlocks0;
	ComputeGridSize(numDistanceConstraints, 256, numBlocks0, numThreads0);

	for (uint i = 0; i < numIteration; i++)
	{
		UpdateDistanceConstraintsD << <numBlocks0, numThreads0 >> > (dPredictPos, dInvMasses, dDistanceConstraint, compressionStiffness, stretchStiffness, dRestLength, numDistanceConstraints);
	}
	cudaDeviceSynchronize(); CUERR
}

//---------------------------------------------------------------------------
// Determine final vel
//---------------------------------------------------------------------------

__global__ void DetermineFinalVelD(double3 *dVel, //output
	double3 *dPredictPos, double3 *dPos, double dt, uint numParticles)
{
	uint index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= numParticles)
		return;

	dVel[index] = (dPredictPos[index] - dPos[index]) / dt;
}

void DetermineFinalVel(double3 *dPredictPos, double3 *dPos, double3 *dVel, double dt, uint numParticles)
{
	uint numThreads, numBlocks;
	ComputeGridSize(numParticles, 256, numBlocks, numThreads);

	DetermineFinalVelD << <numBlocks, numThreads >> > (dVel, dPredictPos, dPos, dt, numParticles);
	cudaDeviceSynchronize(); CUERR
}