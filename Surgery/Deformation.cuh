#include "GLAD\glad.h"
#include "PhysParam.cuh"
#include "helper_cuda.h"

#include <cuda_runtime.h>
#include <cuda_gl_interop.h> //include the appropriate gl headers before including cuda_gl_interop.h


#include <cuda_profiler_api.h >


extern "C" void RegisterCudaBuffer(struct cudaGraphicsResource** dst, GLuint srcvbo);
extern "C" void UnRegisterGLBufferObject(struct cudaGraphicsResource *cuda_vbo_resource);
extern "C" void AllocateArray(void **devPtr, size_t size);
extern "C" void FreeArray(void *devPtr);
extern "C" void CopyArrayToDevice(void *device, const void *host, int offset, int size);
extern "C" void CopyArrayFromDevice(void *host, const void *device, struct cudaGraphicsResource **cuda_vbo_resource, int size);
extern "C" void CopyArrayFromDeviceToDevice(void *dTarg, void *dSrc, struct cudaGraphicsResource **cuda_vbo_resource, int size);
extern "C" void SetValueToDevice(void *device, int value, int size);
extern "C" void* MapGLBufferObject(struct cudaGraphicsResource **cuda_vbo_resource);
extern "C" void UnMapGLBufferObject(struct cudaGraphicsResource *cuda_vbo_resource);
extern "C" void SetParameters(PhysParam *hostParams);

extern "C" void ComputeGridSize(uint numParticles, uint defaultBlockSize, uint &numBlocks, uint &numThreads);
extern "C" void CalcHash(uint *gridParticleHash, uint *gridParticleIndex, double3 *pos, int numParticles);
extern "C" void SortParticles(uint *dGridParticleHash, uint *dGridParticleIndex, uint numParticles);
extern "C" void FindCellStartEnd(uint *cellStart, uint *cellEnd, uint *gridParticleHash, uint numParticles, uint numCells);
extern "C" void ReorderData(double3 *sortedPos, double3 *sortedVel, uint *gridParticleIndex, double3 *oldPos, double3 *oldVel, uint numParticles);
extern "C" void FindNeighborsWithinDst(int* neighborIndex, uint *neighborCount, uint *gridParticleIndex, double3 *sortedPos, double dst, uint *cellStart, uint *cellEnd, uint numParticles);

extern "C" void UpdateDistanceConstraints(int2 *dDistanceConstraints, double *dRestLength, uint *dNumDistanceConstraints, double3 *dPos, int* dNeighborIndex, uint *dNeighborCount, uint numParticles);


//PBD约束操作，updateVelocity弃用，改用semiImplicitEuler；IntegrateExplicitWithDamping弃用
extern "C" unsigned int nextPow2(unsigned int x);
extern "C" void FindParticlesWithinDstFromTarget(int* nearbyIndex, uint *nearbyCount, uint *gridParticleIndex, double3 *sortedPos, double3 targ, double dst, uint *cellStart, uint *cellEnd);
extern "C" void UpdateVelocities(double3 *dVel, double dt, int* nearbyIndex, uint* nearbyCount, double3 hapticForce, double3 environmentForce, double *dInvMasses, uint numParticles);
void SemiImplicitEuler(double3 *dPredictPos, double3 *dPos, double3 *dVel, double dt, int* nearbyIndex, uint *nearbyCount, double3 hapticForce, double3 environmentForce, double *dInvMasses, uint numParticles);
template<class T> void SumFromArray(T *dOut, T *dSrc, uint numParticles);
extern "C" void IntegrateExplicitWithDamping(double3 *dPredictPos, double3 *dPos, double3 *dVel, double damping, double dt, uint numParticles);
void UpdateAllConstraints(double3 *dPredictPos, double *dInvMasses, int2 *dDistanceConstraint, double compressionStiffness, double stretchStiffness, double *dRestLength, int numDistanceConstraints, uint numIteration);
void DetermineFinalVel(double3 *dPredictPos, double3 *dPos, double3 *dVel, double dt, uint numParticles);
