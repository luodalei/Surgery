//#include "math\CVector3d.h"//包含会引起错误

#include "Deformation.cuh"

#include <iostream>

const int MySIZE = 43250;

// error makro
#define CUERR {                                                              \
	cudaError_t err;                                                         \
	if ((err = cudaGetLastError()) != cudaSuccess) {                         \
		std::cout << "CUDA error: " << cudaGetErrorString(err) << " : "      \
					<< __FILE__ << ", line " << __LINE__ << std::endl;       \
		/*exit(1);                                                            */ \
	}                                                                        \
}

struct cudaGraphicsResource *vertsPos_resource;

__device__ float3 d_pVelocities[MySIZE];


void InitVertsVelocities()
{
	float3 *d_pVel;
	
	cudaGetSymbolAddress((void**)&d_pVel, d_pVelocities); CUERR
	cudaMemset(d_pVel, 0, MySIZE * sizeof(float3));CUERR
}

extern "C" void RegisterCudaBuffer(struct cudaGraphicsResource** dst, GLuint srcvbo)
{
	cudaGraphicsGLRegisterBuffer(dst, srcvbo, cudaGraphicsMapFlagsNone); CUERR
}

size_t GetCudaMapPointer(struct cudaGraphicsResource** res, void **ptr)
{
	size_t num_bytes;
	cudaGraphicsMapResources(1, res, 0); CUERR
	cudaGraphicsResourceGetMappedPointer((void **)ptr, &num_bytes, *res); CUERR

	return num_bytes;
}

__global__ void _OutputVertsPos(double* _poses)
{
	for (unsigned i = 0; i < 10; i++)
	{
		printf("%f %f %f\n", _poses[i * 3], _poses[i * 3 + 1], _poses[i * 3 + 2]);
	}
}

extern "C" void OutputVertsPos()
{
	double *poses;
	GetCudaMapPointer(&vertsPos_resource, (void**)&poses);

	_OutputVertsPos << <1, 1 >> > (poses);
	cudaDeviceSynchronize(); CUERR
}