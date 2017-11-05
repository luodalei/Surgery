#include "GLAD\glad.h"

#include <cuda_runtime.h>
#include <cuda_gl_interop.h> //include the appropriate gl headers before including cuda_gl_interop.h
#include <thrust\device_vector.h>

#include <cuda_profiler_api.h >

extern struct cudaGraphicsResource *vertsPos_resource;

extern "C" void RegisterCudaBuffer(struct cudaGraphicsResource** dst, GLuint srcvbo);
extern "C" void OutputVertsPos();