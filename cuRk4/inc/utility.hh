#include <cuda_runtime.h>
#include <cuda.h>
#include <stdio.h>

template<class state_type>
__global__ void sumk(int n, state_type* out, state_type* v, state_type* coef, int n_v) 
{

	int grain = ceil((double)n/(gridDim.x*blockDim.x));
	int tid = blockIdx.x*blockDim.x*grain+threadIdx.x;

	for(int i=0; i<grain; i++) {
		if(tid<n) {
			out[tid] = 0;
			for(int j=0; j<n_v; j++) {
				out[tid]+=coef[j]*v[j*n+tid];
			}
			//printf("out[%d] = %f\n", tid, out[tid]);
		}
		tid+=blockDim.x;
		
	}

}

