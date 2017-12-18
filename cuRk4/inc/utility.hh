#include <cuda_runtime.h>
#include <cuda.h>
#include <math.h>


template<class state_type>
__global__ void sumk(int n, state_type* out, state_type* v, state_type* coef, int n_v) {

	int grain = ceil((double)n/(gridDim.x*blockDim.x));
	int index = blockIdx.x*blockDim.x*grain+threadIdx.x;

	for(int i=0; i<grain; i++) {
		if(index<n) {
			out[index] = 0;
			for(int j=0; j<n_v; j++) {
				out[index]+=coef[j]*v[j*n+index];
			}
		}
		index+=blockDim.x;
	}

}


