#include <cuda_runtime.h>
#include <cuda.h>
#include <stdio.h>

template<class state_type>
__global__ void sumk(int n, state_type* out, state_type* v, state_type* coef, int n_v) 
{
	/*0if(threadIdx.x==0) {
		printf("\n\t-- sumk --\n");
		printf("%d coefs: ", n_v);
		for(int i=0; i<n_v; i++) {
			printf("%f\t", coef[i]);
		}
		printf("\n");
	}
	__syncthreads();*/

	int grain = ceil((double)n/(gridDim.x*blockDim.x));
	int index = blockIdx.x*blockDim.x*grain+threadIdx.x;

	for(int i=0; i<grain; i++) {
		if(index<n) {
			out[index] = 0;
			for(int j=0; j<n_v; j++) {
				out[index]+=coef[j]*v[j*n+index];
			}
			//printf("out[%d] = %f\n", index, out[index]);
		}
		index+=blockDim.x;
		
	}

}

