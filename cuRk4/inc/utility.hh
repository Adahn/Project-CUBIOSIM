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

void load_coef(double* g_coefs)
{
	double coefs[4*5];
	coefs[0]=1;		coefs[1]=step/2;
	coefs[2]=coefs[3]=coefs[4]=0;

	coefs[5]=1;	coefs[6]=0;	coefs[7]=step/2;
	coefs[8]=coefs[9]=0;

	coefs[10]=1;	coefs[11]=0;	coefs[12]=0;	coefs[13]=step;
	coefs[14]=0;

	coefs[15]=1;		coefs[16]=step/6;	coefs[17]=step/3;
	coefs[18]=step/3;	coefs[19]=step/6;

	cudaMemcpy(g_coefs, coefs, 4*5*sizeof(double), cudaMemcpyHostToDevice);

}
