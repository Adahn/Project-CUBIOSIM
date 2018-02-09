#pragma once
#include <cuda_runtime.h>
#include <cuda.h>
#include <stdio.h>
#include <string>


// Computes the sum of "n_v" arrays of size "n" stored at "v", 
// each weighted by its corresponding coefficient at "coef"
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
		}
		tid+=blockDim.x;

	}

}


// Loads Runge-Kutta coefficients arrays for each successive computation of 
// k1, k2, k3, and k4 at g_coefs on the GPU.
void load_coef(double* g_coefs, double step)
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


// simple debug display function
template<class state_type>
void debug_GPU(state_type* gpu_v, int dim, std::string message){

	state_type* cpu_v = new state_type[dim];
	cudaMemcpy(cpu_v, gpu_v, dim*sizeof(state_type), cudaMemcpyDeviceToHost);

	std::cout << message << '\n';
	for (size_t i = 0; i < dim; i++) {
		std::cout << cpu_v[i] << "\t";
	}
	std::cout << '\n';

	delete[] cpu_v;
}
