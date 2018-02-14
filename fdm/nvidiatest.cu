#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <time.h>   
#include <sys/time.h>  


#include <cuda.h>
#include<cuda_runtime.h>
#include<device_launch_parameters.h>

#define BLOCKSIZE_x 32
#define BLOCKSIZE_y 32


#define l 100
#define dt 0.01
#define D 10.0
#define d 0.4


//using namespace std;
double my_gettimeofday(){
  struct timeval tmp_time;
  gettimeofday(&tmp_time, NULL);
  return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}


/*****************/
/* CUDA MEMCHECK */
/*****************/
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, char *file, int line, bool abort = true)
{
    if (code != cudaSuccess)
    {
      fprintf(stderr, "GPUassert: %s %s %dn", cudaGetErrorString(code), file, line);
      if (abort) { exit(code); }
    }
}


int iDivUp(int hostPtr, int b){ return ((hostPtr % b) != 0) ? (hostPtr / b + 1) : (hostPtr / b); }


__global__ void solve(double A[l][l])
{
	int i = threadIdx.x+blockIdx.x * blockDim.x;
	int j = threadIdx.y+blockIdx.y * blockDim.y;

	if( i < l-1 && j < l-1 && (i!=0) && (j!=0) )
	{
		A[i][j] = A[i][j]*(1-d*dt)/l+(A[i-1][j] + A[i+1][j] + A[i][j+1] + A[i][j-1] - 4*A[i][j])*D*dt/l;
	}
}




int main(){
	float phi0=0.4;
	double cpu_mesh[l][l]; 
	double cpu_res[l][l];
	double (*gpu_mesh)[l]; //pointers to arrays of dimension N
	double (*gpu_res)[l];

	/* Initializing cpu_mesh with source at the center*/
	for(int i=0 ; i< l; i++){
		for(int j=0 ; j<l ; ++j){
			cpu_mesh[i][j]=0.0;
		}
	}
	cpu_mesh[l/2-1][l/2 -1]=phi0;
        cpu_mesh[l/2][l/2 -1]=phi0;
        cpu_mesh[l/2][l/2]=phi0;
        cpu_mesh[l/2-1][l/2]=phi0;

	/* Allocation */
	
	cudaMalloc((void**)&gpu_mesh, (l*l)*sizeof(double));
	cudaMalloc((void**)&gpu_res, (l*l)*sizeof(double));

	//copying from host to device
	double debut = my_gettimeofday();
	double debutTransfert = my_gettimeofday();
	gpuErrchk(cudaMemcpy(gpu_mesh, cpu_mesh, (l*l)*sizeof(double), cudaMemcpyHostToDevice));
	double finTransfert = my_gettimeofday();
	std::cout << "Transfert CPU vers GPU :" << finTransfert-debutTransfert << std::endl;
	dim3 gridSize(iDivUp(l, BLOCKSIZE_x), iDivUp(l, BLOCKSIZE_y));
	dim3 blockSize(BLOCKSIZE_y, BLOCKSIZE_x);
	//solve <<<gridSize, blockSize>>> (gpu_mesh, gpu_res, D, dt, d);
	for(int i=0; i<1000; ++i){
		solve<<<gridSize, blockSize>>> (gpu_mesh);
	}
	
	debutTransfert = my_gettimeofday();
	cudaMemcpy(cpu_res, gpu_mesh, (l*l)*sizeof(double), cudaMemcpyDeviceToHost);
	std::cout << "Transfert GPU vers CPU :" << finTransfert-debutTransfert << std::endl;
	finTransfert = my_gettimeofday();
	
	double fin= my_gettimeofday();
	std::cout << "Temps calcul :" << fin-debut << std::endl;
	/*for (int i = 0; i < l; i++){
		for (int j = 0; j < l; j++){
      			std::cout << cpu_res[i][j] << " ";
		}
		std::cout << std::endl; 
	}*/
	return 0;
}

