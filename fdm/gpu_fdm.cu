#include <iostream>
#include <cstdlib>
#include <cstdio>


#include<cuda.h>
#include<cuda_runtime.h>
#include<device_launch_parameters.h>

#define BLOCKSIZE_x 16
#define BLOCKSIZE_y 16

#define Nrows 10
#define Ncols 10


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




__global__ void solve(double *X, size_t pitch, int l, float D, float dt, float d){
        int idx = threadIdx.x+blockIdx.x * blockDim.x;
        int idy = threadIdx.y+blockIdx.y * blockDim.y;


	if(idx<l-1 && idy < l-1 && (idx!=0) && (idy!=0) )
        {
		double *Xres = (double *)((char*)X + idy * pitch);
		//double laplacien = -4*Xres[idx];
		//double laplacien = -4*Xres[idx]+Xres[idx-l]+Xres[idx+l]+Xres[idx-1]+Xres[idx+1];
		//laplacien/=l;
		Xres[idy]=Xres[idy+idx*(l+1)]+Xres[idy+idx*(l-1)]+Xres[(idy+1)+idx*l]+Xres[(idy-1)+idx*l]-4*Xres[idy+idx*l];
		//Xres[idx]=dt*D*laplacien +(1- d *dt)*Xres[idx];
		//Xres[idy]=-4*Xres[idx]+Xres[idx-1]+Xres[idx+1]+X[idy-1]+X[idy+1];
		//Xres[idy]=Xres[idy-1]+Xres[idy+1]-2*Xres[idy];
		//Xres[idx]=(Xres[idx-1]+Xres[idx+1]-2*Xres[idx])/l;
		//Xres[idx]=Xres[idx]*dt*D + (1-d*dt)*Xres[idx];
	}

}

int main(){
	float phi0=0.4;
	float D=0.4;
	float d=0.4;
	float dt=0.01;
	int nx=6;
	int ny=6;
	int hx=1;
	int hy=1;
	int lx=nx*hx;
	int ly=ny*hy;
	double cpu_mesh[lx][ly]; 
	double cpu_res[lx][ly];
	double* gpu_mesh;
	double* gpu_res;
	size_t pitch1;
	//size_t pitch2;

	/* Initializing cpu_mesh with source at the center*/

	for(int i=0 ; i< lx; i++){
		for(int j=0 ; j<ly ; ++j){
			cpu_mesh[i][j]=0.0;
		}
	}
	cpu_mesh[lx/2-1][ly/2 -1]=phi0;
        cpu_mesh[lx/2][ly/2 -1]=phi0;
        cpu_mesh[lx/2][ly/2]=phi0;
        cpu_mesh[lx/2-1][ly/2]=phi0;



	/* Allocation */
	gpuErrchk(cudaMallocPitch(&gpu_mesh, &pitch1, ly * sizeof(double), lx));
	gpuErrchk(cudaMemcpy2D(gpu_mesh, pitch1, cpu_mesh, ly*sizeof(double), ly*sizeof(double), lx, cudaMemcpyHostToDevice));
	dim3 gridSize(iDivUp(ly, BLOCKSIZE_x), iDivUp(lx, BLOCKSIZE_y));
	dim3 blockSize(BLOCKSIZE_y, BLOCKSIZE_x);
	solve <<<gridSize, blockSize>>> (gpu_mesh, pitch1, lx, D, dt, d);
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());

	gpuErrchk(cudaMemcpy2D(cpu_res, ly * sizeof(double), gpu_mesh, pitch1, ly * sizeof(double), lx, cudaMemcpyDeviceToHost));


	for (int i = 0; i < lx; i++){
		for (int j = 0; j < ly; j++){
      			std::cout << cpu_res[i][j] << " ";
		}
		std::cout << std::endl; 
	}


	return 0;
}
