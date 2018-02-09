#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include <cuda_runtime.h>
#include "cublas_v2.h"

#include <utility.hh>

using namespace std;

// This kernel computes reactions rates
__global__ void computeR(double* R, double* Y, double Ktl, double Ktr, double KR, double nR,
						int n) {
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	int nreac = 2*n+1;

	if(tid<n+1) {
		R[tid] = Ktl*Y[tid+n+1];
	}
	else if(tid<2*n) {
		R[tid] = Ktr/(1 + pow(Y[tid-n]/KR, nR));
	}
	else if(tid<nreac) {
		R[tid] = Ktr/(1 + pow(Y[0]/KR, nR));
	}
}

// This kernel computes the decay, a component-wise product
__global__ void computeDecay(double* res, double* Y, double* decay, int n) {
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	int nspec = 2*(n+1);
	
	if(tid<nspec) {
		res[tid] = decay[tid]*Y[tid];
	}
}


class Repressilator_ODE
{
	private:
		// size
		int _n;
		int _nreac;
		int _nspec;

		// system parameters
		double* R;	// reaction rates
		double* decay;	// decay vector
		double* S;	// stochiometric matrix

		// Rate parameters
		double _Ktl;
		double _Ktr;
		double _KR;
		double _nR;
		
		// other variables
		cublasHandle_t _handle;
		cublasStatus_t _status;
		string _filename;
		
	public:	

		Repressilator_ODE(int n, double dprot, double dmRNA, double Ktl, double Ktr, double KR, double nR, string filename):
			_n(n), _nreac( 2*n+1 ), _nspec( 2*(n+1) ),
			_Ktl(Ktl), _Ktr(Ktr), _KR(KR), _nR(nR),
			_filename(filename)
		{
			
			// initialize cuBLAS
			_status = cublasCreate(&_handle);
			
			// gpu memory alloc
			cudaMalloc(&R, _nspec*sizeof(double));
			cudaMalloc(&S, _nspec*_nreac*sizeof(double));
			cudaMalloc(&decay, _nspec*sizeof(double));
			
			// fill system matrices
			double* tmp_S = (double*)calloc( _nspec * _nreac, sizeof(double) );
			double* tmp_decay = (double*)calloc(_nspec, sizeof(double));
	
			// Assemble S and decay on CPU
			// S is stored in column major format (equivalent to S transpose)
			// TODO use sparse format for S
			int i,j;
			for( i=0; i<_n+1; i++ ) {
				tmp_decay[i] = dprot;
				tmp_S[i + i*_nspec] = 1;
			}
			for( j=i; j<2*_n+1; j++ ) {		
				tmp_decay[j] = dmRNA;
				tmp_S[j + j*_nspec] = 1;
			}
			tmp_decay[j] = dmRNA;
			tmp_S[_nspec*_nreac-1] = 1;
			
			// copy H2D
			cudaMemcpy(S, tmp_S, _nspec*_nreac*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(decay, tmp_decay, _nspec*sizeof(double), cudaMemcpyHostToDevice);
			
			free(tmp_S);
			free(tmp_decay);
		}
	
		// Overloading operator () to compute the derivate	
		void operator()( const int dim, double* Y, double* dYdt, const double t ) {

			int thread_per_block = 512;
			int blockSize = ceil((float)dim/thread_per_block);
			computeR<<<blockSize,thread_per_block>>>(R, Y, _Ktl, _Ktr, _KR, _nR, _n);

			// Compute d*Y
			computeDecay<<<blockSize,thread_per_block>>>(dYdt, Y, decay, _n);

			// Model: dYdt = S*R - d*Y
			//TODO use cuSparse 
			double alpha = 1; double beta = -1;
			_status = cublasDgemv(_handle, CUBLAS_OP_N, _nspec, _nreac, &alpha, S, _nspec,
				R, 1, &beta, dYdt, 1);
				
		}

		// This function stores the state given in argument
		void observer( const int dim, double* Y, const double t )
		{
			ofstream file;
			file.open(_filename,ios::app); // write at the end of file
	
			file << t;
			for( int i=0; i<dim; i++) {
				file << ";" << Y[i];	
			}
			file << "\n";
	
			file.close();
		}

	~Repressilator_ODE() {
		cublasDestroy(_handle);
		cudaFree(S);	cudaFree(R);	cudaFree(decay);
	}
};

