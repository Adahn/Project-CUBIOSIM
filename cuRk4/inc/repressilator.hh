#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include <cuda_runtime.h>
#include "cublas_v2.h"

using namespace std;

// tested successfully
__global__ void computeR(double* R, double* Y, double Ktl, double Ktr, double KR, double nR,
						int n) {
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	int nspec = 2*(n+1);

	//for( int i=0; i<_n+1; i++ )
	if(tid<n+1) {
		R[tid] = Ktl*Y[tid+n+1];
	}
	//for( int i=1; i<_n; i++ )
	else if(tid<2*n+1) {
		R[tid] = Ktr/(1 + pow(Y[tid-n]/KR, nR));
	}
	else if(tid<nspec) {
		R[2*n+1] = Ktr/(1 + pow(Y[0]/KR, nR));
	}
}

// tested successfully
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
		double* R;
		double* decay;
		double* S;

		// Rate parameters
		double _Ktl;
		double _Ktr;
		double _KR;
		double _nR;
		
		// other variables
		cublasHandle_t _handle;
		cublasStatus_t _status;
		string _filename;
		
			// Model: dYdt = S*R - decay*Y
			/*for( int i=0; i<_nspec; i++) {
				dYdt[i] = 0;
				for( int j=0; j<_nreac; j++) {
					dYdt[i] += S[i *_nreac + j] * R[j];
				}
				dYdt[i] -= decay[i]*Y[i];
			}*/
			
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
	
			int i,j;
			// S is stored in column major format (equivalent to S transpose)
			for( i=0; i<_n+1; i++ ) {
				tmp_decay[i] = dprot;	// prot
				tmp_S[i + i*_nspec] = 1;	// S looks like Identity
			}
			for( j=i; j<2*_n+1; j++ ) {		
				tmp_decay[j] = dmRNA;	// RNA
				tmp_S[j + j*_nspec] = 1;	// S looks like Identity
			}
			tmp_decay[j] = dmRNA;		// last RNA
			tmp_S[_nspec*_nreac-1] = 1;		// S last line is the same as the previous one
			
			// copy H2D
			cudaMemcpy(S, tmp_S, _nspec*_nreac*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(decay, tmp_decay, _nspec*sizeof(double), cudaMemcpyHostToDevice);
			
			delete[] tmp_S;
			delete[] tmp_decay;
		}
		
		void operator()( const int dim, double* Y, double* dYdt, const double t ) {

			int thread_per_block = 512;
			int blockSize = ceil((float)dim/thread_per_block);
			computeR<<<blockSize,thread_per_block>>>(R, Y, _Ktl, _Ktr, _KR, _nR, _n);

			// Compute d*Y
			computeDecay<<<blockSize,thread_per_block>>>(dYdt, Y, decay, _n);

			// Model: dYdt = S*R - d*Y
			//TODO sparse ?
			double alpha = 1; double beta = -1;
			_status = cublasDgemv(_handle, CUBLAS_OP_N, _nspec, _nreac, &alpha, S, _nspec,
				R, 1, &beta, dYdt, 1);
				
		}

		// utility functions
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

