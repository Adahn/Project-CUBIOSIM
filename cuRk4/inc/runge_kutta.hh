#include <iostream>
#include <stdio.h>
#include <cuda_runtime.h>
#include <cuda.h>
#include <utility.hh>

#define gpuErrorCheck(ans) { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
/* CUDA error checking */
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


template<class state_type, class system>
state_type* rk4(int dim, system f, double t0, state_type* u0, double dt, double* coefs) {

	// Gpu memory alloc and parameters
	// state_matrix is by line :
	//	u0
	//	f0
	//	...
	//	f3
	state_type* state_matrix;
	cudaMalloc(&state_matrix, 5*dim*sizeof(state_type));
	cudaMemcpy(state_matrix, u0, dim*sizeof(state_type), cudaMemcpyHostToDevice);

	state_type *g_u;	// temporary variable
	cudaMalloc(&g_u, dim*sizeof(state_type));

	// TODO flexible dimensions, maybe in arguments
	int thread_per_block = 512;
	int blockSize = ceil((float)dim/thread_per_block);

	// CPU memory alloc
	state_type* u_sol = new state_type[dim];

	double t1 = t0 + dt/2.0;
	double t2 = t0 + dt/2.0;
	double t3 = t0 + dt;

	//  Get four sample values of the derivative.
	// k1 <=> f0
	f(dim, state_matrix, state_matrix+dim, t0);

	// k2 <=> f1
	sumk<state_type><<<blockSize, thread_per_block>>>(dim, g_u, state_matrix, coefs, 2);
	//cudaDeviceSynchronize();
	f(dim, g_u, state_matrix+2*dim, t1);

	// k3 <=> f2
	sumk<state_type><<<blockSize, thread_per_block>>>(dim, g_u, state_matrix, coefs+5, 3);
	//cudaDeviceSynchronize();
	f(dim, g_u, state_matrix+3*dim, t2);

	// k4 <=> f3
	sumk<state_type><<<blockSize, thread_per_block>>>(dim, g_u, state_matrix, coefs+2*5, 4);
	//cudaDeviceSynchronize();
	f(dim, g_u, state_matrix+4*dim, t3);

	//  Combine them to estimate the solution.
	sumk<state_type><<<blockSize, thread_per_block>>>(dim, g_u, state_matrix, coefs+3*5, 5);
	cudaDeviceSynchronize();
	cudaMemcpy(u_sol, g_u, dim*sizeof(state_type), cudaMemcpyDeviceToHost);

	//  Free memory.
	cudaFree(g_u);	cudaFree(state_matrix);

	return u_sol;

}


template<class state_type, class system>
void rk4_wrapper(int dim, system f, state_type* initial_u,
        double t0, double t_max, double step) {

	state_type* u0 = new state_type[dim];
	state_type* u1;
	int c = 0;

	// copy initial point to avoid erasing it
	for(int i=0; i<dim; i++) {
		u0[i] = initial_u[i];
	}
	f.observer(dim, u0, t0);

	// load coefficients array on GPU, to be fed to the sumk primitive
	// each line corresponds to the coef. of one term un the runge kutta method
	// for ex 1 row represents coef [1, step/2, 0, 0, 0] of this term:
	//		f1 = u0 + dt*f0/2
	double* g_coefs;	// coefficients array
	cudaMalloc(&g_coefs, 4*5*sizeof(double));
	load_coef(g_coefs, step);

	// loop over time
	while(t_max > t0) {

		u1 = rk4<state_type, system>(dim, f, t0, u0, step, g_coefs);

		t0 += step;		c += 1;
		delete[] u0;	u0 = u1;

		if( (c%10000)==0 ) {
			f.observer(dim, u0, t0);
			cout << "t=" << t0 << endl;
		}
	}

	delete[] u0;
	cudaFree(g_coefs);
}

