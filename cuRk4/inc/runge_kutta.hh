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

// computes and returns the next state of the system from u0=u(t0)
// with time step dt
// u0 is assumed to be on the CPU, the returned vector is also on CPU
template<class state_type, class system>
state_type* rk4(int dim, system f, double t0, state_type* u0, double dt, double* coefs) {

	// Gpu memory alloc and parameters
	// state_matrix is by line :
	//	u0
	//	k1
	//	...
	//	k4
	state_type* state_matrix;
	cudaMalloc(&state_matrix, 5*dim*sizeof(state_type));
	cudaMemcpy(state_matrix, u0, dim*sizeof(state_type), cudaMemcpyHostToDevice);

	state_type *g_u;	// temporary variable on the GPU
	cudaMalloc(&g_u, dim*sizeof(state_type));

	// TODO flexible dimensions for the kernel call
	int thread_per_block = 512;
	int blockSize = ceil((float)dim/thread_per_block);

	// CPU memory alloc
	state_type* u_sol = new state_type[dim]; // solution to be returned

	double t1 = t0 + dt/2.0;
	double t2 = t0 + dt/2.0;
	double t3 = t0 + dt;

	//  Get four sample values of the derivative.
	// k1
	f(dim, state_matrix, state_matrix+dim, t0);

	// k2
	sumk<state_type><<<blockSize, thread_per_block>>>(dim, g_u, state_matrix, coefs, 2);
	f(dim, g_u, state_matrix+2*dim, t1);

	// k3
	sumk<state_type><<<blockSize, thread_per_block>>>(dim, g_u, state_matrix, coefs+5, 3);
	f(dim, g_u, state_matrix+3*dim, t2);

	// k4
	sumk<state_type><<<blockSize, thread_per_block>>>(dim, g_u, state_matrix, coefs+2*5, 4);
	f(dim, g_u, state_matrix+4*dim, t3);

	//  Combine them to estimate the solution.
	sumk<state_type><<<blockSize, thread_per_block>>>(dim, g_u, state_matrix, coefs+3*5, 5);
	cudaDeviceSynchronize();
	cudaMemcpy(u_sol, g_u, dim*sizeof(state_type), cudaMemcpyDeviceToHost);

	//  Free memory.
	cudaFree(g_u);	cudaFree(state_matrix);

	return u_sol;

}


// This function wraps the rk4 numeric scheme to iterate over the timespan
// [t0,t_max] with step "step"
// the initial_u is assumed to be on the CPU
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
	// each line corresponds to the coefs of one term of the runge kutta method
	// example: row 1 represents coef [1, step/2, 0, 0, 0] of this term:
	//		k2 = u0 + dt*k1/2
	double* g_coefs;
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

