#include "./inc/utility.hh"
#include <iostream>

#define N 16
#define N_V 2

using namespace std;

int main() {

	// cpu memory
	cout << "cpu memory..." << flush;
	int v[N], out[N], coef[N_V];
	
	for(int i=0; i<N; i++) {
		v[i] = i;
	}
	
	for(int i=0; i<N_V; i++) {
		coef[i] = i+1;
	}
	
	cout << "done" << endl;

	// gpu memory
	cout << "gpu memory..." << flush;
	int *gpu_v,  *g_coef, *gpu_out;
	cudaMalloc(&gpu_out, sizeof(int*)*N);
	cudaMalloc(&gpu_v, sizeof(int*)*2*N);
	cudaMalloc(&g_coef, sizeof(int)*N_V);
	cout << "done" << endl;

	// H2D
	cout << "H2D..." << flush;
	for(int i=0; i<N_V; i++ ) {
		cudaMemcpy(gpu_v+i*N, v, sizeof(int)*N, cudaMemcpyHostToDevice);
	}
	cudaMemcpy(g_coef, &coef, sizeof(int)*N_V, cudaMemcpyHostToDevice);
	cout << "done" << endl;

	// kernel call
	cout << "kernel call...\n" << flush;
	sumk<int><<<2, 4>>>(N, gpu_out, gpu_v, g_coef, N_V);
	cudaDeviceSynchronize();
	cout << "done" << endl;

	// D2H
	cout << "D2H..." << flush;
	cudaMemcpy(out, gpu_out, sizeof(int)*N, cudaMemcpyDeviceToHost);
	cout << "done" << endl;

	// display	
	for(int i=0; i<N; i++) {
		cout << out[i] << "\t";
	}
	cout << endl;

	return 0;
	
}
