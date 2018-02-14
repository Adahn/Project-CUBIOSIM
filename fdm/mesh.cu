#include <iostream>
#include <stdio.h>

class fixedFunction{
public:
 __host__ fixedFunction() {}
 __host__ __device__ double operator()(double x) {
    return x*x;
 }
};

__host__ __device__ double f1(double x){
  return x*x;
}

typedef double (*pf) (double var);

__device__ pf f1_d = f1;

class genericFunction{
public:
  __host__ genericFunction(double (*h_infunc)(double), double (*d_infunc)(double)) : h_func(h_infunc),d_func(d_infunc){}
  __host__ __device__ double operator()(double x) {
#ifdef __CUDA_ARCH__
    return d_func(x);
#else
    return h_func(x);
#endif
  }
private:
  pf h_func;
  pf d_func;
};

__global__ void kernel1(fixedFunction* g1){
  unsigned int tid = blockIdx.x *blockDim.x + threadIdx.x;
  printf("Func val is: %f\n", (*g1)(tid));
}

__global__ void kernel2(genericFunction* g1){
  unsigned int tid = blockIdx.x *blockDim.x + threadIdx.x;
  printf("Func val is: %f\n", (*g1)(tid));
}

int main(){

  fixedFunction h_g1;
  fixedFunction* d_g1;
  cudaMallocManaged(&d_g1, sizeof(h_g1));

  //Host call
  std::cout << h_g1(2.0) << "\n";

  //device call
  kernel1<<<1,32>>>(d_g1);
  cudaDeviceSynchronize();
  pf d_f1;
  cudaMemcpyFromSymbol(&d_f1, f1_d, sizeof(void*));
  genericFunction h_g2(f1, d_f1);
  genericFunction* d_g2;
  cudaMallocManaged(&d_g2, sizeof(h_g2));
  cudaMemcpy(d_g2, &h_g2, sizeof(h_g2), cudaMemcpyDefault);
  //Host call
  std::cout << h_g2(3.0) << "\n";

  //device call
  kernel2<<<1,32>>>(d_g2);
  cudaDeviceSynchronize();
}
