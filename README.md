# Project-CUBIOSIM
University project aiming to improve the performance in the calculation of biological systems

## File Structure
    - cuRk4 : Parallelised version of the runge-kutta method (RK4) using CUDA, check out the README.txt in this directory to get started
    - fdm : TODO LILIANE
    - odeint : Makefile2 to execute thrust+odeint (repressilator_thrust.cc in /src (n species in /src/repressilator_n_thrust.cc))
    - read : Read.file.cc (to read resp. the dimensions, the coefficients of the vector of concentrations, the matrix S, and the vector of degradation)
    - rk4 : Sequential version of runge-kutta method (RK4) to solve biological systems, check out the README.txt in this directory to get started
    - sources: Collection of helpful articles, presentations, etc.

## Requirements
CUDA requirements: In order to use the programs, one needs a CUDA compatible GPU with compute capability 3.5 or higher, and the appropriate driver. We are also using the cuBLAS library from the CUDA toolkit available at https://developer.nvidia.com/cuda-toolkit.

Other requirements to use your programs 
odeint: Library boost/odeint

## Getting started
If you want to check the sequential version of the Runge-Kutta method solvers, check out the directory rk4 and read the README.txt

If you want to check out the parallelised version of the Runge-Kutta method solver, check out the directory cuRk4 and read the README.txt

If you want to check out the parallelised version of odeint, check out the directory odeint and the command make -f Makefile2


