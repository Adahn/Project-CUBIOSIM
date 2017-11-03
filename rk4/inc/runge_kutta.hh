#include <cstdlib>
#include <iostream>

using namespace std;

void rk4_wrapper(int dim, void f(int dim, double u[], double dudt[], double t),
	double initial_u[], double t0, double t_max, double step,
	void observer(int dim, double u[], double t)=NULL);

double* rk4(int dim, double t0, double u0[], double dt,
        void f(int dim, double u[], double dudt[], double t) );
