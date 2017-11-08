#include <cstdlib>
#include <iostream>
#include <math.h>


using namespace std;

void rk4_wrapper(int dim, void f(int dim, double u[], double dudt[], double t),
	double initial_u[], double t0, double t_max, double step,
	void observer(int dim, double u[], double t)=NULL);

void rk45_wrapper(int dim, void f(int dim, double u[], double dudt[], double t),
	double initial_u[], double t0, double t_max, double step,
	void observer(int dim, double u[], double t), double eps);

double* rk4(int dim, double t0, double u0[], double dt,
        void f(int dim, double u[], double dudt[], double t) );

double* rk45(int dim, double t0, double u0[], double dt,
	void f(int dim, double u[], double dudt[], double t), double eps, double* R, double* delta );
