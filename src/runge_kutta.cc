#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <new>

using namespace std;

#include "runge_kutta.hh" 


double* rk4(double t0, double u0[], double dt, int dim, void f (double t, int dim, double u[], double dudt[]))
// 	  takes one Runge-Kutta step for a scalar ODE
//    It is assumed that an initial value problem, of the form
//
//      du/dt = f ( t, u )
//      u(t0) = u0
//
//    is being solved.
//
//    If the user can supply current values of t, u, a stepsize dt, and a
//    function to evaluate the derivative, this function can compute the
//    fourth-order Runge Kutta estimate to the solution at time t+dt.
//
//
//  Parameters:
//
//	 Input
//    - double T0, the current time.
//	  - int dim, the spatial dimension.
//	  - double U0[dim], the solution estimate at the current time.
//	  - double dt, the time step.
//	  - double *F ( double T, int dim, double U[] ), a function which evaluates
//                           the derivative, or right hand side of the problem.
//
//   Output, 
//	  - double RK4VEC[dim], the fourth-order Runge-Kutta solution estimate
//                        at time T0+DT.
//
{
  /*double *f0 = new double(dim);
  double *f1 = new double(dim);
  double *f2 = new double(dim);
  double *f3 = new double(dim);*/

  double *f0 = new double[dim];
  double *f1 = new double[dim];
  double *f2 = new double[dim];
  double *f3 = new double[dim];
  double *u  = new double[dim];
  double *u1 = new double[dim];
  double *u2 = new double[dim];
  double *u3 = new double[dim];

  int i;
  double t1 = t0 + dt / 2.0;
  double t2 = t0 + dt / 2.0;
  double t3 = t0 + dt;

//
//  Get four sample values of the derivative.
//
  // k1
  f ( t0, dim, u0, f0 );

  // k2
  for ( i = 0; i < dim; i++ ) {
    u1[i] = u0[i] + dt * f0[i] / 2.0;
  }
  f ( t1, dim, u1, f1 );

  // k3
  for ( i = 0; i < dim; i++ ) {
    u2[i] = u0[i] + dt * f1[i] / 2.0;
  }
  f ( t2, dim, u2, f2 );

  // k4
  for ( i = 0; i < dim; i++ ) {
     u3[i] = u0[i] + dt * f2[i];
  }
  f ( t3, dim, u3, f3 );

//
//  Combine them to estimate the solution.
//
  for ( i = 0; i < dim; i++ )
  {
     u[i] = u0[i] + dt * ( f0[i] + 2.0 * f1[i] + 2.0 * f2[i] + f3[i] ) / 6.0;
  }
//
//  Free memory.
//
  delete [] f0;
  delete [] f1;
  delete [] f2;
  delete [] f3;
  delete [] u1;
  delete [] u2;
  delete [] u3;

  return u;
}
//****************************************************************************80

/*
void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}*/