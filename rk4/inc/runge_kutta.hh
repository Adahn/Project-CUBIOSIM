// Normaly no need to modify this file with the runge kutta solvers

#include <iostream>
#include "utility.hh"

/* Solver Runge-Kutta 4 method */
template<class state_type, class system>
state_type* rk4(int dim, system f, double t0, state_type* u0, double dt) {

	// define the 4 intermediate calculation terms of the RK4 method
	state_type* k1 = new state_type[dim];
	state_type* k2 = new state_type[dim];
	state_type* k3 = new state_type[dim];
	state_type* k4 = new state_type[dim];

	state_type* tmp = new state_type[dim];

	state_type* u = new state_type[dim];

	int i;
	double t1 = t0 + dt/2.0;
	double t2 = t0 + dt/2.0;
	double t3 = t0 + dt;

	//  Get four sample values of the derivative.
	// k1
	f(dim, u0, k1, t0);

	// k2
	for (i=0; i<dim; i++) {
		tmp[i] = u0[i] + dt*k1[i]/2.0;
	}
	f(dim, tmp, k2, t1);

	// k3
	for ( i = 0; i < dim; i++ ) {
		tmp[i] = u0[i] + dt*k2[i]/2.0;
	}
	f(dim, tmp, k3, t2);

	// k4
	for(i=0; i<dim; i++) {
		tmp[i] = u0[i] + dt*k3[i];
	}
	f(dim, tmp, k4, t3);

	//  Combine them to estimate the solution.
	for ( i = 0; i < dim; i++ ) {
		u[i] = u0[i] + dt*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0;
	}

	//  Free memory.
	delete[] k1;	delete[] k2;
	delete[] k3;	delete[] k4;
	delete[] tmp;

	return u;
}

/* Wrapper calling the RK4 method every time step */
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

	// break when end time is reached
	while(t_max > t0) {
		// call rk4 solver
		u1 = rk4<state_type, system>(dim, f, t0, u0, step);

		t0 += step;		c += 1;
		delete[] u0;	u0 = u1;

		// write results every 10000's iteration
		if( (c%10000)==0 ) {
			f.observer(dim, u0, t0);
		}
	}

	delete[] u0;
}


/* adaptive runge-kutta method of order 4 and 5
https://math.okstate.edu/people/yqwang/teaching/math4513_fall11/Notes/rungekutta.pdf
*/
template<class state_type, class system>
state_type* rk45(int dim, system f, double t0, state_type* u0, double dt,
	double eps, state_type *R, double *delta) {

	// define the 6 intermediate calculation terms of the RK4 method
	state_type *k1 = new state_type[dim];
	state_type *k2 = new state_type[dim];
	state_type *k3 = new state_type[dim];
	state_type *k4 = new state_type[dim];
	state_type *k5 = new state_type[dim];
	state_type *k6 = new state_type[dim];

	state_type *tmp = new state_type[dim];

	state_type *w1 = new state_type[dim];
	state_type *w2 = new state_type[dim];

	int i;
	double t1 = t0 + dt/4.0;
	double t2 = t0 + 3*dt/8.0;
	double t3 = t0 + 12*dt/13.0;
	double t4 = t0 + dt;
	double t5 = t0 + dt/2;

	//  Get four sample values of the derivative.
	// k1
	f(dim, u0, k1, t0);

	// k2
	for (i=0; i<dim; i++) {
		tmp[i] = u0[i] + dt*k1[i]/4.0;
	}
	f(dim, tmp, k2, t1);

	// k3
	for ( i = 0; i < dim; i++ ) {
		tmp[i] = u0[i] + 3*dt*k1[i]/32.0 + 9*dt*k2[i]/32.0;
	}
	f(dim, tmp, k3, t2);

	// k4
	for(i=0; i<dim; i++) {
		tmp[i] = u0[i] + dt/2197 * (1932*k1[i] - 7200*k2[i] + 7296*k3[i]);
	}
	f(dim, tmp, k4, t3);

	// k5
	for(i=0; i<dim; i++) {
		tmp[i] = u0[i] + dt * (439.0/216*k1[i] - 8*k2[i] + 3680/513*k3[i] - 845/4104*k4[i]);
	}
	f(dim, tmp, k5, t4);

	// k6
	for(i=0; i<dim; i++) {
		tmp[i] = u0[i] + dt * (-8.0/27*k1[i] + 2*k2[i] - 3544/2565*k3[i] + 1859/4104*k4[i] - 11/40*k5[i]);
	}
	f(dim, tmp, k6, t5);

	//  Combine them to estimate the solution.
	for ( i = 0; i < dim; i++ ) {
		w1[i] = u0[i] + dt*(25.0/216*k1[i] + 1408.0/2565*k3[i] + 2197.0/4104*k4[i] - k5[i]/5);
		w2[i] = u0[i] + dt*(16.0/135*k1[i] + 6656.0/12825*k3[i] + 28561.0/56430*k4[i] - 9.0/50*k5[i] + 2.0/55*k6[i]);
	}

	// calculate |w1 - w2| ; norm L2
	*R = 0;
	for (int i = 0; i < dim; ++i)
	{
		*R += (w1[i] - w2[i])*(w1[i] - w2[i]);
	}

	*R = 1.0/dt * sqrt(*R);

	*delta = 0.84* pow(eps/(*R), 1.0/4.0); // step

	//  Free memory.
	delete[] k1;	delete[] k2;
	delete[] k3;	delete[] k4;
	delete[] k5;	delete[] k6;
	delete[] tmp;	delete[] w2;

	return w1;
}

/* wrapps runge kutta Fehlberg method (RK45) */
template<class state_type, class system>
void rk45_wrapper(int dim, system f, state_type* initial_u,
        double t0, double t_max, double step, double eps) {

	state_type* u0 = new state_type[dim];
	state_type* u1;

	state_type R;
	double delta;

	// copy initial point to avoid erasing it
	for(int i=0; i<dim; i++) {
		u0[i] = initial_u[i];
	}
	f.observer(dim, u0, t0);

	// breaks when end time is reached
	while(t_max > t0) {
		// calling rk45 solver
		u1 = rk45<state_type, system>(dim, f, t0, u0, step, eps, &R, &delta );

		if ( R <= eps ) {
			t0 += step;		step=delta*step;
			delete[] u0;	u0 = u1;
		}
		else {
			step = delta*step;
		}

		f.observer(dim, u0, t0);
	}

	delete[] u0;
}
