#include <iostream>

template<class state_type, class system>
state_type* rk4(int dim, system f, double t0, state_type* u0, double dt) {

	state_type* f0 = new state_type[dim];
	state_type* f1 = new state_type[dim];
	state_type* f2 = new state_type[dim];
	state_type* f3 = new state_type[dim];

	state_type* u1 = new state_type[dim];
	state_type* u2 = new state_type[dim];
	state_type* u3 = new state_type[dim];

	state_type* u = new state_type[dim];

	int i;
	double t1 = t0 + dt/2.0;
	double t2 = t0 + dt/2.0;
	double t3 = t0 + dt;

	//  Get four sample values of the derivative.
	// k1
	f(dim, u0, f0, t0);

	// k2
	for (i=0; i<dim; i++) {
		u1[i] = u0[i] + dt*f0[i]/2.0;
	}
	f(dim, u1, f1, t1);

	// k3
	for ( i = 0; i < dim; i++ ) {
		u2[i] = u0[i] + dt*f1[i]/2.0;
	}
	f(dim, u2, f2, t2);

	// k4
	for(i=0; i<dim; i++) {
		u3[i] = u0[i] + dt*f2[i];
	}
	f(dim, u3, f3, t3);

	//  Combine them to estimate the solution.
	for ( i = 0; i < dim; i++ ) {
		u[i] = u0[i] + dt*(f0[i] + 2.0*f1[i] + 2.0*f2[i] + f3[i])/6.0;
	}

	//  Free memory.
	delete[] f0;
	delete[] f1;
	delete[] f2;
	delete[] f3;
	delete[] u1;
	delete[] u2;
	delete[] u3;

	return u;
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

	// loop over time
	while(t_max > t0) {
		u1 = rk4<state_type>(dim, f, t0, u0, step);

		t0 += step;	c += 1;
		delete[] u0;	u0 = u1;

		if( (c%10000)==0 ) {
			f.observer(dim, u0, t0);
		}
	}
	
	delete[] u0;
}

