#include <runge_kutta.hh>
//#include <complex> // std::norm

/* wrapps runge kutta 4 method */
void rk4_wrapper(int dim, void f(int dim, double u[], double dudt[], double t),
	double initial_u[], double t0, double t_max, double step,
	void observer(int dim, double u[], double t)) { 

	double* u0 = new double[dim];
	double* u1;
	int c = 0;

	cout << "Runge Kutta Order 4..." << flush;

	// copy initial point to avoid erasing it
	for(int i=0; i<dim; i++) {
		u0[i] = initial_u[i];
	}
	observer(dim, u0, t0);

	// loop over time
	while(t_max > t0) {
		u1 = rk4(dim, t0, u0, step, f );

		t0 += step;	
		c += 1;
		delete[] u0;	u0 = u1;

		if( observer!=NULL && (c%10000)==0 ) {
			observer(dim, u0, t0);
		}
	}
	
	delete[] u0;
	cout << "done" << endl;
}


/* Rung-kutta method order 4 */
double* rk4(int dim, double t0, double u0[], double dt,
	void f(int dim, double u[], double dudt[], double t) ) {

	double *f0 = new double[dim];
	double *f1 = new double[dim];
	double *f2 = new double[dim];
	double *f3 = new double[dim];

	double *u1 = new double[dim];
	double *u2 = new double[dim];
	double *u3 = new double[dim];

	double *u = new double[dim];


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


/* wrapps runge kutta Fehlberg method */
void rk45_wrapper(int dim, void f(int dim, double u[], double dudt[], double t),
	double initial_u[], double t0, double t_max, double step,
	void observer(int dim, double u[], double t), double eps) { 

	double* u0 = new double[dim];
	double* u1;
	int c = 0;
	double* R;
	double* delta;

	cout << "Runge Kutta Fehlberg Order 45..." << flush;

	// copy initial point to avoid erasing it
	for(int i=0; i<dim; i++) {
		u0[i] = initial_u[i];
	}
	observer(dim, u0, t0);



	// loop over time
	while(t_max > t0) {
		u1 = rk45(dim, t0, u0, step, f, eps, R, delta );

		if ( *R <= eps )
		{
			t0 += step;	
			delete[] u0;
			u0 = u1; 
			c += 1;
			step = *delta*step;
		}
		else
		{
			step = *delta*step;
		}

		if( observer!=NULL && (c%10000)==0 ) {
			observer(dim, u0, t0);
		}
	}
	
	delete[] u0;
	cout << "done" << endl;
}



/* adaptive runge-kutta method of order 4 and 5
https://math.okstate.edu/people/yqwang/teaching/math4513_fall11/Notes/rungekutta.pdf
*/
double* rk45(int dim, double t0, double u0[], double dt,
	void f(int dim, double u[], double dudt[], double t), double eps, double* R, double* delta ) {

	double *f0 = new double[dim];
	double *f1 = new double[dim];
	double *f2 = new double[dim];
	double *f3 = new double[dim];
	double *f4 = new double[dim];
	double *f5 = new double[dim];


	double *u1 = new double[dim];
	double *u2 = new double[dim];
	double *u3 = new double[dim];
	double *u4 = new double[dim];
	double *u5 = new double[dim];

	double *w1 = new double[dim];
	double *w2 = new double[dim];


	int i;
	double t1 = t0 + dt/4.0;
	double t2 = t0 + 3*dt/8.0;
	double t3 = t0 + 12*dt/13.0;
	double t4 = t0 + dt;
	double t5 = t0 + dt/2;

	//  Get four sample values of the derivative.
	// k1
	f(dim, u0, f0, t0);

	// k2
	for (i=0; i<dim; i++) {
		u1[i] = u0[i] + dt*f0[i]/4.0;
	}
	f(dim, u1, f1, t1);

	// k3
	for ( i = 0; i < dim; i++ ) {
		u2[i] = u0[i] + 3*dt*f0[i]/32.0 + 9*dt*f1[i]/32.0;
	}
	f(dim, u2, f2, t2);

	// k4
	for(i=0; i<dim; i++) {
		u3[i] = u0[i] + dt/2197 * (1932*f0[i] - 7200*f1[i] + 7296*f2[i]);
	}
	f(dim, u3, f3, t3);

	// k5 
	for(i=0; i<dim; i++) {
		u4[i] = u0[i] + dt * (439.0/216*f0[i] - 8*f1[i] + 3680/513*f2[i] - 845/4104*f3[i]);
	}
	f(dim, u4, f4, t4);

	// k6
	for(i=0; i<dim; i++) {
		u5[i] = u0[i] + dt * (-8.0/27*f0[i] + 2*f1[i] - 3544/2565*f2[i] + 1859/4104*f3[i] - 11/40*f4[i]);
	}
	f(dim, u5, f5, t5);

	//  Combine them to estimate the solution.
	for ( i = 0; i < dim; i++ ) {
		w1[i] = u0[i] + dt*(25.0/216*f0[i] + 1408.0/2565*f2[i] + 2197.0/4104*f3[i] - f4[i]/5);
		w2[i] = u0[i] + dt*(16.0/135*f0[i] + 6656.0/12825*f2[i] + 28561.0/56430*f3[i] - 9.0/50*f4[i] + 2.0/55*f5[i]);
	}

	// calculate |w1 - w2| ; norm L2
	for (int i = 0; i < dim; ++i)
	{
		*R += (w1[i] - w2[i])*(w1[i] - w2[i]);
	}
	
	*R = 1.0/dt * sqrt(*R);

	delta = 0.84* pow(eps/R, 1.0/4.0); // step

	//  Free memory.
	delete[] f0;
	delete[] f1;
	delete[] f2;
	delete[] f3;
	delete[] f4;
	delete[] f5;
	delete[] u1;
	delete[] u2;
	delete[] u3;
	delete[] u4;
	delete[] u5;
	delete[] w2;

	return w1;
}
