#include <iostream>

#include <boost/numeric/odeint.hpp>

#include <fstream>

#include <runge_kutta.hh>
//
#define NSPECIES 8
#define NREACTIONS 7

using namespace std;
using namespace boost::numeric::odeint;

// Parameter definition
const double Ktl = 1e-2;
const double Ktr = 1e-3;
const double KR = 1e-2;
const double nR = 2;
const double dprot = 1e-3;
const double dmRNA = 1e-2;

// decay vector
const double decay_v[NSPECIES] = { dprot,	// Species 1
						  dprot,	// Species 2 
						  dprot,	// Species 3
						  dprot,	// Species 4
						  dmRNA,	// Species 5
						  dmRNA,	// Species 6
						  dmRNA,	// Species 7
						  dmRNA };	// Species 8


// Stochiometric matrix
int S[NSPECIES*NREACTIONS] = {	1,  0,  0,  0,  0,  0,  0,     // Species 1
				0,  1,  0,  0,  0,  0,  0,     // Species 2
				0,  0,  1,  0,  0,  0,  0,     // Species 3
				0,  0,  0,  1,  0,  0,  0,     // Species 4
				0,  0,  0,  0,  1,  0,  0,     // Species 5
				0,  0,  0,  0,  0,  1,  0,     // Species 6
				0,  0,  0,  0,  0,  0,  1,     // Species 7
				0,  0,  0,  0,  0,  0,  1 };   // Species 8
		// Reaction     A   B   C   D   E   F   G

void Repressilator_ODE( double t, int dim, double* Y, double* dYdt ) {

	// Reaction Rates
	double R[NREACTIONS] = {Ktl*Y[4] , 					// Reaction A
					Ktl*Y[5] ,				// Reaction B
					Ktl*Y[6] ,				// Reaction C
					Ktl*Y[7] ,				// Reaction D
					Ktr/(1 + pow(Y[1]/KR, nR)) ,	// Reaction E
					Ktr/(1 + pow(Y[2]/KR, nR)) ,	// Reaction F
					Ktr/(1 + pow(Y[0]/KR, nR)) };	// Reaction G

	// Model: dYdt = S * R - d * Y
	for (int i = 0; i < NSPECIES; ++i)
	{
		dYdt[i] = 0;
		for (int j = 0; j < NREACTIONS; ++j)
		{
			dYdt[i] += S[j+i*NREACTIONS] * R[j];

		}
	
		dYdt[i] -= decay_v[i]*Y[i];
	}

}

int main(int argc, char **argv) {
	
	double X0[] = { 1e-6, 0, 0, 0, 0, 0, 0, 0 }; // initial condition

	// result matrix
	rk4_wrapper( Repressilator_ODE , X0 , 0.0 , 100000.0 , 1e-2, NSPECIES );

}

