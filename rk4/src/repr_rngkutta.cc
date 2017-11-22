#include <iostream>
#include <fstream>
#include <cmath>

#include <chTimer.hpp>

#include <runge_kutta.hh>

#define NSPECIES 8
#define NREACTIONS 7

using namespace std;

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

void Repressilator_ODE( int dim, double* Y, double* dYdt, double t ) {

	// Reaction Rates
	double R[NREACTIONS] = {Ktl*Y[4],			// Reaction A
				Ktl*Y[5] ,			// Reaction B
				Ktl*Y[6] ,			// Reaction C
				Ktl*Y[7] ,			// Reaction D
				Ktr/(1 + pow(Y[1]/KR, nR)) ,	// Reaction E
				Ktr/(1 + pow(Y[2]/KR, nR)) ,	// Reaction F
				Ktr/(1 + pow(Y[0]/KR, nR)) };	// Reaction G

	// Model: dYdt = S * R - d * Y
	for (int i = 0; i < NSPECIES; ++i) {

		dYdt[i] = 0;
		for (int j = 0; j < NREACTIONS; ++j) {
			dYdt[i] += S[j+i*NREACTIONS] * R[j];
		}
		dYdt[i] -= decay_v[i]*Y[i];

	}

}

void test( int dim, double* Y, double* dYdt, double t ) {
	int i;
	for (i = 0; i < NSPECIES/2; ++i) {
		dYdt[i] = 1;
	}

	for(; i < NSPECIES; ++i) {
		dYdt[i] = 1-Y[i];
	}

}

void write_ODE_result_rk4(int dim, double Y[], double t ) {
	ofstream file;
	file.open ("bin/result_rk4.csv",ios::app); // write at the end of file
	
	file << t;
	for(int i=0; i<NSPECIES; i++) {
		file << ";" << Y[i];	
	}
	file << "\n";
	
	file.close();
}

void write_ODE_result_rk45(int dim, double Y[], double t ) {
	ofstream file;
	file.open ("bin/result_rk45.csv",ios::app); // write at the end of file
	
	file << t;
	for(int i=0; i<NSPECIES; i++) {
		file << ";" << Y[i];	
	}
	file << "\n";
	
	file.close();
}


int main(int argc, char **argv) {
	
	ChTimer rk4, rk45;

	double X0[] = { 1, 0, 0, 0, 0, 0, 0, 0 }; // initial condition
	double eps = 1e-10;

	// Runge-Kutta method
	/*rk4.start();
	rk4_wrapper( NSPECIES, Repressilator_ODE , X0 , 0.0 , 100000 , 1e-2,
			write_ODE_result_rk4 );
	rk4.stop();*/
			
	// adaptive step size Runge-Kutta-Fehlberg method
	rk45.start();
	rk45_wrapper( NSPECIES, Repressilator_ODE , X0 , 0.0 , 100000 , 1e-2,
			write_ODE_result_rk45, eps );
	rk45.stop();

	cout << "rk4:\t" << rk4.getTime() << "s" << endl
		<< "rk45:\t" << rk45.getTime() << "s" << endl;
}

