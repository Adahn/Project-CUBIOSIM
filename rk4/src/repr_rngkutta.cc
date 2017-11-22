#include <iostream>
#include <fstream>
#include <cmath>

#include <runge_kutta.hh>

#define NSPECIES 8
#define NREACTIONS 7

//using namespace std;

class Repressilator_ODE {

  private:
	// Parameter definition
	const double Ktl = 1e-2;
	const double Ktr = 1e-3;
	const double KR = 1e-2;
	const double nR = 2;
	const double dprot = 1e-3;
	const double dmRNA = 1e-2;

	// decay vector
	const double decay_v[NSPECIES] =
			{ dprot,	// Species 1
			  dprot,	// Species 2 
			  dprot,	// Species 3
			  dprot,	// Species 4
			  dmRNA,	// Species 5
			  dmRNA,	// Species 6
			  dmRNA,	// Species 7
			  dmRNA };	// Species 8

	// Stochiometric matrix
	int S[NSPECIES*NREACTIONS] =
		{	1,  0,  0,  0,  0,  0,  0,     // Species 1
			0,  1,  0,  0,  0,  0,  0,     // Species 2
			0,  0,  1,  0,  0,  0,  0,     // Species 3
			0,  0,  0,  1,  0,  0,  0,     // Species 4
			0,  0,  0,  0,  1,  0,  0,     // Species 5
			0,  0,  0,  0,  0,  1,  0,     // Species 6
			0,  0,  0,  0,  0,  0,  1,     // Species 7
			0,  0,  0,  0,  0,  0,  1 };   // Species 8
	// Reaction     A   B   C   D   E   F   G

  public:
	Repressilator_ODE() {}
	~Repressilator_ODE() {}

	void operator()( int dim, double* Y, double* dYdt, double t ) {
	
		// Reaction Rates
		double R[NREACTIONS] = {
			Ktl*Y[4],			// Reaction A
			Ktl*Y[5],			// Reaction B
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

	void observer(int dim, double Y[], double t ) {
		std::ofstream file;
		file.open ("bin/Result.csv",std::ios::app);
	
		file << t;
		for(int i=0; i<NSPECIES; i++) {
			file << ";" << Y[i];	
		}
		file << "\n";
	
		file.close();
	}

};

int main(int argc, char **argv) {
	
	Repressilator_ODE repr;
	double X0[] = { 1, 0, 0, 0, 0, 0, 0, 0 }; // initial condition

	// result matrix
	std::cout << "Runge Kutta Fehlberg Order 4...\t" << std::flush;
	rk4_wrapper<double, Repressilator_ODE>
		( NSPECIES, repr, X0 , 0.0 , 10000 , 1e-2);
	std::cout << "done" << std::endl;


	std::cout << "Runge Kutta Fehlberg Order 45...\t" << std::flush;
	rk45_wrapper<double, Repressilator_ODE>
		( NSPECIES, repr, X0 , 0.0 , 10000 , 1e-2, 1e-6);
		
	std::cout << "done" << std::endl;

}

