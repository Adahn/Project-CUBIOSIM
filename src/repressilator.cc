#include <iostream>
//#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>

#include <fstream>

#define NSPECIES 8
#define NREACTIONS 7

using namespace std;
using namespace boost::numeric::odeint;

typedef boost::array< double , NSPECIES > state_type;


// Parameter definition
const double Ktl = 1e-2;
const double Ktr = 1e-9;
const double KR = 1e-8;
const double nR = 2;
const double dprot = 1e-3;
const double dmRNA = 1e-2;

// decay vector
double decay[NSPECIES] = { dprot, // Species 1
						  dprot, // Species 2 
						  dprot, // Species 3
						  dprot, // Species 4
						  dmRNA, // Species 5
						  dmRNA, // Species 6
						  dmRNA, // Species 7
						  dmRNA }; // Species 8


// Stoichiometric matrix
int S[NSPECIES*NREACTIONS] = {	1,  0,  0,  0,  0,  0,  0,     // Species 1
							    0,  1,  0,  0,  0,  0,  0,     // Species 2
								0,  0,  1,  0,  0,  0,  0,     // Species 3
								0,  0,  0,  1,  0,  0,  0,     // Species 4
								0,  0,  0,  0,  1,  0,  0,     // Species 5
								0,  0,  0,  0,  0,  1,  0,     // Species 6
								0,  0,  0,  0,  0,  0,  1,     // Species 7
								0,  0,  0,  0,  0,  0,  1 };    // Species 8
					// Reaction A   B   C   D   E   F   G


// const double sigma = 10.0;
// const double R = 28.0;
// const double b = 8.0 / 3.0;

// typedef boost::array< double , 3 > state_type;

// void lorenz( const state_type &x , state_type &dxdt , double t )
// {
//     dxdt[0] = sigma * ( x[1] - x[0] );
//     dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
//     dxdt[2] = -b * x[2] + x[0] * x[1];
// }

// void write_lorenz( const state_type &x , const double t )
// {
//     cout << t << '\t' << x[0] << '\t' << x[1] << '\t' << x[2] << endl;
// }

// int main(int argc, char **argv)
// {
//     state_type x = { 10.0 , 1.0 , 1.0 }; // initial conditions
//     integrate( lorenz , x , 0.0 , 25.0 , 0.1 , write_lorenz );
//}

void Repressilator_ODE( const state_type &Y, state_type &dYdt, const double t ) {
	
	printf("t = %f\n", t);

	// Reaction Rates
	double R[NREACTIONS] = {(Y[4]>0)? Ktl*Y[4]:0 ,  					 // Reaction A
							(Y[5]>0)? Ktl*Y[5]:0 ,				     // Reaction B
							(Y[6]>0)? Ktl*Y[6]:0 ,				     // Reaction C
							(Y[7]>0)? Ktl*Y[7]:0 ,				     // Reaction D
							(Y[1]>0)? Ktr/(1 + pow(Y[1]/KR, nR)):0 ,       // Reaction E
							(Y[2]>0)? Ktr/(1 + pow(Y[2]/KR, nR)):0 ,       // Reaction F
							(Y[0]>0)? Ktr/(1 + pow(Y[0]/KR, nR)):0 };      // Reaction G

	// Model: dYdt = S * R - d * Y
	for (int i = 0; i < NSPECIES; ++i)
	{
		dYdt[i] = 0;
		for (int j = 0; j < NREACTIONS; ++j)
		{
			dYdt[i] += S[j+i*NREACTIONS] * R[j];

		}
	
		dYdt[i] -= decay[i]*Y[i];
	}

}


void write_ODE_result( const state_type &Y, const double t )
{
	ofstream file;
	file.open ("bin/Result.csv",ios::app); // write at the end of file
	
	file << Y[0];
	for(int i=1; i<NSPECIES; i++) {
		file << ";" << Y[i];	
	}
	file << "\n";
	
	file.close();
}

int main(int argc, char **argv) {

	/*
		List of species ( gives correspondance between index in the species vector and actual species names)
		1 - tetR
		2 - LacI
		3 - cI
		4 - GFP
		5 - mRNA_tetR
		6 - mRNA_LacI
		7 - mRNA_ cI
		8 - mRNA_GFP
	*/

	/*
		List of reactions
		A - tetR production
		B - LacI production
		C - cI production
		D - GFP production
		E - Repression of promoter Lac
		F - Repression of promoter lambda
		G - Repression of promoter tet
	*/


	// Simulation Control
	// double tstart = 0.0;
	// double tfinal = 100000;
	// double step = 0.1;

	state_type X0 = { 1e-6, 0, 0, 0, 0, 0, 0, 0 }; // initial condition

	// result matrix
	integrate( Repressilator_ODE , X0 , 0.0 , 100000.0 , 0.1 , write_ODE_result );

}

