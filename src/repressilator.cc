#include <iostream>
#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>

#define NSPECIES 8
#define NREACTIONS 7

using namespace std;
using namespace boost::numeric::odeint;

// Parameter definition
float Ktl = 1e-2;
float Ktr = 1e-9;
float KR = 1e-8;
float nR = 2;
float dprot = 1e-3;
float dmRNA = 1e-2;

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


template<typename T>void compute_reaction_rates(const T* X, T* dxdt, int n ) {
	// Reaction Rates
	T R[NREACTIONS] = {	Ktl * Y(5) ,  					 // Reaction A
						Ktl * Y(6) ,				     // Reaction B
						Ktl * Y(7) ,				     // Reaction C
						Ktl * Y(8) ,				     // Reaction D
						Ktr / (1 + (Y(2)/KR)^nR) ,       // Reaction E
						Ktr / (1 + (Y(3)/KR)^nR) ,       // Reaction F
						Ktr / (1 + (Y(1)/KR)^nR) };      // Reaction G



}


int main(int argc, char **argv) {

	// Simulation Control
	int tfinal = 100000;
	double X0[NSPECIES] = { 1e-6, 0, 0, 0, 0, 0, 0, 0 }; // initial condition

	// Stoichiometric matrix
	int S[NSPECIES*NREACTIONS] = {	1,  0,  0,  0,  0,  0,  0,     // Species 1
									0,  1,  0,  0,  0,  0,  0,     // Species 2
									0,  0,  1,  0,  0,  0,  0,     // Species 3
									0,  0,  0,  1,  0,  0,  0,     // Species 4
									0,  0,  0,  0,  1,  0,  0,     // Species 5
									0,  0,  0,  0,  0,  1,  0,     // Species 6
									0,  0,  0,  0,  0,  0,  1,     // Species 7
									0,  0,  0,  0,  0,  0,  1 };    // Species 8
}