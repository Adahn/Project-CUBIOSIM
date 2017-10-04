#include <iostream>

#include <boost/numeric/odeint.hpp>

#include <fstream>
#include <repressilator.hh>

using namespace std;
using namespace boost::numeric::odeint;

//typedef boost::array< double , NSPECIES > state_type;

// Parameter definition
const double Ktl = 1e-2;
const double Ktr = 1e-3;
const double KR = 1e-2;
const double nR = 2;
const double dprot = 1e-3;
const double dmRNA = 1e-2;


void write_ODE_result( const double* Y, const double t )
{
	ofstream file;
	file.open ("bin/Result.csv",ios::app); // write at the end of file
	
	file << Y[0];
	for(int i=1; i<8; i++) {
		file << ";" << Y[i];	
	}
	file << "\n";
	
	file.close();
}

int main(int argc, char **argv) {

	double X0[8] = { 1e-6, 0, 0, 0, 0, 0, 0, 0 }; // initial condition
	Repressilator_ODE repr(3, dprot, dmRNA, Ktl, Ktr, KR, nR);

	// result matrix
	integrate( repr , X0 , 0.0 , 100000.0 , 1e-2 , write_ODE_result );

	return 0;

}

