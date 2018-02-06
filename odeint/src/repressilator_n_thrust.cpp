/***************

Project cuBIOSIM
Polytech Paris UPMC

authors:
	Liliane Kissita
	Elise Grojean
	Adrian Ahne
	Lucas Gaudelet

****************/

#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <algorithm>
#include <cstdlib>
#include <boost/numeric/odeint/external/thrust.hpp>
#include <fstream>
#include <vector> 

// command line parser and timer
#include <chCommandLine.h>
#include <chTimer.hpp>


// Constants
const static int DEFAULT_N = 3;

const static double Ktl = 1e-2;	// Translation rate
const static double Ktr = 1e-3;	// Transcription rate
const static double KR = 1e-2;  // Strength of repressors
const static double nR = 2;	// Hills coefficient of repressors
const static double dprot = 1e-3;	// Protein degradation rate
const static double dmRNA = 1e-2;	// mRNA degradation rate

int NSPECIES;
int NREACTIONS;

// function prototypes
void write_ODE_result( const vector<double> &Y, const double t );
void print_help(char* argv);




struct repressilator
{

	struct R_vect {
		R_vect(vector <double> Y) : Y_(Y) {}
		vector <double> operator()(vector<double> Y) const { 		
			for(size_t i=0;i<NSPECIES;i++){//Declaration de R
				if(i<NSPECIES/2)
					R_.push_back(Ktl*Y[i+NSPECIES/2]);
				else if(i<NSPECIES-2)
					R_.push_back(Ktr/(1 + pow(Y[i-(NSPECIES/2-1)]/KR, nR)));
				else
					R_.push_back(Ktr/(1 + pow(Y[0]/KR, nR)));
			}
		}
		private:
  	vector <float> Y_;
	};

    struct repressilator_functor
    {
        template< class T >
        __host__ __device__
        void operator()( T t ) const
        {
            // unpack the parameter we want to vary and the variables
            value_type Y = thrust::get< 0 >( t );
            value_type R = thrust::get< 1 >( t );
            thrust::get< 2 >( t ) = R-d*Y;
        }
    };

    repressilator()
    : x( const state_type &y ),R(R_vect(y)) { }

    void operator()(  const state_type &x , state_type &dxdt , value_type t ) const
    {
				R_vect R(x);
        thrust::for_each(
            thrust::make_zip_iterator( thrust::make_tuple( x.begin() , R.begin() , dxdt.begin() ) ) ,
            thrust::make_zip_iterator( thrust::make_tuple( x.end() , R.end() , dxdt.end() ) ) ,
            repressilator_functor() );
    }

		//state_type d_;
   // const state_type &x;
		state_type R;
};





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

const value_type dt = 0.01;

// main
int main(int argc, char **argv) {

	// print help
	bool      help = chCommandLineGetBool("h", argc, argv);
	if(!help) help = chCommandLineGetBool("help", argc, argv);
	if(help)  print_help(argv[0]);

	// repressilator size
	int n = -1;
	chCommandLineGet<int>(&n, "n", argc, argv);
	chCommandLineGet<int>(&n, "repressor-number", argc, argv);
	n = (n!=-1)? n:DEFAULT_N;
	
	if( n%2==0 ) {
		cout << "n must be an odd interger" << endl;
		exit(-1);
	}
	else {
		cout << "repressors : " << n 
			<< "\tspecies: " << 2*(n+1)
			<< "\treactions: " << 2*n+1 << endl;
			NSPECIES=2*(n+1);
			NREACTIONS=2*n+1;
	}

	state_type decay(NSPECIES); // device
	thrust::fill(decay.begin(), decay.begin() + NSPECIES/2, dprot);//copie des 4 premiers elts à dprot
	thrust::fill(decay.begin()+NSPECIES/2, decay.begin() + NSPECIES, dmRNA);//copie des

    // initialize the initial state X
    vector< value_type > X0_host( NSPECIES );
   	X0_host.push_back(1e-6);
		for(size_t i=1;i<NSPECIES;i++)
			X0_host.push_back(0);
    state_type X0 = X0_host;


    // integrate
		repressilator rep();
    integrate_const( stepper_type() , boost::ref( rep ) , X0 , 0.0 , 10.0 , dt );

    return 0;
}



// utility functions
void write_ODE_result( const vector<double> &Y, const double t )
{
	ofstream file;
	file.open ("bin/Result.csv",ios::app); // write at the end of file
	
	file << t;
	for( double d : Y) {
		file << ";" << d;	
	}
	file << "\n";
	
	file.close();
}

