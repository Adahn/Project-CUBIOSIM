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

#define NSPECIES 8
#define NREACTIONS 7

using namespace std;
using namespace boost::numeric::odeint;

typedef double value_type;
typedef thrust::device_vector< double > state_type;
typedef thrust::device_vector< device_vector< double > > matrix;
typedef runge_kutta4< state_type , value_type , state_type , value_type , thrust_algebra , thrust_operations > stepper_type;


// Parameter definition
const double Ktl = 1e-2;
const double Ktr = 1e-3;
const double KR = 1e-2;
const double nR = 2;
const double dprot = 1e-3;
const double dmRNA = 1e-2;


state_type decay(NSPECIES); // device
thrust::fill(decay.begin(), decay.begin() + 4, dprot);//copie des 4 premiers elts à dprot
thrust::fill(decay.begin()+4, decay.begin() + NSPECIES, dmRNA);//copie des


struct repressilator
{

	struct R_vect {
		R_vect(vector <double> Y) : Y_(Y) {}
		vector <double> operator()(vector<double> Y) const { 		
			for(size_t i=0;i<8;i++){//Declaration de R
				if(i<4)
					R_.push_back(Ktl*Y[i+4]);
				else if(i<6)
					R_.push_back(Ktr/(1 + pow(Y[i-3]/KR, nR)));
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

int main( int arc , char* argv[] )
{
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

	// result matrix
	//integrate( Repressilator_ODE , X0 , 0.0 , 100000.0 , 1e-2 , write_ODE_result );

