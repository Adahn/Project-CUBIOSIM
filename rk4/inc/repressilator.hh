//
// Comments of this forme:
// !! this is a sample comment !!
//show where you can modify parameters, files, etc.

#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

/* Repressilator for this model: dYdt = S*R - d*Y */
class Repressilator_ODE
{
	private:
		// sizes
		int _n;
		int _nreac;
		int _nspec;

		// system parameters
		// !! Add/Modify vectors depending on the system you study !!
		double* R;	// reaction rates
		double* decay;	// decay vector
		int* S;		// stochiometric matrix

		// Rate parameters
		// !! Modify vectors depending on the system you study!!
		double _Ktl;
		double _Ktr;
		double _KR;
		double _nR;

		// utility variables
		string _filename;

	public:

		Repressilator_ODE(int n, double dprot, double dmRNA, double Ktl, double Ktr, double KR, double nR, string filename):
			_n(n), _nreac( 2*n+1 ), _nspec( 2*(n+1) ),
			_Ktl(Ktl), _Ktr(Ktr), _KR(KR), _nR(nR),
			_filename(filename)
		{
			S = new int[_nspec * _nreac];
			R = new double[_nspec];
			decay = new double[_nspec];

			// create stochiometric matrix S and decay vector
			int i,j;
			for( i=0; i<_n+1; i++ ) {
				decay[i] = dprot;	// prot
				S[i *_nreac+ i] = 1;	// S looks like Identity
			}
			for( j=i; j<2*_n+1; j++ ) {
				decay[j] = dmRNA;	// RNA
				S[j *_nreac+ j] = 1;	// S looks like Identity
			}
			decay[j] = dmRNA;		// last RNA
			S[j *_nreac+ j-1] = 1;		// S last line is the same as the previous one
		}

		// Overloading operator () to compute the derivate
		void operator()( const int dim, double* Y, double* dYdt, const double t ) {

			// Reaction Rates
			for( int i=0; i<_n+1; i++ ) {
				R[i] = _Ktl*Y[i+_n+1];
			}
			for( int i=1; i<_n; i++ ) {
				R[i+_n] = _Ktr/(1 + pow(Y[i]/_KR, _nR));
			}
			R[2*_n] = _Ktr/(1 + pow(Y[0]/_KR, _nR));

			// Model: dYdt = S*R - d*Y
			for( int i=0; i<_nspec; i++) {
				dYdt[i] = 0;
				for( int j=0; j<_nreac; j++) {
					dYdt[i] += S[i *_nreac + j] * R[j];
				}
				dYdt[i] -= decay[i]*Y[i];
			}

		}

		// This function writes the state given in argument to a file
		void observer( const int dim, double* Y, const double t ) {
			ofstream file;
			file.open(_filename,ios::app); // write at the end of file

			file << t;
			for( int i=0; i<dim; i++) {
				file << ";" << Y[i];
			}
			file << "\n";

			file.close();
		}

	~Repressilator_ODE() {
		delete[] S;	delete[] R;	delete[] decay;
	}
};
