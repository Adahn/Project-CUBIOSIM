#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;


class Repressilator_ODE
{
	private:
		// size
		int _n;
		int _nreac;
		int _nspec;

		// system parameters
		double* R;
		double* decay;
		int* S;

		// Rate parameters
		double _Ktl;
		double _Ktr;
		double _KR;
		double _nR;
		
		// other variables
		string _filename;
	
	public:	

		Repressilator_ODE(int n, double dprot, double dmRNA, double Ktl, double Ktr, double KR, double nR, string filename):
			_n(n), _nreac( 2*n+1 ), _nspec( 2*(n+1) ),
			_Ktl(Ktl), _Ktr(Ktr), _KR(KR), _nR(nR),
			_filename(filename)
		{
			S = new int[_nspec * _nreac]; //(int*)calloc( _nspec * _nreac, sizeof(int) );
			R = new double[_nspec]; //(double*)calloc(_nspec, sizeof(double));
			decay = new double[_nspec]; //(double*)calloc(_nspec, sizeof(double));
	
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

		// utility functions
		void observer( const int dim, double* Y, const double t )
		{
			ofstream file;
			file.open(_filename,ios::app); // write at the end of file
	
			file << t;
			for( int i=0; i<dim; i++) {
				file << ";" << Y[i];	
			}
			file << "\n";
	
			file.close();
		}

		/*void check() {
			for(int i=0; i<_nspec; i++) {
				cout << decay[i] << "\t";
			}
			cout << endl;
			cout << endl;
			for(int i=0; i<_nspec; i++) {
				for(int j=0; j<_nreac; j++) {
					cout << S[i *_nreac+ j] << "\t";
				}
				cout << endl;
			}
		}*/
		
	~Repressilator_ODE() {
		//free(S);	free(R);	free(decay);
		delete[] S;	delete[] R;	delete[] decay;
	}
};

