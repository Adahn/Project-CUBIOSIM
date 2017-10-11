#pragma once
#include <iostream>

using namespace std;


class Repressilator_ODE
{
	private:
		// size
		int _n;
		int _nreac;
		int _nspec;

		// system parameters
		vector<double> decay;
		vector<int> S;

		// Rate parameters
		double _Ktl;
		double _Ktr;
		double _KR;
		double _nR;
		
		vector<double> R;
	
	public:	

		Repressilator_ODE(int n, double dprot, double dmRNA, double Ktl, double Ktr, double KR, double nR):
			_n(n),
			_nreac( 2*n+1 ),
			_nspec( 2*(n+1) ),
			_Ktl(Ktl),
			_Ktr(Ktr),
			_KR(KR),
			_nR(nR)
		{
			S = vector<int>(_nspec * _nreac);
			R = vector<double>(_nspec);
			decay = vector<double>(_nspec);
	
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


		void operator()( const vector<double>& Y, vector<double>& dYdt, const double t ) {

			// Reaction Rates
			for( int i=0; i<_n+1; i++ ) {
				R[i] = _Ktl*Y[i+_n+1];
			}
			for( int i=0; i<_n+1; i++ ) {
				R[i+_n+1] = _Ktr/(1 + pow(Y[i]/_KR, _nR));
			}

			// Model: dYdt = S*R - d*Y
			for( int i=0; i<_nspec; i++) {
				dYdt[i] = 0;
				for( int j=0; j<_nreac; j++) {
					dYdt[i] += S[i *_nreac + j] * R[j];
				}
				dYdt[i] -= decay[i]*Y[i];
			}

		}


/*		void check() {
			for(int i=0; i<i_nspec; i++) {
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
};
