#pragma once
//#include <boost/numeric/odeint.hpp>
#include <iostream>

struct Repressilator_ODE
{
	private:
		int nreactions;
		int nspecies;
		int __n;
		double __Ktl;
		double __Ktr;
		double __KR;
		double __nR;
		double* decay;
		int* S;
		double* R;
	
	public:	

		Repressilator_ODE(int n, double dprot, double dmRNA, double Ktl, double Ktr, double KR, double nR):
			__n(n),
			nspecies( 2*(n+1) ), // n repressors, n+1 reactors, 1 indicator
			nreactions( 2*n+1 ),
			__Ktl(Ktl),
			__Ktr(Ktr),
			__KR(KR),
			__nR(nR)
		{
			S = new int[nspecies*nreactions]();
			R = new double[nspecies]();
			decay = new double[nspecies]();
	
			int i,j;
			
			for( i=0; i<n+1; i++ ) {
				decay[i] = dprot;			// prot
				S[i*nreactions + i] = 1;	// S looks like Identity
			}
	
			for( j=i; j<2*n+1; j++ ) {		
				decay[j] = dmRNA;			// 
				S[j*nreactions + j] = 1;	// S looks like Identity
			}
			
			decay[j] = dmRNA;				// last RNA
			S[j*nreactions + j-1] = 1;	// S last line is the same as the previous one
		}

		~Repressilator_ODE() {	}


		void operator()( const double* Y, double* dYdt, const double t ) {

			// Reaction Rates
			for( int i=0; i<__n+1; i++ ) {
				R[i] = __Ktl*Y[i+__n+1];
			}
			for( int i=0; i<__n+1; i++ ) {
				R[i+__n+1] = __Ktr/(1 + pow(Y[i]/__KR, __nR));
			}

			// Model: dYdt = S*R - d*Y
			for (int i = 0; i < nspecies; ++i) {
			
				dYdt[i] = 0;
				
				for (int j = 0; j < nreactions; ++j) {
					dYdt[i] += S[i*nreactions + j] * R[j];
				}
	
				dYdt[i] -= decay[i]*Y[i];
				
			}

		}


/*		void check() {
			for(int i=0; i<nspecies; i++) {
				std::cout << decay[i] << "\t";
			}
			std::cout << std::endl;
			std::cout << std::endl;
			for(int i=0; i<nspecies; i++) {
				for(int j=0; j<nreactions; j++) {
					std::cout << S[i*nreactions+j] << "\t";
				}
				std::cout << std::endl;
			}
		}*/
};
