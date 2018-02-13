#ifndef FDM_H
#define FDM_H

#include "writeout.hh"

//we use dirichlet condition, boundary nodes = 0
class FDM
{
	public:
		FDM(Mesh mesh, float D=10.0, float d=0.4, int Nt=10000);
		double xx_discretization(int i, int j);
		double yy_discretization(int i, int j);
		double laplacien(int i, int j);
		void calculation(); //calculate node value temporal discretization
		void simulation();
		void results(int i);

	private:
		Mesh mesh_;
		float D_; //diffusion coefficient
		float d_; //decay coefficient
		int Nt_; //number of time steps
		//Solution out_; better create it in a function
};


#endif
