#include "fdm.hh"
//#include "writeout.cc"

FDM::FDM(Mesh mesh, float D, float d, int Nt) :	mesh_(mesh), D_(D), d_(d), Nt_(Nt)
{}

double FDM::xx_discretization(int i, int j)
{
	double quotient= (double)(mesh_.getNx()/mesh_.getHx());
	quotient=quotient*quotient;
	return ( mesh_(i+1,j)-2*mesh_(i,j)+mesh_(i-1,j) )/quotient  ;
}

double FDM::yy_discretization(int i, int j)
{
	double quotient= (double)(mesh_.getNy()/mesh_.getHy());
	quotient=quotient*quotient;
	return ( mesh_(i,j+1)-2*mesh_(i,j)+mesh_(i,j-1) ) / quotient;
}

double FDM::laplacien(int i, int j)
{
	return xx_discretization(i,j) +yy_discretization(i,j);
}

void FDM::calculation()
{
	mesh_.setXold(mesh_.getX());
	double dt= (double) 1/Nt_;
	for(size_t i=1 ; i<mesh_.getX().size()-1 ; i++ )
	{
		for(size_t j=1 ; j<mesh_.getX()[i].size()-1 ; j++)
		{
			mesh_(i,j) += dt*D_*laplacien(i,j) - dt*d_*mesh_(i,j);
		}
	}
}

void FDM::simulation()
{
	for(int i=0 ; i<Nt_ ; ++i)
	{
		calculation();
		results(i);
	}
}

void FDM::results(int i)
{
	Solution sol(mesh_);
	sol.writeSol(i);
}


