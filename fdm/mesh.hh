#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <vector>

using namespace std;

class Mesh
{
	public :
		Mesh(int nx, int ny, int hx, int hy, float phi0=0.4);
		Mesh(const Mesh& mesh);
		int getNx() const;
		int getNy() const;
		int getHx() const ;
		int getHy() const;
		vector< vector<double> > getX() const;
		vector< vector<double> > getXold() const;
		ostream& display(ostream& os) const; //print X_. This function is called into another one, to display through a mesh function
		friend ostream &operator << (ostream & os, Mesh const &mesh);
		double operator()(int r,int c)const; //get value of X_[r][c] through a mesh instance
      double& operator()(int r,int c);  //set value X_[r][c] through a mesh instance
		Mesh& operator=(Mesh const& mesh);
		void setXold(vector< vector<double> > const & v);
		//plus gmsh function

	private :
		int nx_; //grid dimension x axis
		int ny_; //grid dimension y axis
		int hx_; //mesh step x axis
		int hy_; //mesh step y axis
		vector< vector<double> > X_; //vector that stores the current X results
		vector< vector<double> > Xold_; //vector that stores the old X results, to calculate the error
		float phi0_; //value at the center of the grid X
};

#endif
