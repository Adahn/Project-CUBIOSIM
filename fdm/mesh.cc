#include "mesh.hh"

Mesh :: Mesh(int nx, int ny, int hx, int hy, float phi0_)
	:nx_(nx), ny_(ny), hx_(hx), hy_(hy)
	{
		int ax = (int)nx_*hx_;
		int ay = (int)ny_*hy_;
		X_.resize(ax);
		Xold_.resize(ax);

		for (int i = 0; i < ax ; ++i){
			X_[i].resize(ay);
			Xold_[i].resize(ay);
		}

		X_[ax/2 -1][ay/2 -1]=phi0_; //source placed at the center of the grid
		X_[ax/2][ay/2 -1]=phi0_;
		X_[ax/2][ay/2]=phi0_;
		X_[ax/2 -1][ay/2]=phi0_;

		Xold_[ax/2 -1][ay/2 -1]=phi0_;
		Xold_[ax/2][ay/2 -1]=phi0_;
		Xold_[ax/2][ay/2]=phi0_;
		Xold_[ax/2 -1][ay/2]=phi0_;
	}

Mesh::Mesh(const Mesh& mesh) :
nx_(mesh.nx_), ny_(mesh.ny_), hx_(mesh.hx_), hy_(mesh.hy_), X_(mesh.X_), Xold_(mesh.Xold_), phi0_(mesh.phi0_)
{}


int Mesh::getNx() const
{
	return nx_;
}

int Mesh::getNy() const
{
	return ny_;
}
int Mesh::getHx() const
{
	return hx_;
}

int Mesh::getHy() const
{
	return hy_;
}

vector< vector<double> > Mesh::getX() const
{
	return X_;
}
vector< vector<double> > Mesh::getXold() const
{
	return Xold_;
}

void Mesh::setXold(vector< vector<double> > const & v)
{
	Xold_=v;
}

ostream& Mesh:: display(ostream& os) const
{
	for (auto & vec : X_)
	{
		for (auto & val : vec)
		{
			os << val << " ";
		}
		os << '\n';
	}
	os << '\n';
	/*for (unsigned int i=0; i< X_.size() ; ++i)
	{
		for (unsigned int j=0; j< X_[i].size() ; ++j)
		{
			os << X_[i][j] - Xold_[i][j] << " ";
		}
		os << '\n';
	}*/
	
	return os;
}

ostream &operator << (ostream & os, Mesh const &mesh)
{
   return mesh.display(os);
}

double Mesh::operator()(int r,int c) const
{
	return X_[r][c];
}

double& Mesh::operator()(int r,int c)
{
	return X_[r][c];
};

Mesh& Mesh::operator=(Mesh const& mesh)
{
	nx_=mesh.nx_;
	ny_=mesh.ny_;
	hx_=mesh.hx_;
	hy_=mesh.hy_;
	X_=mesh.X_;
	//Xold_=mesh.Xold_; //to calculate the error 
	return *this;
}
/*
int main()
{
	Mesh mesh(10, 10, 1, 1);
	vector< vector<double> > v= mesh.getX();
	mesh(0,0)=5;
	cout << mesh<<endl;
	return 0;
}*/
