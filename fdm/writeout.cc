#include "writeout.hh"
//#include "mesh.cc"

Solution::Solution(Mesh mesh):mesh_(mesh){}

Solution::Solution(const Solution& sol):mesh_(sol.mesh_){}

void Solution::writeSol(int i)
{
	ofstream myfile;
	myfile.open("solution"+to_string(i)+".txt");
	myfile << mesh_;
	myfile.close();
}

/*
int main()
{
	Mesh mesh(10,10,1,1);
	Solution sol(mesh);
	sol.writeSol();
	return 0;
}
*/
