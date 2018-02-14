#ifndef WRITEOUT_H
#define WRITEOUT_H

#include "mesh.hh"
#include <fstream>
#include <string>

using namespace std;

class Solution
{
	public:
		Solution(Mesh mesh);
		Solution(const Solution& sol);
		void writeSol(int i); //write current solution

	private:
		Mesh mesh_;
};

#endif

