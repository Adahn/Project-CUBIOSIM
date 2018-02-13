#include "mesh.hh"
#include "writeout.hh"
#include "fdm.hh"

int main()
{
	Mesh mesh(10, 10, 1, 1);
	FDM fdm(mesh);
	fdm.simulation();
	//fdm.results();
	return 0;
}
