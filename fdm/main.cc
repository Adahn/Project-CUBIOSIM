#include "mesh.hh"
#include "writeout.hh"
#include "fdm.hh"

#include <time.h>   
#include <sys/time.h>  


double my_gettimeofday(){
  struct timeval tmp_time;
  gettimeofday(&tmp_time, NULL);
  return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}

int main()
{
	Mesh mesh(100, 100);
	FDM fdm(mesh);
	double debutCPU = my_gettimeofday();
	fdm.simulation();
	double finCPU = my_gettimeofday();
	cout << "temps cpu : " <<  finCPU-debutCPU << endl;
	//fdm.results();
	return 0;
}
