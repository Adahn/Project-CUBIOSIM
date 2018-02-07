#include <stdio.h>
#include <string>

// prints out message followed by vector v
template<class state_type>
void print_array(state_type* v, int dim, std::string message){

	std::cout << message << endl;
	for (int i = 0; i < dim; i++) {
		std::cout << v[i] << "\t";
	}
	std::cout << endl;

}
