#include <stdio.h>
#include <string>

template<class state_type>
void debug_GPU(state_type* v, int dim, std::string message){

	std::cout << message << '\n';
	for (size_t i = 0; i < dim; i++) {
		std::cout << v[i] << "\t";
	}
	std::cout << '\n';

}
