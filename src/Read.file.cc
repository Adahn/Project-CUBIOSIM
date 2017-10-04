#include <iostream>
#include "Read.file.hh"

using namespace std;

int main(){
    Espece a("test.txt");
	int nbesp=a.get_nb_esp();
	int* degrad=new int[nbesp];
	degrad=a.degrad();
	for(int i=0;i<nbesp;i++){
		cout<<degrad[i];
	}
	cout<<endl;
  return 0;
 }
  
