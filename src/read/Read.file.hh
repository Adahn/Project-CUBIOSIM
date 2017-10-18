#ifndef READFILE_H_
#define READFILE_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;


class Espece{
 	public: 
 		Espece(string nom_fichier);
 		~Espece();
		string get_name();
		int get_nb_react();
		int get_nb_esp();
		int * get_speed();
		int * get_Y();
		int ** get_stoch();
		int * degrad();
 	private:
 	  string file_name_;
 		int nb_esp_;
		int nb_react_;
 		string name_;
 		int * Y_; //concentration
 		int * speed_; //vitesse
   	int ** stoch_; //matrice stochastique
    int * degrad_; //vecteur de degradation
 };


#endif


