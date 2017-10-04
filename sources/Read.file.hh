#ifndef READFILE_H_
#define READFILE_H_

#include <iostream>
#include "Read.file.hh"
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

Espece::Espece(string nom_fichier):
file_name_(nom_fichier)
{
	//ifstream fichier(file_name_, ios::in);  // on ouvre le fichier en lecture
	ifstream fichier(file_name_.c_str());
	//fichier.open(file_name_);
 
        if(fichier)  // si l'ouverture a réussi
	{       
		getline(fichier, name_); //nom de l'espèce
		cout << "name: " << name_ << endl;  // on affiche la ligne
		fichier >> nb_esp_;
		cout << "nb_esp: " << nb_esp_ << endl;  // on affiche la ligne
		fichier >> nb_react_;
		cout << "nb_react: " << nb_react_ << endl;  // on affiche la ligne
		Y_=new int[nb_esp_];
		for (int i=0;i<nb_esp_;i++){
			fichier >> Y_[i];
		}
		//cout << "dernier Y: " << Y_[nb_esp_-1] << endl;  // on affiche la ligne
		speed_=new int[nb_react_];
		for (int i=0;i<nb_react_;i++){
			fichier >> speed_[i];
		}
		//cout << "dernier speed: " << speed_[nb_esp_-1] << endl;  // on affiche la ligne
		stoch_=new int*[nb_esp_];
		for(int i=0;i<nb_esp_;i++){
		       stoch_[i]=new int[nb_react_];
		}
		
		for (int i=0;i<nb_esp_;i++){
			for(int j=0;j<nb_react_;j++){
			        fichier >> stoch_[i][j];
				//cout << "stoch: " << stoch_[i][j];  // on affiche la ligne
			}
		}
		degrad_=new int[nb_esp_];
		for (int i=0;i<nb_esp_;i++){
			fichier >> degrad_[i];
		}
		//cout << "dernier degrad: " << degrad_[nb_esp_-1] <<endl;  // on affiche la ligne
	        fichier.close();  // on ferme le fichier
	}
        else  // sinon
                cerr << "Impossible d'ouvrir le fichier !" << endl;

}

Espece::~Espece(){
{
	for (int i=0;i<nb_esp_;i++){
		delete stoch_[i];
	}
	delete stoch_;
	}
	delete Y_;
	delete speed_;
	delete degrad_;
} //désallouer

string Espece::get_name(){
	return name_;
}
int Espece::get_nb_react(){
	return nb_react_;
}
int Espece::get_nb_esp(){
	return nb_esp_;
}
int * Espece::get_speed(){
	return speed_;
}
int * Espece::get_Y(){
	return Y_;
}
int ** Espece::get_stoch(){
	return stoch_;
}
int * Espece::degrad(){
	return degrad_;
}

#endif


