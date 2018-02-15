CODE SEQUENTIEL : 

classe Maillage dans (mesh.*) : 
construit le maillage

classe Solution (dans writeout.*) :
donne la solution pour une itération i

classe fdm (dans fdm.*) :
implémente la méthode des différences finies pour le maillage
attention : le code n'est pas optimal. Il fonctionne pour une grille avec un pas de 1.
Pour un pas suivant x et y, il faut réecrire la formule à chaque point en fonction de hx et hy



CODE PARALLELE

dans nvidia.cu

