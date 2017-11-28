#include "Individual.hpp"
#include <vector>
using namespace std;

/**
 * Implementacion de las funciones de la clase Indiividuo
 * sirve para encapsular todos los datos que necesitamos 
 * conservar por cada elemento en la poblacion
 */

/**
 * Funcion que inicializa las variables de un individuo
 * @param variables - vector de variables 
 */
void Individual::init(vector<double>variables)
{
	this->variables = variables;
	this->rank = -1;
	this->fitness = 1e10;
}

/**
 * Se sobreescribio el constructor de la clase
 */
Individual::Individual()
{
	this->rank = -1;
	this->fitness = 1e10;
}

