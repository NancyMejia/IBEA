#ifndef __INDIVIDUAL__
#define __INDIVIDUAL__

#include <vector>
#include <string>
using namespace std;

class Individual
{
public:
	Individual();
	// Metodo que permite inicializar el inviduo dado el vector de variables
	void init(vector<double>variables);

	// Parametros de la clase
	vector<double> variables;
	vector<double> objectives;
	int rank;
	double fitness;
};


#endif