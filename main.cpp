#include "IBEA.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cstring>
#include <ctime>

using namespace std;

/**
 * Mejia Juarez Nancy Arlette
 * Tarea 03 larga
 * 
 * Archivo de pruebas para resolver cada una de las funciones de prueba
 */

/**
 * Metodo princiipal que nos permite realizar las pruebas del algoritmo
 */
int main(int argc, char **argv)
{
	if(argc < 4)
	{
		cout << "ERROR: Se debe de pasar como parametro el nombre de la funcion, numero de generaciones y factor de escala" << endl;
		return -1;
	}
	srand(time(NULL));
	//cout << "IBEA " << endl;
	int trials = 5;
	int pop_size = 100;
	int obj_num = 2;
	double epsilon = 0.1;
	string function = argv[1];
	if(function != "WFG2" && function != "WFG3" && function != "WFG8")
	{
		cout  << "ERROR: EL nombre de la funcion no es valida. (WFG2, WFG3, WFG8)" << endl;
		return -1; 
	}
	cout << "function " << function << endl;
	int generations = atoi(argv[2]);
	double scale_factor = atof(argv[3]);

	double mut_prob = 0.01;
	double cross_prob = 1.0;
	IBEA *ibea = new IBEA(obj_num,pop_size, scale_factor, epsilon, function, mut_prob, cross_prob);

	for(int trial = 0; trial < trials; trial++)
	{
		char filename_solution[180];
		vector<Individual*> solution = ibea->solve(generations);
		//Archivo en el que se escribiran las soluciones
		ofstream outputfile;
		if (function == "WFG2")
		{
			strcpy(filename_solution, "results/WFG2_solution");
		}
		else if(function == "WFG8")
		{
			strcpy(filename_solution, "results/WFG8_solution");
		}
		else if(function == "WFG3")
		{
			strcpy(filename_solution, "results/WFG3_solution");
		}

		strcat(filename_solution, to_string(trial).c_str());
		strcat(filename_solution, "_scale_");
		strcat(filename_solution, to_string(scale_factor).c_str());
		strcat(filename_solution, ".txt");
		cout << filename_solution << endl;
		
		outputfile.open(filename_solution);
		
		for(int i = 0; i < solution.size();i++)
		{
			for(int j = 0; j < obj_num; j++)
			{
				outputfile << solution.at(i)->objectives.at(j)  << " ";
			}
			outputfile << endl;
		}		
		// Cerramos el archivo
		outputfile.close();
		
	}	
	return 0;
}