#include "Individual.hpp"
#include "WFG2.hpp"
#include "WFG3.hpp"
#include "WFG8.hpp"
#include <string>
#include <vector>
using namespace std;

class IBEA
{
public:
	// COnstructor de la clase
	IBEA(int obj_num, int pop_size, double scale_factor, double epsilon, string function, double mut_prob, double cross_prob);
	// Funcion prinicipal para resolver el problema
	vector<Individual*> solve(int generations);
	//variables de la clase
	int pop_size;
	int obj_num;
	double scale_factor;
	double epsilon;
	double mutation_prob;
	double cross_prob;
	string function;	
	WFG2 *wfg2;
	WFG3 *wfg3;
	WFG8 *wfg8;

private:
	// FUnciones que son usadas dentro de la funcion principal que se encarga de resolver el problema
	vector<Individual*> generatePopulation(string function, int pop_size);
	void evaluatePopulation(vector<Individual*> *population);
	double epsilonIndicator(Individual *ind1, Individual *ind2);
	void appendOffspring(vector<Individual*>*population, vector<Individual*> new_pop);
	vector<Individual*> selection(vector<Individual*> population);
	void crossoverInd(Individual *parent1, Individual *parent2, Individual *ind1, Individual*ind2);
	Individual* tournament(Individual*ind1, Individual*ind2);
	void mutateInd(Individual* ind);	
	void mutation(vector<Individual*> *population);
	Individual* getBest(vector<Individual*> population, int *position);
	void updateFitness(Individual* best, vector<Individual*> *population);
	vector<Individual*> getNonDominated(vector<Individual*> population);
	bool dominated(Individual* ind1, Individual* ind2);
	double next_double( const double bound);
};