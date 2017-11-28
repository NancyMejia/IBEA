#include "Individual.hpp"
#include "IBEA.hpp"
#include <string>
#include <vector>
#include <cstdlib>	
#include <cassert>
using namespace std;

/**
 * Tarea 03 Larga
 * Mejia Juarez Nancy Arlette
 * Implementacion del algoritmo IBEA
 * Indicator Based Evolutionary Algorithm, resolver las funciones
 * WFG2, WFG3 y WFG8.
 */


//** Using a uniform random distribution, generate a number in [0,bound]. ***
double IBEA::next_double( const double bound = 1.0 )
{
  assert( bound > 0.0 );

  return bound * rand() / static_cast< double >( RAND_MAX );
}

/**
 * Constructor de la clase
 * @param obj_num  - numero de objetivos del problema 
 * @param pop_size - tamanio de la poblacion
 * @param scale    - factor de escala
 * @param eps      -  valor de epsilon para el indicador 
 * @param funct    - cadena con el nombre de la funcion que se quiere resolver
 */
IBEA::IBEA(int obj_num, int pop_size, double scale, double eps, string funct, double mut_prob, double cross_prob)
{
	this->obj_num = obj_num;
	this->pop_size = pop_size;
	this->scale_factor = scale;
	this->epsilon = eps;
	this->mutation_prob = mut_prob;
	this->cross_prob = cross_prob;
	this->function = funct;
	if(funct == "WFG2")
	{
		this->wfg2 = new WFG2();
		this->wfg2->init(obj_num);
		this->wfg3 = NULL;
		this->wfg8 = NULL;
	}
	else if(funct == "WFG3")
	{
		this->wfg3 = new WFG3();
		this->wfg3->init(obj_num);
		this->wfg2 = NULL;
		this->wfg8 = NULL;
	}
	else if(funct == "WFG8")
	{
		this->wfg8 = new WFG8();
		this->wfg8->init(obj_num);
		this->wfg2 = NULL;
		this->wfg3 = NULL;
	}
}


/**
 * Funcion auxiliar que genera la poblacion inicial. 
 * Esta se genera de acuerdo a la funcion que se quiere resolver.
 * Recibe como parametros
 * Nombre de la funcion a resolver 
 * Tamanio de la poblacion que se quiere obtener
 */
vector<Individual*> IBEA::generatePopulation(string function, int pop_size)
{
	vector<Individual*> population;
	for (int i = 0; i < pop_size; i++)
	{
		Individual *new_individual = new Individual();
		vector<double> variables;
		if(function == "WFG2")
		{
			variables = this->wfg2->WFG_2_thru_7_random_soln(this->wfg2->k, this->wfg2->l);
		}
		else if(function == "WFG3")
		{
			variables = this->wfg3->WFG_2_thru_7_random_soln(this->wfg3->k, this->wfg3->l);
		}
		else if(function == "WFG8")
		{
			variables = this->wfg8->WFG_8_random_soln(this->wfg8->k, this->wfg8->l);
		}
		new_individual->init(variables);
		population.push_back(new_individual);
	}
	return population;
}

/**
 * FUncion que realiza el calculo del indicador epsilon para una pareja
 * de individuos, en este caso como se desea minimizar un problema se debe 
 * tomar el maximo de estas diferencias
 * @param  ind1 - Invididuo uno que se quiere comparar
 * @param  ind2 - Individuo dos que se quiere comparar
 * @return      - Valor del indicador epsilon
 */
double IBEA::epsilonIndicator(Individual *ind1, Individual *ind2)
{
	vector<double> values; 
	for(int i=0; i < this->obj_num; i++)
	{
		// queremos ver cuanto falta para pasar el ind1 al ind2 
		double epsilon = ind1->objectives.at(i) - ind2->objectives.at(i);
		values.push_back(epsilon);
	}	
	double min_value = values.at(0);
	for(int i = 0; i < this-> obj_num; i++)
	{
		if (values.at(i) >= min_value)
		{
			min_value = values.at(i);
		}
	}
	return min_value;
}

/**
 * FUncion que calcula el valor del fitness de la poblacion
 * en este caso como solo se cuenta con el conjunto de variables de cada
 * individuo se debe calcular primero el conjunto de objetivos y apartir
 * de estos valores se puede obtener el fitness usando el indicador epsilon
 * @param population - poblacion que se desea evaluar
 */
void IBEA::evaluatePopulation(vector<Individual*> *population)
{
	vector<double> vector_aux;
	for(int i = 0; i < this->obj_num; i++)
	{
		vector_aux.push_back(0.0);
	}
	// Getting objectives values 
	for(int i = 0; i < population->size(); i++)
	{
		if (this->function == "WFG2")
		{
			this->wfg2->evaluate(population->at(i)->variables, vector_aux);		
		}
		else if(this->function == "WFG3")
		{
			this->wfg3->evaluate(population->at(i)->variables, vector_aux);		
		}
		else if(this->function == "WFG8")
		{
			this->wfg8->evaluate(population->at(i)->variables, vector_aux);		
		}
		population->at(i)->objectives = vector_aux;		
	}
	// Evaluating with epsilon indicator  
	// Se obtiene el valor de fitness para cada individuo usando el
	// indicador epsilon
	for(int i = 0; i < population->size(); i++)
	{
		double fitness_val = 0.0;
		for(int j = 0; j < population->size(); j++)
		{
			if(i != j)
			{
		
				fitness_val += -exp(-epsilonIndicator(population->at(j), population->at(i))/ this->scale_factor);	
			}
		}
		population->at(i)->fitness = fitness_val;
	}
}

/**
 * Funcion que realiza la cruza entre dos individuos y los almacena en las variables
 * child1, child2
 * @param parent1 - Primer padre de los individuos
 * @param parent2 - Segundo padre de los individuos
 * @param ind1    - Primer hijo generado
 * @param ind2    - Segundo hijo generado
 */
void IBEA::crossoverInd(Individual *parent1, Individual *parent2, Individual *ind1, Individual*ind2)
{
	double pcross = this->cross_prob;
	double EPS = 1.0e-14;
	double y1, y2, ylow, yupper, rnd, alpha, beta, betaq;
	double c1, c2;
	double eta_c = 10;
	int nvar;
	if(this->function == "WFG2")
	{
		nvar = this->wfg2->k + this->wfg2->l;
	}
	else if(this->function == "WFG8")
	{
		nvar = this->wfg8->k + this->wfg8->l;
	}
	if(next_double() <= pcross)
	{
		for(int i = 0; i < nvar; i++)
		{
			if(next_double() <= 0.5)
			{
				if(fabs(parent1->variables.at(i) - parent2->variables.at(i)) > EPS)
				{	
					if(parent1->variables.at(i) < parent2->variables.at(i))
					{
						y1 = parent1->variables.at(i);
						y2 = parent2->variables.at(i);
					}
					else
					{
						y1 = parent2->variables.at(i);
						y2 = parent1->variables.at(i);
					}

					ylow = 0.0;
					yupper = 2.0 * (i+1);
					rnd = next_double();
					beta = 1.0 + (2.0*(y1- ylow)/(y2-y1));
					alpha = 2.0 - pow(beta, -(eta_c+1.0));
					if(rnd <= (1.0/alpha))
					{
						betaq = pow((rnd*alpha), (1.0/(eta_c+1.0)));
					}
					else
					{
						betaq = pow((1.0/(2.0 - rnd*alpha)), (1.0/(eta_c+1.0))); 
					}
					c1 = 0.5 * ((y1+y2) - betaq*(y2-y1));
					beta = 1.0 + (2.0*(yupper -y2)/(y2-y1));
					alpha = 2.0 - pow(beta, -(eta_c+1.0));
					if(rnd <= (1.0/alpha))
					{
						betaq = pow((rnd*alpha), (1.0/(eta_c+1.0)));
					}
					else
					{
						betaq = pow((1.0/(2.0 - rnd*alpha)), (1.0/(eta_c+1.0))); 
					}
					c2 = 0.5 * ((y1+y2)+betaq*(y2-y1));
					if(c1 < ylow)
						c1 = ylow;
					if (c2 < ylow)
						c2 = ylow;
					if (c1 > yupper)
						c1 = yupper;
					if (c2 > yupper)
						c2 = yupper;
					if(next_double() <= 0.5)
					{
						ind1->variables.at(i) = c2;
						ind2->variables.at(i) = c1;
					}
					else
					{
						ind1->variables.at(i) = c1;
						ind2->variables.at(i) = c2;
					}
				}
				else
				{
					ind1->variables.at(i) = parent1->variables.at(i);
					ind2->variables.at(i) = parent2->variables.at(i);
				}
			}
			else
			{
				ind1->variables.at(i) = parent1->variables.at(i);
				ind2->variables.at(i) = parent2->variables.at(i);
			}
		}
	}
	else
	{
		for(int i = 0; i < nvar; i++)
		{
			ind1->variables.at(i) = parent1->variables.at(i);
			ind2->variables.at(i) = parent2->variables.at(i);
		}
	}
}

/**
 * Funcion auxiliar que recibe como parametro dos individyos y se verifica la dominancia
 * @return bool si es que lo domina o no   
 */
bool IBEA::dominated(Individual* ind1, Individual* ind2)
{
	bool dominated = false;
	if((ind1->objectives.at(0) <= ind2->objectives.at(0) && ind1->objectives.at(1) < ind2->objectives.at(1))  \
		|| (ind1->objectives.at(0) < ind2->objectives.at(0) && ind1->objectives.at(1) <= ind2->objectives.at(1)))
	{
		dominated = true;
	}
	return dominated;
}

/**
 * Funcion que implementa un torneo binario entre dos individuos de la poblacion
 * @param  ind1 - Primer individuo que se quiere comparar
 * @param  ind2 - Segundo individuo que se quiere comparar 
 * @return      - individuo ganador 
 */
Individual* IBEA::tournament(Individual*ind1, Individual*ind2)
{
	//bool flag = dominated(ind1, ind2);
	if (dominated(ind1, ind2))
	{
		return ind1;
	}
	else if(dominated(ind2, ind1))
	{
		return ind2;
	}
	if(ind1->fitness < ind2->fitness)
	{
		return ind1;
	}
	if(ind2->fitness < ind1->fitness)
	{
		return ind2;
	}
	if(next_double() <= 0.5)
	{
		return ind1;
	}
	else
	{
		return ind2;
	}
}

/**
 * Funcion que se encarga de hacer la seleccion de los individuos 
 * dada la poblacion de los padres. 
 * Para ello utiliza un torneo y ademas aplica la cruza
 */
vector<Individual*> IBEA::selection(vector<Individual*> population)
{

	vector<Individual*> selected = generatePopulation(this->function, this->pop_size);
	int a1[this->pop_size];
	int a2[this->pop_size];
	int rnd, rnd2, temp;
	for(int i = 0; i < this->pop_size; i++)
	{
		a1[i] = i;
		a2[i] = i;
	}

	for(int i = 0; i < this->pop_size; i++)
	{
		rnd = rand() % (this->pop_size - i) + i;
		temp = a1[rnd];
		a1[rnd] = a1[i];
		a1[i] = temp;
		rnd2 = rand() % (this->pop_size - i) + i;
		temp = a2[rnd2];
		a2[rnd2] = a2[i];
		a2[i] = temp;
	}
	for(int i = 0; i < this->pop_size; i+=4)
	{
		Individual *parent1 = tournament(population.at(a1[i]), population.at(a1[i+1]));
		Individual *parent2 = tournament(population.at(a1[i+2]), population.at(a1[i+3]));
		crossoverInd(parent1, parent2, selected.at(i), selected.at(i+1));
		Individual *parent3 = tournament(population.at(a2[i]), population.at(a2[i+1]));
		Individual *parent4 = tournament(population.at(a2[i+2]), population.at(a2[i+3]));
		crossoverInd(parent3, parent4, selected.at(i+2), selected.at(i+3));
	}
	return selected;
}

/**
 * Funcion que se encarga de mutar un individuo de acuerdo a la probabilidad
 * de mutacion que se establecio.
 * @param ind - individuo que se desea mutar.
 */
void IBEA::mutateInd(Individual* ind)
{
	double pmut_r = this->mutation_prob;
	double eta_m = 50;
	double y, ylow, yupper, delta1, delta2, deltaq, val, xy, rnd, mut_pow;
	int nvar;
	if(this->function == "WFG2")
	{
		nvar = this->wfg2->k + this->wfg2->l;
	}
	else if(this->function == "WFG8")
	{
		nvar = this->wfg8->k + this->wfg8->l;
	}
	for(int j = 0; j < nvar; j++)
	{
		if(next_double() <= pmut_r)
		{
			y = ind->variables.at(j);
			ylow = 0.0;
			yupper = 2.0 * (j+1);
			delta1 = (y-ylow)/(yupper-ylow);
			delta2 = (yupper-y)/(yupper-ylow);
			rnd = next_double();
			mut_pow = 1.0/(eta_m+1.0);
			if(rnd <= 0.5)
			{
				xy = 1.0 - delta1;
				val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(eta_m+1.0)));
				deltaq = pow(val, mut_pow)- 1.0;
			}
			else
			{
				xy = 1.0 - delta2;
				val = 2.0 *(1.0-rnd)+ 2.0*(rnd-0.5)*(pow(xy,(eta_m+1.0)));
				deltaq = 1.0 - pow(val, mut_pow);
			}
			y = y + deltaq*(yupper-ylow);
			if (y < ylow)
				y = ylow;
			if (y > yupper)
				y = yupper;
			ind->variables.at(j) = y;
		}
	}
}

/**
 * Funcion que se encarga de aplicarle la mutacion a cada uno de los individuos de la poblacion
 * @param population - poblacion que se desea mutar
 */
void IBEA::mutation(vector<Individual*> *population)
{
	for(int i = 0; i < population->size(); i++)
	{
		mutateInd(population->at(i));
	}
}

/**
 * Funcion auxiliar que permite hacer la concatenacion de dos poblacion
 * la poblacion actual y la poblacion de hijos, por lo tanto la poblacion 
 * resultante sera de tamanio 2 veces el tamanio de la poblacion
 * @param population - poblacion a la que se quieren concatenar los elementos
 * @param new_pop    - poblacion que se quiere concatenar
 */
void IBEA::appendOffspring(vector<Individual*>*population, vector<Individual*> new_pop)
{
	for(int i = 0; i < new_pop.size(); i++)
	{
		population->push_back(new_pop.at(i));
	}
}

/**
 * Funcion que obtiene el mejor elemento de la poblacion, es decir, aquel elemento cuyo 
 * valor de fitness es menor 
 * @param  population - poblacion actual
 * @param  position   - indice del mejor elemento
 * @return            - regresa el elemento con el mejor fitness
 */
Individual* IBEA::getBest(vector<Individual*> population, int*position)
{
	double min_fitness = population.at(0)->fitness;
	Individual* best = population.at(0);
	for(int i = 0; i < population.size(); i++)
	{
		if (population.at(i)->fitness < min_fitness)
		{
			min_fitness = population.at(i)->fitness;
			best = population.at(i);
			*position = i;
		}
	}
	return best;
}

/**
 * FUncion que nos permite hacer la actualizacion de fitness de cada uno de los 
 * individuos de la poblacion, esto es que dado el elemento que eliminaremos
 * le debemos de sumar el valor de fitness que este aportaba al resto de los
 * individuos en la poblacion
 * @param best       - mejor elemento que se elimino
 * @param population - poblacion actual
 */
void IBEA::updateFitness(Individual* best, vector<Individual*> *population)
{
	for(int i = 0; i < population->size(); i++)
	{
		population->at(i)->fitness += exp(-epsilonIndicator(best, population->at(i))/ this->scale_factor);	
	}
}

/**
 * Funcion que permite obtener el conjunto de elementos no dominados de un conjunto
 * en este caso este calculo de forma cuadratica y se compara un elemento contra 
 * el resto de los elementos de la poblacion
 */
vector<Individual*> IBEA::getNonDominated(vector<Individual*> population)
{
	vector<Individual*> non_dominated;
	for(int i = 0; i < population.size(); i++)
	{
		for(int j = 0; j < population.size(); j++)
		{
			if(i != j)
			{	
				if((population.at(j)->objectives.at(0) < population.at(i)->objectives.at(0) && population.at(j)->objectives.at(1) <= population.at(i)->objectives.at(1))\
				   || (population.at(j)->objectives.at(0) <= population.at(i)->objectives.at(0) && population.at(j)->objectives.at(1) < population.at(i)->objectives.at(1)))
				{
					break;				
				}
			}
			if (j == population.size() - 1)
			{
				non_dominated.push_back(population.at(i));
			}
		}
	}
	return non_dominated;
}

/**
 * FUncion principal que se encarga de encontrar el conjunto de optimo de pareto para
 * intentar resolver una funcion
 * @para generations - maximo numero de generaciones, criterio de paro
 */
vector<Individual*> IBEA::solve(int generations)
{
	vector<Individual*> solution;
	vector<Individual*> P =  generatePopulation(this->function, this->pop_size);
	int generation = 0;
	while (true)
	{
		evaluatePopulation(&P);
		while (P.size() > this->pop_size)
		{
			int position = 0;
			Individual* best = getBest(P, &position);
			P.erase(P.begin() + position);
			updateFitness(best, &P);
		}
		if(generation >= generations)
		{
			solution = getNonDominated(P);
			break;	
		}
		vector <Individual*> P_aux = selection(P);
		mutation(&P_aux);
		appendOffspring(&P, P_aux);
		generation += 1;
	}
	return solution;
} 