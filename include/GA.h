#pragma once
#include <numeric>
#include <fstream>
#include <utility>
#include <functional>
#include <unordered_set>
#include <algorithm>
#include <random>
#include <set>
#include <cassert>
#include <iomanip>
#include <limits>


#include "interactionRewards.h"
#include "model.h"
#include "chromosome.h"
#include "definitions.h"
#include "randnumgen.h"

// #include "tbb/tbb.h"
// #include "tbb/parallel_for.h"
// #include <tbb/blocked_range.h>
// #include <tbb/concurrent_vector.h>
// #include <tbb/spin_mutex.h>
// // // #include "tbb/compat/thread"
// #include "tbb/mutex.h"
// #include "tbb/queuing_mutex.h"


// Use to choose between serial and parallel execution
#define SERIAL

// choose between ID (intra-schedule dependencies) and XD (cross-schedule dependencies)
//#define XD

using namespace std;


class GA {

	friend class CPScheduler;
	friend class LOG;
	friend class vtScheduler;
	friend class REWARD;

public:
	/* Default Constructor */	
	GA( MODEL* _model, int maxGens);

	/* Constructor with objLimit*/
	GA(MODEL *_model, int maxGens, double objLimit);
	
	/* Constractor with provided values */
	GA( MODEL* _model, int _POPULATION_SIZE, double _pMUTATION, double _pXOVER, double _pELITE, int _XOVER_TYPE, int _OBJECTIVE);

	int POPULATION_SIZE;
	vector<CHROMOSOME> population;

	/// Debug variables
	int d_one, d_two, d_three, d_total; // debugging variables
	
	vector<pair<int, int> > _bias;

	void printIndividual(int idx);
	void evolveNextGeneration();
	void initializePopulation();
	void createSelectionPool();
	void evaluate(CHROMOSOME & chromo);
	void evaluateXD(CHROMOSOME & chromo);
	static void createGantt(CHROMOSOME & chromo, MODEL *model);
	double getTimeSteps(CHROMOSOME & chromo, MODEL *model,vector<vector<pair<double, double> > > &timeSteps, bool logFlag);
	void pcFix(CHROMOSOME & chromo);

	void getTimeline();
	int getNumEvolvedGenerations();
	int getBestIndividualsIndex();
	
	double getBestCost();
	double getMaxSteps();
	vector<pair<double, double> > ecdf;

	int getMutationNumber();

	/// Debugging Methods
	void printOrderedPlan(vector<int> &vector, int ag);
	void countTasks(CHROMOSOME & chromo);

	~GA() { std::cout << "GA object destroyed!" << std::endl; };

private:
	MODEL* model;

	int MUTATION_NUMBER, XOVER_NUMBER, ELITE;
	int CURRENT_GENERATION = 0;
	int _maxGens;


	double pMUTATION, pXOVER, pELITE, pVirtual, pParallel, _objLimit = 0;
	int XOVER_TYPE, OBJECTIVE, PENALTY = 1000000;
	
	/// Helper Variables
	int _halfSelectionPool, _parent1, _parent2, _mutationsDone;
	vector<int> _selectionPool, _splitPoints;
	
	std::vector<gantt> timeline;
	vector<CHROMOSOME> _childrenPopulation, _elitePopulation;
	std::vector<int> depots;
	vector<int> depotStart;
	int _bestIndividualsIndex;
	int INDIVIDUAL_SIZE;
	double _bestCost = std::numeric_limits<double>::max();

	void initializeChromosomeSRST(CHROMOSOME & chromo);
	void initializeAnders(CHROMOSOME &chromo);

	void xover();
	void mutate(CHROMOSOME & chromo);
	void elitismSelection();
	void localRefinement(CHROMOSOME & chromo);
	

	/// Helper Methods
	void updateSigma(vector<vector<int> > &_sigma, int rndT);
	void selectParentsFromPool();
	void selectSplitPoints();
	void withinAgentSwap(CHROMOSOME &chromo);
	void withinAgentJump(CHROMOSOME &chromo);
	void crossAgentSwap(CHROMOSOME &chromo);
	void crossAgentJump(CHROMOSOME &chromo);
	//void geneSwap(vector<int> &individual, int g1, int g2);
	void rankPop(vector<pair<int, double> > & ranking);
	void updatePopElite();
	void updatePopBreed();
	void pickDestinationDepot(vector<int> & chromo, int tmp, int i);

	void twoOptLocalSearch(CHROMOSOME & chromo);
	void GreedySearch(CHROMOSOME &chromo);
	void twoOptSwap(vector<int>& inGenes, vector<int>& outGenes, int iGene1, int iGene2);
	void extractPath(CHROMOSOME & chromo, vector<vector<int> > &path);
	void insertPath(vector<int> & tmp, vector<int> &path);
	void performPCfix(CHROMOSOME & chromo, int _agent1, int _agent2, int &_start, int _end, int _idxTask2Move, int _count);
	void ERX(CHROMOSOME & _tmpChromo, vector<unordered_set<int> > & _aMatrix);
	void genAdjMatrix(vector<unordered_set<int> > & _return);
	void allocateRest(vector<unordered_set<int> > & aMatrix, CHROMOSOME & chromo);
	void insertTask(CHROMOSOME &chromo, vector<unordered_set<int> > & aMatrix, int agent, int insertIdx, int task2Insert);

	int pickNextTask(int task, int agent, vector<unordered_set<int> > &);
	int getRandomFromSet(unordered_set<int> &);
	int checkColor(int Agent, int Task);
	int findTask(vector<int> &_vector, int ag, int start, int end, int taskIdx, int & tmpOld);
	int insertionIdx(vector<int> &_vector,int _start, int _end, int _tasksLeft);
	static pair<double, double> calculateAgentsPosition(MODEL *model, int t1, int t2, int a, int timeStep, bool &nextTask);

	pair<int, int> crossAgentTaskSearch(CHROMOSOME & chromo, int _agent, int _task, bool flag);

	double evaluateLight(vector<int> & individual, vector<int> & _tmpPC, int i, bool flag);
	double evaluateLightXD(vector<int> & individual, vector<gantt> &timeline, int i, bool flag);

};