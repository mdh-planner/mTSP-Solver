#pragma once
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include "definitions.h"

using namespace std;

class MODEL {

friend class GA;
friend class LOG;
friend class CPScheduler;
friend class vtScheduler;

public:

	MODEL(Problem &model);
	//MODEL(Problem &model, bool rModel );
	/* Get Source Depot Number */
	int getSourceDepotNumber() { return sigma; }

	/* Get Task Number */
	int getTaskNumber() { return V; }

	/* Get Destination Depot Number */
	int getDestinationDepotNumber() { return delta; }

	/* Get Weight Matrix */
	std::vector<std::vector<std::vector<int> > > getCostMatrix() const { return w; }
	//std::vector<std::vector<std::vector<double>>> getCostMatrix() const { return w; }
	
	/* Get Edge Weight */
	/*int getEdgeCost(int i, int j, int s) const { return w[i][j][s]; }*/
	double getEdgeCost(int i, int j, int s) const { return w[i][j][s]; }

	/* Get Task Duration */
	int getTaskDuration(int i);

	/* Get a vector of Agents that can performe a certain task */
	vector<int> findAgents(int Task);

	/* Get Agents : Tasks Matrix */
	std::vector<std::vector<int> > *getTasksPerAgents() { return &_tasksPerAgents; }

	/* Get vector of Agents : Tasks Matrix */
	std::vector<std::vector<int> > getListTasksPerAgents() { return _listTasksPerAgents; }

	/* Get vector of Agents : Tasks Matrix */
	std::vector<std::vector<int> > *getVectorTasksPerAgents() { return &_vectorTasksPerAgents; }

	/* Get Virtual Tasks */
	std::vector<int>  getVirtualTasks() const { return H; }

	/* Get Parallel Tasks */
	std::vector<std::vector<int> > getParallelMatrix() const { return R; }

	/* Get Precedence Constraints */
	std::vector<std::vector<int> > getPrecedenceMatrix() const { return P; }

	/* Get Agents per Task */
	std::vector<int>  getAgentsPerTask() const { return K; }

	/* Get the number of agents to use in the mission, this is added to use with Anders' work */
	int getNumberOfAgentsToUse() const { return nAU; }

	/* set range, for UVA work */
	void setRange(double _range)  { range = _range; }

	/* set velocity, for UVA work */
	void setAgentVelocityForAll(double _v)  { 
		for (int i = 0; i < _model.A.size(); i++) {
			_model.A[i].v = _v; 
		}
	}

	
	~MODEL() { std::cout << "MODEL object destroyed!" << std::endl; };
private:

	struct PC {
		vector<int> before;
		vector<int> after;
	};

	int nAU;
	double avgWeight, range;
	int sigma, V, delta, Vtilde;
	double _pVirtual, _pParallel;

	Problem _model;
	std::vector<int> K, H;
	std::vector<std::vector<int> > _PV, P, R, aPpT, _tasksPerAgents, _vectorTasksPerAgents;
	std::vector<std::vector<int> >_listTasksPerAgents;
	std::vector<std::vector<std::vector<int> > > w;
	//std::vector<std::vector<std::vector<double>>> w;
	std::vector<int> taskDuration;
	std::vector<PC> vpc;
	

	std::vector<std::pair<double, double> > Vtilde2;
};

