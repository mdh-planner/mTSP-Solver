#pragma once
#include <vector>
#include "../include/definitions.h"
struct penalties {
	std::vector<int> colorList;
	int color;
	int pc;
	int mt;
};


class CHROMOSOME {
public:
//	CHROMOSOME(MODEL &A);
	
	double rFitness;
	//size_t cost;
	double cost;
	bool feasibility;
	penalties penal;
	//std::vector<size_t> indCost;
	std::vector<double> indCost;
	std::vector<std::vector<int> > individual;
	std::vector<int> numTasks;
	std::vector < std::vector<int> > indTask;
	std::vector<gantt> timeline;
	std::vector<int> interPenalties;

	// variables for the work with Lauren
	int maxSteps;
	double weight;
	double distWeight;
	int interNum;
	int weightedInteractions;
	std::vector<std::vector<double> > rewardsPerAgent;
	//std::vector<pair<int,int>> depots;
	friend class GA;
};