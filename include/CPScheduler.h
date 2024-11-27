#pragma once
#include "include/GA.h"

// Magic tricks to have CPLEX behave well:
#ifndef IL_STD
#define IL_STD
#endif

#include <cstring>
#include "ilcp/cp.h"

ILOSTLBEGIN
// End magic tricks

class CPScheduler 
{

public:
	CPScheduler(shared_ptr<GA> Genetic, MODEL *_model,string &path, int run, double runTime, double _timeLimit);
	CPScheduler(MODEL *_model, string &path, double _timeLimit);
	void solverFull();
	void solverPartial();
	void solverFull_warmStart();
private:
	shared_ptr<GA> ga;
	MODEL* fullModel;

	vector<vector<int>> _fixedAgents2Tasks;
	string path;
	IloInt weight, run;
	double timeLimit, timeLimitPartial;

};

