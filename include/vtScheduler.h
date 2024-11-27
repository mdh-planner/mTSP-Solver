#pragma once
#include "include/GA.h"
#include "include/randnumgen.h"

struct triplet {
	int agentID, taskID, diff;
};

struct startTime {
	int earliest;
	int latest;
};

class vtScheduler
{
public:
	vtScheduler(shared_ptr<GA> Genetic, MODEL *_fullModel);

	void solver();

private:
	shared_ptr<GA> ga, ga2;
	MODEL* fullModel;
	CHROMOSOME chromo, originalChromo;
	vector<int> _possibleParallelTasks, _depots;
	vector<vector<int>> transitionGaps;
	std::vector<gantt> fullTimeline, originalFullTimeline;
	int mem_task;
	int count = 1;

	void findParallel(vector<int> &_parallelVector);
	void updatePlan(int _task, triplet &_tri, int ii, vector<int> &a2t_ii);
	void createTransitionGaps();
	void scheduleTask(int task);
	int pickTask(int _taskD);
	int checkColor(int ag, int _task);
	startTime checkPC(int ag, int _task);
	bool findTask(int _task, triplet &_tri, int ii);
	triplet findGap(vector<int> &a2t, int _task);
	int checkXD();
	/// Helper Methods
	pair<int, int> findPClimits(int _task);
	void updateTimeline(int _task, triplet &_tri, int _tmpVal);
	int insertTask(int _task, triplet &_tri);
	int findEarliestSchedulingTime(vector<int> &tasks2check, int _tmp);
	bool findLimits(triplet &_tri,vector<pair<int, int>> &startTimes, int task, int insertTask);
	void updateTimes(vector<pair<int, int>> & startTimes, int start, int end);
	/// Debug methods
	void printIndividual(vector<int> &_vector, int ag);
	void ganttLog(int count);
	
};



