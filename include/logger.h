// #pragma once
#include <iostream>
#include <string>
#include <cstdio>
#include <ctime>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <filesystem>

#include "GA.h"

using namespace std;
namespace filesys = std::__fs::filesystem;

class LOG
{
public:
	LOG(shared_ptr<GA> Genetic, int instance, int nRun);
	LOG(int instance);

	void writeOutput( double runTime, int run);
	void planLog(int instance, int run, vector<Problem> & model);
	void orderedPlanLog(int instance, int run, int id, vector<Problem> & model);
	void ganttLog(int instance, int run, vector<Problem> & model);
	void debug();
	void timeStepPlan(Problem &model, vector<vector<pair<double,double>>> &timeSteps , int nameID);
	void allocationOutput(Problem &model, int id, int nameID);
	void stepsPerAgent(Problem &model, int id, int nameID);
	string getPath() { return path; }

private:
	shared_ptr<GA> ga;
	string path;
	int _nRun;
	char dir_path[14] = "Test Results\\";

	void createDir(const char *name, char *dir_path);
	void createDirS(string bla);
	bool checkIfDirectory(std::string &filePath);
	string pathCreator(int instance);
	void instanceSettings();

};

