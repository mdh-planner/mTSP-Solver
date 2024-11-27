#include "include/logger.h"

// Constructor
LOG::LOG(shared_ptr<GA> Genetic, int instance, int nRun)
{
	ga = Genetic;
	_nRun = nRun;

	string directory_path(dir_path);
	bool resDir = checkIfDirectory(directory_path);

	if (!resDir)
	{
		createDir("", dir_path);
	}

	path = pathCreator(instance);
	instanceSettings();
}

LOG::LOG(int instance)
{
	string directory_path(dir_path);
	bool resDir = checkIfDirectory(directory_path);

	if (!resDir)
	{
		createDir("", dir_path);
	}

	path = pathCreator(instance);
}

// Public:

void LOG::writeOutput(double runTime, int run)
{
	string path3 = path + to_string(run) + ".txt";
	ofstream myfile(path3, ios_base::app | ios_base::out);

	// for (int i = 0; i < ga->population[ga->getBestIndividualsIndex()].indCost.size(); i++) {
	// myfile << "Cost " << ga->population[ga->getBestIndividualsIndex()].indCost[i] << endl;
	//}

	myfile << "Cost " << ga->population[ga->getBestIndividualsIndex()].cost << endl;
	myfile << "Run Time " << runTime << endl;
	myfile.close();

	// LOG eCDF, i.e., log best cost on every generation

	string path2 = path + to_string(run) + "_eCDF.txt";
	ofstream myfile2(path2, ios_base::app | ios_base::out);

	for (int ecdfRun = 0; ecdfRun < ga->ecdf.size(); ecdfRun++)
	{
		myfile2 << ga->ecdf[ecdfRun].first << "\t" << ga->ecdf[ecdfRun].second << endl;
	}

	myfile2.close();
}

void LOG::orderedPlanLog(int instance, int run, int id, vector<Problem> &model)
{

	string path3 = path + to_string(run) + "_plan.txt";
	ofstream myfile(path3, ios_base::app | ios_base::out);

	auto sigma = model[instance].A.size();
	auto V = model[instance].T.size();
	auto Vtilde = sigma + V + model[instance].dest.size();

	for (int i = 0; i < ga->population[id].individual.size(); i++)
	{

		auto vector = ga->population[id].individual[i];

		int _k = i;
		int count = 0;
		// myfile << _k << "\t" << 0 << "\t" << flush;

		while (true)
		{

			if (vector[_k] >= sigma + V)
			{
				myfile << vector[_k] << "\t" << flush;
				break;
			}
			else
			{
				myfile << vector[_k] + 1 << "," << flush;
			}

			_k = vector[_k];

			if (count > Vtilde)
			{
				break;
			}
			count++;
		}
		myfile << endl;
	}
}

void LOG::planLog(int instance, int run, vector<Problem> &model)
{

	string path2 = path + to_string(run) + "_drawTasks.txt";
	ofstream myfile(path2, ios_base::app | ios_base::out);

	for (int i = 0; i < ga->population[ga->getBestIndividualsIndex()].individual.size(); ++i)
	{
		size_t k = i;

		myfile << i << "\t" << model[instance].A[i].x << "\t" << model[instance].A[i].y << endl;

		while (true)
		{

			int task = ga->population[ga->getBestIndividualsIndex()].individual[i][k] - model[instance].A.size();

			if (task >= model[instance].T.size())
			{
				myfile << task + model[instance].A.size() << "\t" << model[instance].dest[task - model[instance].T.size()].x << "\t" << model[instance].dest[task - model[instance].T.size()].y << endl;
				break;
			}
			else
			{
				myfile << task + model[instance].A.size() << "\t" << model[instance].T[task].x << "\t" << model[instance].T[task].y << endl;
			}

			k = task + model[instance].A.size();
		}
		myfile << "END" << endl;
	}
	myfile << "--" << endl;
	myfile.close();
}

void LOG::ganttLog(int instance, int run, vector<Problem> &model)
{
	// write tasks
	string path2 = path + to_string(run) + "_drawGantt.txt";
	// remove("c:\\Users\\bmc01\\Documents\\Visual Studio 2017\\Projects\\MR_MT_Solver\\MR_MT_Solver\\drawGantt.txt");
	// string path2 = ("c:\\Users\\bmc01\\Documents\\Visual Studio 2017\\Projects\\MR_MT_Solver\\MR_MT_Solver\\drawGantt.txt");
	ofstream myfile(path2, ios_base::app | ios_base::out);

	for (int i = 0; i < ga->timeline.size(); i++)
	{
		myfile << ga->timeline[i].startTime << "\t" << ga->timeline[i].endTime << "\t" << flush;
		myfile << ga->timeline[i].agentIdx << "\t" << model[instance].T[i].virt << "\t" << model[instance].T[i].index << endl;
	}

	// write destination depots
	string pathD = path + to_string(run) + "_drawGantt_depot.txt";
	// remove("c:\\Users\\bmc01\\Documents\\Visual Studio 2017\\Projects\\MR_MT_Solver\\MR_MT_Solver\\drawGantt_depot.txt");
	// string pathD = ("c:\\Users\\bmc01\\Documents\\Visual Studio 2017\\Projects\\MR_MT_Solver\\MR_MT_Solver\\drawGantt_depot.txt");
	ofstream myfile2(pathD, ios_base::app | ios_base::out);

	for (int i = 0; i < ga->depotStart.size(); i++)
	{
		myfile2 << ga->depotStart[i] << endl;
	}

	myfile2.close();
}

void LOG::debug()
{
	remove("c:\\Users\\bmc01\\Documents\\Visual Studio 2017\\Projects\\MR_MT_Solver\\MR_MT_Solver\\debug.txt");
	ofstream myfile2("debug.txt", ios_base::app | ios_base::out);

	for (int i = 0; i < ga->_bias.size(); i++)
	{
		myfile2 << ga->_bias[i].first << "\t" << ga->_bias[i].second << endl;
	}
	myfile2.close();
}

void LOG::timeStepPlan(Problem &model, vector<vector<pair<double,double>>> & timeSteps, int nameID )
{
	string path = "timesteps" + to_string(nameID) + ".txt";
	remove(path.c_str());
	
	ofstream myfile(path, ios_base::app | ios_base::out);
	for (int column = 0; column < timeSteps[0].size(); column++)
	{
		for (int row = 0; row < model.A.size(); row++)
		{

			myfile << timeSteps[row][column].first << "\t" << timeSteps[row][column].second << "\t" << flush;
		}
		myfile << endl;
	}
}

void LOG::allocationOutput (Problem &model, int id, int nameID) {

	string path = "allocation" + to_string(nameID) + ".txt";
	remove(path.c_str());
	ofstream myfile(path, ios_base::app | ios_base::out);

	auto sigma = model.A.size();
	auto V = model.T.size();
	auto Vtilde = sigma + V + model.dest.size();

	for (int i = 0; i < ga->population[id].individual.size(); i++)
	{

		auto vector = ga->population[id].individual[i];

		int _k = i;
		int count = 0;
		// myfile << _k << "\t" << 0 << "\t" << flush;
		std::vector<int> tmp(V+1);
		
		while (true)
		{

			if (vector[_k] >= sigma + V)
			{
				//myfile << vector[_k] << "\t" << flush;
				tmp[vector[_k]-sigma] = count+1;
				break;
			}
			else
			{
				//myfile << vector[_k] + 1 << "," << flush;
				tmp[vector[_k]-sigma] = count+1;
			}

			_k = vector[_k];

			if (count > Vtilde)
			{
				break;
			}
			count++;
		}
		for (int ii = 0; ii < tmp.size(); ii++){
			myfile << tmp[ii] << "\t" << flush;
		}
		myfile << endl;
	}

}

void LOG::stepsPerAgent (Problem &model, int id, int nameID) {
	
	string path = "stepspa" + to_string(nameID) + ".txt";
	remove(path.c_str());

	ofstream myfile(path, ios_base::app | ios_base::out);
	for (int i = 0; i < ga->population[id].rewardsPerAgent[0].size(); i++)
	{
		for (int j = 0; j < ga->population[id].rewardsPerAgent.size(); j++) {
				myfile << ga->population[id].rewardsPerAgent[j][i] << "\t"  << flush;	
		}
		myfile << endl;
	}
	//myfile << ga->model->avgWeight << endl;
}
// Private:

void LOG::createDir(const char *name, char *dir_path)
{

	// strcat(dir_path, path);
	strcat(dir_path, name);
	filesys::path dir(dir_path);

	if (filesys::create_directory(dir))
	{
		cout << endl;
		std::cout << "Directory has been Successfully created"
				  << "\n";
		cout << endl;
	}
}

void LOG::createDirS(string bla)
{
	// strcat(dir_path, path);
	// strcat(dir_path, name);
	filesys::path dir(bla);

	if (filesys::create_directory(dir))
	{
		cout << endl;
		std::cout << "Directory has been Successfully created"
				  << "\n";
		cout << endl;
	}
}

bool LOG::checkIfDirectory(std::string &filePath)
{
	try
	{
		// Create a Path object from given path string
		filesys::path pathObj(filePath);
		// Check if path exists and is of a directory file
		if (filesys::exists(pathObj) && filesys::is_directory(pathObj))
			return true;
	}
	catch (filesys::filesystem_error &e)
	{
		std::cerr << e.what() << std::endl;
	}
	return false;
}

string LOG::pathCreator(int instance)
{

	int count;
	ifstream infile("counter.txt");
	infile >> count;

	string path = dir_path + to_string(count) + "\\";

	if (!checkIfDirectory(path))
	{
		createDirS(path);
	}

	path += "Instance_" + to_string(instance);

	if (!checkIfDirectory(path))
	{
		createDirS(path);
	}

	return path + "\\";
};

void LOG::instanceSettings()
{

	ofstream myfile(path + "instanceSettings.txt", ios_base::app | ios_base::out);

	myfile << "Number_of_tasks " << ga->model->V << endl;
	myfile << "Number_of_runs " << _nRun << endl;
	myfile << "Population_size " << ga->POPULATION_SIZE << endl;
	myfile << "Number_of_generations " << ga->_maxGens << endl;
	myfile << "Chromosome_length " << ga->INDIVIDUAL_SIZE << endl;
	myfile << "Crossover_probability " << ga->pXOVER << endl;
	myfile << "Mutation_probability " << ga->pMUTATION << endl;
	myfile << "Kept_individuals " << ga->pELITE << endl;
	myfile << "Crossover_type " << ga->XOVER_TYPE << endl;
	// myfile << "Benchmark_Seed " << ga.seed << endl;
	// myfile << "Instance_Seed " << instanceSeed << endl;
	// myfile << "Run_Seed	" << flush;
	myfile.close();
}
