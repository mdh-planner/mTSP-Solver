#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <iomanip>
#include <stdio.h>

#include "include/definitions.h"
#include "include/GA.h"
#include "include/logger.h"
#include "include/importProblem.h"


using namespace std;

vector<Problem> problemGenerator();

int main(int argc, const char *argv[])
{
	int GENERATIONS = 500;
	int maxRun = 1;
	int run = 0;
	int iStart = 0, iEnd = 1;

	double tmpCost = __INT_MAX__;
	double objLimit = 0;
	double runTime;
	double coeff = 1.2;
	// cout << "Enter the upper bound value of the solution quality compared to the best one (e.g., 1.2 - this equals to 20% worse): " << flush;
	// cin >> coeff;

	bool log = false;
	
	/*Update counter file*/
	int count = 0;
	ifstream infile("../counter.txt");
	infile >> count;
	ofstream count_file;
	count_file.open("../counter.txt");
	count_file << count + 1;
	count_file.close();

	// auto model = problemGenerator();

	vector<Problem> model = importProblem::importInstanceUVA();

	// Addition for UVA
	std::ifstream inputFile("../parameters.txt");
	double v, range;
	inputFile >> v >> range >> coeff;

	//
	/*=============== Loop Over Problem Instances ================*/

	for (int instance = iStart; instance < iEnd; instance++)
	{
		MODEL context(model[instance]);
		// for UVA
		
		context.setRange(range);
		context.setAgentVelocityForAll(v);
		//
		auto Genetic = shared_ptr<GA>(new GA(&context, GENERATIONS));

		LOG logger(Genetic, instance, run);

		for (; run < maxRun; run++)
		{
			clock_t tStart = clock();
			Genetic->initializePopulation();

			for (int i = 0; i < Genetic->POPULATION_SIZE; i++)
			{
				Genetic->pcFix(Genetic->population[i]);
				Genetic->evaluate(Genetic->population[i]); // for different problems like MRMT you will have to pick here which evaluation to use. If no virtual tasks - ST, if no parallel tasks - SR;
			}

			Genetic->createSelectionPool();

			// Start the optimization process
			for (int i = 0; i < GENERATIONS; i++)
			{
				Genetic->evolveNextGeneration();
				// vector<vector<pair<double, double>>> timeSteps(model[instance].A.size());
				// Genetic->getTimeSteps(Genetic->population[Genetic->getBestIndividualsIndex()], &context, timeSteps, true);
				if (i % static_cast<int>(std::round(static_cast<float>(GENERATIONS * 0.1))) == 0)
				{
					if (Genetic->population[Genetic->getBestIndividualsIndex()].feasibility)
					{
						cout << "Generation: " << left << setw(5) << Genetic->getNumEvolvedGenerations()
							 << "   Best Cost: " << left << setw(10) << fixed
							 << setprecision(4) << Genetic->getBestCost() << " #interactions: " << left
							 << setw(10) << fixed << setprecision(4) << Genetic->population[Genetic->getBestIndividualsIndex()].interNum 
							 << " Weight: " << left << setw(10) << fixed << setprecision(4) << Genetic->population[Genetic->getBestIndividualsIndex()].weight 
							 << " weightedInterNum: " << left << setw(10) << fixed << setprecision(4) << Genetic->population[Genetic->getBestIndividualsIndex()].weightedInteractions << "\n";
					}
					else
					{
						cout << endl;
						cout << "No feasible solution found. " << endl;
						cout << endl;
					}
				}
			}

			Genetic->printIndividual(Genetic->getBestIndividualsIndex());

			for (int i = 0; i < model[instance].A.size(); i++)
			{
				cout << Genetic->population[Genetic->getBestIndividualsIndex()].indTask[i].size() << endl;

				for (int j = 0; j < Genetic->population[Genetic->getBestIndividualsIndex()].indTask[i].size(); j++)
				{
					cout << Genetic->population[Genetic->getBestIndividualsIndex()].indTask[i][j] << "," << flush;
				}
				cout << endl;
			}

			cout << "Best Cost: " << Genetic->getBestCost() << endl;
			cout << endl;

			objLimit = Genetic->getBestCost();

			/// Logging stuff
			if (log)
			{
				// Plan Log
				logger.planLog(instance, run, model);
				// Cost, Time, and eCDF
				logger.writeOutput(runTime, run);
			}

			vector<vector<pair<double, double>>> timeSteps(model[instance].A.size());
			int bestIndividualIdx = Genetic->getBestIndividualsIndex();
			// int bestIndividualIdx = 1;
			double num = Genetic->getTimeSteps(Genetic->population[bestIndividualIdx], &context, timeSteps, true);

			logger.timeStepPlan(model[instance], timeSteps, 1);
			logger.stepsPerAgent(model[instance], bestIndividualIdx, 1);
			logger.allocationOutput(model[instance], bestIndividualIdx, 1);
		}

		run = 0;

		//--=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-==--=-=-=-==--=-=-==--=-=-=-=-==-=--=
		//--=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-==--=-=-=-==--=-=-==--=-=-=-=-==-=--=
		//--=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-==--=-=-=-==--=-=-==--=-=-=-=-==-=--=

		if (objLimit != -1)
		{

			auto Genetic2 = shared_ptr<GA>(new GA(&context, GENERATIONS, Genetic->getMaxSteps() * coeff));
			LOG logger2(Genetic2, instance, run);

			for (; run < maxRun; run++)
			{
				clock_t tStart = clock();
				Genetic2->initializePopulation();

				// Genetic2->population[0] = Genetic->population[Genetic->getBestIndividualsIndex()];
				Genetic2->population = Genetic->population;

				for (int i = 0; i < Genetic2->POPULATION_SIZE; i++)
				{
					Genetic2->pcFix(Genetic2->population[i]);
					Genetic2->evaluate(Genetic2->population[i]); // for different problems like MRMT you will have to pick here which evaluation to use. If no virtual tasks - ST, if no parallel tasks - SR;
				}

				Genetic2->createSelectionPool();

				// Start the optimization process
				for (int i = 0; i < GENERATIONS; i++)
				{
					Genetic2->evolveNextGeneration();

					if (i % static_cast<int>(std::round(static_cast<float>(GENERATIONS * 0.1))) == 0)
					{
						if (Genetic2->population[Genetic2->getBestIndividualsIndex()].feasibility)
						{
							cout << "Generation: " << left << setw(5) << Genetic2->getNumEvolvedGenerations()
								 << "   Best Cost: " << left << setw(10) << fixed
								 << setprecision(4) << Genetic2->getBestCost() << " #interactions: " << left
								 << setw(10) << fixed << setprecision(4) << Genetic2->population[Genetic2->getBestIndividualsIndex()].interNum << " Weight: " << left << setw(10) << fixed << setprecision(4) << Genetic2->population[Genetic2->getBestIndividualsIndex()].weight
								 << " weightedInterNum: " << left << setw(10) << fixed << setprecision(4) << Genetic2->population[Genetic2->getBestIndividualsIndex()].weightedInteractions << "\n";
						}
						else
						{
							cout << endl;
							cout << "No feasible solution found. " << endl;
							cout << endl;
						}
					}
				}
			}

			vector<vector<pair<double, double>>> timeSteps(model[instance].A.size());
			int bestIndividualIdx = Genetic2->getBestIndividualsIndex();
			// int bestIndividualIdx = 1;
			double num = Genetic2->getTimeSteps(Genetic2->population[bestIndividualIdx], &context, timeSteps, true);

			logger2.timeStepPlan(model[instance], timeSteps, 2);
			logger2.stepsPerAgent(model[instance], bestIndividualIdx, 2);
			logger2.allocationOutput(model[instance], bestIndividualIdx, 2);

		}
		return 0;
	}
 }
// Put this out of the loops.
// remove("c:\\Users\\bmc01\\Documents\\Visual Studio 2017\\Projects\\MR_MT_Solver\\MR_MT_Solver\\resultsgif.txt");

// if (Genetic->getBestCost() < tmpCost) {
//	//remove("c:\\Users\\bmc01\\Documents\\Visual Studio 2017\\Projects\\MR_MT_Solver\\MR_MT_Solver\\resultsgif.txt");
//	ofstream myfile("resultsgif.txt", ios_base::app | ios_base::out);

//	for (int i = 0; i < Genetic->population[Genetic->getBestIndividualsIndex()].individual.size(); ++i) {
//		int k = i;

//		myfile << i << "\t" << model[instance].A[i].x << "\t" << model[instance].A[i].y << endl;
//		while (true) {

//			int task = Genetic->population[Genetic->getBestIndividualsIndex()].individual[i][k] - model[instance].A.size();

//			if (task == model[instance].T.size()) {

//				myfile << task + model[instance].A.size() << "\t" << model[instance].dest[task - model[instance].T.size()].x << "\t" << model[instance].T[task - model[instance].T.size()].y << endl;
//				break;
//			}
//			else {
//				myfile << task + model[instance].A.size() << "\t" << model[instance].T[task].x << "\t" << model[instance].T[task].y << endl;
//			}

//			k = task + model[instance].A.size();

//			//cout << Genetic->population[Genetic->getBestIndividualsIndex()].individual[i][k] << model[instance].A.size() + model[instance].T.size() << endl;

//		}
//		myfile << "END" << endl;
//	}
//	myfile << "--" << endl;
//	myfile.close();
//	tmpCost = Genetic->getBestCost();
//}

// Genetic->population[0].individual[0] = { 31,    27,    12   , 35 ,   54 ,    0 ,   60,    57,    24 ,   56 ,   36  ,  26  ,  59  ,  50  ,  49  ,   4  ,  10  ,   2  ,  70    ,25  ,  58  ,  29  ,  45  ,  47  ,  53  ,   9};
// Genetic->population[0].individual[1] = { 0   ,9   ,0   ,0   ,0   ,0   ,0   ,0   ,10  ,21  ,13  ,14  ,0   ,51  ,18  ,0   ,0   ,20  ,27  ,0   ,52  ,8   ,0   ,0   ,0   ,32  ,0   ,39  ,0   ,0   ,0   ,0   ,11  ,0   ,17  ,0   ,0   ,0   ,0   ,54  ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,34  ,25 };
// Genetic->population[0].individual[2] = { 0   ,0   ,28  ,0   ,0   ,36  ,49  ,6   ,0   ,0   ,0   ,0   ,43  ,0   ,0   ,40  ,7   ,0   ,0   ,0   ,0   ,0   ,48  ,45  ,15  ,0   ,5   ,0   ,24  ,0   ,47  ,0   ,0   ,0   ,0   ,0   ,23  ,12  ,0   ,0   ,16  ,0   ,0   ,22  ,0   ,54  ,0   ,37  ,26  ,30  ,0   ,0   ,0 };

/*	Genetic->population[0].individual[0] = { 80,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	Genetic->population[0].individual[1] = { 0,35,0,0,64,0,74,0,14,0,0,0,0,80,6,0,39,0,0,0,0,0,0,0,0,0,0,0,60,62,63,16,0,49,0,31,0,0,0,58,30,0,0,0,0,0,0,0,0,51,0,57,0,8,53,0,0,13,4,0,40,28,33,54,61,0,0,0,0,0,0,0,0,0,29,0,0,0,0 };
	Genetic->population[0].individual[3] = { 0,0,0,19,0,55,0,0,0,0,0,0,67,0,0,0,0,0,0,5,41,44,12,78,0,20,0,0,0,0,0,0,0,0,56,0,0,0,75,0,0,21,0,0,70,80,0,0,38,0,0,0,0,0,0,22,48,0,0,0,0,0,0,0,0,0,0,34,0,0,45,0,25,0,0,23,0,0,72 };
	Genetic->population[0].individual[2] = { 0,0,43,0,0,0,0,9,0,52,27,76,0,0,0,71,0,18,79,0,0,0,0,0,66,0,50,37,0,0,0,0,59,0,0,0,69,26,0,0,0,0,77,47,0,0,17,73,0,0,42,0,68,0,0,0,0,0,0,11,0,0,0,0,0,36,65,0,15,46,0,10,0,7,0,0,24,32,0 };
*/

// Genetic->printIndividual(0);