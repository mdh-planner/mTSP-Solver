#include "include/GA.h"

#pragma region == Constructor ==

GA::GA(MODEL *_model, int maxGens)
{
	model = _model;
	POPULATION_SIZE = 500;
	pMUTATION = 0.15;
	pXOVER = 0.7;
	pVirtual = model->_pVirtual;
	MUTATION_NUMBER = (int)(POPULATION_SIZE * pMUTATION);
	XOVER_NUMBER = (int)(POPULATION_SIZE * pXOVER);
	pELITE = 0.05;
	ELITE = (int)(1 + POPULATION_SIZE * pELITE);
	XOVER_TYPE = 0;
	OBJECTIVE = 1;
	_elitePopulation.resize(ELITE);
	INDIVIDUAL_SIZE = model->sigma + model->V;
	_maxGens = maxGens;
	// ecdf.resize(maxGens);
}

GA::GA(MODEL *_model, int maxGens, double objLimit)
{
	model = _model;
	POPULATION_SIZE = 200;
	pMUTATION = 0.15;
	pXOVER = 0.7;
	pVirtual = model->_pVirtual;
	MUTATION_NUMBER = (int)(POPULATION_SIZE * pMUTATION);
	XOVER_NUMBER = (int)(POPULATION_SIZE * pXOVER);
	pELITE = 0.05;
	ELITE = (int)(1 + POPULATION_SIZE * pELITE);
	XOVER_TYPE = 0;
	OBJECTIVE = 1;
	_elitePopulation.resize(ELITE);
	INDIVIDUAL_SIZE = model->sigma + model->V;
	_maxGens = maxGens;
	_objLimit = objLimit;
	// ecdf.resize(maxGens);
}

GA::GA(MODEL *_model, int _POPULATION_SIZE, double _pMUTATION, double _pXOVER, double _pELITE, int _XOVER_TYPE, int _OBJECTIVE)
{

	model = _model;
	POPULATION_SIZE = _POPULATION_SIZE;
	pMUTATION = _pMUTATION;
	pXOVER = _pXOVER;
	MUTATION_NUMBER = (int)(POPULATION_SIZE * pMUTATION);
	XOVER_NUMBER = (int)(POPULATION_SIZE * pXOVER);
	pELITE = _pELITE;
	ELITE = (int)(1 + POPULATION_SIZE * pELITE);
	XOVER_TYPE = _XOVER_TYPE;
	OBJECTIVE = _OBJECTIVE;
	_elitePopulation.resize(ELITE);
}

#pragma endregion

#pragma region == Public Methods ==

void GA::evolveNextGeneration()
{

	// cout << "CURRENT_GENERATION: " << CURRENT_GENERATION << endl;

	xover();

	for (int i = 0; i < _childrenPopulation.size(); i++)
	{

		countTasks(_childrenPopulation[i]);
	}

#ifdef SERIAL
	///* Serial Execution*/
	for (int i = 0; i < _childrenPopulation.size(); i++)
	{

		// cout << i << endl;
#ifdef XD
		createGantt(_childrenPopulation[i]);
		evaluateXD(_childrenPopulation[i]);
#else
		pcFix(_childrenPopulation[i]);
		evaluate(_childrenPopulation[i]); // for different problems like MRMT you will have to pick here which evaluation to use. If no virtual tasks - ST, if no parallel tasks - SR;
#endif
	}
	updatePopBreed();

	for (int i = 0; i < POPULATION_SIZE; i++)
	{
		// cout << "i: " << i << endl;
		mutate(population[i]);
		localRefinement(population[i]);

#ifdef XD
		createGantt(population[i]);
		evaluateXD(population[i]);
#else
		pcFix(population[i]);
		evaluate(population[i]); // for different problems like MRMT you will have to pick here which evaluation to use. If no virtual tasks - ST, if no parallel tasks - SR;
		countTasks(population[i]);
#endif
	}
#else
	///* Parallel Execution*/
	tbb::parallel_for(tbb::blocked_range<int>(0, _childrenPopulation.size()), [&](const tbb::blocked_range<int> &x)
					  {
		for (int i = x.begin(); i < x.end(); i++) {
#ifdef XD
			createGantt(_childrenPopulation[i]);
			evaluateXD(_childrenPopulation[i]);
#else
			pcFix(_childrenPopulation[i]);
			evaluate(_childrenPopulation[i]); // for different problems like MRMT you will have to pick here which evaluation to use. If no virtual tasks - ST, if no parallel tasks - SR;
#endif
			//if (_childrenPopulation[i].feasibility == false) {
			//	cout << "_childrenPopulation[i].feasibility == false" << endl;
			//}
		} });
	updatePopBreed();

	tbb::parallel_for(tbb::blocked_range<int>(0, POPULATION_SIZE), [&](const tbb::blocked_range<int> &x)
					  {
		for (int i = x.begin(); i < x.end(); i++) {

			//cout << "i " << i << endl;

			mutate(population[i]);

			//if (i % 2 == 0) {
			localRefinement(population[i]);
			//}

#ifdef XD
			createGantt(population[i]);
			evaluateXD(population[i]);
#else
			pcFix(population[i]);
			evaluate(population[i]); // for different problems like MRMT you will have to pick here which evaluation to use. If no virtual tasks - ST, if no parallel tasks - SR;
#endif
		} });

#endif // EXE

	updatePopElite();

	createSelectionPool();

	/*ecdf[CURRENT_GENERATION].first = _bestCost;*/
	CURRENT_GENERATION++;
}

void GA::getTimeline()
{

	timeline.resize(model->V);
	int cost = 0, memCost = 0, memTask;
	vector<int> tmpV;
	// cout << "Timeline formation." << endl;
	for (int i = 0; i < model->sigma; i++)
	{
		// printOrderedPlan(population[_bestIndividualsIndex].individual[i], i);
		int _tmp = i;
		int cost = 0;
		int endT = 0;
		while (_tmp < population[_bestIndividualsIndex].individual[i].size())
		{
			/* Calculate the cost of the chromosome */
			gantt tmpGantt;
			if (population[_bestIndividualsIndex].individual[i][_tmp] >= population[_bestIndividualsIndex].individual[i].size())
			{
				break;
			}
			// cout << _tmp << "\t" << population[_bestIndividualsIndex].individual[i][_tmp] << endl;

			if (model->H[_tmp] == 0 && model->H[population[_bestIndividualsIndex].individual[i][_tmp]] == 1)
			{
				cost = 0;
				memCost = 0;
				memTask = _tmp;
			}
			else if (model->H[_tmp] == 1 && model->H[population[_bestIndividualsIndex].individual[i][_tmp]] == 0)
			{
				memCost += model->getTaskDuration(_tmp);
				if (memCost > model->getEdgeCost(memTask, population[_bestIndividualsIndex].individual[i][_tmp], i))
				{
					cost = 0;
				}
				else
				{
					cost = model->getEdgeCost(memTask, population[_bestIndividualsIndex].individual[i][_tmp], i) - memCost;
				}
				memCost = 0;
			}
			else if (model->H[_tmp] == 0 && model->H[population[_bestIndividualsIndex].individual[i][_tmp]] == 0)
			{
				cost = model->getEdgeCost(_tmp, population[_bestIndividualsIndex].individual[i][_tmp], i);
			}
			else if (model->H[_tmp] == 1 && model->H[population[_bestIndividualsIndex].individual[i][_tmp]] == 1)
			{
				cost = 0;
				memCost += model->getTaskDuration(_tmp);
			}

			cout << " |" << cost << "| |" << model->getTaskDuration(population[_bestIndividualsIndex].individual[i][_tmp]) << "| " << flush;
			// tmpGantt.startTime = endT + model->getEdgeCost(_tmp, population[_bestIndividualsIndex].individual[i][_tmp], i);
			tmpGantt.startTime = endT + cost;
			tmpGantt.endTime = tmpGantt.startTime + model->getTaskDuration(population[_bestIndividualsIndex].individual[i][_tmp]);

			// cost += model->getEdgeCost(_tmp, population[_bestIndividualsIndex].individual[i][_tmp], i) + model->getTaskDuration(_tmp);

			tmpGantt.agentIdx = i;
			tmpGantt.virtualTask = model->H[population[_bestIndividualsIndex].individual[i][_tmp]];

			tmpV.emplace_back(population[_bestIndividualsIndex].individual[i][_tmp]);
			timeline[population[_bestIndividualsIndex].individual[i][_tmp] - model->sigma] = tmpGantt;

			// Update tmp
			_tmp = population[_bestIndividualsIndex].individual[i][_tmp];
			endT = tmpGantt.endTime;
		}
		if (!tmpV.empty())
		{

			for (int j = 0; j < tmpV.size(); j++)
			{
				// cout << "j: " << j << "tmpV[j] " << tmpV[j] << endl;
				timeline[tmpV[j] - model->sigma].taskIdx = tmpV[j] - model->sigma;
			}
		}

		depots.emplace_back(population[_bestIndividualsIndex].individual[i][_tmp] - model->V - model->sigma);
		depotStart.emplace_back(population[_bestIndividualsIndex].indCost[i]);
		// population[_bestIndividualsIndex].depots.emplace_back(make_pair(_tmp, population[_bestIndividualsIndex].indCost[i]));
	}
}

int GA::getNumEvolvedGenerations()
{
	return CURRENT_GENERATION;
}

int GA::getMutationNumber()
{
	return _mutationsDone;
}

double GA::getBestCost()
{
	return population[_bestIndividualsIndex].cost;
	// return accumulate(population[_bestIndividualsIndex].indCost.begin(), population[_bestIndividualsIndex].indCost.end(), 0);
	//  return *max_element(population[_bestIndividualsIndex].indCost.begin(), population[_bestIndividualsIndex].indCost.end());
	// return  *max_element(population[_bestIndividualsIndex].indCost.begin(), population[_bestIndividualsIndex].indCost.end()) + 0.1 * accumulate(population[_bestIndividualsIndex].indCost.begin(), population[_bestIndividualsIndex].indCost.end(), 0);
}

double GA::getMaxSteps()
{
	return population[_bestIndividualsIndex].maxSteps;
}

int GA::getBestIndividualsIndex()
{
	return _bestIndividualsIndex;
}

void GA::printIndividual(int idx)
{

	for (int i = 0; i < population[idx].individual.size(); i++)
	{
		cout << endl;
		cout << "||Agent " << i << "||" << endl;
		for (auto a = 0; a < population[idx].individual[0].size(); a++)
		{
			cout << setw(4) << a << flush;
		}

		cout << endl;

		for (int j = 0; j < population[idx].individual[i].size(); j++)
		{
			cout << setw(4) << population[idx].individual[i][j] << "," << flush;
		}

		cout << endl;
		cout << endl;
	}

	cout << endl;
}

double GA::getTimeSteps(CHROMOSOME &chromo, MODEL *model, vector<vector<pair<double, double>>> &timeSteps, bool logFlag)
{

	//int range = 225; // meters squared, 15m
	double range = model->range;
	double coeff = 1;
	// if (_objLimit > 0) {
	// 	coeff = 0.6;
	// } 
	
	double startVal = 0.15;
	double endVal = 0.8;
	int maxSteps = 0;
	chromo.interNum = 0;
	chromo.weight = 0;
	chromo.weightedInteractions = 0;
	vector<int> stepsPA(model->sigma);
	chromo.interPenalties.resize(model->sigma);

	for (int tickAg = 0; tickAg < model->sigma; tickAg++)
	{
		timeSteps[tickAg].resize(chromo.indCost[tickAg]);
		int _tmp = tickAg;
		int localStep = 1; // resets every time a task is reached
		bool nextTask = false;

		for (int tick = 0; tick < chromo.indCost[tickAg]; tick++)
		{
			if (_tmp >= chromo.individual[tickAg].size())
			{
				timeSteps[tickAg].resize(tick - 1);
				stepsPA[tickAg] = tick - 1;
				break;
			}

			pair<double, double> coord = calculateAgentsPosition(model, _tmp, chromo.individual[tickAg][_tmp], tickAg, localStep, nextTask);
			timeSteps[tickAg][tick] = coord;
			localStep++;

			if (nextTask)
			{
				// Update tmp
				_tmp = chromo.individual[tickAg][_tmp];
				localStep = 1;
				nextTask = false;
			}
		}

		if (timeSteps[tickAg].size() > maxSteps)
		{
			maxSteps = timeSteps[tickAg].size();
		}
	}

	vector<vector<int>> change(timeSteps.size(), vector<int>(timeSteps.size())); // Defaults to zero initial value
	vector<shared_ptr<REWARD>> _vRew(model->sigma);
	vector<vector<int>> inContact(maxSteps, vector<int>(3 * model->sigma));

	// making it a matrix for matlab
	for (int i = 0; i < model->sigma; i++)
	{
		timeSteps[i].resize(maxSteps);
		// defining start and end of the reward distribution
		int start = startVal * stepsPA[i]; //
		int end = endVal * stepsPA[i];	   //
		// inst reward objects for each agent
		_vRew[i] = shared_ptr<REWARD>(new REWARD(stepsPA[i], start, end, model->avgWeight));
	}

	vector<bool> halftimeFlag(model->sigma);
	for (auto a : halftimeFlag)
	{
		a = false;
	}

	for (int a = 0; a < timeSteps.size() - 1; a++)
	{
		for (int b = a + 1; b < timeSteps.size(); b++)
		{
			for (int t = 0; t < maxSteps; t++)
			{
				if (t < stepsPA[a] && t < stepsPA[b])
				{
					double dist = (timeSteps[a][t].first - timeSteps[b][t].first) * (timeSteps[a][t].first - timeSteps[b][t].first)
								 + (timeSteps[a][t].second - timeSteps[b][t].second) * (timeSteps[a][t].second - timeSteps[b][t].second);
					//cout << t << " " << (stepsPA[a] / 2) + chromo.interPenalties[a] << "  " << halftimeFlag[a] << endl;
					if (t > ((stepsPA[a] * coeff) + chromo.interPenalties[a]) && halftimeFlag[a] == false)
					{
						chromo.interPenalties[a]++;
					}
					else if (t > ((stepsPA[b] * coeff) + chromo.interPenalties[b]) && halftimeFlag[b] == false)
					{
						chromo.interPenalties[b]++;
					}

					if (dist <= range)
					{
						inContact[t][a] = a + 1;
						inContact[t][b] = b + 1;
						inContact[t][2 * model->sigma + a] = 1;
						inContact[t][2 * model->sigma + b] = 1;

						if (change[a][b] == 0)
						{
							double reward = _vRew[a]->interactionReward(t);
							reward += _vRew[b]->interactionReward(t);

							if (reward > 0)
							{
								chromo.weightedInteractions++;

								if (t < stepsPA[a] * coeff)
								{
									halftimeFlag[a] = true;
								}
								else if (t < stepsPA[b] * coeff)
								{
									halftimeFlag[b] = true;
								}
							}

							chromo.weight += reward;

							change[a][b] = 1;
							inContact[t][model->sigma + a] = 1;
							inContact[t][model->sigma + b] = 1;
							chromo.interNum++;
						}
					}
					else
					{
						change[a][b] = 0;
					}
				}
			}
		}
	}

	if (logFlag)
	{
		chromo.rewardsPerAgent.resize(4 * model->sigma); // Defaults to zero initial value
		for (int i = 0; i < model->sigma; i++)
		{
			// defining start and end of the reward distribution
			int start = startVal * stepsPA[i]; //
			int end = endVal * stepsPA[i];	   //
			// inst reward objects for each agent
			_vRew[i] = shared_ptr<REWARD>(new REWARD(stepsPA[i], start, end, model->avgWeight));
		}
		for (int i = 0; i < model->sigma; i++)
		{
			chromo.rewardsPerAgent[i].resize(maxSteps);
			for (int j = 0; j < maxSteps; j++)
			{
				if (j <= stepsPA[i])
				{
					chromo.rewardsPerAgent[i][j] = _vRew[i]->interactionReward2(j);
					if (inContact[j][model->sigma + i] == 1) // if there are agents interacting 1
					{

						if (inContact[j][i] > 0) //
						{
							chromo.rewardsPerAgent[i][j] = _vRew[i]->interactionReward(j);
						}
					}
				}
			}
		}

		for (int i = model->sigma; i < 4 * model->sigma; i++)
		{
			chromo.rewardsPerAgent[i].resize(maxSteps);
			for (int j = 0; j < maxSteps; j++)
			{
				chromo.rewardsPerAgent[i][j] = inContact[j][i - model->sigma];
			}
		}
	}

	return maxSteps;
}
#pragma endregion

#pragma region == GA Core Methods ==

void GA::initializePopulation()
{
	// Reset generation counter
	CURRENT_GENERATION = 0;

	// Initialize the population (memory allocation)
	population = vector<CHROMOSOME>(POPULATION_SIZE);

	// creation of initial plans - candidate solutions
	for (int i = 0; i < POPULATION_SIZE; i++)
	{
		// cout << i << endl;

		initializeChromosomeSRST(population[i]); // for different problems like MRMT you will have to pick here which inichrom to use. If no virtual tasks - ST, if no parallel tasks - SR;
												 //	initializeAnders(population[i]);
	}
}

void GA::createSelectionPool()
{

	// Usual methods for SELECTION in a Genetic Algorithm are:
	//
	// - Fitness Proportionate Selection
	// - Roulette Wheel Selection
	// - Stochastic Universal Sampling
	// - Rank Selection
	// - Elitist Selection
	// - Tournament Selection
	// - Truncation Selection
	//
	// This implementation combines some of the ideas.
	//
	// When a new generation will be evolved, not every individual in the population
	// will get a chance to reproduce. A given percentage of the best individualss will
	// be placed in a pool from which parents will be randomly selected
	// for reproduction. An individual with a better fitness value will be inserted
	// more often into the selection pool than an individual with a low fitness value.
	// An individual which is represented more often in the pool has a higher
	// probability to be selected as a parent.

	// The individuals will not be copied to the ranking list or the selection pool.
	// Both are vectors of indices pointing to the individuals in the population.

	vector<pair<int, double>> ranking(POPULATION_SIZE);

	// Rank the population based on the cost
	rankPop(ranking);

	// Take some percentage of the best DNAs from the ranking as candidates
	// for the selection pool. Better ranked candidates have a higher
	// probability to get inserted into the pool.

	_selectionPool.clear();
	_bestIndividualsIndex = ranking[0].first;
	_bestCost = ranking[0].second;

	const int candidateCount = static_cast<int>(ranking.size() * pXOVER);

	const double candidateQuotient = 1.0 / candidateCount;

	double insertionProbability = 0.0;

	for (int rank = 0; rank < candidateCount; rank++)
	{

		/// Elitism
		if (rank < ELITE)
		{
			_elitePopulation[rank] = population[ranking[rank].first];
		}

		insertionProbability = 1.0 - rank * candidateQuotient;

		// Optional: Strong (exponential) falloff
		// insertionProbability = std::pow(insertionProbability, 2.0);

		if (getRandomTrueWithProbability(insertionProbability)) // Bernoulli probability distribution
			_selectionPool.emplace_back(ranking[rank].first);
	}

	_halfSelectionPool = _selectionPool.size() / 2;
}

void GA::evaluate(CHROMOSOME &chromo)
{
	// TO-DO:
	// evaluate cost - checked
	// check for MR constraint violation
	// check for Precedence constraint violation
	// can we have color violation? - checked
	// how are we gonna solve parallel tasks? how to check for parallel task validity?

	/********************************************/
	vector<vector<pair<double, double>>> timeSteps(model->sigma);
	chromo.indCost.clear();
	chromo.indCost.resize(model->sigma);
	chromo.penal.pc = 0;
	chromo.penal.color = 0;
	chromo.cost = 0;

	for (int i = 0; i < model->sigma; i++)
	{
		// printOrderedPlan(population[_bestIndividualsIndex].individual[i], i);
		int _tmp = i;
		vector<int> _tmpPC(model->Vtilde); // mozda moze i van ove petlje?

		while (_tmp < chromo.individual[i].size())
		{
			/* Calculate the cost of the chromosome */
			// chromo.indCost[i] += model->getEdgeCost(_tmp, chromo.individual[i][_tmp], i) + model->getTaskDuration(_tmp);
			chromo.indCost[i] += model->getEdgeCost(_tmp, chromo.individual[i][_tmp], i) + model->getTaskDuration(chromo.individual[i][_tmp]);
			// cout << "Edge: (" << _tmp << "," << chromo.individual[i][_tmp] << ")" << " Edge Cost: " << model->getEdgeCost(_tmp, chromo.individual[i][_tmp], i) << " Task Dur: " << model->getTaskDuration(_tmp) << endl;
			/* Check color constraints */
			if ((chromo.individual[i][_tmp] > model->sigma) && (chromo.individual[i][_tmp] < chromo.individual[i].size()))
			{
				chromo.penal.color += checkColor(i, chromo.individual[i][_tmp]);
			}

			if (chromo.penal.color > 0)
			{
				cout << "Color Penal: " << chromo.individual[i][_tmp] << " Agent: " << i << endl;
			}

			/* Check Precedence Constraints */
			if (_tmpPC[chromo.individual[i][_tmp]] != 0)
			{
				_tmpPC[chromo.individual[i][_tmp]] = 0;
			}

			if (!model->_PV[chromo.individual[i][_tmp]].empty() && chromo.individual[i][_tmp] < chromo.individual[i].size())
			{
				for (int j = 0; j < model->_PV[chromo.individual[i][_tmp]].size(); j++)
				{
					if (model->_PV[chromo.individual[i][_tmp]][j] < (model->V + model->sigma))
					{
						_tmpPC[model->_PV[chromo.individual[i][_tmp]][j]] = 1;
					}
				}
			}

			// Update tmp
			_tmp = chromo.individual[i][_tmp];
		}
		chromo.penal.pc += accumulate(_tmpPC.begin(), _tmpPC.end(), 0);
	}

	///* addition for anders' paper */
	// int aSum = 0;
	// for (int a = 0; a < chromo.indTask.size(); a++) {
	//	if (chromo.indTask[a].empty()) {
	//		aSum++;
	//	}
	// }

	// createGantt(chromo, model);

	chromo.maxSteps = getTimeSteps(chromo, model, timeSteps, false);

	// int andersPenalty = model->sigma - (aSum + model->getNumberOfAgentsToUse());
	int interactionPen = accumulate(chromo.interPenalties.begin(), chromo.interPenalties.end(), 0);
	int _totalPenalties = chromo.penal.pc + chromo.penal.color + interactionPen; // +abs(andersPenalty);
	/* Sum over All */
	// chromo.cost = accumulate(chromo.indCost.begin(), chromo.indCost.end(), 0);
	/* MinMax + Sum over all*/
	// float target = 8;
	// float weight = abs(target - (inRangeCount * (-1)));

	if (_objLimit > 0)
	{

		// double cost = *max_element(chromo.indCost.begin(), chromo.indCost.end());
		double cost = chromo.maxSteps;
		double diff = cost - _objLimit;
		if (diff > 0)
		{
			_totalPenalties += diff;
		}
		chromo.cost = chromo.weightedInteractions * diff;
	}
	else if (_objLimit == -1)
	{
		// chromo.cost = *max_element(chromo.indCost.begin(), chromo.indCost.end()) + 0.1 * accumulate(chromo.indCost.begin(), chromo.indCost.end(), 0) - chromo.weight;
		//chromo.cost = chromo.maxSteps + 0.1 * accumulate(chromo.indCost.begin(), chromo.indCost.end(), 0) - chromo.weight;
		chromo.cost = chromo.maxSteps - chromo.weight;
	} else {
		//chromo.cost = *max_element(chromo.indCost.begin(), chromo.indCost.end()) + 0.1 * accumulate(chromo.indCost.begin(), chromo.indCost.end(), 0);
		chromo.cost = chromo.maxSteps;
	}

	// chromo.cost = *max_element(chromo.indCost.begin(), chromo.indCost.end()) + 0.1 * accumulate(chromo.indCost.begin(), chromo.indCost.end(), 0) - chromo.weight + 10*chromo.distWeight;

	/* MinMax */
	// chromo.cost = *max_element(chromo.indCost.begin(), chromo.indCost.end());

	if (_totalPenalties == 0)
	{
		chromo.feasibility = true;
		return;
	}
	else
	{
		chromo.feasibility = false;
		/* Sum all costs and penalties */
		chromo.cost += _totalPenalties * PENALTY;
		return;
	}
}

void GA::evaluateXD(CHROMOSOME &chromo)
{
	chromo.indCost.clear();
	chromo.indCost.resize(model->sigma);
	chromo.penal.pc = 0;
	chromo.penal.color = 0;
	chromo.cost = 0;
	chromo.timeline.resize(model->V);
	int try3 = 0;
	for (int i = 0; i < model->sigma; i++)
	{
		// printOrderedPlan(population[_bestIndividualsIndex].individual[i], i);
		int _tmp = i;
		int cpen = 0;
		int cost = 0;
		int memTask;
		int memCost = 0;
		try3 = 0;
		while (_tmp < chromo.individual[i].size())
		{
			/* Calculate the cost of the chromosome */

			if (model->H[_tmp] == 0 && model->H[chromo.individual[i][_tmp]] == 1)
			{
				cost = model->getTaskDuration(_tmp);
				memCost = 0;
				memTask = _tmp;
			}
			else if (model->H[_tmp] == 1 && model->H[chromo.individual[i][_tmp]] == 0)
			{
				memCost += model->getTaskDuration(_tmp);
				if (memCost > model->getEdgeCost(memTask, chromo.individual[i][_tmp], i))
				{
					cost = memCost;
				}
				else
				{
					cost = model->getEdgeCost(memTask, chromo.individual[i][_tmp], i);
				}
				memCost = 0;
			}
			else if (model->H[_tmp] == 0 && model->H[chromo.individual[i][_tmp]] == 0)
			{
				cost = model->getEdgeCost(_tmp, chromo.individual[i][_tmp], i) + model->getTaskDuration(_tmp);
			}
			else if (model->H[_tmp] == 1 && model->H[chromo.individual[i][_tmp]] == 1)
			{
				cost = 0;
				memCost += model->getTaskDuration(_tmp);
			}

			// chromo.indCost[i] += model->getEdgeCost(_tmp, chromo.individual[i][_tmp], i) + model->getTaskDuration(_tmp);
			try3 += model->getEdgeCost(_tmp, chromo.individual[i][_tmp], i) + model->getTaskDuration(_tmp);
			chromo.indCost[i] += cost;

			/* Check color constraints */
			if ((chromo.individual[i][_tmp] > model->sigma) && (chromo.individual[i][_tmp] < chromo.individual[i].size()))
			{
				cpen = checkColor(i, chromo.individual[i][_tmp]);
				chromo.penal.color += cpen;
			}

			if (cpen > 0)
			{
				cout << "Color Penal: " << chromo.individual[i][_tmp] << " Agent: " << i << endl;
			}

			/* Check for XD dependencies */
			if (!model->_PV[chromo.individual[i][_tmp]].empty() && chromo.individual[i][_tmp] < chromo.individual[i].size())
			{
				for (int j = 0; j < model->_PV[chromo.individual[i][_tmp]].size(); j++)
				{
					if (model->_PV[chromo.individual[i][_tmp]][j] < (model->V + model->sigma))
					{
						if (chromo.timeline[model->_PV[chromo.individual[i][_tmp]][j] - model->sigma].startTime < chromo.timeline[chromo.individual[i][_tmp] - model->sigma].endTime)
						{
							chromo.penal.pc++;
							// cout << "PC penal: " << chromo.individual[i][_tmp] << " <- " << model->_PV[chromo.individual[i][_tmp]][j] << endl;
						}
					}
				}
			}

			// Update tmp
			_tmp = chromo.individual[i][_tmp];
		}
	}

	int _totalPenalties = chromo.penal.pc + chromo.penal.color;
	/* Sum over All */
	// chromo.cost = accumulate(chromo.indCost.begin(), chromo.indCost.end(), 0);
	/* MinMax + Sum over all*/
	// chromo.cost = *max_element(chromo.indCost.begin(), chromo.indCost.end()) + 0.01 * accumulate(chromo.indCost.begin(), chromo.indCost.end(), 0);
	/* MinMax */
	chromo.cost = *max_element(chromo.indCost.begin(), chromo.indCost.end());

	if (_totalPenalties == 0)
	{
		chromo.feasibility = true;
		return;
	}
	else
	{
		chromo.feasibility = false;
		/* Sum all costs and penalties */
		chromo.cost += _totalPenalties * PENALTY;
		return;
	}
}

void GA::pcFix(CHROMOSOME &chromo)
{
	vector<int> tempPC;

	for (int i = 0; i < model->sigma; i++)
	{
		int _tmp = i, _tmpOld = _tmp, _tmpOld2 = _tmp, _taskCounter = 1;
		bool flag = false;

		while (_tmp < chromo.individual[i].size() && chromo.individual[i][_tmp] < chromo.individual[i].size())
		{

			/* Kill the while loop when the destination deopt is reached */
			// if ((model->Vtilde - model->delta) == chromo.individual[i][_tmp]) { break; }

			/*if (!model->_PV[chromo.individual[i][_tmp]].empty()) {*/

			/* It works only in case of 1 prec constraint?!? */

			if (!model->_PV[chromo.individual[i][_tmp]].empty() && model->_PV[chromo.individual[i][_tmp]][0] < (model->V + model->sigma))
			{
				// printOrderedPlan(chromo.individual[i], i);
				//  The case where tasks are on the same agent but the ordering is wrong
				auto _val = findTask(chromo.individual[i], i, chromo.individual[i][_tmp], model->sigma + model->V - 1, model->_PV[chromo.individual[i][_tmp]][0], _tmpOld2);

				if (_val >= 0)
				{
					// Fix within agent plan
					int _start, _end, _prevT, _idxTask2Move = 0;
					// printOrderedPlan(chromo.individual[i], i);
					if (getRandomTrueWithProbability(0.5))
					{
						_idxTask2Move = _val;
						_start = chromo.individual[i][_tmp];
						_end = model->sigma + model->V - 1;
						_prevT = _tmpOld2;
					}
					else
					{
						_idxTask2Move = _tmp;
						_start = i;
						_end = _val;
						_prevT = _tmp;
					}

					// cout << " PC FIX 1 " << endl;
					int _insertIdx = insertionIdx(chromo.individual[i], _start, _end, _taskCounter);

					int _tmpVal = chromo.individual[i][_insertIdx];
					int _tmpVal2 = chromo.individual[i][chromo.individual[i][_idxTask2Move]];

					chromo.individual[i][_insertIdx] = chromo.individual[i][_idxTask2Move];
					chromo.individual[i][chromo.individual[i][_insertIdx]] = _tmpVal;
					chromo.individual[i][_idxTask2Move] = _tmpVal2;

					// printOrderedPlan(chromo.individual[i], i);
					flag = true;
				}
				else if (_val == -1)
				{

					// bool _flag = true;

					// fix cross agent plan
					int _start, _end, _prevT, _idxTask2Move = 0, _agent1, _agent2; // a task is removed from agent 1 and inserted in agent's 2 plan

					// printOrderedPlan(chromo.individual[i], i);
					auto _tmpPair = crossAgentTaskSearch(chromo, i, model->_PV[chromo.individual[i][_tmp]][0], true);

					auto it = find(model->_vectorTasksPerAgents[chromo.individual[_tmpPair.first][_tmpPair.second] - model->sigma].begin(), model->_vectorTasksPerAgents[chromo.individual[_tmpPair.first][_tmpPair.second] - model->sigma].end(), i);
					auto it2 = find(model->_vectorTasksPerAgents[chromo.individual[i][_tmp] - model->sigma].begin(), model->_vectorTasksPerAgents[chromo.individual[i][_tmp] - model->sigma].end(), _tmpPair.first);

					if (it == model->_vectorTasksPerAgents[chromo.individual[_tmpPair.first][_tmpPair.second] - model->sigma].end() && it2 == model->_vectorTasksPerAgents[chromo.individual[i][_tmp] - model->sigma].end())
					{
						// ako ni jedan od ova dva ne podrzava oba taska //REMEMBER TO UPDATE numTasks and indTasks;
						vector<int> _intersection;
						std::set_intersection(model->_vectorTasksPerAgents[chromo.individual[_tmpPair.first][_tmpPair.second] - model->sigma].begin(), model->_vectorTasksPerAgents[chromo.individual[_tmpPair.first][_tmpPair.second] - model->sigma].end(),
											  model->_vectorTasksPerAgents[chromo.individual[i][_tmp] - model->sigma].begin(), model->_vectorTasksPerAgents[chromo.individual[i][_tmp] - model->sigma].end(),
											  back_inserter(_intersection));

						if (_intersection.empty())
						{
							cout << "Precedence Constraint problem, bad model: " << _tmp << " < " << model->_PV[chromo.individual[i][_tmp]][0] << " not feasible" << endl;
						}

						_agent1 = getRandomIntegerInRange(0, (int)_intersection.size() - 1);
						_start = _intersection[_agent1];
						_end = model->sigma + model->V - 1;

						_idxTask2Move = _tmp;
						_agent2 = i; //_flag = false;
						performPCfix(chromo, _intersection[_agent1], _agent2, _start, _end, _idxTask2Move, _taskCounter);

						_idxTask2Move = _tmpPair.second;
						_agent2 = _tmpPair.first;
						performPCfix(chromo, _intersection[_agent1], _agent2, _start, _end, _idxTask2Move, _taskCounter);

						flag = true;
					}
					else
					{
						if (it == model->_vectorTasksPerAgents[chromo.individual[_tmpPair.first][_tmpPair.second] - model->sigma].end() && it2 != model->_vectorTasksPerAgents[chromo.individual[i][_tmp] - model->sigma].end())
						{
							// if only the second agent can perform both tasks //REMEMBER TO UPDATE numTasks and indTasks;
							_idxTask2Move = _tmp;
							_end = chromo.individual[_tmpPair.first][_tmpPair.second];
							_agent1 = _tmpPair.first;
							_agent2 = i;
							_start = _agent1; //_flag = false;
						}
						else if (it != model->_vectorTasksPerAgents[chromo.individual[_tmpPair.first][_tmpPair.second] - model->sigma].end() && it2 == model->_vectorTasksPerAgents[chromo.individual[i][_tmp] - model->sigma].end())
						{
							// if only the first agent can perform both tasks //REMEMBER TO UPDATE numTasks and indTasks;
							_idxTask2Move = _tmpPair.second;
							_end = model->sigma + model->V - 1;
							_agent1 = i;
							_agent2 = _tmpPair.first;
							_start = chromo.individual[i][_tmp]; //_start = _agent1;// _flag = true;
						}
						else
						{
							// if both agents can perform both tasks //REMEMBER TO UPDATE numTasks and indTasks;
							if (getRandomTrueWithProbability(0.5))
							{
								_idxTask2Move = _tmp;
								_end = chromo.individual[_tmpPair.first][_tmpPair.second];
								_agent1 = _tmpPair.first;
								_agent2 = i;
								_start = _agent1; //_flag = false;
							}
							else
							{
								_idxTask2Move = _tmpPair.second;
								_end = INDIVIDUAL_SIZE;
								_agent1 = i;
								_agent2 = _tmpPair.first;
								_start = chromo.individual[i][_tmp]; //_start = _agent1; //_flag = true;
							}
						}

						performPCfix(chromo, _agent1, _agent2, _start, _end, _idxTask2Move, _taskCounter);

						flag = true;
					}
				}
			}

			// Update tmp
			_tmpOld = _tmp;

			if (flag)
			{
				_tmp = i;
				flag = false;
			}
			else
			{
				_tmp = chromo.individual[i][_tmp];
			}
			_taskCounter++;
		}
	}
}

void GA::xover()
{
	// Select parents from the pool
	_splitPoints.resize(model->sigma);
	_childrenPopulation.clear();
	_childrenPopulation.resize(_halfSelectionPool);

	vector<vector<unordered_set<int>>> _aMatrix(_halfSelectionPool, vector<unordered_set<int>>(INDIVIDUAL_SIZE));

#ifdef SERIAL
	///* Serial Execution*/
	// Generate Adjacency Matrix
	for (int i = 0; i < _halfSelectionPool; i++)
	{
		genAdjMatrix(_aMatrix[i]);
	}

	for (int i = 0; i < _halfSelectionPool; i++)
	{
		// cout << i << endl;
		selectParentsFromPool();
		selectSplitPoints();
		ERX(_childrenPopulation[i], _aMatrix[i]);
	}

#else
	/* Parallel Execution */
	/// Generate Adjacency Matrix
	tbb::parallel_for(tbb::blocked_range<int>(0, _halfSelectionPool), [&](const tbb::blocked_range<int> &x)
					  {
		for (int i = x.begin(); i < x.end(); i++) {
			genAdjMatrix(_aMatrix[i]);
		} });

	tbb::parallel_for(tbb::blocked_range<int>(0, _halfSelectionPool), [&](const tbb::blocked_range<int> &x)
					  {
		for (int i = x.begin(); i < x.end(); i++) {
			selectParentsFromPool();
			selectSplitPoints();
			ERX(_childrenPopulation[i], _aMatrix[i]);
		} });

#endif
}

void GA::mutate(CHROMOSOME &chromo)
{
	// for (int i = 0; i < POPULATION_SIZE; i++) {

	//	cout << "AAAAAA MUTATION: " << i << endl;

	if (getRandomTrueWithProbability(pMUTATION))
	{
		withinAgentSwap(chromo);
	}
	if (getRandomTrueWithProbability(pMUTATION))
	{
		withinAgentJump(chromo);
	}

	if (model->sigma > 1)
	{
		if (getRandomTrueWithProbability(pMUTATION))
		{
			crossAgentJump(chromo);
		}
		if (getRandomTrueWithProbability(pMUTATION))
		{
			crossAgentSwap(chromo);
		}
	}

	//}
}

void GA::elitismSelection()
{
}

void GA::localRefinement(CHROMOSOME &chromo)
{
	GreedySearch(chromo);
	twoOptLocalSearch(chromo);
}

#pragma endregion

#pragma region == Helper Methods ==

void GA::initializeChromosomeSRST(CHROMOSOME &chromo)
{
	// allocate mem
	chromo.numTasks.resize(model->sigma);
	// allocate mem
	chromo.indTask.resize(model->sigma);
	// Initialize the Individual (allocate memory)
	chromo.individual = vector<vector<int>>(model->sigma, vector<int>(model->V + model->sigma));
	// Get the list of tasks per Agents
	vector<vector<int>> _sigma = model->getListTasksPerAgents();
	// Initialize temp variable containing the last assigned task in the assignment process for each agent
	vector<int> lastAssigned(model->sigma);
	std::iota(lastAssigned.begin(), lastAssigned.end(), 0);

	/// Add the first task in the plan
	for (int i = 0; i < model->sigma; i++)
	{

		// Randomly pick a task
		int rndT = getRandomIntegerInRange(0, static_cast<int>(_sigma[i].size() - 1));

		if (rndT == -1)
		{
			continue;
		}
		// Add the task to the chromosome
		chromo.individual[i][i] = _sigma[i][rndT] + model->sigma;
		// Update the last assigned tasks
		lastAssigned[i] = _sigma[i][rndT] + model->sigma;
		// Remove assigned task from other lists
		updateSigma(_sigma, _sigma[i][rndT]);
		// update number of tasks assigned to this specific agent
		chromo.numTasks[i] += 1;

		chromo.indTask[i].emplace_back(lastAssigned[i]);
	}

	int count = 0;
	for (int i = 0; i < model->sigma; i++)
	{
		if (_sigma[i].empty())
		{
			count++;
		}
	}
	/// Allocate the rest of the Tasks
	if (count < model->sigma)
	{
		for (int i = 0; i < (model->V - model->sigma); i++)
		{

			// cout << i << "---------" << endl;
			int rndA, rndT;

			// Randomly pick an agent
			do
			{
				rndA = getRandomIntegerInRange(0, static_cast<int>(model->sigma - 1)); // FIX THIS WITHOUT WHILE LOOP SOMEHOW!!!
			} while (_sigma[rndA].size() == 0);

			// Randomly pick a task
			rndT = getRandomIntegerInRange(0, static_cast<int>(_sigma[rndA].size() - 1));

			// Add the task to the chromosome
			// cout << rndA  << " " << rndT << " " << model->sigma << " " << lastAssigned[rndA] << endl;
			if (rndT >= 0)
			{
				chromo.individual[rndA][lastAssigned[rndA]] = _sigma[rndA][rndT] + model->sigma;

				// Update the last assigned tasks
				lastAssigned[rndA] = _sigma[rndA][rndT] + model->sigma;
				// Remove assigned task from other lists
				updateSigma(_sigma, _sigma[rndA][rndT]);
				chromo.numTasks[rndA] += 1;

				// update indTask
				chromo.indTask[rndA].emplace_back(lastAssigned[rndA]);
			}
		}
	}
	/// Add final destination :P i.e., destination depots
	for (int i = 0; i < model->sigma; i++)
	{

		// Randomly pick a destination depot
		// int rndD = (model->V + model->sigma) + getRandomIntegerInRange(0, static_cast<int>(model->delta - 1));
		// chromo.individual[i][lastAssigned[i]] = rndD;

		// Pick closet destination depot
		pickDestinationDepot(chromo.individual[i], lastAssigned[i], i);
		// Add the task to the chromosome

		chromo.numTasks[i] += 1;
	}
}

void GA::initializeAnders(CHROMOSOME &chromo)
{
	// allocate mem
	chromo.numTasks.resize(model->sigma);
	// allocate mem
	chromo.indTask.resize(model->sigma);
	// Initialize the Individual (allocate memory)
	chromo.individual = vector<vector<int>>(model->sigma, vector<int>(model->V + model->sigma));
	// Get the list of tasks per Agents
	vector<vector<int>> _sigma = model->getListTasksPerAgents();
	// Initialize temp variable containing the last assigned task in the assignment process for each agent
	vector<int> lastAssigned(model->sigma);
	std::iota(lastAssigned.begin(), lastAssigned.end(), 0);

	vector<int> agentsToPick(model->sigma);
	std::iota(agentsToPick.begin(), agentsToPick.end(), 0);

	vector<int> agList;
	for (int nn = 0; nn < model->getNumberOfAgentsToUse(); nn++)
	{
		int a = getRandomIntegerInRange(0, static_cast<int>(agentsToPick.size() - 1));
		agList.push_back(agentsToPick[a]);
		agentsToPick.erase(agentsToPick.begin() + a);
	}

	/// Add the first task in the plan
	for (int i = 0; i < model->getNumberOfAgentsToUse(); i++)
	{

		// Randomly pick a task
		int rndT = getRandomIntegerInRange(0, static_cast<int>(_sigma[agList[i]].size() - 1));

		if (rndT == -1)
		{
			continue;
		}
		// Add the task to the chromosome
		chromo.individual[agList[i]][agList[i]] = _sigma[agList[i]][rndT] + model->sigma;
		// Update the last assigned tasks
		lastAssigned[agList[i]] = _sigma[agList[i]][rndT] + model->sigma;
		// Remove assigned task from other lists
		updateSigma(_sigma, _sigma[agList[i]][rndT]);
		// update number of tasks assigned to this specific agent
		chromo.numTasks[agList[i]] += 1;

		chromo.indTask[agList[i]].emplace_back(lastAssigned[agList[i]]);
	}

	int count = 0;
	for (int i = 0; i < model->getNumberOfAgentsToUse(); i++)
	{
		if (_sigma[agList[i]].empty())
		{
			count++;
		}
	}
	/// Allocate the rest of the Tasks
	if (count < model->getNumberOfAgentsToUse())
	{
		for (int i = 0; i < (model->V - model->getNumberOfAgentsToUse()); i++)
		{

			// cout << i << "---------" << endl;
			int rndA, rndT;

			// Randomly pick an agent
			do
			{
				rndA = getRandomIntegerInRange(0, static_cast<int>(model->getNumberOfAgentsToUse() - 1)); // FIX THIS WITHOUT WHILE LOOP SOMEHOW!!!
			} while (_sigma[rndA].size() == 0);
			rndA = agList[rndA];

			// Randomly pick a task
			rndT = getRandomIntegerInRange(0, static_cast<int>(_sigma[rndA].size() - 1));

			// Add the task to the chromosome
			// cout << rndA  << " " << rndT << " " << model->sigma << " " << lastAssigned[rndA] << endl;
			if (rndT >= 0)
			{
				chromo.individual[rndA][lastAssigned[rndA]] = _sigma[rndA][rndT] + model->sigma;

				// Update the last assigned tasks
				lastAssigned[rndA] = _sigma[rndA][rndT] + model->sigma;
				// Remove assigned task from other lists
				updateSigma(_sigma, _sigma[rndA][rndT]);
				chromo.numTasks[rndA] += 1;

				// update indTask
				chromo.indTask[rndA].emplace_back(lastAssigned[rndA]);
			}
		}
	}
	/// Add final destination :P i.e., destination depots
	for (int i = 0; i < model->sigma; i++)
	{

		// Randomly pick a destination depot
		// int rndD = (model->V + model->sigma) + getRandomIntegerInRange(0, static_cast<int>(model->delta - 1));
		// chromo.individual[i][lastAssigned[i]] = rndD;

		// Pick closet destination depot
		pickDestinationDepot(chromo.individual[i], lastAssigned[i], i);
		// Add the task to the chromosome

		chromo.numTasks[i] += 1;
	}
}

pair<double, double> GA::calculateAgentsPosition(MODEL *model, int t1, int t2, int a, int timeStep, bool &nextTask)
{
	double dx = 0, dy = 0, x = 0, y = 0;
	double t1x, t1y, t2x, t2y;

	if (t1 < model->sigma)
	{
		t1x = model->_model.A[t1].x;
		t1y = model->_model.A[t1].y;
		t2x = model->_model.T[t2 - model->sigma].x;
		t2y = model->_model.T[t2 - model->sigma].y;
	}
	else if (t2 >= model->sigma + model->V)
	{
		t1x = model->_model.T[t1 - model->sigma].x;
		t1y = model->_model.T[t1 - model->sigma].y;
		t2x = model->_model.dest[t2 - (model->sigma + model->V)].x;
		t2y = model->_model.dest[t2 - (model->sigma + model->V)].y;
	}
	else
	{
		t1x = model->_model.T[t1 - model->sigma].x;
		t1y = model->_model.T[t1 - model->sigma].y;
		t2x = model->_model.T[t2 - model->sigma].x;
		t2y = model->_model.T[t2 - model->sigma].y;
	}

	dx = t2x - t1x;
	dy = t2y - t1y;

	double distance = model->getEdgeCost(t1, t2, a);
	if (distance == 0.0)
	{
		// The vehicle has reached the destination
		nextTask = true;
		return {t2x, t2y};
	}

	double ratio = model->_model.A[a].v * timeStep / distance;

	if (ratio >= 1.0)
	{
		// The vehicle has reached or surpassed the destination
		nextTask = true;
		return {t2x, t2y};
	}

	if (t1 < model->sigma)
	{
		x = model->_model.A[t1].x + ratio * dx;
		y = model->_model.A[t1].y + ratio * dy;
	}
	else
	{
		x = model->_model.T[t1 - model->sigma].x + ratio * dx;
		y = model->_model.T[t1 - model->sigma].y + ratio * dy;
	}
	return {x, y};
}

void GA::updateSigma(vector<vector<int>> &_sigma, int rndT)
{

	for (int rndA = 0; rndA < model->sigma; rndA++)
	{

		std::vector<int>::iterator itr = std::find(_sigma[rndA].begin(), _sigma[rndA].end(), rndT);

		if (_sigma[rndA].size() > 1 && itr != _sigma[rndA].cend())
		{

			auto index = std::distance(_sigma[rndA].begin(), itr);
			std::iter_swap(_sigma[rndA].begin() + index, _sigma[rndA].end() - 1);
			_sigma[rndA].pop_back();
		}
		else if (_sigma[rndA].size() <= 1 && itr != _sigma[rndA].cend())
		{
			_sigma[rndA].clear();
		}
	}
}

int GA::pickNextTask(int task, int agent, vector<unordered_set<int>> &_aMatrix)
{
	int _tmp = 5; // I just need some big number here. 5 is good because 4 is the maximum it can happen.
	vector<int> _v;
	_v.reserve(4);

	if (task >= INDIVIDUAL_SIZE)
	{
		return -1;
	}

	// Erase added task from the whole adjacency Matrix
	for (int i = 0; i < _aMatrix.size(); i++)
	{
		auto it = find(_aMatrix[i].begin(), _aMatrix[i].end(), task);
		if (it != _aMatrix[i].end())
		{
			_aMatrix[i].erase(it);
			if (_aMatrix[i].empty())
			{
				_aMatrix[i].insert(-1);
			}
		}
	}

	// Pick a new task to be added
	for (auto it = _aMatrix[task].begin(); it != _aMatrix[task].end(); it++)
	{
		if (*it > -1 && _aMatrix[*it].size() >= 1 && *_aMatrix[*it].begin() != -1)
		{
			if (_aMatrix[*it].size() <= _tmp)
			{
				if (_aMatrix[*it].size() < _tmp)
				{
					_v.clear();
				}

				if (!checkColor(agent, *it))
				{
					_v.emplace_back(*it);
					_tmp = _aMatrix[*it].size();
				}
			}
		}
		else
		{
			return -1;
		}
	}

	if (_aMatrix[task].empty())
	{
		return -1;
	}
	_aMatrix[task].clear();

	if (_v.empty())
	{
		return -1;
	}
	else
	{
		return _v[getRandomIntegerInRange(static_cast<int>(0), static_cast<int>(_v.size() - 1))];
	}
}

int GA::getRandomFromSet(unordered_set<int> &_set)
{
	// if (!_set.empty()) {
	int _rnd = getRandomIntegerInRange(static_cast<int>(0), static_cast<int>(INDIVIDUAL_SIZE - 1));

	auto _ins = _set.insert(_rnd);

	if (_ins.second == false)
	{
		return _rnd;
	}

	auto _ptr = _ins.first++;

	if (_ins.first == _set.end())
	{
		_set.erase(_ptr);
		return *_set.begin();
	}
	else
	{
		_set.erase(_ptr);
		return *_ins.first;
	}
	/*}
	else { return -1; }*/
}

int GA::checkColor(int Agent, int Task)
{
	auto *tmp = model->getTasksPerAgents();

	if ((*tmp)[Agent][Task - model->sigma] == 1)
	{
		return 0;
	}

	return 1;
}

int GA::findTask(vector<int> &vector, int _ag, int _start, int _end, int _taskIdx, int &tmp)
{
	// No Prec violation -> return -2;
	// Prec violation, both tasks are assigned to the same agent return task index;
	// Prec violation, tasks are not assigned to the same agent return -1;

	tmp = _start;

	while (tmp < vector.size())
	{
		if (vector[tmp] == _taskIdx)
		{
			return -2;
		}
		tmp = vector[tmp];
	}

	int _iter = 0;
	tmp = _ag;

	while (tmp != vector[_start])
	{
		if (vector[tmp] == _taskIdx)
		{
			return tmp;
		}

		tmp = vector[tmp];
		_iter++;
	}

	return -1;
}

int GA::insertionIdx(vector<int> &vector, int _start, int _end, int _tasksLeft)
{

	int count = 0, start = _start;
	/*	cout << "start: " << _start << " end: " << _end << endl;

		auto it = find_if(vector.begin(), vector.begin() + model->sigma, [](int n) { return  n > 0; });

		printOrderedPlan(vector, *it);*/

	while (_start != _end && _start < INDIVIDUAL_SIZE)
	{
		// cout << _start << "  " << flush;
		_start = vector[_start];

		count++;
	}
	// cout << endl;

	if (count == 0)
	{
		return _start;
	}

	int rnd = getRandomIntegerInRange(0, count - 1);

	//_bias.emplace_back(make_pair(rnd, count-1));

	for (int i = 0; i < rnd; i++)
	{

		start = vector[start];
	}

	return start;

	// the function should not reach this point
	cout << "Something went horribly wrong in insertionIdx function!" << endl;
}

void GA::selectParentsFromPool()
{
	/*	  Randomly select two individuals from the selection pool.
		  Individuals with better fitness values have a higher chance to get selected.
		  It is tolerable that both parents might be the same. */

	const auto poolIdx1 =
		getRandomIntegerInRange<int>(0, _selectionPool.size() - 1);

	const auto poolIdx2 =
		getRandomIntegerInRange<int>(0, _selectionPool.size() - 1);

	_parent1 = _selectionPool[poolIdx1];
	_parent2 = _selectionPool[poolIdx2];
}

void GA::selectSplitPoints()
{
	int limit;

	for (int i = 0; i < model->sigma; i++)
	{

		if (population[_parent1].indTask[i].size() > population[_parent2].indTask[i].size())
		{
			limit = population[_parent2].indTask[i].size();
		}
		else
		{
			limit = population[_parent1].indTask[i].size();
		}

		_splitPoints[i] = getRandomIntegerInRange<int>(0, limit - 1);
	}
}

void GA::withinAgentSwap(CHROMOSOME &chromo)
{

	int rndA = getRandomIntegerInRange(0, static_cast<int>(model->sigma - 1)); // what if i pick agent without tasks?
	int rndT1 = getRandomIntegerInRange(0, static_cast<int>(chromo.indTask[rndA].size() - 1));
	if (rndT1 < 0)
	{
		return;
	}
	int rndT2 = getRandomIntegerInRangeExcluding(0, static_cast<int>(chromo.indTask[rndA].size() - 1), rndT1);

	if (rndT2 > -1)
	{

		auto it = find(chromo.individual[rndA].begin(), chromo.individual[rndA].end(), chromo.indTask[rndA][rndT1]);
		int index1 = it - chromo.individual[rndA].begin();
		int tmp1 = chromo.individual[rndA][chromo.indTask[rndA][rndT1]];

		it = find(chromo.individual[rndA].begin(), chromo.individual[rndA].end(), chromo.indTask[rndA][rndT2]);
		int index2 = it - chromo.individual[rndA].begin();

		if (chromo.indTask[rndA][rndT1] == index2)
		{
			chromo.individual[rndA][chromo.indTask[rndA][rndT1]] = chromo.individual[rndA][chromo.indTask[rndA][rndT2]];

			chromo.individual[rndA][index1] = chromo.indTask[rndA][rndT2];
			chromo.individual[rndA][chromo.indTask[rndA][rndT2]] = chromo.indTask[rndA][rndT1];
		}
		else if (chromo.indTask[rndA][rndT2] == index1)
		{

			chromo.individual[rndA][chromo.indTask[rndA][rndT1]] = chromo.indTask[rndA][rndT2];
			chromo.individual[rndA][index1] = tmp1;
			chromo.individual[rndA][index2] = chromo.indTask[rndA][rndT1];
		}
		else
		{
			chromo.individual[rndA][index2] = chromo.indTask[rndA][rndT1];
			chromo.individual[rndA][chromo.indTask[rndA][rndT1]] = chromo.individual[rndA][chromo.indTask[rndA][rndT2]];

			chromo.individual[rndA][index1] = chromo.indTask[rndA][rndT2];
			chromo.individual[rndA][chromo.indTask[rndA][rndT2]] = tmp1;
		}

		_mutationsDone++;
	}

	/// Sanity check

	// auto a = chromo.individual[rndA];
	// sort(a.begin(), a.end());
	// auto last = unique(a.begin(), a.end());
	// a.erase(last, a.end());

	// if (a.size() != chromo.indTask[rndA].size() + 2) {
	//	cout << "Govno1" << endl;
	// }

	// cout << "stop" << endl;
}

void GA::withinAgentJump(CHROMOSOME &chromo)
{
	int rndA = getRandomIntegerInRange(0, static_cast<int>(model->sigma - 1)); // what if i pick agent without tasks?
	int rndT = getRandomIntegerInRange(0, static_cast<int>(chromo.indTask[rndA].size() - 1));

	if (rndT < 0)
	{
		return;
	}

	int rndP = getRandomIntegerInRangeExcluding(0, static_cast<int>(chromo.indTask[rndA].size() - 1), rndT);

	if (rndP > -1)
	{
		auto it = find(chromo.individual[rndA].begin(), chromo.individual[rndA].end(), chromo.indTask[rndA][rndP]);
		int index2 = it - chromo.individual[rndA].begin();
		int tmp = chromo.individual[rndA][index2];
		int tmp3 = chromo.individual[rndA][tmp];

		it = find(chromo.individual[rndA].begin(), chromo.individual[rndA].end(), chromo.indTask[rndA][rndT]);
		int index1 = it - chromo.individual[rndA].begin();
		int tmp2 = chromo.individual[rndA][index1];

		if (index2 != chromo.indTask[rndA][rndT])
		{
			chromo.individual[rndA][index2] = chromo.indTask[rndA][rndT];
			int tmp4 = chromo.individual[rndA][chromo.indTask[rndA][rndT]];
			chromo.individual[rndA][chromo.indTask[rndA][rndT]] = tmp;
			chromo.individual[rndA][index1] = tmp4;
		}
		else if (index1 == chromo.indTask[rndA][rndP])
		{

			cout << "index1 = rndP" << endl;
			return;
		}
		else
		{
			// geneSwap(chromo.individual[rndA], chromo.indTask[rndA][rndT], chromo.indTask[rndA][rndP]);
			chromo.individual[rndA][index1] = tmp;
			chromo.individual[rndA][tmp] = chromo.indTask[rndA][rndT];
			chromo.individual[rndA][tmp2] = tmp3;

			_mutationsDone++;
		}

		/// Sanity check

		// auto a = chromo.individual[rndA];
		// sort(a.begin(), a.end());
		// auto last = unique(a.begin(), a.end());
		// a.erase(last, a.end());

		// if (a.size() != chromo.indTask[rndA].size() + 2) {
		//	cout << "Govno1" << endl;
		// }
	}
	// cout << "stop" << endl;
}

void GA::crossAgentSwap(CHROMOSOME &chromo)
{
	// auto *tmpTPA = model->getVectorTasksPerAgents();
	int index1, index2, tmp, tmp1, tmp2, tmpback;

	int rndA1 = getRandomIntegerInRange(0, static_cast<int>(model->sigma - 1)); // what if i pick agent without tasks?

	int rndT1 = getRandomIntegerInRange(0, static_cast<int>(chromo.indTask[rndA1].size() - 1));

	if (rndT1 == -1)
	{ // how should I solve this differently?
		return;
	}
	int rnd2A1 = find(model->_vectorTasksPerAgents[chromo.indTask[rndA1][rndT1] - model->sigma].begin(), model->_vectorTasksPerAgents[chromo.indTask[rndA1][rndT1] - model->sigma].end(), rndA1) - model->_vectorTasksPerAgents[chromo.indTask[rndA1][rndT1] - model->sigma].begin();

	int _rndA2 = getRandomIntegerInRangeExcluding(0, static_cast<int>(model->_vectorTasksPerAgents[chromo.indTask[rndA1][rndT1] - model->sigma].size() - 1), rnd2A1); // what if i pick agent without tasks?

	if (_rndA2 == -1)
	{ // how should I solve this differently?
		return;
	}

	int rndA2 = model->_vectorTasksPerAgents[chromo.indTask[rndA1][rndT1] - model->sigma][_rndA2];

	if (rndA2 == -1)
	{ // how should I solve this differently?
		return;
	}

	int rndT2 = getRandomIntegerInRange(0, static_cast<int>(chromo.indTask[rndA2].size() - 1));

	if (rndT2 == -1)
	{ // how should I solve this differently?
		return;
	}

	auto _it = find(model->_vectorTasksPerAgents[chromo.indTask[rndA2][rndT2] - model->sigma].begin(), model->_vectorTasksPerAgents[chromo.indTask[rndA2][rndT2] - model->sigma].end(), rndA1);

	if (_it == model->_vectorTasksPerAgents[chromo.indTask[rndA2][rndT2] - model->sigma].end())
	{
		return;
	}

	if (chromo.indTask[rndA2].size() > 0)
	{

		/*for (auto a : chromo.individual[rndA1]) { cout << a << " " << flush; }
		cout << endl;
		cout << endl;
		cout << chromo.indTask[rndA1][rndT] << endl;;*/
		auto it = find(chromo.individual[rndA1].begin(), chromo.individual[rndA1].end(), chromo.indTask[rndA1][rndT1]); // there is some bug here...
		if (it != chromo.individual[rndA1].end())
		{
			index1 = it - chromo.individual[rndA1].begin();
			tmp = chromo.individual[rndA1][chromo.indTask[rndA1][rndT1]];
			tmp1 = chromo.individual[rndA1][index1];
		}
		else
		{
			return;
		}

		it = find(chromo.individual[rndA2].begin(), chromo.individual[rndA2].end(), chromo.indTask[rndA2][rndT2]);
		if (it != chromo.individual[rndA2].end())
		{
			index2 = it - chromo.individual[rndA2].begin();
			tmp2 = chromo.individual[rndA2][index2];
		}
		else
		{
			return;
		}

		chromo.individual[rndA1][index1] = chromo.indTask[rndA2][rndT2];
		chromo.individual[rndA1][chromo.indTask[rndA2][rndT2]] = tmp;
		chromo.individual[rndA1][chromo.indTask[rndA1][rndT1]] = 0;

		chromo.individual[rndA2][index2] = tmp1;
		chromo.individual[rndA2][chromo.indTask[rndA1][rndT1]] = chromo.individual[rndA2][tmp2];
		chromo.individual[rndA2][chromo.indTask[rndA2][rndT2]] = 0;

		if (chromo.indTask[rndA1].size() > 0)
		{
			std::iter_swap(chromo.indTask[rndA1].begin() + rndT1, chromo.indTask[rndA1].end() - 1);
			tmpback = chromo.indTask[rndA1].back();
			chromo.indTask[rndA1].pop_back();
			chromo.indTask[rndA1].emplace_back(chromo.indTask[rndA2][rndT2]);

			std::iter_swap(chromo.indTask[rndA2].begin() + rndT2, chromo.indTask[rndA2].end() - 1);
			chromo.indTask[rndA2].pop_back();
			chromo.indTask[rndA2].emplace_back(tmpback);
		}
		else if (chromo.indTask[rndA1].size() <= 1)
		{
			chromo.indTask[rndA1].clear();
		}
		// chromo.indTask[rndA2].emplace_back(chromo.individual[rndA2][index2]);

		_mutationsDone++;
	}

	//////check consistency
	// auto a = chromo.individual[rndA1]; auto b = chromo.individual[rndA2];
	// sort(a.begin(), a.end());
	// auto last = unique(a.begin(), a.end());
	// a.erase(last, a.end());

	// if (a.size() != chromo.indTask[rndA1].size() + 2) {
	//	sort(chromo.indTask[rndA1].begin(), chromo.indTask[rndA1].end());
	//	cout << "Govno1" << endl;
	// }

	// sort(b.begin(), b.end());
	// last = unique(b.begin(), b.end());
	// b.erase(last, b.end());

	// if (b.size() != chromo.indTask[rndA2].size() + 2) {
	//	cout << "Govno2" << endl;
	// }
}

void GA::crossAgentJump(CHROMOSOME &chromo)
{
	auto *tmpTPA = model->getVectorTasksPerAgents();
	int index1, index2, rndP, rnd2A1;

	int rndA1 = getRandomIntegerInRange(0, static_cast<int>(model->sigma - 1)); // what if i pick agent without tasks?

	int rndT = getRandomIntegerInRange(0, static_cast<int>(chromo.indTask[rndA1].size() - 1));
	if (rndT < 0)
	{
		return;
	}

	if (model->_vectorTasksPerAgents[chromo.indTask[rndA1][rndT] - model->sigma].size() > 1)
	{
		rnd2A1 = find(model->_vectorTasksPerAgents[chromo.indTask[rndA1][rndT] - model->sigma].begin(), model->_vectorTasksPerAgents[chromo.indTask[rndA1][rndT] - model->sigma].end(), rndA1) - model->_vectorTasksPerAgents[chromo.indTask[rndA1][rndT] - model->sigma].begin();
	}
	else
	{
		return;
	}

	int rndA2 = (*tmpTPA)[chromo.indTask[rndA1][rndT] - model->sigma][getRandomIntegerInRangeExcluding(0, static_cast<int>((*tmpTPA)[chromo.indTask[rndA1][rndT] - model->sigma].size() - 1), rnd2A1)]; // what if i pick agent without tasks?

	if (rndA2 >= 0 && chromo.indTask[rndA2].size() > 1)
	{
		rndP = getRandomIntegerInRange(0, static_cast<int>(chromo.indTask[rndA2].size() - 1));

		/*for (auto a : chromo.individual[rndA1]) { cout << a << " " << flush; }
		cout << endl;
		cout << endl;
		cout << chromo.indTask[rndA1][rndT] << endl;;*/
		auto it = find(chromo.individual[rndA1].begin(), chromo.individual[rndA1].end(), chromo.indTask[rndA1][rndT]); // there is some bug here...
		if (it != chromo.individual[rndA1].end())
		{
			index1 = it - chromo.individual[rndA1].begin();
			chromo.individual[rndA1][index1] = chromo.individual[rndA1][chromo.indTask[rndA1][rndT]];
			chromo.individual[rndA1][chromo.indTask[rndA1][rndT]] = 0;
		}
		else
		{
			return;
		}

		it = find(chromo.individual[rndA2].begin(), chromo.individual[rndA2].end(), chromo.indTask[rndA2][rndP]);
		if (it != chromo.individual[rndA2].end())
		{
			index2 = it - chromo.individual[rndA2].begin();
			int tmp = chromo.individual[rndA2][index2];
			chromo.individual[rndA2][index2] = chromo.indTask[rndA1][rndT];
			chromo.individual[rndA2][chromo.individual[rndA2][index2]] = tmp;
		}
		else
		{
			return;
		}

		if (chromo.indTask[rndA1].size() > 1)
		{
			std::iter_swap(chromo.indTask[rndA1].begin() + rndT, chromo.indTask[rndA1].end() - 1);
			chromo.indTask[rndA1].pop_back();
		}
		else if (chromo.indTask[rndA1].size() <= 1)
		{
			chromo.indTask[rndA1].clear();
		}
		chromo.indTask[rndA2].emplace_back(chromo.individual[rndA2][index2]);
		_mutationsDone++;

		////check consistency
		// auto a = chromo.individual[rndA1]; auto b = chromo.individual[rndA2];
		// sort(a.begin(), a.end());
		// auto last = unique(a.begin(), a.end());
		// a.erase(last, a.end());

		// if (a.size() != chromo.indTask[rndA1].size() + 2) {
		//	cout << "Govno1" << endl;
		// }

		// sort(b.begin(), b.end());
		// last = unique(b.begin(), b.end());
		// b.erase(last, b.end());

		// if (b.size() != chromo.indTask[rndA2].size() + 2) {
		//	cout << "Govno2" << endl;
		// }
	}
}

void GA::rankPop(vector<pair<int, double>> &ranking)
{

	for (int iDNA = 0; iDNA < POPULATION_SIZE; iDNA++)
		ranking[iDNA] = make_pair(iDNA, population[iDNA].cost);

	if (OBJECTIVE == 1)
	{
		std::sort(
			ranking.begin(),
			ranking.end(),
			[](pair<int, double> &left, pair<int, double> &right)
			{
				return left.second < right.second;
			});
	}
	else if (OBJECTIVE == -1)
	{
		std::sort(
			ranking.begin(),
			ranking.end(),
			[](pair<int, double> &left, pair<int, double> &right)
			{
				return left.second > right.second;
			});
	}
}

void GA::updatePopElite()
{
	vector<pair<int, double>> ranking(POPULATION_SIZE);

	rankPop(ranking);
	int k = 0;

	for (auto i = ranking.size() - 1; i > ranking.size() - ELITE - 1; --i)
	{
		population[ranking[i].first] = _elitePopulation[k];
		k++;
	}
}

void GA::updatePopBreed()
{
	vector<pair<int, double>> ranking(POPULATION_SIZE);

	rankPop(ranking);
	int k = 0;

	for (auto i = ranking.size() - 1; i > ranking.size() - _childrenPopulation.size() - 1; --i)
	{
		population[ranking[i].first] = _childrenPopulation[k];
		k++;
	}
}

void GA::pickDestinationDepot(vector<int> &_tmpChromo, int tmp, int i)
{
	double cost = 999999999999;
	int k = INDIVIDUAL_SIZE;

	for (; k < model->Vtilde; k++)
	{
		if (model->w[tmp][k][i] < cost)
		{
			cost = model->w[tmp][k][i];
		}
	}
	_tmpChromo[tmp] = k - 1;
}

void GA::createGantt(CHROMOSOME &chromo, MODEL *model)
{
	chromo.timeline.clear();
	chromo.timeline.resize(model->V);

	for (int i = 0; i < model->sigma; i++)
	{
		// printOrderedPlan(population[_bestIndividualsIndex].individual[i], i);
		// cout << i << endl;
		int _tmp = i;
		int endT = 0, cost = 0, memCost = 0, memTask;
		while (_tmp < chromo.individual[i].size())
		{
			/* Calculate the cost of the chromosome */

			if (chromo.individual[i][_tmp] >= chromo.individual[i].size())
			{
				break;
			}
			// cout << _tmp << "\t" << chromo.individual[i][_tmp] << endl;

			if (model->H[_tmp] == 0 && model->H[chromo.individual[i][_tmp]] == 1)
			{
				cost = 0;
				memCost = 0;
				memTask = _tmp;
			}
			else if (model->H[_tmp] == 1 && model->H[chromo.individual[i][_tmp]] == 0)
			{
				memCost += model->getTaskDuration(_tmp);
				if (memCost > model->getEdgeCost(memTask, chromo.individual[i][_tmp], i))
				{
					cost = 0;
				}
				else
				{
					cost = model->getEdgeCost(memTask, chromo.individual[i][_tmp], i) - memCost;
				}
				memCost = 0;
			}
			else if (model->H[_tmp] == 0 && model->H[chromo.individual[i][_tmp]] == 0)
			{
				cost = model->getEdgeCost(_tmp, chromo.individual[i][_tmp], i);
			}
			else if (model->H[_tmp] == 1 && model->H[chromo.individual[i][_tmp]] == 1)
			{
				cost = 0;
				memCost += model->getTaskDuration(_tmp);
			}

			chromo.timeline[chromo.individual[i][_tmp] - model->sigma].startTime = endT + cost;
			chromo.timeline[chromo.individual[i][_tmp] - model->sigma].endTime = chromo.timeline[chromo.individual[i][_tmp] - model->sigma].startTime + model->getTaskDuration(chromo.individual[i][_tmp]);
			chromo.timeline[chromo.individual[i][_tmp] - model->sigma].agentIdx = i;

			/*chromo.timeline[chromo.individual[i][_tmp] - model->sigma].startTime = endT + model->getEdgeCost(_tmp, chromo.individual[i][_tmp], i);
			chromo.timeline[chromo.individual[i][_tmp] - model->sigma].endTime = chromo.timeline[chromo.individual[i][_tmp] - model->sigma].startTime + model->getTaskDuration(chromo.individual[i][_tmp]);
			chromo.timeline[chromo.individual[i][_tmp] - model->sigma].agentIdx = i;*/

			// Update tmp

			endT = chromo.timeline[chromo.individual[i][_tmp] - model->sigma].endTime;
			_tmp = chromo.individual[i][_tmp];
		}
	}
}

void GA::twoOptLocalSearch(CHROMOSOME &chromo)
{
	// '2-opt' is a local search algorithm for optimizing a Travelling Salesman
	// Problem. It takes a route that crosses over itself and reorders the
	// locations to eliminate the cross over.

	// The reordering will be applied to all possible location pairs as long as
	// it can improve the total length of the route.

	// bool hasImproved = true;
	CHROMOSOME tmp = chromo, tChromo;
	vector<vector<int>> tempChromo, twoOptOutcome;
	vector<int> _tmpPC(model->Vtilde);
	tempChromo.resize(model->sigma);

	extractPath(chromo, twoOptOutcome);

	for (int i = 0; i < model->sigma; ++i)
	{

		int improve = 1;
		// bool flag = false;
		while (improve > 0)
		{

			improve = 0;

			for (int iGene1 = 1; iGene1 < twoOptOutcome[i].size() - 2; iGene1++)
			{
				for (int iGene2 = iGene1 + 1; iGene2 < twoOptOutcome[i].size() - 1; iGene2++)
				{

					tempChromo[i].clear();
					tempChromo[i].resize(twoOptOutcome[i].size());
					twoOptSwap(twoOptOutcome[i], tempChromo[i], iGene1, iGene2);
					double _cost = 0;
					/*	printOrderedPlan(chromo.individual[i], i);
						cout << endl;
						for (auto a : tempChromo[i]) { cout << a << "   " << flush; }
						cout << endl;*/
					// cout << iGene1 << "\t" << iGene2 << "\t" << i << endl;
					//  very slow, think of a better way to do it;
#ifdef XD
					vector<gantt> tmpTimeline = chromo.timeline;
					_cost = evaluateLightXD(tempChromo[i], tmpTimeline, i, false);
					if (_cost < chromo.indCost[i])
					{
						chromo.timeline.swap(tmpTimeline);
					}

#else
					_cost = evaluateLight(tempChromo[i], _tmpPC, i, false);
#endif
					if (_cost < chromo.indCost[i])
					{

						improve++;
						insertPath(tmp.individual[i], tempChromo[i]);

						chromo.individual[i].swap(tmp.individual[i]);

						twoOptOutcome[i].swap(tempChromo[i]);
						chromo.indCost[i] = _cost;
						/*flag = true;
						break;*/
					}
					// for (int j = 0; j < model->sigma; ++j) {
					//	chromo.individual[j].swap(tmp.individual[j]);
					//	twoOptOutcome[i].swap(tempChromo[i]);
					//	chromo.cost = tmp.cost;
					// }

					// ofstream myfile("resultsgif.txt", ios_base::app | ios_base::out);

					// for (int i = 0; i < chromo.individual.size(); ++i) {
					//	int k = i;

					//	myfile << i << "\t" << model->_model.A[i].x << "\t" << model->_model.A[i].y << endl;
					//	while (true) {

					//		int task = chromo.individual[i][k] - model->_model.A.size();

					//		if (task == model->_model.T.size()) {

					//			myfile << task + model->_model.A.size() << "\t" << model->_model.dest[task - model->_model.T.size()].x << "\t" << model->_model.T[task - model->_model.T.size()].y << endl;
					//			break;
					//		}
					//		else {
					//			myfile << task + model->_model.A.size() << "\t" << model->_model.T[task].x << "\t" << model->_model.T[task].y << endl;
					//		}

					//		k = task + model->_model.A.size();

					//		//cout << Genetic->population[Genetic->getBestIndividualsIndex()].individual[i][k] << model[instance].A.size() + model[instance].T.size() << endl;

					//	}
					//	myfile << "END" << endl;
					//}
					// myfile << "--" << endl;
					// myfile.close();

					//}
				}
				// if (flag) { break; }
			}
		}
	}
}

void GA::GreedySearch(CHROMOSOME &chromo)
{

	vector<int> tmpIndividual(INDIVIDUAL_SIZE), tmpPC(INDIVIDUAL_SIZE);
	double d;
	int index, task, idx = 0;
	vector<gantt> tmpTimeline = chromo.timeline;

	for (int i = 0; i < model->sigma; i++)
	{

		task = i;

		if (!chromo.indTask[i].empty())
		{

			while (idx < chromo.indTask[i].size())
			{

				double dMin = numeric_limits<double>::infinity();

				for (int j = idx; j < chromo.indTask[i].size(); j++)
				{

					d = model->getEdgeCost(task, chromo.indTask[i][j], i);

					if (d < dMin)
					{
						index = j;
						dMin = d;
					}
				}

				tmpIndividual[task] = chromo.indTask[i][index];
				swap(chromo.indTask[i][idx], chromo.indTask[i][index]);
				task = tmpIndividual[task];
				idx++;
			}

			pickDestinationDepot(tmpIndividual, task, i);
			double cost = 0;

#ifdef XD
			cost = evaluateLightXD(tmpIndividual, tmpTimeline, i, true);
#else
			cost = evaluateLight(tmpIndividual, tmpPC, i, true);
#endif
			if (cost < chromo.indCost[i])
			{
				chromo.individual[i].swap(tmpIndividual);
				chromo.indCost[i] = cost;
			}

			tmpIndividual.clear();
			tmpIndividual.resize(INDIVIDUAL_SIZE);
			idx = 0;
		}
	}

	return;
}

void GA::twoOptSwap(vector<int> &inGenes, vector<int> &outGenes, int iGene1, int iGene2)
{
	// Take inGenes[0] to inGenes[iGene1 - 1]
	// and add them in order to outGenes

	for (int iGene = 0; iGene <= iGene1 - 1; iGene++)
		outGenes[iGene] = inGenes[iGene];

	// Take inGenes[iGene1] to inGenes[iGene2] and
	// add them in reverse order to outGenes

	int iter = 0;
	for (int iGene = iGene1; iGene <= iGene2; iGene++)
	{
		outGenes[iGene] = inGenes[iGene2 - iter];
		iter++;
	}

	// Take inGenes[iGene2 + 1] to end of inGenes
	// and add them in order to outGenes

	for (int iGene = iGene2 + 1; iGene < inGenes.size(); iGene++)
		outGenes[iGene] = inGenes[iGene];
}

void GA::extractPath(CHROMOSOME &chromo, vector<vector<int>> &path)
{

	// Extract the ordered path from the existing crhomosome design, e.g., 0 4 5 6 7 2  -> 0-4-5-6-7-2-1
	//																	   4 5 6 7 2 1

	path.resize(model->sigma);

	for (int i = 0; i < model->sigma; ++i)
	{

		path.reserve(chromo.indTask[i].size() + 2);
		path[i].emplace_back(i);

		int k = i;

		while (true)
		{

			path[i].emplace_back(chromo.individual[i][k]);

			if (chromo.individual[i][k] >= model->sigma + model->V)
			{
				break;
			}

			k = chromo.individual[i][k];
		}
	}
}

void GA::insertPath(vector<int> &tmp, vector<int> &path)
{
	// Insert the ordered path, e.g., 3-4-5-6-7-2-1 in the existing chromosome design, e.g., 0 4 5 6 7 2
	//																				         4 5 6 7 2 1
	int s = tmp.size();
	tmp.clear();
	tmp.resize(s);

	for (int i = 0; i < path.size() - 1; ++i)
	{
		tmp[path[i]] = path[i + 1];
	}
}

void GA::performPCfix(CHROMOSOME &chromo, int _agent1, int _agent2, int &_start, int _end, int _idxTask2Move, int _count)
{
	int _insertIdx, _tmpVal, _tmpVal2;

	//	printOrderedPlan(chromo.individual[_agent1], _agent1);
	// printOrderedPlan(chromo.individual[_agent2], _agent2);
	// cout << " Perform PC FIX 1 " << endl;
	_insertIdx = insertionIdx(chromo.individual[_agent1], _start, _end, _count);
	_tmpVal = chromo.individual[_agent1][_insertIdx];

	_tmpVal2 = chromo.individual[_agent2][_idxTask2Move];

	chromo.individual[_agent2][_idxTask2Move] = chromo.individual[_agent2][chromo.individual[_agent2][_idxTask2Move]];
	chromo.individual[_agent2][_tmpVal2] = 0;

	// insert task
	chromo.individual[_agent1][_insertIdx] = _tmpVal2;
	chromo.individual[_agent1][_tmpVal2] = _tmpVal;

	// printOrderedPlan(chromo.individual[_agent1], _agent1);
	// printOrderedPlan(chromo.individual[_agent2], _agent2);

	// update the vector holding task indices
	chromo.indTask[_agent1].emplace_back(_tmpVal2);

	if (chromo.indTask[_agent2].size() > 1)
	{
		auto _it = find(chromo.indTask[_agent2].begin(), chromo.indTask[_agent2].end(), _tmpVal2) - chromo.indTask[_agent2].begin();
		std::iter_swap(chromo.indTask[_agent2].begin() + _it, chromo.indTask[_agent2].end() - 1);
		chromo.indTask[_agent2].pop_back();
	}
	else if (chromo.indTask[_agent2].size() <= 1)
	{
		chromo.indTask[_agent2].clear();
	}

	_start = _tmpVal2;
}

void GA::ERX(CHROMOSOME &_tmpChromo, vector<unordered_set<int>> &_aMatrix)
{
	_tmpChromo.individual.resize(model->sigma);
	_tmpChromo.indTask.resize(model->sigma);

	vector<int> _lastAssigned(model->sigma);
	int min, max, _rnd, _rnd2, _idA;

	int assignedTasks = 0;

#ifdef NO_FIX_SEED

	std::random_device rd;
	std::mt19937 g(rd());
	std::vector<int> v(model->sigma);
	std::iota(v.begin(), v.end(), 0);
	std::shuffle(begin(v), end(v), g);

#endif // NO_FIX_SEED

	for (int iii = 0; iii < model->sigma; iii++)
	{

#ifdef NO_FIX_SEED
		int i = v[iii];
#else
		int i = iii;
#endif
		if (_tmpChromo.indTask[i].size() == 0)
		{
			_lastAssigned[i] = i;
		}

		_tmpChromo.individual[i].resize(INDIVIDUAL_SIZE);

		if (population[_parent1].indTask[i].size() < population[_parent2].indTask[i].size())
		{
			min = population[_parent1].indTask[i].size();
			max = population[_parent2].indTask[i].size();
		}
		else
		{
			max = population[_parent1].indTask[i].size();
			min = population[_parent2].indTask[i].size();
		}

		if (!_aMatrix[i].empty() && *_aMatrix[i].begin() != -1)
		{
			_idA = i;
		}
		else
		{
			vector<int> _tmpTaskVector;
			int ii = 0;

			for (auto a : _aMatrix)
			{
				if (!a.empty() && *a.begin() != -1)
				{
					_tmpTaskVector.emplace_back(ii);
				}
				ii++;
			}
			if (_tmpTaskVector.empty())
			{
				for (int k = model->sigma; k < INDIVIDUAL_SIZE; k++)
				{
					if (!_aMatrix[k].empty())
					{
						_tmpTaskVector.emplace_back(k);
					}
				}
			}

			if (_tmpTaskVector.empty())
			{
				//_lastAssigned[i] = i;
				pickDestinationDepot(_tmpChromo.individual[i], _lastAssigned[i], i);
				// break;
				continue;
			}
			else
			{
				_idA = _tmpTaskVector[getRandomIntegerInRange(0, (int)_tmpTaskVector.size() - 1)];
			}
		}

		_rnd = getRandomFromSet(_aMatrix[_idA]);

		if (_rnd == -1)
		{
			_rnd = _idA;
		}

		if (checkColor(i, _rnd))
		{
			int _rnd3 = model->_vectorTasksPerAgents[_rnd - model->sigma][getRandomIntegerInRange(0, (int)model->_vectorTasksPerAgents[_rnd - model->sigma].size() - 1)];

			if (_tmpChromo.individual[_rnd3].size() == 0)
			{
				_tmpChromo.individual[_rnd3].resize(INDIVIDUAL_SIZE);
				_lastAssigned[_rnd3] = _rnd3;
			}

			_tmpChromo.indTask[_rnd3].emplace_back(_rnd);
			_tmpChromo.individual[_rnd3][_lastAssigned[_rnd3]] = _rnd;

			_lastAssigned[_rnd3] = _rnd;
			_aMatrix[_lastAssigned[_rnd3]].clear();
			_rnd = pickNextTask(_rnd, i, _aMatrix);

			if (_rnd == -1)
			{
				continue;
			}
		}
		else if (_aMatrix[_rnd].empty())
		{
			_aMatrix[_rnd].insert(-1);
		}

		if (min == max)
		{
			_rnd2 = min;
		}
		else
		{
			_rnd2 = getRandomIntegerInRange(min, max - 1);
		}

		for (int j = 0; j < _rnd2; j++)
		{
			// cout << j << endl;

			_tmpChromo.indTask[i].emplace_back(_rnd);
			_tmpChromo.individual[i][_lastAssigned[i]] = _rnd;
			if (_rnd > model->V + model->sigma)
			{
				cout << "bla bla" << endl;
			}
			_lastAssigned[i] = _rnd;
			/*}*/

			_rnd = pickNextTask(_rnd, i, _aMatrix);
			_aMatrix[_lastAssigned[i]].clear();

			//_lastAssigned[i] = _tmpChromo.individual[i][_lastAssigned[i]];

			if (_rnd == -1 && j == 0 && (_lastAssigned[i] != i) && _lastAssigned[i] == 0)
			{
				_lastAssigned[i] = i;
				break;
			}
			else if (_rnd == -1)
			{
				break;
			}
		}
	}

	for (int i = 0; i < model->sigma; i++)
	{

		pickDestinationDepot(_tmpChromo.individual[i], _lastAssigned[i], i);
	}
	allocateRest(_aMatrix, _tmpChromo);

	return;
}

pair<int, int> GA::crossAgentTaskSearch(CHROMOSOME &chromo, int _agent, int _task, bool flag)
{
	for (int i = 0; i < model->sigma; i++)
	{
		int _tmp = i;

		if (flag && i != _agent)
		{ // true if I want to skip search through agent [i], flase if I want to search for a task over all agents
			while (_tmp < chromo.individual[i].size())
			{

				if (chromo.individual[i][_tmp] == _task)
				{
					return pair<int, int>(i, _tmp);
				} // return pair, first agent index, second task index

				// Update tmp
				_tmp = chromo.individual[i][_tmp];
			}
		}
	}
	cout << "Something went horribly wrong!" << endl;
	// the function should not reach this point
	printOrderedPlan(chromo.individual[0], 0);
	pair<int, int> dummy;
	return dummy;
}

void GA::genAdjMatrix(vector<unordered_set<int>> &_return)
{
	// vector<vector<unordered_set<int>>> _aMatrix(INDIVIDUAL_SIZE);

	for (int i = 0; i < model->sigma; i++)
	{
		int _tmp = i;

		// Parent 1
		while (population[_parent1].individual[i][_tmp] < INDIVIDUAL_SIZE)
		{

			_return[_tmp].insert(population[_parent1].individual[i][_tmp]);
			if (_tmp >= model->sigma)
			{
				_return[population[_parent1].individual[i][_tmp]].insert(_tmp);
			}
			_tmp = population[_parent1].individual[i][_tmp];
		}

		// Parent 2
		_tmp = i;
		while (population[_parent2].individual[i][_tmp] < INDIVIDUAL_SIZE)
		{

			_return[_tmp].insert(population[_parent2].individual[i][_tmp]);
			if (_tmp >= model->sigma)
			{
				_return[population[_parent2].individual[i][_tmp]].insert(_tmp);
			}
			_tmp = population[_parent2].individual[i][_tmp];
		}
	}

	return;
}

void GA::allocateRest(vector<unordered_set<int>> &aMatrix, CHROMOSOME &chromo)
{
	for (int i = model->sigma; i < INDIVIDUAL_SIZE; i++)
	{
		if (!aMatrix[i].empty())
		{

			// I have 2 options, to assign it randomly, or next to the task it has in its set.... first one is easier, second one might be better?

			// if (aMatrix[i].size() == 1 && *aMatrix[i].begin() == -1) {
			int _rndAidx = getRandomIntegerInRange(0, (int)model->_vectorTasksPerAgents[i - model->sigma].size() - 1);
			int _rndA = model->_vectorTasksPerAgents[i - model->sigma][_rndAidx];
			// cout << " allocate Rest 1 " << endl;
			int _iidx2 = insertionIdx(chromo.individual[_rndA], _rndA, INDIVIDUAL_SIZE, chromo.indTask[_rndA].size());
			insertTask(chromo, aMatrix, _rndA, _iidx2, i);
			aMatrix[i].clear();
			//}
			// else {
			//	int _rnd = getRandomFromSet(aMatrix[i]);

			//	if (aMatrix[_rnd].empty()) {
			//		pair<int, int> _tmpPair = crossAgentTaskSearch(chromo, -1, _rnd, false);
			//	}
			//	else {
			//		vector<int> _intersection;
			//		std::set_intersection(model->_vectorTasksPerAgents[i - model->sigma].begin(), model->_vectorTasksPerAgents[i - model->sigma].end(),
			//			model->_vectorTasksPerAgents[_rnd - model->sigma].begin(), model->_vectorTasksPerAgents[_rnd - model->sigma].end(),
			//			back_inserter(_intersection));

			//		int _rnd2 = _intersection[getRandomIntegerInRange(0, (int)_intersection.size() - 1)];
			//		//cout << " allocate Rest 2 " << endl;
			//		int _iidx = insertionIdx(chromo.individual[_rnd2], _rnd2, INDIVIDUAL_SIZE, chromo.indTask[_rnd2].size());

			//		insertTask(chromo, aMatrix, _rnd2, _iidx, i);
			//		aMatrix[i].clear();

			//		if (getRandomTrueWithProbability(0.5)) { _iidx = i; } // 50% chance to insert the other task before the task i and 50% to insert it after;

			//		insertTask(chromo, aMatrix, _rnd2, _iidx, _rnd);
			//		aMatrix[_rnd].clear();
			//	}
			//}
			// cout << endl;
		}
	}
}

void GA::insertTask(CHROMOSOME &chromo, vector<unordered_set<int>> &aMatrix, int agent, int insertIdx, int task2Insert)
{

	int _tmpVal = chromo.individual[agent][insertIdx];

	// insert task
	chromo.individual[agent][insertIdx] = task2Insert;
	chromo.individual[agent][task2Insert] = _tmpVal;

	// printOrderedPlan(chromo.individual[_agent1], _agent1);
	// printOrderedPlan(chromo.individual[_agent2], _agent2);

	// update the vector holding task indices
	chromo.indTask[agent].emplace_back(task2Insert);

	// update adjacency matrix
	for (int i = 0; i < aMatrix.size(); i++)
	{
		auto it = find(aMatrix[i].begin(), aMatrix[i].end(), task2Insert);
		if (it != aMatrix[i].end())
		{
			aMatrix[i].erase(it);
			if (aMatrix[i].empty())
			{
				aMatrix[i].insert(-1);
			}
		}
	}
}

double GA::evaluateLight(vector<int> &individual, vector<int> &_tmpPC, int i, bool flag)
{
	double _cost = 0;
	int _penal = 0;

	_tmpPC.clear();
	_tmpPC.resize(model->Vtilde);
	if (!flag)
	{
		for (int _tmp = 0; _tmp < individual.size() - 1; _tmp++)
		{

			/* Calculate the cost of the chromosome */
			_cost += model->getEdgeCost(individual[_tmp], individual[_tmp + 1], i) + model->getTaskDuration(individual[_tmp]);

			/* Check Precedence Constraints */
			if (_tmpPC[individual[_tmp]] != 0)
			{
				_tmpPC[individual[_tmp]] = 0;
			}

			if (!model->_PV[individual[_tmp]].empty() && individual[_tmp] < individual.size() && model->_PV[individual[_tmp]][0] < (model->V + model->sigma))
			{
				_tmpPC[model->_PV[individual[_tmp]][0]] = 1;
			}
		}
	}
	else
	{
		int _tmp = i;
		while (_tmp < INDIVIDUAL_SIZE)
		{
			/* Calculate the cost of the chromosome */
			_cost += model->getEdgeCost(_tmp, individual[_tmp], i) + model->getTaskDuration(_tmp);

			/* Check Precedence Constraints */
			if (_tmpPC[individual[_tmp]] != 0)
			{
				_tmpPC[individual[_tmp]] = 0;
			}

			if (!model->_PV[individual[_tmp]].empty() && individual[_tmp] < INDIVIDUAL_SIZE)
			{
				for (int j = 0; j < model->_PV[individual[_tmp]].size(); j++)
				{
					if (model->_PV[individual[_tmp]][j] < (model->V + model->sigma))
					{
						_tmpPC[model->_PV[individual[_tmp]][j]] = 1;
					}
				}
			}

			// Update tmp
			_tmp = individual[_tmp];
		}
	}

	_penal = accumulate(_tmpPC.begin(), _tmpPC.end(), 0);

	/* Sum all costs and penalties */
	return _cost + _penal * PENALTY;
}

double GA::evaluateLightXD(vector<int> &individual, vector<gantt> &tmpTimeline, int i, bool flag)
{

	double _cost = 0, _cost2 = 0;
	int _penal = 0;
	double cost = 0;
	int memCost = 0;
	int memTask;

	// 2opt
	if (!flag)
	{
		int endT = 0;
		for (int _tmp = 0; _tmp < individual.size() - 1; _tmp++)
		{
			/* Calculate the cost of the chromosome */

			if (individual[_tmp + 1] >= INDIVIDUAL_SIZE)
			{
				break;
			}
			// cout << _tmp << "\t" << individual[_tmp] << endl;
			tmpTimeline[individual[_tmp + 1] - model->sigma].startTime = endT + model->getEdgeCost(individual[_tmp], individual[_tmp + 1], i);
			tmpTimeline[individual[_tmp + 1] - model->sigma].endTime = tmpTimeline[individual[_tmp + 1] - model->sigma].startTime + model->getTaskDuration(individual[_tmp + 1]);
			tmpTimeline[individual[_tmp + 1] - model->sigma].agentIdx = i;

			// Update tmp

			endT = tmpTimeline[individual[_tmp + 1] - model->sigma].endTime;
			//_tmp = individual[_tmp];
		}

		for (int _tmp = 0; _tmp < individual.size() - 1; _tmp++)
		{

			if (model->H[_tmp] == 0 && model->H[individual[_tmp + 1]] == 1)
			{
				cost = model->getTaskDuration(_tmp);
				memCost = 0;
				memTask = _tmp;
			}
			else if (model->H[_tmp] == 1 && model->H[individual[_tmp + 1]] == 0)
			{
				memCost += model->getTaskDuration(_tmp);
				if (memCost > model->getEdgeCost(memTask, individual[_tmp + 1], i))
				{
					cost = memCost;
				}
				else
				{
					cost = model->getEdgeCost(memTask, individual[_tmp + 1], i);
				}
				memCost = 0;
			}
			else if (model->H[_tmp] == 0 && model->H[individual[_tmp + 1]] == 0)
			{
				cost = model->getEdgeCost(_tmp, individual[_tmp + 1], i) + model->getTaskDuration(_tmp);
			}
			else if (model->H[_tmp] == 1 && model->H[individual[_tmp + 1]] == 1)
			{
				cost = 0;
				memCost += model->getTaskDuration(_tmp);
			}
			_cost += cost;

			/* Calculate the cost of the chromosome */
			_cost2 += model->getEdgeCost(individual[_tmp], individual[_tmp + 1], i) + model->getTaskDuration(individual[_tmp]);

			/* Check Precedence Constraints */
			if (!model->_PV[individual[_tmp]].empty() && individual[_tmp] < INDIVIDUAL_SIZE)
			{
				for (int j = 0; j < model->_PV[individual[_tmp]].size(); j++)
				{
					if (model->_PV[individual[_tmp]][j] < (model->V + model->sigma))
					{
						if (tmpTimeline[model->_PV[individual[_tmp]][j] - model->sigma].startTime < tmpTimeline[individual[_tmp] - model->sigma].endTime)
						{
							_penal++;
						}
					}
				}
			}
		}
	}
	// greedy search
	else
	{
		int _tmp = i;
		int endT = 0;
		while (_tmp < individual.size())
		{
			/* Calculate the cost of the chromosome */

			if (individual[_tmp] >= individual.size())
			{
				break;
			}
			// cout << _tmp << "\t" << chromo.individual[i][_tmp] << endl;
			tmpTimeline[individual[_tmp] - model->sigma].startTime = endT + model->getEdgeCost(_tmp, individual[_tmp], i);
			tmpTimeline[individual[_tmp] - model->sigma].endTime = tmpTimeline[individual[_tmp] - model->sigma].startTime + model->getTaskDuration(individual[_tmp]);
			tmpTimeline[individual[_tmp] - model->sigma].agentIdx = i;

			// Update tmp

			endT = tmpTimeline[individual[_tmp] - model->sigma].endTime;
			_tmp = individual[_tmp];
		}

		_tmp = i;
		while (_tmp < INDIVIDUAL_SIZE)
		{

			if (model->H[_tmp] == 0 && model->H[individual[_tmp]] == 1)
			{
				cost = model->getTaskDuration(_tmp);
				memCost = 0;
				memTask = _tmp;
			}
			else if (model->H[_tmp] == 1 && model->H[individual[_tmp]] == 0)
			{
				memCost += model->getTaskDuration(_tmp);
				if (memCost > model->getEdgeCost(memTask, individual[_tmp], i))
				{
					cost = memCost;
				}
				else
				{
					cost = model->getEdgeCost(memTask, individual[_tmp], i);
				}
				memCost = 0;
			}
			else if (model->H[_tmp] == 0 && model->H[individual[_tmp]] == 0)
			{
				cost = model->getEdgeCost(_tmp, individual[_tmp], i) + model->getTaskDuration(_tmp);
			}
			else if (model->H[_tmp] == 1 && model->H[individual[_tmp]] == 1)
			{
				cost = 0;
				memCost += model->getTaskDuration(_tmp);
			}

			/* Calculate the cost of the chromosome */
			_cost2 += model->getEdgeCost(_tmp, individual[_tmp], i) + model->getTaskDuration(_tmp);
			_cost += cost;

			// cout << _tmp << endl;
			/* Check Precedence Constraints */
			if (!model->_PV[individual[_tmp]].empty() && individual[_tmp] < INDIVIDUAL_SIZE)
			{
				for (int j = 0; j < model->_PV[individual[_tmp]].size(); j++)
				{
					if (model->_PV[individual[_tmp]][j] < (model->V + model->sigma))
					{
						if (tmpTimeline[model->_PV[individual[_tmp]][j] - model->sigma].startTime < tmpTimeline[individual[_tmp] - model->sigma].endTime)
						{
							_penal++;
						}
					}
				}
			}

			/*if (!model->_PV[individual[_tmp]].empty() && individual[_tmp] < INDIVIDUAL_SIZE) {
				for (int j = 0; j < model->_PV[individual[_tmp]].size(); j++) {
					if (model->_PV[individual[_tmp]][j] < (model->V + model->sigma)) {
						_tmpPC[model->_PV[individual[_tmp]][j]] = 1;
					}
				}
			}*/

			// Update tmp
			_tmp = individual[_tmp];
		}
	}

	/* Sum all costs and penalties */
	return _cost + _penal * PENALTY;
}

#pragma endregion

#pragma region == Debugging Methods ==

void GA::countTasks(CHROMOSOME &chromo)
{

	vector<int> path;
	int count = 0;
	for (int i = 0; i < model->sigma; ++i)
	{

		int k = i;

		while (true)
		{

			if (chromo.individual[i][k] >= model->sigma + model->V)
			{
				break;
			}
			path.emplace_back(chromo.individual[i][k]);
			count++;
			k = chromo.individual[i][k];
		}
	}

	sort(path.begin(), path.end());
	auto it = std::unique(path.begin(), path.end()); // 10 20 30 20 10 ?  ?  ?  ?
													 //                ^
	path.resize(std::distance(path.begin(), it));	 // 10 20 30 20 10

	if (count != model->V || path.size() != model->V)
	{
		cout << "ALARM!" << endl;
	}
}

void GA::printOrderedPlan(vector<int> &vector, int _ag)
{
	string path = "c:\\Users\\bmc01\\Documents\\Visual Studio 2017\\Projects\\MR_MT\\MR_MT\\MIPstart.txt";

	// ofstream myfile(path, ios_base::app | ios_base::out);
	cout << endl;
	int _k = _ag;
	int count = 0;
	// myfile << _k << "\t" << 0 << "\t" << flush;
	// cout << _k << "\t" << 0 << "\t" << flush;

	while (true)
	{

		if (vector[_k] >= model->sigma + model->V)
		{
			// myfile << vector[_k] << "\t" << depotStart[_ag] << "\t" << flush;
			//			cout << vector[_k] << "\t" << depotStart[_ag] << "\t" << flush;
			cout << vector[_k] << "\t" << flush;
			break;
		}
		else
		{
			// cout << _k << endl;
			// myfile << vector[_k] << "\t" << timeline[vector[_k] - model->sigma].startTime << "\t" << flush;
			// cout << vector[_k] << "\t" << timeline[vector[_k] - model->sigma].startTime << "\t" << flush;*/

			cout << vector[_k] + 1 << "," << flush;
		}

		_k = vector[_k];

		if (count > model->Vtilde)
		{
			break;
		}
		count++;
	}
	// myfile << endl;
}

#pragma endregion

// double GA::getTimeSteps(CHROMOSOME &chromo, MODEL *model, vector<vector<pair<double, double>>> &timeSteps)
// {

// 	// int maxDur = *max_element(chromo.indCost.begin(), chromo.indCost.end());
// 	int realDur = 0;
// 	for (int tickAg = 0; tickAg < model->sigma; tickAg++)
// 	{
// 		// int maxDur = chromo.indCost[tickAg];
// 		timeSteps[tickAg].resize(maxDur);
// 		int _tmp = tickAg;
// 		int localStep = 1;
// 		for (int tick = 0; tick < maxDur; tick++)
// 		{

// 			if (_tmp >= chromo.individual[tickAg].size())
// 			{
// 				if (realDur < tick) {
// 					realDur =  tick;
// 				}

// 				break;
// 			}

// 			pair<double, double> coord = calculateAgentsPosition(model, _tmp, chromo.individual[tickAg][_tmp], tickAg, localStep);
// 			if (coord.first > -1)
// 			{
// 				timeSteps[tickAg][tick] = coord;
// 				localStep++;
// 			}
// 			else
// 			{
// 				// Update tmp
// 				_tmp = chromo.individual[tickAg][_tmp];
// 				localStep = 1;
// 			}
// 		}
// 	}

// 	for (int i = 0; i < model->sigma; i++) {
// 		timeSteps[i].resize(realDur);
// 	}

// 	int durInRange = 0;
// 	int _old = 0;
// 	int range = 225; // meters squared, 30m
// 	double countInRange = 0;
// 	vector<vector<int>> change(timeSteps.size(), vector<int>(timeSteps.size())); // Defaults to zero initial value
// 	vector<int> inContact(realDur);

// 	int totalSteps = realDur;
// 	int start = 0.30 * totalSteps; // 20% of total steps
// 	int end = 0.80 * totalSteps;	  // 70% of total steps
// 	int someNumber = model->avgWeight;
// 	chromo.interNum = 0;

// 	REWARD rw(totalSteps,start,end,model->avgWeight);
// 	rw.interactionReward(t);
// 	for (int t = 0; t < realDur; t++)
// 	{
// 		for (int a = 0; a < timeSteps.size() - 1; a++)
// 		{
// 			for (int b = a + 1; b < timeSteps.size(); b++)
// 			{
// 				if (timeSteps[a][t].first != 0 && timeSteps[b][t].first != 0 && timeSteps[a][t].second != 0 && timeSteps[b][t].second != 0)
// 				{
// 					double dist = (timeSteps[a][t].first - timeSteps[b][t].first) * (timeSteps[a][t].first - timeSteps[b][t].first) + (timeSteps[a][t].second - timeSteps[b][t].second) * (timeSteps[a][t].second - timeSteps[b][t].second);
// 					if (dist < range)
// 					{
// 						durInRange++;
// 						if (change[a][b] == 0)
// 						{
// 							change[a][b] = 1;
// 							inContact[t] = 1;
// 							// countInRange -= (float)t/(float)maxDur;
// 							// countInRange -= float(maxDur) - abs((float)t - float(maxDur/2));
// 							// countInRange--;
// 								chromo.interNum++;

// 							// if (t < triangleStart || t > triangleEnd)
// 							// {
// 							// 	countInRange = 0.0; // First and last 20% are 0
// 							// }

// 							// else if (t <= ((totalSteps) / 2))
// 							// {
// 							// 	double x = static_cast<double>(t - triangleStart) / (totalSteps / 2 - triangleStart);
// 							// 	countInRange += x * someNumber; // Increase to 1000
// 							// }
// 							// else
// 							// {
// 							// 	double x = static_cast<double>(t - (totalSteps / 2)) / (triangleEnd - (totalSteps / 2));
// 							// 	countInRange += (1.0 - x) * someNumber; // Decrease back to 0
// 							// }
// 						}
// 						else if (change[a][b] == 1)
// 						{
// 							// countInRange++;

// 						}
// 					}
// 					else
// 					{
// 						change[a][b] = 0;
// 					}
// 				}
// 			}
// 		}
// 		if (_old < durInRange)
// 		{
// 			// inContact[t] = 1;
// 			_old = durInRange;
// 		}
// 	}

// 	if (chromo.interNum > 0 ){
// 	std::vector<int> indices;
// 	indices.push_back(start);
//     for (int i = 0; i < inContact.size(); ++i) {
//         if (inContact[i] == 1) {
//             indices.push_back(i);
//         }
//     }
// 	indices.push_back(end);
// 	double meanDistanceCrossings = 0;
// 	double numInd = indices.size();
// 	if (!empty(indices)) {

// 		// meanDistanceCrossings = accumulate(indices.begin(),indices.end(),0) / numInd;
// 		meanDistanceCrossings = (end - start) / numInd;
// 	}

// 	double distWeight = 0;
// 	for (int d = 0; d < numInd-1; d++) {
// 		double diff = abs((double)(indices[d+1] - indices[d]) - meanDistanceCrossings);
// 		if (diff > model->avgWeight/4){
// 		distWeight += diff;
// 		}
// 	}

// 	chromo.distWeight = distWeight;
// 	} else {
// 		chromo.distWeight = 0;
// 	}
// 	// remove("timesteps.txt");
// 	// string path = "timesteps.txt";
// 	// ofstream myfile(path, ios_base::app | ios_base::out);
// 	// for (int column = 0; column < timeSteps[0].size(); column++)
// 	// {
// 	// 	for (int row = 0; row < model->sigma; row++)
// 	// 	{

// 	// 		myfile << timeSteps[row][column].first << "\t" << timeSteps[row][column].second << "\t" << flush;
// 	// 	}
// 	// 	myfile << endl;
// 	// }

// 	return countInRange;
// }