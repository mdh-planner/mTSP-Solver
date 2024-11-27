#include "../include/definitions.h"
#include "../include/randnumgenFixed.h"
#include <numeric>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <string_view>
#include <vector>
#include <fstream>

using namespace std;

void pushVirtual(Problem &object);
void pushParallel(Problem &object);
int selectNumTask(vector<double> & probV);
int findMaxAgentNumber(int color, Problem & object);
void write2file(Problem &object);

vector<double> parallelTasksProbability(vector<int> & paratemp);

vector<Problem> problemGenerator() {

	Problem object;
	
	/* Number of Test Instances*/
	int N = 1;
	vector<Problem> _return(N);

	/* Max number of Colors */
	constexpr auto COLOR_NUMBER = 1;
	std::vector<int> colors(COLOR_NUMBER), tempPara;
	std::iota(colors.begin(), colors.end(), 0);

	/* Probabilities for tasks having precedence constraint "pcProb", for task being eligable for parallel execution "paraProb", and for task being virtual "virtProb" */
	constexpr auto pcProb = 0;
	constexpr auto paraProb = 0;
	double virtProb = 0;
	/*cout << "Enter % of Virtual Tasks (0.0 - 1.0):  " << flush;
	cin >> virtProb;*/
	
//	double paraProb = 0.5;
	/*cout << "Enter % of Parallel Tasks (0.1 - 0.9):  " << flush;
	cin >> paraProb;*/

	object.pParallel = paraProb;	object.pVirtual = virtProb;

	/* Data Set; This can be read from a file as well */

	/* used in JINT paper */
	
	vector<int> Agents = { 4,1,2,2,3,3,4,4,5,5 };
	vector<int> Tasks = { 23,8,10,12,10,16,18,20,25, 30 };
	std::vector<int> SourceDepots = {3,1,2,2,3,3,4,4,5,5};
	vector<int> DestinationDepots = { 1,1,2,2,3,3,4,4,5,5 };
	//vector<int> SourceDepots = { 1,1,2,2,3,3,4,4,5,5 };
	//vector<int> DestinationDepots = { 1,1,2,2,3,3,4,4,5,5 };
	//vector<int> Agents = { 1,1,2,2,3,3,4,4,5,5 };
	//vector<int> Tasks = { 5,8,10,12,10,16,18,20,25, 30 };

	/* Navigation Area x,y [meters] */
	double xMin = 0, xMax = 300, yMin = 0, yMax = 200;

	/*Task duration min max [seconds] */
	double taskMin = 1, taskMax = 1;

	/* Max velocity of an agent */
	double velMax = 5.0;

	double commRadius = 225; // this is squared
	/* Instance Generation Loop*/
	for (int inst = 0; inst < N; inst++) {
		int pccount = 0;

		//Initialize vectors
		object.src.resize(SourceDepots[inst]);
		object.dest.resize(DestinationDepots[inst]);
		object.A.resize(Agents[inst]);
		object.T.resize(Tasks[inst]);

		/* Create Source Depots with randomized X and Y coordiantes */
		for (int sdepot = 0; sdepot < SourceDepots[inst]; sdepot++) {
			object.src[sdepot].x = getRandomRealInRangeF(xMin, xMax);
			object.src[sdepot].y = getRandomRealInRangeF(yMin, yMax);
		}

		/* Create Destination Depots with randomized X and Y coordiantes */
		for (int ddepot = 0; ddepot < DestinationDepots[inst]; ddepot++) {
			object.dest[ddepot].x = getRandomRealInRangeF(xMin, xMax);
			object.dest[ddepot].y = getRandomRealInRangeF(yMin, yMax);
		}

		/* Create Agents with randomized X and Y coordiantes and randomized colors */
		for (int agent = 0; agent < Agents[inst]; agent++) {
			auto tempColors = colors;

			auto rndSource = getRandomIntegerInRangeF(0, static_cast<int>(SourceDepots[inst] - 1));
			object.A[agent].x = object.src[rndSource].x;
			object.A[agent].y = object.src[rndSource].y;
			object.A[agent].v = 2.0;
			object.A[agent].commRadius = commRadius;
			
			int rnd = getRandomIntegerInRangeF(1, COLOR_NUMBER);
			object.A[agent].color.resize(rnd);

			for (int i = 0; i < rnd; i++) {
				int tmp = tempColors.size() - 1;
				auto rnd2 = getRandomIntegerInRangeF(0, tmp);
				object.A[agent].color[i] = tempColors[rnd2];
				tempColors.erase(tempColors.begin() + rnd2);
			}
		}

		/* Create Agents with randomized X and Y coordiantes and randomized colors */
		int cc = 0;
		for (int task = 0; task < Tasks[inst]; task++) {
			object.T[task].x = getRandomRealInRangeF(xMin, xMax);
			object.T[task].y = getRandomRealInRangeF(yMin, yMax);
			object.T[task].duration = getRandomRealInRangeF(taskMin, taskMax);
			object.T[task].color = getRandomIntegerInRangeF(0, COLOR_NUMBER - 1);

			/* Find the maximum number of agents with certain color. If no agent is found, the color is added to one of the agents randomly */
			int MAX = findMaxAgentNumber(object.T[task].color, object);

			/* Randomly set the maximum number of agents required per task */
			object.T[task].reqA = getRandomIntegerInRangeF(1, MAX);

			// Probability of task being Virtual
			auto limit = getRandomRealInRangeF(0.0, 1.0);

			if (limit < virtProb) {
				object.T[task].virt = 1;
				//tempPara.emplace_back(inst);
				//tempPara.emplace_back(task);
				tempPara.emplace_back(Tasks[inst] - cc - 1);
				cc++;
			}
			//
			//// Probability of having Precedence Constraint
			//limit = getRandomRealInRangeF(0.0, 1.0);

			//if (limit < pcProb) {
			//	object.T[task].prec = getRandomIntegerInRangeExcludingF(0, Tasks[inst] - 1, task);
			//}

		}

		/* Push Virtual Tasks to the end of the Task array */
		pushVirtual(object);

		
		//vector<double> probV = parallelTasksProbability(tempPara);

		//int start = Tasks[inst] - round(Tasks[inst] * paraProb);

		//for (int task = start; task < Tasks[inst]; task++) {

		//	if (tempPara.size() > 1) {

		//		auto tempPara2 = tempPara;
		//		int rnd = selectNumTask(probV);

		//		// remove task -> it cannot be parallel with itself
		//		auto it = find(tempPara2.begin(), tempPara2.end(), task);
		//		if (it != tempPara2.end()) {
		//			tempPara2.erase(it);
		//		}

		//		//remove tasks that have PC relations, they cannot be parallel at the same time while one precedes the other
		//		it = find(tempPara2.begin(), tempPara2.end(), object.T[task].prec);
		//		if (it != tempPara2.end()) {
		//			tempPara2.erase(it);
		//		}

		//		for (int t = 0; t < object.T.size(); t++) {
		//			if (object.T[t].prec == task) {
		//				it = find(tempPara2.begin(), tempPara2.end(), t);
		//				if (it != tempPara2.end()) {
		//					tempPara2.erase(it);
		//				}
		//			}
		//		}

		//		for (int ij = 0; ij < rnd; ij++) {
		//			//cout << ij << endl;
		//			auto rnd2 = getRandomIntegerInRangeF(0, (int)tempPara2.size() - 1);

		//			object.T[task].para.emplace_back(tempPara2[rnd2]);

		//			auto it = find(tempPara2.begin(), tempPara2.end(), tempPara2[rnd2]);
		//			if (it != tempPara2.end()) {
		//				tempPara2.erase(it);
		//			}

		//		}
		//	}
		//}
		vector<int> taskList(Tasks[inst]);
		iota(taskList.begin(), taskList.end(), 0);
//		cout << inst << endl;
		vector<double> probV = parallelTasksProbability(taskList);

		for (int task = 0; task < Tasks[inst]; task++) {

			// Probability of parallel Tasks
			auto limit = getRandomRealInRangeF(0.0, 1.0);

			if (limit < paraProb) {

				auto tempPara2 = tempPara;

				// remove task -> it cannot be parallel with itself
				auto it = find(tempPara2.begin(), tempPara2.end(), task);
				if (it != tempPara2.end()) {
					tempPara2.erase(it);
				}

				//remove tasks that have PC relations, they cannot be parallel at the same time while one precedes the other
				it = find(tempPara2.begin(), tempPara2.end(), object.T[task].prec);
				if (it != tempPara2.end()) {
					tempPara2.erase(it);
				}
				for (int t = 0; t < object.T.size(); t++) {
					if (object.T[t].prec == task) {
						it = find(tempPara2.begin(), tempPara2.end(), t);
						if (it != tempPara2.end()) {
							tempPara2.erase(it);
						}
					}
				}

				if (!tempPara2.empty()) {
					int rnd = getRandomIntegerInRangeF(0, (int)tempPara2.size() - 1); //bilo je pocinjalo od 1 ne od 0. ne znam zasto.
					object.T[task].para.emplace_back(tempPara2[rnd]);
				}
			}

		}

		pushParallel(object);

		std::vector<int> l(Tasks[inst]);
		std::iota(l.begin(), l.end(), 0);

		int N = object.T.size();

		/* Read Parallel Tasks */
		std::vector<std::vector<int>> R = vector<vector<int>>(N, vector<int>(N));
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < object.T[i].para.size(); j++) {
				//cout << i << " " << object.T[i].para[j] << endl;
				R[i][object.T[i].para[j]] = 1;
				R[object.T[i].para[j]][i] = 1;
			}

		}

		for (int task = 0; task < Tasks[inst]; task++) {
			// Probability of having Precedence Constraint
			auto limit = getRandomRealInRangeF(0.0, 1.0);

			if (limit < pcProb) {
				vector<int> tmpTasksForPrec;
				pccount++;
				if (!object.T[task].para.empty()) {
					vector<int> jediGovna;
					for (int govno = 0; govno < R[task].size(); govno++) {
						if (R[task][govno] == 1) {
							jediGovna.emplace_back(govno - Agents[inst]);
						}
					}

					sort(jediGovna.begin(), jediGovna.end());
					set_difference(l.begin(), l.end(), jediGovna.begin(), jediGovna.end(),  inserter(tmpTasksForPrec, tmpTasksForPrec.begin()));
				}
				else {
					tmpTasksForPrec = l;
				}

				object.T[task].prec = getRandomIntegerInRangeExcludingF(0, (int)(tmpTasksForPrec.size() - 1), task);
			}
			else {
				object.T[task].prec = -1;
			}

		}

		//cout << "Instance " << inst << flush;
		//for (auto a : object.T) { cout << a.color << "\t" << flush; }
		//cout << endl;

		// object.A[0].x = 50; object.A[0].y = 20;
		// object.A[1].x = 150; object.A[1].y = 20;

		// object.T[0].x = 50; object.T[0].y = 40;
		// object.T[1].x = 50; object.T[1].y = 70;
		// object.T[2].x = 50; object.T[2].y = 100;
		// object.T[3].x = 50; object.T[3].y = 130;
		// object.T[4].x = 50; object.T[4].y = 160;
		// object.T[5].x = 150; object.T[5].y = 40;
		// object.T[6].x = 150; object.T[6].y = 70;
		// object.T[7].x = 150; object.T[7].y = 100;
		// object.T[8].x = 150; object.T[8].y = 130;
		// object.T[9].x = 150; object.T[9].y = 160;

		// object.dest[0].x = 100; object.dest[0].y = 180;

		// object.A[0].x = 50; object.A[0].y = 20;
		// object.A[1].x = 150; object.A[1].y = 20;

		// object.T[0].x = 20; object.T[0].y = 40;
		// object.T[1].x = 30; object.T[1].y = 70;
		// object.T[2].x = 40; object.T[2].y = 100;
		// object.T[3].x = 50; object.T[3].y = 130;
		// object.T[4].x = 30; object.T[4].y = 160;
		// object.T[5].x = 140; object.T[5].y = 40;
		// object.T[6].x = 150; object.T[6].y = 70;
		// object.T[7].x = 160; object.T[7].y = 100;
		// object.T[8].x = 150; object.T[8].y = 130;
		// object.T[9].x = 140; object.T[9].y = 160;

		// object.dest[0].x = 100; object.dest[0].y = 180;

		write2file(object);
		_return[inst] = object;
		object.src.clear();;
		object.dest.clear();
		object.A.clear();
		object.T.clear();
		tempPara.clear();
	}

	return _return;
}

int findMaxAgentNumber(int color, Problem &object) {
	int _return = 0;

	/* Count how many agents have the certain color */
	for (int i = 0; i < object.A.size(); i++) {
		for (int j = 0; j < object.A[i].color.size(); j++) {
			if (color == object.A[i].color[j]) {
				_return++;
				break;
			}
		}
	}

	/* A fix if no agent can perform the task */
	if (_return == 0) {
		int tmpA = object.A.size() - 1;
		int agentIdx = getRandomIntegerInRangeF(0, tmpA);
		object.A[agentIdx].color.emplace_back(color);
		return 1;
	}

	return _return;
}

void pushVirtual(Problem &object) {
	vector<Task> T(object.T.size());
	int _physical = 0;
	int _virtual = T.size() - 1;

	for (int i = 0; i < T.size(); i++) {
		if (object.T[i].virt == 0) {
			T[_physical] = object.T[i];
			_physical++;
		}
		else {
			T[_virtual] = object.T[i];
			_virtual--;
		}
	}

	object.T = T;

	return;
}

void pushParallel(Problem &object) {
	int N = object.T.size();

	/* Read Parallel Tasks */
	std::vector<std::vector<int>> R = vector<vector<int>>(N, vector<int>(N));
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < object.T[i].para.size(); j++) {
		//	cout << i << " " << object.T[i].para[j] << endl;
			R[i][object.T[i].para[j]] = 1;
			R[object.T[i].para[j]][i] = 1;
		}

	}

	vector<Task> T(N);
	int _notParallel = 0;
	int _Parallel = T.size() - 1;

	for (int i = 0; i < T.size(); i++) {
		int sum = accumulate(R[i].begin(), R[i].end(), 0);
		if (sum == 0 || object.T[i].virt == 0) {
			T[_notParallel] = object.T[i];
			_notParallel++;
		}
		else {
			T[_Parallel] = object.T[i];
			_Parallel--;
		}
	}

	object.T = T;

	return;
}

vector<double> parallelTasksProbability(vector<int> & paratemp) {

	vector<double> tmpP(paratemp.size());
	//double sum = (paratemp.size() + 1) * paratemp.size() / 2;
	//tmpP[0] = paratemp.size() / sum;
	//for (int i = paratemp.size() - 1; i > 0; i--) {
	//	tmpP[paratemp.size() - i] = tmpP[paratemp.size() - i - 1] + i / sum;
	//}

	
	
	for (int i = 0; i < paratemp.size(); i++) {
		tmpP[i] = 1.0 / paratemp.size();
	}

	return tmpP;

}

int selectNumTask(vector<double> & probV) {
	int binNum = 0;
	double p = getRandomRealInRangeF(0.0, 1.0);
	if (0 <= p && p < probV[binNum])
		return binNum + 1;

	for (; binNum < probV.size() - 1; ++binNum) {
		if (probV[binNum] <= p && p < probV[binNum + 1]) {
			return binNum + 1;
		}
	}

	return probV.size() - 1;

}

void write2file(Problem &object) {
	remove("instance.txt");
	string path = "instance.txt";
	ofstream myfile(path, ios_base::app | ios_base::out);

	for (int i = 0; i < object.A.size(); i++) {
		myfile << object.A[i].x << "\t" << object.A[i].y << endl;
	}

	myfile << -1 << "\t" << -1 << endl;

	for (int i = 0; i < object.T.size(); i++) {
		myfile << object.T[i].x << "\t" << object.T[i].y << endl;
	}

	myfile << -1 << "\t" << -1 << endl;

	for (int i = 0; i < object.dest.size(); i++) {
		myfile << object.dest[i].x << "\t" << object.dest[i].y << endl;
	}
}