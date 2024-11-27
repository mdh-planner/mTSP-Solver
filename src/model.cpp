#include "include/model.h"

MODEL::MODEL(Problem & model)
{
	// Constructor 
	//model.T[1].prec = 9;
	//model.T[4].prec = 8;
	///*model.T[6].prec = 8;
	//model.T[0].prec = 3;*/
	//model.T[5].prec = 6;
	////model.T[8].prec = 2;
	//model.T[8].prec = 2;
	//model.T[3].prec = 7;
	//model.T[7].prec = 5;
	//model.T[0].prec = 4;
	//model.T[8].duration = 50;
	//model.T[9].duration = 60;
	//model.T[7].duration = 100;
	//model.T[4].duration = 80;
	//model.T[5].duration = 80;
	//model.T[1].duration = 100;
	//model.T[6].virt = 0;
	//model.T[7].virt = 1;
	//model.T[8].virt = 1;
	//model.T[9].virt = 1;
	//model.T[4].virt = 1;
	//model.T[5].virt = 1;
	//model.T[7].para.emplace_back(1);
	//model.T[2].para.emplace_back(5);
	///*model.T[1].para.emplace_back(4);*/
	//model.T[3].x -= 150;
	//model.T[3].y -= 100;

	//model.dest[0].x += 100;
	//model.dest[0].y += 50;

	//model.T[2].color = 4;
	////model.T[6].color = 4;
	//model.T[8].color = 4;

	//model.T[0].color = 1;
	//model.T[4].color = 1;

	//model.T[6].color = 4;
	//model.T[5].color = 4;
	//model.T[3].color = 0;
	//model.T[1].color = 0;
	//model.T[7].color = 0;
	//model.T[9].color = 0;
	//model.T[6].duration = 100;
	//model.A[2].color.clear(); model.A[2].color.emplace_back(4);
	//model.A[1].color.clear(); model.A[1].color.emplace_back(1);
	//model.A[0].color.clear(); model.A[0].color.emplace_back(0);
	sigma = model.A.size();
	V = model.T.size();
	delta = model.dest.size();
	Vtilde = sigma + V + delta;
	double minX = std::numeric_limits<double>::max(), maxX = 0, minY = std::numeric_limits<double>::max(), maxY = 0;
	

	/* Read Parallel Tasks */
	R = vector<vector<int> >(Vtilde, vector<int>(Vtilde));
	for (int i = 0; i < Vtilde; i++) {
		if (i >= sigma && i < V + sigma) {
			model.T[i - sigma].index = i - sigma;
			for (int j = 0; j < model.T[i - sigma].para.size(); j++) {
				//cout << i << " " << j << endl;
				R[i][model.T[i - sigma].para[j] + sigma] = 1;
				R[model.T[i - sigma].para[j] + sigma][i] = 1;
			}
		}
	}
	
	int id = 0;
	//for (auto a : R) {
	//	cout << id << "\t" << flush;
	//	for (auto b : a) { 
	//		cout << b << " " << flush; 
	//	} 
	//	cout << endl; 
	//	id++;
	//}
	//cout << "============" << endl;
	/* This is done in a very shitty way */
	if (true) {
		int ii = 0;
		for (int i = 0; i < V; i++) {
			if (model.T[i].virt == 1) {
				//for (auto a : R[ii + sigma]) { cout << a << " " << flush; } cout << endl;
				std::vector<int>::iterator it = find(R[ii + sigma].begin(), R[ii + sigma].end(), 1);
				if (it != R[ii + sigma].end()) {
					model.T.erase(model.T.begin() + i);
					i--; V--; Vtilde--;
				}
			}
			ii++;
		}
	}

	_model = model;

	H.resize(Vtilde); K.resize(Vtilde);

	for (int i = sigma; i < Vtilde - delta; i++) {
		/* Read Virtual Tasks */
		H[i] = model.T[i - sigma].virt;

		/* Read Agents per Task */
		K[i] = model.T[i - sigma].reqA;
	}

	/* Create Agents : Tasks feasibility matrix */
	_tasksPerAgents = vector<vector<int> >(sigma, vector<int>(V));
	_vectorTasksPerAgents.resize(V);
	_listTasksPerAgents.resize(sigma);

	for (int i = 0; i < V; i++) {
		for (int j = 0; j < sigma; j++) {
			for (int k = 0; k < model.A[j].color.size(); k++) {
				if (model.T[i].color == model.A[j].color[k]) {
					_tasksPerAgents[j][i] = 1;
					_listTasksPerAgents[j].emplace_back(i);
					_vectorTasksPerAgents[i].emplace_back(j);
					//_vtpaID[i].emplace_back()
					break;
				}
			}
			sort(_listTasksPerAgents[j].begin(), _listTasksPerAgents[j].end());
			sort(_vectorTasksPerAgents[i].begin(), _vectorTasksPerAgents[i].end());
		}
	}

	/* Read Precedence Constraints */
	P = vector<vector<int> >(V, vector<int>(V));
	_PV.resize(Vtilde);
	vpc.resize(V);
	
	for (int i = 0; i < V; i++) {
		for (int j = 0; j < V; j++) {
			if (model.T[i].prec == j) {
				P[i][j] = 1;
				//cout << i << "\t" << j << endl;
				vpc[i].before.emplace_back(j);
				vpc[j].after.emplace_back(i);
			}
		}
	}

	for (int i = sigma; i < V + sigma; i++) {
		if (model.T[i - sigma].prec >= 0) {
			_PV[i].emplace_back(model.T[i - sigma].prec + sigma);
		}
	}


	taskDuration = vector<int>(Vtilde);

	for (int i = 0; i < Vtilde; i++) {
		if (i >= sigma && i < V + sigma) {
			taskDuration[i] = model.T[i - sigma].duration;
		}
	}

	if (model.T[0].x > -1) {
		/* Weight matrix */
		Vtilde2.resize(Vtilde);

		for (int i = 0; i < Vtilde2.size(); i++) {
			if (i < sigma) {
				Vtilde2[i].first = model.A[i].x;
				Vtilde2[i].second = model.A[i].y;
			}
			else if (i >= sigma && i < sigma + V) {
				Vtilde2[i].first = model.T[i - sigma].x;
				Vtilde2[i].second = model.T[i - sigma].y;
			}
			else {
				Vtilde2[i].first = model.dest[i - sigma - V].x;
				Vtilde2[i].second = model.dest[i - sigma - V].y;
			}
		}



		/*w = vector<vector<vector<int>>>(Vtilde2.size(), vector<vector<int>>(Vtilde2.size(), vector<int>(sigma)));*/
		w = vector<vector<vector<int> > >(Vtilde2.size(), vector<vector<int> >(Vtilde2.size(), vector<int>(sigma)));

		for (int i = 0; i < Vtilde2.size(); i++) {
			for (int j = i; j < Vtilde2.size(); j++) {
				for (int s = 0; s < sigma; s++) {

					if (j >= sigma && j < sigma + V && H[j] == 1) { // weight to virtual tasks is 0
						w[i][j][s] = 0;
						w[j][i][s] = 0;
					}
					else if ((i < sigma || i >= sigma + V) && (j < sigma || j >= sigma + V)) {
						w[i][j][s] = 0;
						w[j][i][s] = 0;
					}
					else {
						w[i][j][s] = sqrt((Vtilde2[i].first - Vtilde2[j].first) * (Vtilde2[i].first - Vtilde2[j].first) + (Vtilde2[i].second - Vtilde2[j].second) * (Vtilde2[i].second - Vtilde2[j].second)) / model.A[s].v;
						w[j][i][s] = w[i][j][s];

						if (Vtilde2[i].first > maxX) {
							maxX = Vtilde2[i].first;
						}
						 if (Vtilde2[i].first < minX) {
							minX = Vtilde2[i].first;
						 }

							if (Vtilde2[i].second > maxY) {
							maxY = Vtilde2[i].second;
						}
						 if (Vtilde2[i].second < minY) {
							minY = Vtilde2[i].second;
						 }


					}


				}

			}
		}
	}
	else {
		//nAU = model.numberOfUsedAgents;

		//w = vector<vector<vector<double>>>(Vtilde, vector<vector<double>>(Vtilde, vector<double>(sigma)));
		//for (int i = 0; i < Vtilde; i++) {
		//	for (int j = 0; j < Vtilde; j++) {
		//		for (int s = 0; s < sigma; s++) {
		//			if ((i < sigma || i >= sigma + V) && (j < sigma || j >= sigma + V)) {
		//				w[i][j][s] = 0;
		//			
		//			}
		//			else {
		//				w[i][j][s] = model.transitionMatrix[i][j];
		//			}
		//			
		//			
		//		}
		//		//cout << w[i][j][0] << "\t" << flush;
		//	}
		//	//cout << endl;
		//}
	}
	double totalDistance = 0.0;
    int numPairs = 0;
	for (int i = 0; i < Vtilde; ++i) {
        for (int j = i + 1; j < Vtilde; ++j) {
            totalDistance += w[i][j][0];
            numPairs++;
        }
    }
	double area = (maxX-minX) * (maxY - minY) / model.A[0].commRadius; // it is assumed that the comm radius is the same for every agent.
	// cout << maxX << " " << maxY << "  " << minX << "  " << minY << "  " << area << endl;
    // avgWeight = totalDistance / numPairs + area/(sigma-1);
	avgWeight = totalDistance / numPairs;
	//for (int i = 0; i < Vtilde2.size(); i++) {
	//	for (int j = 0; j < Vtilde2.size(); j++) {
	//		cout << "  " << w[i][j][0] << flush;
	//	}
	//	cout << endl;
	//}


	//for (int i = 0; i < w.size(); i++) {
	//	for (int j = 0; j < w[i].size(); j++) {
	//		cout << w[i][j][0] << "\t" << flush;
	//	}
	//	cout << endl;
	//}
}

int MODEL::getTaskDuration(int i) 
{
	return taskDuration[i];
}

vector<int> MODEL::findAgents(int Task) {
	vector<int> _return;

	for (int i = 0; i < _tasksPerAgents.size(); i++) {
		for (int j = 0; j < _tasksPerAgents[i].size(); j++) {
			if (Task == _tasksPerAgents[i][j]) {
				_return.emplace_back(i);
				break;
			}
		}
	}

	return _return;
}



//
//MODEL::MODEL(Problem &model)
//{
	//// Constructor 
	//sigma = model.A.size();
	//V = model.T.size();
	//delta = model.dest.size();
	//Vtilde = sigma + V + delta;
	//_model = model;
	//
	///* Weight matrix */
	//Vtilde2.resize(Vtilde);
	//
	//for (int i = 0; i < Vtilde2.size(); i++) {
	//	if (i < sigma) {
	//		Vtilde2[i].first = model.A[i].x;
	//		Vtilde2[i].second = model.A[i].y;
	//	}
	//	else if (i >= sigma && i < sigma + V) {
	//		Vtilde2[i].first = model.T[i - sigma].x;
	//		Vtilde2[i].second = model.T[i - sigma].y;
	//	}
	//	else {
	//		Vtilde2[i].first = model.dest[i - sigma - V].x;
	//		Vtilde2[i].second = model.dest[i - sigma - V].y;
	//	}
	//}

	//H.resize(Vtilde); K.resize(Vtilde);

	//for (int i = sigma; i < Vtilde - delta; i++) {
	//	/* Read Virtual Tasks */
	//	H[i] = model.T[i - sigma].virt;

	//	/* Read Agents per Task */
	//	K[i] = model.T[i - sigma].reqA;
	//}

	//w = vector<vector<vector<int>>>(Vtilde2.size(), vector<vector<int>>(Vtilde2.size(), vector<int>(sigma)));

	//for (int i = 0; i < Vtilde2.size(); i++) {
	//	for (int j = i; j < Vtilde2.size(); j++) {
	//		for (int s = 0; s < sigma; s++) {

	//			if (j >= sigma && j < sigma + V && H[j] == 1) {
	//				w[i][j][s] = 0;
	//				w[j][i][s] = 0;
	//			} else {
	//				w[i][j][s] = sqrt((Vtilde2[i].first - Vtilde2[j].first) * (Vtilde2[i].first - Vtilde2[j].first) + (Vtilde2[i].second - Vtilde2[j].second) * (Vtilde2[i].second - Vtilde2[j].second)) / model.A[s].v;
	//				w[j][i][s] = w[i][j][s];
	//			}
	//		}
	//	}
	//}

	///* Create Agents : Tasks feasibility matrix */
	//_tasksPerAgents = vector<vector<int>>(sigma, vector<int>(V ));
	//_vectorTasksPerAgents.resize(V);
	//_listTasksPerAgents.resize(sigma);

	//for (int i = 0; i < V; i++) {
	//	for (int j = 0; j < sigma; j++) {
	//		for (int k = 0; k < model.A[j].color.size(); k++) {
	//			if (model.T[i].color == model.A[j].color[k]) {
	//				_tasksPerAgents[j][i] = 1;
	//				_listTasksPerAgents[j].emplace_back(i);
	//				_vectorTasksPerAgents[i].emplace_back(j);
	//				break;
	//			}
	//		}
	//		sort(_listTasksPerAgents[j].begin(), _listTasksPerAgents[j].end());
	//		sort(_vectorTasksPerAgents[i].begin(), _vectorTasksPerAgents[i].end());
	//	}
	//}

	///* Read Precedence Constraints */
	//P = vector<vector<int>>(V, vector<int>(V));
	//_PV.resize(Vtilde);

	//for (int i = 0; i < V; i++) {
	//	for (int j = 0; j < V; j++) {	
	//		if (model.T[i].prec == j) {
	//			P[i][j] = 1;
	//		}
	//	}
	//}
	//
	//for (int i = sigma; i < V + sigma; i++) {
	//	if (model.T[i-sigma].prec != 0) {
	//		_PV[i].emplace_back(model.T[i-sigma].prec+sigma);
	//	}
	//}

	///* Read Parallel Tasks */
	//R = vector<vector<int>>(Vtilde, vector<int>(Vtilde));
	//for (int i = 0; i < Vtilde; i++) {
	//	if (i >= sigma && i < V + sigma) {
	//		for (int j = 0; j < model.T[i - sigma].para.size(); j++) {
	//			R[i][model.T[i - sigma].para[j] + sigma] = 1;
	//			R[model.T[i - sigma].para[j] + sigma][i] = 1;
	//		}
	//	}
	//}
	//
	//taskDuration = vector<int>(Vtilde);

	//for (int i = 0; i < Vtilde; i++) {
	//	if (i >= sigma && i < V + sigma) {
	//		taskDuration[i] = model.T[i - sigma].duration;
	//	}
	//}
//}