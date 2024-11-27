#include "include/CPScheduler.h"

void findParallel(vector<int> & _vector, vector<int> & _tasksPerAgents, vector<int> & _return, int sigma, int delta, int task);

class IncumbentSolutionCallback : public IloCP::Callback {
private:
	ILOSTD(ostream&) _out;
	string pathbb, pathobj;
	IloTimer t;
	IloIntervalVarArray D;
	MODEL * _model;

public:
	IncumbentSolutionCallback(ostream& out, IloIntervalVarArray &D, MODEL * _model, string &_pathbb, string &_pathobj, IloTimer &timer) : _out(out), D(D), _model(_model), pathbb(_pathbb), pathobj(_pathobj), t(timer) { }

	void invoke(IloCP cp, IloCP::Callback::Reason reason) {

		if (reason == IloCP::Callback::Solution) {
			ofstream myfile(pathobj, ios_base::app | ios_base::out);
			
			auto obj = cp.getObjValue();

			//cout << cp.getObjValue() << endl;
			//cout << endl;
			//for (int i = 0; i < D.getSize(); i++) {
			//	cout << "Depot Value: " << cp.getStart(D[i]) << endl;
			//	for (int j = 0; j < _model->getDestinationDepotNumber(); j++) {
			//		if (cp.getStart(D[i]) == _model->getEdgeCost(i, j + _model->getSourceDepotNumber() + _model->getTaskNumber(), i)) {
			//			obj -= cp.getStart(D[i]);
			//			break;
			//		}
			//	}
			//}

			//cout << obj << endl;
			myfile << obj << "\t" << cp.getTime() << endl;
			myfile.close();
		}

		if (reason == IloCP::Callback::ObjBound) {
			ofstream myfile2(pathbb, ios_base::app | ios_base::out);
			myfile2 << cp.getObjBound() << "\t" << cp.getTime() << endl;
			myfile2.close();
		}
	}
};

CPScheduler::CPScheduler(shared_ptr<GA> Genetic, MODEL *_model, string &_path, int _run, double runTime, double _timeLimit)
{
	ga = Genetic;
	fullModel = _model;
	timeLimit = _timeLimit;
	timeLimitPartial = timeLimit - runTime;
	weight = 0.1;
	path = _path;
	run = _run;
	_fixedAgents2Tasks = fullModel->_vectorTasksPerAgents;

	for (int i = 0; i < ga->model->sigma; i++) {
		int _tmp = i;

		while (_tmp < ga->population[ga->_bestIndividualsIndex].individual[i].size()) {

			if (ga->population[ga->_bestIndividualsIndex].individual[i][_tmp] >= ga->population[ga->_bestIndividualsIndex].individual[i].size()) {
				break;
			}

			_fixedAgents2Tasks[ga->population[ga->_bestIndividualsIndex].individual[i][_tmp] - ga->model->sigma].clear();
			_fixedAgents2Tasks[ga->population[ga->_bestIndividualsIndex].individual[i][_tmp] - ga->model->sigma].emplace_back(i);

			// Update tmp
			_tmp = ga->population[ga->_bestIndividualsIndex].individual[i][_tmp];
		}
	}
}

CPScheduler::CPScheduler(MODEL * _model, string & _path, double _timeLimit)
{

	fullModel = _model;
	timeLimit = _timeLimit;
	weight = 0.1;
	path = _path;
	run = 0;

}

/* Solvers */

void CPScheduler::solverPartial()
{
	/* CP MODEL */
	int sigma = fullModel->sigma;
	int delta = fullModel->delta;
	int Vp = ga->model->V;
	int V = fullModel->V;
	int Vtilde = V + sigma + delta;

	// create a list of interelated tasks
	set<int> pcc;
	for (int i = 0; i < fullModel->_PV.size(); i++) {
		if (!fullModel->_PV[i].empty()) {
			pcc.insert(i);
			pcc.insert(fullModel->_PV[i][0]);
		}
	}


	IloEnv env;
	IloModel model(env);
	string name;

	IloIntervalVar tmp(env);
	IloIntervalVarArray Tasks(env, V), D(env, sigma), SrcDepot(env, sigma);
	IloIntervalVarArray2 A(env), TasksV(env, V), noOverlapArray(env), DestDepot2(env, sigma);
	IloIntArray POS(env);

	vector<int> index;
	int count = 0;

	/* Create array of Physical Tasks */
	for (int i = 0; i < Vtilde; i++) {
		if (i < sigma || i >= V + sigma) {
			POS.add(count++);
			index.emplace_back(i);
		}
		else if (fullModel->H[i] != 1) {
			POS.add(count++);
			index.emplace_back(i);
		}
	}

	/* Intervals representing the source depot of agents */
	for (int i = 0; i < sigma; i++) {
		IloIntervalVar tmp(env);
		SrcDepot[i] = tmp;
		SrcDepot[i].setSizeMin(1); SrcDepot[i].setSizeMax(1);
		SrcDepot[i].setStartMin(-1); SrcDepot[i].setStartMax(-1);

		name = "S_" + to_string(i);
		SrcDepot[i].setName(name.c_str());
		name = "";
	}

	/* Intervals representing the destination depot of agents */

	vector<int> uDepot = ga->depots;
	sort(uDepot.begin(), uDepot.end());
	unique(uDepot.begin(), uDepot.end());

	//for (int i = 0; i < delta; i++) {
	//	IloIntervalVar tmp(env);
	//	DestDepot[i] = tmp;
	//	DestDepot[i].setSizeMin(1); DestDepot[i].setSizeMax(1);
	//	
	//	DestDepot[i].setStartMin(ga->depotStart[i]);
	//	name = "D_" + to_string(i);
	//	DestDepot[i].setName(name.c_str());
	//	name = "";
	//}

	/* Intervals representing the destination depot of agents */
	for (int i = 0; i < sigma; i++) {
		IloIntervalVar tmp(env);
		D[i] = tmp;
		D[i].setSizeMin(1);
		D[i].setSizeMax(1);
		//D[i].setStartMax(ga->depotStart[i]);

		name = "D_" + to_string(i);
		D[i].setName(name.c_str());
		name = "";
	}


	//for (int i = 0; i < sigma; i++) {
	//	DestDepot2[i] = IloIntervalVarArray(env, delta);
	//	for (int j = 0; j < delta; j++) {
	//		IloIntervalVar tmp(env);
	//		DestDepot2[i][j] = tmp;
	//		DestDepot2[i][j].setSizeMin(1); DestDepot2[i][j].setSizeMax(1);
	//		//DestDepot2[i][j].setOptional();
	//		name = "DestDepot_" + to_string(i) + to_string(j);
	//		DestDepot2[i][j].setName(name.c_str());
	//		name = "";
	//	}
	//}

	//for (int i = 0; i < sigma; i++) {	
	//		DestDepot2[i][ga->depots[i]].setPresent();
	//		DestDepot2[i][ga->depots[i]].setStartMin(ga->depotStart[i]);	
	//}

	/* Create Tasks and set their duration */
	for (int i = 0; i < V; i++) {
		IloIntervalVar tmp(env);
		Tasks[i] = tmp;
	}
	for (int i = 0; i < V; i++) {

		Tasks[i].setSizeMin(fullModel->getTaskDuration(i + sigma));
		Tasks[i].setSizeMax(fullModel->getTaskDuration(i + sigma));
		//cout << fullModel->getTaskDuration(i + sigma) << "\t" << endl;
		if (i < Vp) {
			auto it = find(pcc.begin(), pcc.end(), i + fullModel->sigma);
			if (it == pcc.end()) {
				int ind = fullModel->_model.T[i].index;
				//Tasks[ind].setStartMin(ga->timeline[i].startTime);
			}
		}

	}

	for (int i = 0; i < V; i++) {
		TasksV[i] = IloIntervalVarArray(env, _fixedAgents2Tasks[i].size());

		for (int j = 0; j < _fixedAgents2Tasks[i].size(); j++) {
			IloIntervalVar tmp(env);
			TasksV[i][j] = tmp;
			TasksV[i][j].setPresent(); // OVO NE TREBA OVAKO, valjda treba isPresent!?

			name = "T" + to_string(i) + "A" + to_string(_fixedAgents2Tasks[i][j]);
			TasksV[i][j].setName(name.c_str());
			name = "";
		}
	}

	//for (int i = 0; i < V; i++) {
	//	TasksV[i] = IloIntervalVarArray(env, _fixedAgents2Tasks[i].size());

	//	for (int j = 0; j < _fixedAgents2Tasks[i].size(); j++) {
	//		IloIntervalVar tmp(env);
	//		TasksV[i][j] = tmp;
	//		TasksV[i][j].setOptional(); // OVO NE TREBA OVAKO, valjda treba isPresent!?

	//		name = "T" + to_string(i) + "A" + to_string(_fixedAgents2Tasks[i][j]);
	//		TasksV[i][j].setName(name.c_str());
	//		name = "";
	//	}
	//}

	//for (int i = 0; i < V; i++) {
	//	TasksV[i] = IloIntervalVarArray(env, delta);

	//	for (int j = 0; j < delta; j++) {
	//		//cout << "TaskV" << i - sigma << "_" << j << endl;
	//		IloIntervalVar tmp(env);
	//		TasksV[i][j] = tmp;
	//		
	//		if (j == _fixedAgents2Tasks[i][0]) {
	//			TasksV[i][j].setPresent();
	//		}
	//		else {
	//			TasksV[i][j].setAbsent();
	//		}
	//		

	//		
	//		//TasksV[i][j].setOptional();
	//		//TasksV[i][j].setSizeMin(Tasks[i].getSizeMax()); TasksV[i][j].setSizeMax(Tasks[i].getSizeMax());
	//		name = "T" + to_string(i) + "A" + to_string(_fixedAgents2Tasks[i][j]);
	//		TasksV[i][j].setName(name.c_str());
	//		name = "";
	//	}
	//}

	/* Transition Distances - source tasks destination. It is symmetric*/

	vector<IloTransitionDistance> vM(sigma);
	for (int s = 0; s < sigma; s++) {
		int ii = 0, jj = 0;
		IloTransitionDistance M(env, POS.getSize());
		vM[s] = M;
		for (int i = 0; i < Vtilde; i++) {
			if (fullModel->H[i] == 0) {
				for (int j = 0; j < Vtilde; j++) {
					if (fullModel->H[j] == 0) {
						if (ii != jj) { M.setValue(ii, jj, fullModel->getEdgeCost(i, j, s)); }
						else { M.setValue(ii, jj, 0); }
						jj++;
					}
				}
				jj = 0; ii++;
			}
		}
	}

	/*IloStateFunction - State function for each agent based on transition distances */
	IloStateFunctionArray POSV(env);

	for (int i = 0; i < sigma; i++) {
		IloStateFunction tmpPOSV(env, vM[i]);
		POSV.add(tmpPOSV);
	}

	/*IloEndBeforeStart - Precedence Constraints */
	for (int i = sigma; i < fullModel->_PV.size() - delta; i++) {
		for (int j = 0; j < fullModel->_PV[i].size(); j++) {
			if (fullModel->_PV[i][j] != 0) {
				//cout << i - sigma << "\t " << fullModel->_PV[i][j] - sigma << endl;
				model.add(IloEndBeforeStart(env, Tasks[i - sigma], Tasks[fullModel->_PV[i][j] - sigma]));
			}
		}

	}

	IloIntervalVarArray tmpNoOverlap(env);
	vector<int> _v;

	for (int i = sigma; i < V + sigma; i++) {
		for (int j = i + 1; j < V + sigma; j++) {

			auto it = set_intersection(_fixedAgents2Tasks[i - sigma].begin(), _fixedAgents2Tasks[i - sigma].end(), _fixedAgents2Tasks[j - sigma].begin(), _fixedAgents2Tasks[j - sigma].end(), back_inserter(_v));

			for (int s = 0; s < _v.size(); s++) {

				if (fullModel->R[i][j] != 1) {
					auto it = find(_fixedAgents2Tasks[i - sigma].begin(), _fixedAgents2Tasks[i - sigma].end(), _v[s]);
					int a = distance(_fixedAgents2Tasks[i - sigma].begin(), it);

					tmpNoOverlap.add(TasksV[i - sigma][a]);

					it = find(_fixedAgents2Tasks[j - sigma].begin(), _fixedAgents2Tasks[j - sigma].end(), _v[s]);
					int a2 = distance(_fixedAgents2Tasks[j - sigma].begin(), it);

					tmpNoOverlap.add(TasksV[j - sigma][a2]);

					//cout << "(TasksV[" << i - sigma << "][" << a << "], TasksV[" << j - sigma << "][" << a2 << "])" << endl;
				}
				model.add(IloNoOverlap(env, tmpNoOverlap));
				tmpNoOverlap.clear();
			}
			_v.clear();
		}
	}

	for (int i = 0; i < Tasks.getSize(); i++) {
		IloIntervalVarArray tmpA(env);
		for (int j = 0; j < _fixedAgents2Tasks[i].size(); j++) {
			tmpA.add(TasksV[i][j]);
			name = "A" + to_string(i);
			tmpA[j].setName(name.c_str());
			name = "";
		}
		A.add(tmpA);
	}

	//for (int i = 0; i < Tasks.getSize(); i++) {
	//	IloIntervalVarArray tmpA(env);
	//	for (int j = 0; j < delta; j++) {
	//		tmpA.add(TasksV[i][j]);
	//		name = "A" + to_string(i);
	//		tmpA[j].setName(name.c_str());
	//		name = "";
	//	}
	//	A.add(tmpA);
	//}

	/// Enforce constraints based on A
	for (int i = 0; i < Tasks.getSize(); i++) {
		model.add(IloAlternative(env, Tasks[i], A[i]));
	}

	/// Enforce constraints based on DestDepot2
	//for (int i = 0; i < sigma; i++) {
	//	model.add(IloAlternative(env, D[i], DestDepot2[i]));
	//}

	/* IloEndBeforeStart - All tasks must finish before the agent reach destination depot */
	for (int i = 0; i < TasksV.getSize(); i++) {
		for (int j = 0; j < TasksV[i].getSize(); j++) {

			//cout << "TasksV[" << i << "][" << j << "]  <---  D[" << fullModel->_vectorTasksPerAgents[i][j] << "]" << "\t" << fullModel->_vectorTasksPerAgents[i][j] << endl;

			//model.add(IloEndBeforeStart(env, TasksV[i][j], D[fullModel->_vectorTasksPerAgents[i][j]]));
			model.add(IloEndBeforeStart(env, TasksV[i][j], D[_fixedAgents2Tasks[i][0]]));

		}
	}

	//for (int i = 0; i < TasksV.getSize(); i++) {
	//	for (int j = 0; j < TasksV[i].getSize(); j++) {
	//		for (int h = 0; h < delta; h++) {
	//			//cout << "TasksV[" << i << "][" << j << "]  <---  D[" << j << "][" << h << "]" << "\t" << graph._agentsPerTasks[i + sigma][j] << endl;
	//			model.add(IloEndBeforeStart(env, TasksV[i][j], DestDepot2[j][h]));

	//		}
	//	}
	//}

	/* IloAlwaysEqual - constraints position of agents */
	for (int i = 0; i < POSV.getSize(); i++) {
		/* Source and Destination Depot */
		model.add(IloAlwaysEqual(env, POSV[i], SrcDepot[i], POS[i]));

		/* Physical Tasks */
		for (int j = sigma; j < index.size() - delta; j++) {
			for (int k = 0; k < _fixedAgents2Tasks[index[j] - sigma].size(); k++) {
				if (_fixedAgents2Tasks[index[j] - sigma][k] == i) {
					//cout << i << " " << index[j] - sigma << " " << k << " " << j << endl;
					model.add(IloAlwaysEqual(env, POSV[i], TasksV[index[j] - sigma][k], j));
				}
			}
		}
		//	for (int j = sigma; j < index.size() - delta; j++) {
	//		for (int k = 0; k < _fixedAgents2Tasks[index[j] - sigma].size(); k++) {
	//			if (_fixedAgents2Tasks[index[j] - sigma][k] == i) {
	//				//cout << i << " " << index[j] - sigma << " " << k << " " << j << endl;
	//				model.add(IloAlwaysEqual(env, POSV[i], TasksV[index[j] - sigma][k], j));
	//			}

	//		}
	//	}
	}

	for (int i = 0; i < sigma; i++) {

		model.add(IloAlwaysEqual(env, POSV[i], D[i], POS.getSize() - delta + ga->depots[i]));

	}

	/*IloExpr expr(env);
	for (int i = 0; i < D.getSize(); i++) {
		expr += IloStartOf(D[i]);
	}
	model.add(IloMinimize(env, expr));
*/
	IloIntExprArray endTimes(env);

	for (int i = 0; i < sigma; i++) {
		endTimes.add(IloStartOf(D[i]));

	}

	model.add(IloMinimize(env, IloMax(endTimes)));



	///* IloEndBeforeStart - All tasks must finish before the agent reach destination depot */
	//for (int i = 0; i < TasksV.getSize(); i++) {
	//	for (int j = 0; j < TasksV[i].getSize(); j++) {
	//		for (int k = 0; k < delta; k++) {
	//			model.add(IloEndBeforeStart(env, TasksV[i][j], DestDepot[k]));
	//		}
	//	}
	//}

	//for (int i = 0; i < POSV.getSize(); i++) {
	//	/* Source and Destination Depot */
	//	model.add(IloAlwaysEqual(env, POSV[i], SrcDepot[i], POS[i]));
	//	//cout << "i: " << i << " POSV.size() " << POSV.getSize() << " POS.size() " << POS.getSize() << " delta " << delta << endl;
	//	model.add(IloAlwaysEqual(env, POSV[i], DestDepot[i], POS.getSize() - delta + i));

	//	/* Physical Tasks */

	//	for (int j = sigma; j < index.size() - delta; j++) {
	//		//if (j >= sigma && j < V + sigma) {
	//		for (int k = 0; k < _fixedAgents2Tasks[index[j] - sigma].size(); k++) {
	//			if (_fixedAgents2Tasks[index[j] - sigma][k] == i) {
	//				//cout << i << " " << index[j] - sigma << " " << k << " " << j << endl;
	//				model.add(IloAlwaysEqual(env, POSV[i], TasksV[index[j] - sigma][k], j));
	//			}
	//			//}
	//		}
	//	}
	//}

	////model.add(IloMinimize(env, IloMax(IloStartOf(DestDepot[0]), IloStartOf(DestDepot[1]))));
	//IloExpr expr(env);
	//for (int i = 0; i < DestDepot.getSize(); i++) {
	//	expr += IloStartOf(DestDepot[i]);
	//}


	IloCP cp(model);
	cp.exportModel("model.cpo");
	cp.setParameter(IloCP::LogSearchTags, IloCP::On);
	cp.setParameter(IloCP::TimeLimit, timeLimitPartial);
	cp.setParameter(IloCP::Presolve, IloCP::On);
	cp.setParameter(IloCP::LogVerbosity, IloCP::Terse);
	IloTimer timer(env);
	timer.start();
	string pathbb = path + to_string(run) + "_incumbentbb_hybrid.txt";
	string pathobj = path + to_string(run) + "_incumbentobj_hybrid.txt";
	IncumbentSolutionCallback cb(cp.out(), D, fullModel, pathbb, pathobj, timer);
	cp.addCallback(&cb);


	bool sol = false;

	try {
		sol = cp.solve();
	}
	catch (const IloException& e) {
		std::cerr << std::endl << std::endl;
		std::cerr << "CP Raised an exception:" << std::endl;
		std::cerr << e << std::endl;
		env.end();
		throw;
	}

	if (sol) {
		cp.out() << "Makespan \t: " << cp.getObjValue() << std::endl;


		for (int i = 0; i < Tasks.getSize(); i++) {
			cp.out() << setw(5) << "Task " << setw(3) << i << ": start time: " << setw(5) << cp.getStart(Tasks[i]) << " | end time: " << setw(5) << cp.getEnd(Tasks[i]) << endl;
		}
		cout << endl;

		set<int> usedAgents;

		//remove("c:\\Users\\bmc01\\Documents\\Visual Studio 2017\\Projects\\MR_MT_Solver\\MR_MT_Solver\\output_hybrid.txt");
		//ofstream myfile3("c:\\Users\\bmc01\\Documents\\Visual Studio 2017\\Projects\\MR_MT_Solver\\MR_MT_Solver\\output_hybrid.txt", ios_base::app | ios_base::out);
		string path2 = path + to_string(run) + "_output_hybrid.txt";
		ofstream myfile3(path2, ios_base::app | ios_base::out);

		for (int i = 0; i < Tasks.getSize(); i++) {
			myfile3 << cp.getStart(Tasks[i]) << "\t" << cp.getEnd(Tasks[i]) << "\t" << flush;
			//cout << cp.getStart(Tasks[i]) << "\t" << cp.getEnd(Tasks[i]) << "\t" << flush;
			for (int j = 0; j < _fixedAgents2Tasks[i].size(); j++) {

				if (cp.isPresent(TasksV[i][j])) {
					usedAgents.insert(_fixedAgents2Tasks[i][j]);
					myfile3 << _fixedAgents2Tasks[i][j] << "\t" << fullModel->H[i + sigma] << endl;
					//cout << _fixedAgents2Tasks[i][j] << "\t V" << fullModel->H[i + sigma] << endl;
				}


			}
		}
		myfile3.close();


		//remove("c:\\Users\\bmc01\\Documents\\Visual Studio 2017\\Projects\\MR_MT_Solver\\MR_MT_Solver\\depot2_hybrid.txt");
		//ofstream myfile4("c:\\Users\\bmc01\\Documents\\Visual Studio 2017\\Projects\\MR_MT_Solver\\MR_MT_Solver\\depot2_hybrid.txt", ios_base::app | ios_base::out);

		string path3 = path + to_string(run) + "_depot2_hybrid.txt";
		ofstream myfile4(path3, ios_base::app | ios_base::out);

		int deductable = 0;

		for (int i = 0; i < sigma; i++) {

			if (cp.isPresent(D[i])) {
				int dt = cp.getStart(D[i]);
				auto it = find(usedAgents.begin(), usedAgents.end(), i);
				if (it == usedAgents.end()) {
					deductable += dt;
					dt = 0;
				}

				myfile4 << dt << endl;
				cout << "D[" << to_string(i) << "] " << dt << "\t" << flush;
			}

		}
		myfile4 << cp.getStatus() << endl;
		myfile4 << cp.getObjGap() << endl;


		cout << endl; cout << endl;
		myfile4.close();

		//cout << "Final obj. value: " << cp.getObjValue() - deductable << endl;
		cout << "Final obj. value: " << cp.getObjValue() << endl;
	}
	else { cp.out() << "No solution found." << std::endl; }
	env.end();
}

void CPScheduler::solverFull()
{

	/* CP MODEL */
	int sigma = fullModel->sigma;
	int delta = fullModel->delta;
	int V = fullModel->V;
	int Vtilde = V + sigma + delta;

	IloEnv env;
	IloModel model(env);
	string name;

	IloIntervalVar tmp(env);
	IloIntervalVarArray Tasks(env, V), D(env, sigma), SrcDepot(env, sigma);
	IloIntervalVarArray2 A(env), TasksV(env, V), noOverlapArray(env), DestDepot2(env, sigma);
	IloIntArray POS(env);

	vector<int> index;
	int count = 0;

	/* Create array of Physical Tasks */
	for (int i = 0; i < Vtilde; i++) {
		if (i < sigma || i >= V + sigma) {
			POS.add(count++);
			index.emplace_back(i);
		}
		else if (fullModel->H[i] != 1) {
			POS.add(count++);
			index.emplace_back(i);
		}
	}

	/* Intervals representing the source depot of agents */
	for (int i = 0; i < sigma; i++) {
		IloIntervalVar tmp(env);
		SrcDepot[i] = tmp;
		SrcDepot[i].setSizeMin(1); SrcDepot[i].setSizeMax(1);
		//SrcDepot[i].setLengthMin(1); SrcDepot[i].setLengthMax(1);
		SrcDepot[i].setStartMin(-1); SrcDepot[i].setStartMax(-1);
		//SrcDepot[i].setEndMin(0); SrcDepot[i].setEndMax(0);
		name = "S_" + to_string(i);
		SrcDepot[i].setName(name.c_str());
		name = "";
	}

	for (int i = 0; i < sigma; i++) {
		DestDepot2[i] = IloIntervalVarArray(env, delta);
		for (int j = 0; j < delta; j++) {
			IloIntervalVar tmp(env);
			DestDepot2[i][j] = tmp;
			DestDepot2[i][j].setSizeMin(1); DestDepot2[i][j].setSizeMax(1);
			//DestDepot2[i][j].setLengthMin(1); DestDepot2[i][j].setLengthMax(1);
			DestDepot2[i][j].setOptional();
			name = "DestDepot_" + to_string(i) + to_string(j);
			DestDepot2[i][j].setName(name.c_str());
			name = "";
		}
	}
	/* Intervals representing the destination depot of agents */
	for (int i = 0; i < sigma; i++) {
		IloIntervalVar tmp(env);
		D[i] = tmp;
		D[i].setSizeMin(1);
		D[i].setSizeMax(1);
		//D[i].setLengthMin(1);
		//D[i].setLengthMax(1);

		name = "D_" + to_string(i);
		D[i].setName(name.c_str());
		name = "";
	}

	/* Create Tasks and set their duration */
	for (int i = 0; i < V; i++) {
		IloIntervalVar tmp(env);
		Tasks[i] = tmp;
		Tasks[i].setSizeMin(fullModel->getTaskDuration(i + sigma)); Tasks[i].setSizeMax(fullModel->getTaskDuration(i + sigma));
		//Tasks[i].setLengthMin(fullModel->getTaskDuration(i + sigma)); Tasks[i].setLengthMax(fullModel->getTaskDuration(i + sigma));
	}

	/* Create double array of possible allocation of tasks to agents */
	for (int i = 0; i < V; i++) {
		TasksV[i] = IloIntervalVarArray(env, fullModel->_vectorTasksPerAgents[i].size());

		for (int j = 0; j < fullModel->_vectorTasksPerAgents[i].size(); j++) {
			//cout << "TaskV" << i - sigma << "_" << j << endl;
			IloIntervalVar tmp(env);
			TasksV[i][j] = tmp;
			TasksV[i][j].setOptional();
			//TasksV[i][j].setSizeMin(Tasks[i].getSizeMax()); TasksV[i][j].setSizeMax(Tasks[i].getSizeMax());
			name = "T" + to_string(i) + "A" + to_string(fullModel->_vectorTasksPerAgents[i][j]);
			TasksV[i][j].setName(name.c_str());
			name = "";
		}
	}

	/* Transition Distances - source tasks destination. It is symmetric*/
	vector<IloTransitionDistance> vM(sigma);
	for (int s = 0; s < sigma; s++) {
		int ii = 0, jj = 0;
		IloTransitionDistance M(env, POS.getSize());
		vM[s] = M;

		for (int i = 0; i < Vtilde; i++) {
			if (fullModel->H[i] == 0) {
				for (int j = 0; j < Vtilde; j++) {
					if (fullModel->H[j] == 0) {
						if (ii == jj) { M.setValue(ii, jj, 0); }
						else { M.setValue(ii, jj, fullModel->getEdgeCost(i, j, s)); }
						jj++;
					}
				}
				jj = 0; ii++;
			}
		}

		//IloTransitionDistance M(env, POS.getSize() - delta);
		//vM[s] = M;
		//for (int i = 0; i < V + sigma; i++) {
		//	if (fullModel->H[i] == 0) {
		//		for (int j = 0; j < V + sigma; j++) {
		//			if (fullModel->H[j] == 0) {
		//				if (ii == jj) { M.setValue(ii, jj, 0); }
		//				else { M.setValue(ii, jj, fullModel->getEdgeCost(i, j, s)); }
		//				jj++;
		//			}
		//		}
		//		jj = 0; ii++;
		//	}
		//}
	}

	//for (int s = 0; s < vM.size(); s++) {
	//	for (int i = 0; i < vM[s].getSize(); i++) {
	//		for (int j = 0; j < vM[s].getSize(); j++) {

	//			cout << vM[s].getValue(i, j) << "\t" << flush;


	//		}
	//		cout << endl;
	//	}
	//	cout << "----------------------------  " << s << endl;
	//}
	//IloTransitionDistance mm(env, 12);
	//int ii = 0, jj = 0;
	//for (int i = 0; i < 1; i++) {

	//	for (int j = 0; j < Vtilde; j++) {

	//		cout << i << "  " << j << endl;
	//		if (ii == jj) { mm.setValue(ii, jj, 0); }
	//		else { mm.setValue(ii, jj, fullModel->getEdgeCost(i, j, 0)); }
	//		jj++;
	//	}
	//}
	//jj = 0; ii++;

	//IloStateFunctionArray POSd(env);
	//IloStateFunction tmpPOSd(env, mm);
	//POSd.add(tmpPOSd);


	/*IloStateFunction - State function for each agent based on transition distances */
	IloStateFunctionArray POSV(env);

	for (int i = 0; i < sigma; i++) {
		IloStateFunction tmpPOSV(env, vM[i]);
		POSV.add(tmpPOSV);
	}

	/*IloEndBeforeStart - Precedence Constraints */
	for (int i = sigma; i < fullModel->_PV.size() - delta; i++) {
		for (int j = 0; j < fullModel->_PV[i].size(); j++) {
			if (fullModel->_PV[i][j] != 0) {
				model.add(IloEndBeforeStart(env, Tasks[i - sigma], Tasks[fullModel->_PV[i][j] - sigma]));
			}
		}

	}

	/* This is neccesary for ECTSP, otherwise comment it out */
	/*vector<int> _v;

	for (int i = sigma; i < fullModel->_PV.size() - delta; i++) {
		for (int j = 0; j < fullModel->_PV[i].size(); j++) {
			if (fullModel->_PV[i][j] != 0) {
				auto it = set_intersection(fullModel->_vectorTasksPerAgents[i - sigma].begin(), fullModel->_vectorTasksPerAgents[i - sigma].end(), fullModel->_vectorTasksPerAgents[fullModel->_PV[i][j] - sigma].begin(), fullModel->_vectorTasksPerAgents[fullModel->_PV[i][j] - sigma].end(), back_inserter(_v));
				if (_v.empty()) {
					cout << " Error: there is no agent that can do both tasks " << i - sigma << " and " << fullModel->_PV[i][j] - sigma << endl;
				}
				else {
					for (int s = 0; s < _v.size(); s++) {
						auto it2 = find(fullModel->_vectorTasksPerAgents[i - sigma].begin(), fullModel->_vectorTasksPerAgents[i - sigma].end(), _v[s]);
						int s1 = distance(fullModel->_vectorTasksPerAgents[i - sigma].begin(), it2);

						auto it3 = find(fullModel->_vectorTasksPerAgents[fullModel->_PV[i][j] - sigma].begin(), fullModel->_vectorTasksPerAgents[fullModel->_PV[i][j] - sigma].end(), _v[s]);
						int s2 = distance(fullModel->_vectorTasksPerAgents[fullModel->_PV[i][j] - sigma].begin(), it3);

						model.add(IloPresenceOf(env, TasksV[i - sigma][s1]) == IloPresenceOf(env, TasksV[fullModel->_PV[i][j] - sigma][s2]));
					}
					_v.clear();
				}
			}
		}
	}*/


	/*-------------------------------------------------------*/

	///*IloNoOverlap - Prevents tasks from overlapping, except the ones that are allowed by matrix R */
	IloIntervalVarArray tmpNoOverlap(env);
	vector<int> _v;
	_v.clear();
	for (int i = sigma; i < V + sigma; i++) {
		for (int j = i + 1; j < V + sigma; j++) {

			auto it = set_intersection(fullModel->_vectorTasksPerAgents[i - sigma].begin(), fullModel->_vectorTasksPerAgents[i - sigma].end(), fullModel->_vectorTasksPerAgents[j - sigma].begin(), fullModel->_vectorTasksPerAgents[j - sigma].end(), back_inserter(_v));

			for (int s = 0; s < _v.size(); s++) {

				if (fullModel->R[i][j] != 1) {
					auto it = find(fullModel->_vectorTasksPerAgents[i - sigma].begin(), fullModel->_vectorTasksPerAgents[i - sigma].end(), _v[s]);
					int a = distance(fullModel->_vectorTasksPerAgents[i - sigma].begin(), it);

					tmpNoOverlap.add(TasksV[i - sigma][a]);

					it = find(fullModel->_vectorTasksPerAgents[j - sigma].begin(), fullModel->_vectorTasksPerAgents[j - sigma].end(), _v[s]);
					int a2 = distance(fullModel->_vectorTasksPerAgents[j - sigma].begin(), it);

					tmpNoOverlap.add(TasksV[j - sigma][a2]);

					//cout << "(TasksV[" << i - sigma << "][" << a << "], TasksV[" << j - sigma << "][" << a2 << "])" << endl;
				}
				model.add(IloNoOverlap(env, tmpNoOverlap));
				tmpNoOverlap.clear();
			}
			_v.clear();
		}
	}

	///* IloAlternative - Allocation of tasks to possible agents */
		/// Create IloIntervalVarArray2 A consisting of all possible allocations
	for (int i = 0; i < Tasks.getSize(); i++) {
		IloIntervalVarArray tmpA(env);
		for (int j = 0; j < fullModel->_vectorTasksPerAgents[i].size(); j++) {
			tmpA.add(TasksV[i][j]);
			name = "A" + to_string(i);
			tmpA[j].setName(name.c_str());
			name = "";
		}
		A.add(tmpA);
	}

	/// Enforce constraints based on A
	for (int i = 0; i < Tasks.getSize(); i++) {
		model.add(IloAlternative(env, Tasks[i], A[i]));
	}

	/// Enforce constraints based on DestDepot2
	for (int i = 0; i < sigma; i++) {
		model.add(IloAlternative(env, D[i], DestDepot2[i]));
	}

	/* IloEndBeforeStart - All tasks must finish before the agent reach destination depot */
	for (int i = 0; i < TasksV.getSize(); i++) {
		for (int j = 0; j < TasksV[i].getSize(); j++) {
			for (int h = 0; h < delta; h++) {
				//cout << "TasksV[" << i << "][" << j << "]  <---  D[" << fullModel->_vectorTasksPerAgents[i][j] << "][" << h << "]" << "\t" << fullModel->_vectorTasksPerAgents[i][j] << endl;
				model.add(IloEndBeforeStart(env, TasksV[i][j], DestDepot2[fullModel->_vectorTasksPerAgents[i][j]][h]));

			}
		}
	}

	/* IloAlwaysEqual - constraints position of agents */
	//for (int i = 0; i < POSV.getSize(); i++) {
	//	/* Source and Destination Depot */
	//	model.add(IloAlwaysEqual(env, POSV[i], SrcDepot[i], POS[i]));

	//	/* Physical Tasks */
	//	for (int j = sigma; j < index.size() - delta; j++) {
	//		for (int k = 0; k < fullModel->_vectorTasksPerAgents[index[j] - sigma].size(); k++) {
	//			if (fullModel->_vectorTasksPerAgents[index[j] - sigma][k] == i) {
	//				cout << i << " " << index[j] - sigma << " " << k << " " << j << endl;
	//				model.add(IloAlwaysEqual(env, POSV[i], TasksV[index[j] - sigma][k], j));
	//			}
	//		}
	//	}
	//}

		/* IloAlwaysEqual - constraints position of agents */
	for (int i = 0; i < POSV.getSize(); i++) {
		/* Source Depot */
		model.add(IloAlwaysEqual(env, POSV[i], SrcDepot[i], POS[i]));
	}

	/* Physical Tasks */
	for (int j = sigma; j < index.size() - delta; j++) {
		for (int k = 0; k < fullModel->_vectorTasksPerAgents[index[j] - sigma].size(); k++) {
			//if (fullModel->_vectorTasksPerAgents[index[j] - sigma][k] == i) {
				//cout << k << " " << index[j] - sigma << " " << k << " " << j << endl;
			model.add(IloAlwaysEqual(env, POSV[k], TasksV[index[j] - sigma][k], j));
		}
	}
	// Destination Depot
	for (int i = 0; i < sigma; i++) {
		for (int j = 0; j < delta; j++) {
			model.add(IloAlwaysEqual(env, POSV[i], DestDepot2[i][j], POS.getSize() - delta + j));
		}
	}

	//model.add(IloAlwaysEqual(env,POSd[0], DestDepot2[0][0], 0));

	IloExpr expr(env);
	for (int i = 0; i < D.getSize(); i++) {
		expr += IloStartOf(D[i]);
	}


	IloIntExprArray endTimes(env);

	for (int i = 0; i < sigma; i++) {
		endTimes.add(IloStartOf(D[i]));
	}

	/* sum */
	model.add(IloMinimize(env, expr));

	/* minMax */
	//model.add(IloMinimize(env, IloMax(endTimes)));

	/* Weighted sum + minMax*/
	//model.add(IloMinimize(env, IloMax(endTimes) + 0.1 * expr));


	IloCP cp(model);
	cp.exportModel("model.cpo");
	cp.setParameter(IloCP::LogSearchTags, IloCP::On);
	cp.setParameter(IloCP::TimeLimit, timeLimit);
	cp.setParameter(IloCP::Presolve, IloCP::On);
	cp.setParameter(IloCP::LogVerbosity, IloCP::Terse);
	cp.setParameter(IloCP::TimeMode, IloCP::ElapsedTime);
	//cp.setParameter(IloCP::DefaultInferenceLevel, IloCP::Extended);
	//cp.setParameter(IloCP::NoOverlapInferenceLevel, IloCP::Extended);
	//cp.setParameter(IloCP::PrecedenceInferenceLevel, IloCP::Extended);
	//cp.setParameter(IloCP::IntervalSequenceInferenceLevel, IloCP::Extended);
	//cp.setParameter(IloCP::StateFunctionInferenceLevel, IloCP::Extended);
	IloTimer timer(env);

	string pathbb = path + to_string(run) + "_incumbentbb.txt";
	string pathobj = path + to_string(run) + "_incumbentobj.txt";
	IncumbentSolutionCallback cb(cp.out(), D, fullModel, pathbb, pathobj, timer);
	cp.addCallback(&cb);

	bool sol = false;

	try {
		timer.start();
		sol = cp.solve();
		cout << timer.getTime() << endl;
	}
	catch (const IloException& e) {
		std::cerr << std::endl << std::endl;
		std::cerr << "CP Raised an exception:" << std::endl;
		std::cerr << e << std::endl;
		env.end();
		throw;
	}

	if (sol) {
		cp.out() << "Makespan \t: " << cp.getObjValue() << std::endl;

		for (int i = 0; i < Tasks.getSize(); i++) {
			cp.out() << "Task " << i + 1 << ": start time: " << cp.getStart(Tasks[i]) << " | end time: " << cp.getEnd(Tasks[i]) << endl;
		}
		set<int> usedAgents;

		//remove("c:\\Users\\bmc01\\Documents\\Visual Studio 2017\\Projects\\MR_MT_Solver\\MR_MT_Solver\\output.txt");
		//ofstream myfile2("output.txt", ios_base::app | ios_base::out);

		string path2 = path + to_string(run) + "_output.txt";
		ofstream myfile2(path2, ios_base::app | ios_base::out);
		vector<vector<int>> outputPlan(sigma);

		for (int i = 0; i < Tasks.getSize(); i++) {
			myfile2 << cp.getStart(Tasks[i]) << "\t" << cp.getEnd(Tasks[i]) << "\t" << flush;
			//cout << cp.getStart(Tasks[i]) << "\t" << cp.getEnd(Tasks[i]) << "\t" << flush;
			for (int j = 0; j < fullModel->_vectorTasksPerAgents[i].size(); j++) {

				if (cp.isPresent(TasksV[i][j])) {
					usedAgents.insert(fullModel->_vectorTasksPerAgents[i][j]);
					outputPlan[fullModel->_vectorTasksPerAgents[i][j]].push_back(i);
					myfile2 << fullModel->_vectorTasksPerAgents[i][j] << "\t" << fullModel->H[i + sigma] << endl;
					cout << fullModel->_vectorTasksPerAgents[i][j] << "\t" << fullModel->H[i + sigma] << endl;
				}
			}
		}
		myfile2.close();


		//remove("c:\\Users\\bmc01\\Documents\\Visual Studio 2017\\Projects\\MR_MT_Solver\\MR_MT_Solver\\depot2.txt");
		//ofstream myfile3("depot2.txt", ios_base::app | ios_base::out);

		path2 = path + to_string(run) + "_depot2.txt";
		ofstream myfile3(path2, ios_base::app | ios_base::out);

		int deductable = 0;
		for (int i = 0; i < sigma; i++) {
			for (int j = 0; j < delta; j++) {
				if (cp.isPresent(DestDepot2[i][j])) {

					int dt = cp.getStart(DestDepot2[i][j]);
					auto it = find(usedAgents.begin(), usedAgents.end(), i);
					if (it == usedAgents.end()) {
						deductable += dt;
						dt = 0;
					}

					myfile3 << dt << endl;
					cout << "Agent " << i << " -> Depot " << j << " D[" << to_string(i) << "] " << dt << "\t" << flush;
				}
			}
		}
		myfile3 << cp.getStatus() << endl;
		myfile3 << cp.getObjGap() << endl;

		myfile3.close();
		cout << endl; cout << endl;
		//cout << "Final obj. value: " << cp.getObjValue() - deductable << endl;
		cout << "Final obj. value: " << cp.getObjValue() << endl;
	}
	else { cp.out() << "No solution found." << std::endl; }
	cout << endl;
	env.end();
}

void CPScheduler::solverFull_warmStart()
{

	/* CP MODEL */
	int sigma = fullModel->sigma;
	int delta = fullModel->delta;
	int V = fullModel->V;
	int Vtilde = V + sigma + delta;

	IloEnv env;
	IloModel model(env);
	string name;

	IloIntervalVar tmp(env);
	IloIntervalVarArray Tasks(env, V), D(env, sigma), SrcDepot(env, sigma);
	IloIntervalVarArray2 A(env, V), TasksV(env, V), noOverlapArray(env), DestDepot2(env, sigma);
	IloIntArray POS(env);

	vector<int> index;
	int count = 0;

	/* Create array of Physical Tasks */
	for (int i = 0; i < Vtilde; i++) {
		if (i < sigma || i >= V + sigma) {
			POS.add(count++);
			index.emplace_back(i);
		}
		else if (fullModel->H[i] != 1) {
			POS.add(count++);
			index.emplace_back(i);
		}
	}

	/* Intervals representing the source depot of agents */
	for (int i = 0; i < sigma; i++) {
		IloIntervalVar tmp(env);
		SrcDepot[i] = tmp;
		SrcDepot[i].setSizeMin(1); SrcDepot[i].setSizeMax(1);
		//SrcDepot[i].setLengthMin(1); SrcDepot[i].setLengthMax(1);
		SrcDepot[i].setStartMin(-1); SrcDepot[i].setStartMax(-1);
		//SrcDepot[i].setEndMin(0); SrcDepot[i].setEndMax(0);
		name = "S_" + to_string(i);
		SrcDepot[i].setName(name.c_str());
		name = "";
	}

	for (int i = 0; i < sigma; i++) {
		DestDepot2[i] = IloIntervalVarArray(env, delta);
		for (int j = 0; j < delta; j++) {
			IloIntervalVar tmp(env);
			DestDepot2[i][j] = tmp;
			DestDepot2[i][j].setSizeMin(1); DestDepot2[i][j].setSizeMax(1);
			//DestDepot2[i][j].setLengthMin(1); DestDepot2[i][j].setLengthMax(1);
			DestDepot2[i][j].setOptional();
			name = "DestDepot_" + to_string(i) + to_string(j);
			DestDepot2[i][j].setName(name.c_str());
			name = "";
		}
	}
	/* Intervals representing the destination depot of agents */
	for (int i = 0; i < sigma; i++) {
		IloIntervalVar tmp(env);
		D[i] = tmp;
		D[i].setSizeMin(1);
		D[i].setSizeMax(1);
		//D[i].setLengthMin(1);
		//D[i].setLengthMax(1);

		name = "D_" + to_string(i);
		D[i].setName(name.c_str());
		name = "";
	}

	/* Create Tasks and set their duration */
	for (int i = 0; i < V; i++) {
		IloIntervalVar tmp(env);
		Tasks[i] = tmp;
		Tasks[i].setSizeMin(fullModel->getTaskDuration(i + sigma)); Tasks[i].setSizeMax(fullModel->getTaskDuration(i + sigma));
		name = "T" + to_string(i);
		Tasks[i].setName(name.c_str());
		name = "";
		//Tasks[i].setStartMin(ga->timeline[i].startTime);
		//Tasks[i].setStartMax(ga->timeline[i].startTime);
		//Tasks[i].setLengthMin(fullModel->getTaskDuration(i + sigma)); Tasks[i].setLengthMax(fullModel->getTaskDuration(i + sigma));
	}

	/* Create double array of possible allocation of tasks to agents */
	for (int i = 0; i < V; i++) {
		TasksV[i] = IloIntervalVarArray(env, fullModel->_vectorTasksPerAgents[i].size());

		for (int j = 0; j < fullModel->_vectorTasksPerAgents[i].size(); j++) {
			//cout << "TaskV" << i - sigma << "_" << j << endl;
			IloIntervalVar tmp(env);
			TasksV[i][j] = tmp;
			TasksV[i][j].setOptional();
			//TasksV[i][j].setSizeMin(Tasks[i].getSizeMax()); TasksV[i][j].setSizeMax(Tasks[i].getSizeMax());
			name = "T" + to_string(i) + "A" + to_string(j);
			TasksV[i][j].setName(name.c_str());
			name = "";
		}
	}

	/* Transition Distances - source tasks destination. It is symmetric*/
	vector<IloTransitionDistance> vM(sigma);
	for (int s = 0; s < sigma; s++) {
		int ii = 0, jj = 0;
		IloTransitionDistance M(env, POS.getSize());
		vM[s] = M;

		for (int i = 0; i < Vtilde; i++) {
			if (fullModel->H[i] == 0) {
				for (int j = 0; j < Vtilde; j++) {
					if (fullModel->H[j] == 0) {
						if (ii == jj) { M.setValue(ii, jj, 0); }
						else { M.setValue(ii, jj, fullModel->getEdgeCost(i, j, s)); }
						jj++;
					}
				}
				jj = 0; ii++;
			}
		}
	}

	//for (int s = 0; s < vM.size(); s++) {
	//	for (int i = 0; i < vM[s].getSize(); i++) {
	//		for (int j = 0; j < vM[s].getSize(); j++) {

	//			cout << vM[s].getValue(i, j) << "\t" << flush;


	//		}
	//		cout << endl;
	//	}
	//	cout << "----------------------------  " << s << endl;
	//}


	/*IloStateFunction - State function for each agent based on transition distances */
	IloStateFunctionArray POSV(env);

	for (int i = 0; i < sigma; i++) {
		IloStateFunction tmpPOSV(env, vM[i]);
		POSV.add(tmpPOSV);
	}

	/*IloEndBeforeStart - Precedence Constraints */
	for (int i = sigma; i < fullModel->_PV.size() - delta; i++) {
		for (int j = 0; j < fullModel->_PV[i].size(); j++) {
			if (fullModel->_PV[i][j] != 0) {
				model.add(IloEndBeforeStart(env, Tasks[i - sigma], Tasks[fullModel->_PV[i][j] - sigma]));
			}
		}

	}

	///*IloNoOverlap - Prevents tasks from overlapping, except the ones that are allowed by matrix R */
	IloIntervalVarArray tmpNoOverlap(env);
	vector<int> _v;

	for (int i = sigma; i < V + sigma; i++) {
		for (int j = i + 1; j < V + sigma; j++) {

			auto it = set_intersection(fullModel->_vectorTasksPerAgents[i - sigma].begin(), fullModel->_vectorTasksPerAgents[i - sigma].end(), fullModel->_vectorTasksPerAgents[j - sigma].begin(), fullModel->_vectorTasksPerAgents[j - sigma].end(), back_inserter(_v));

			for (int s = 0; s < _v.size(); s++) {

				if (fullModel->R[i][j] != 1) {
					auto it = find(fullModel->_vectorTasksPerAgents[i - sigma].begin(), fullModel->_vectorTasksPerAgents[i - sigma].end(), _v[s]);
					int a = distance(fullModel->_vectorTasksPerAgents[i - sigma].begin(), it);

					tmpNoOverlap.add(TasksV[i - sigma][a]);

					it = find(fullModel->_vectorTasksPerAgents[j - sigma].begin(), fullModel->_vectorTasksPerAgents[j - sigma].end(), _v[s]);
					int a2 = distance(fullModel->_vectorTasksPerAgents[j - sigma].begin(), it);

					tmpNoOverlap.add(TasksV[j - sigma][a2]);

					//cout << "(TasksV[" << i - sigma << "][" << a << "], TasksV[" << j - sigma << "][" << a2 << "])" << endl;
				}
				model.add(IloNoOverlap(env, tmpNoOverlap));
				tmpNoOverlap.clear();
			}
			_v.clear();
		}
	}

	///* IloAlternative - Allocation of tasks to possible agents */
		/// Create IloIntervalVarArray2 A consisting of all possible allocations
	for (int i = 0; i < Tasks.getSize(); i++) {
		IloIntervalVarArray tmpA(env);
		for (int j = 0; j < fullModel->_vectorTasksPerAgents[i].size(); j++) {
			tmpA.add(TasksV[i][j]);
			name = "A_" + to_string(i) + "_" + to_string(j);
			tmpA[j].setName(name.c_str());
			name = "";
		}
		A[i] = tmpA;
	}

	/// Enforce constraints based on A
	for (int i = 0; i < Tasks.getSize(); i++) {
		model.add(IloAlternative(env, Tasks[i], A[i]));
	}

	/// Enforce constraints based on DestDepot2
	for (int i = 0; i < sigma; i++) {
		model.add(IloAlternative(env, D[i], DestDepot2[i]));
	}

	/* IloEndBeforeStart - All tasks must finish before the agent reach destination depot */
	for (int i = 0; i < TasksV.getSize(); i++) {
		for (int j = 0; j < TasksV[i].getSize(); j++) {
			for (int h = 0; h < delta; h++) {
				//cout << "TasksV[" << i << "][" << j << "]  <---  D[" << fullModel->_vectorTasksPerAgents[i][j] << "][" << h << "]" << "\t" << fullModel->_vectorTasksPerAgents[i][j] << endl;
				model.add(IloEndBeforeStart(env, TasksV[i][j], DestDepot2[fullModel->_vectorTasksPerAgents[i][j]][h]));

			}
		}
	}

	/* IloAlwaysEqual - constraints position of agents */
	for (int i = 0; i < POSV.getSize(); i++) {
		/* Source and Destination Depot */
		model.add(IloAlwaysEqual(env, POSV[i], SrcDepot[i], POS[i]));

		/* Physical Tasks */
		for (int j = sigma; j < index.size() - delta; j++) {
			for (int k = 0; k < fullModel->_vectorTasksPerAgents[index[j] - sigma].size(); k++) {
				if (fullModel->_vectorTasksPerAgents[index[j] - sigma][k] == i) {
					//cout << i << " " << index[j] - sigma << " " << k << " " << j << endl;
					model.add(IloAlwaysEqual(env, POSV[i], TasksV[index[j] - sigma][k], j));
				}
			}
		}
	}

	for (int i = 0; i < sigma; i++) {
		for (int j = 0; j < delta; j++) {
			model.add(IloAlwaysEqual(env, POSV[i], DestDepot2[i][j], POS.getSize() - delta + j));
		}
	}

	//IloExpr expr(env);
	//for (int i = 0; i < D.getSize(); i++) {
	//	expr += IloStartOf(D[i]);
	//}
	//model.add(IloMinimize(env, expr));

	IloIntExprArray endTimes(env);

	for (int i = 0; i < sigma; i++) {
		endTimes.add(IloStartOf(D[i]));
	}

	model.add(IloMinimize(env, IloMax(endTimes)));
	IloTimer timer(env);

	IloCP cp(model);

	IloSolution solution(env);

	for (int i = 0; i < V; i++) {
		for (int j = 0; j < fullModel->_vectorTasksPerAgents[i].size(); j++) {
			if (fullModel->_vectorTasksPerAgents[i][j] == ga->timeline[i].agentIdx) {
				//cout << "TaskV" << i << "_" << fullModel->_vectorTasksPerAgents[i][j] << endl;
				solution.setStart(TasksV[i][j], ga->timeline[i].startTime);
				solution.setPresent(TasksV[i][j]);
			}
		}
	}

	//for (int d = 0; d < sigma; d++) {
	//	solution.setStartMin(D[d], ga->depotStart[d]);
	//	solution.setStartMax(D[d], ga->depotStart[d]);
	//}
	cp.setStartingPoint(solution);


	cp.exportModel("model.cpo");
	cp.setParameter(IloCP::LogSearchTags, IloCP::On);
	cp.setParameter(IloCP::TimeLimit, timeLimitPartial);
	cp.setParameter(IloCP::Presolve, IloCP::On);
	cp.setParameter(IloCP::LogVerbosity, IloCP::Terse);


	//cp.setParameter(IloCP::TimeMode, IloCP::CPUTime);
	//cp.setParameter(IloCP::LogVerbosity, IloCP::Verbose);
	string pathbb = path + to_string(run) + "_incumbentbb.txt";
	string pathobj = path + to_string(run) + "_incumbentobj.txt";
	IncumbentSolutionCallback cb(cp.out(), D, fullModel, pathbb, pathobj, timer);
	cp.addCallback(&cb);

	bool sol = false;

	try {
		timer.start();
		sol = cp.solve();
	}
	catch (const IloException& e) {
		std::cerr << std::endl << std::endl;
		std::cerr << "CP Raised an exception:" << std::endl;
		std::cerr << e << std::endl;
		env.end();
		throw;
	}

	if (sol) {
		cp.out() << "Makespan \t: " << cp.getObjValue() << std::endl;

		for (int i = 0; i < Tasks.getSize(); i++) {
			cp.out() << "Task " << i + 1 << ": start time: " << cp.getStart(Tasks[i]) << " | end time: " << cp.getEnd(Tasks[i]) << endl;
		}
		set<int> usedAgents;

		//remove("c:\\Users\\bmc01\\Documents\\Visual Studio 2017\\Projects\\MR_MT_Solver\\MR_MT_Solver\\output.txt");
		//ofstream myfile2("output.txt", ios_base::app | ios_base::out);

		string path2 = path + to_string(run) + "_output.txt";
		ofstream myfile2(path2, ios_base::app | ios_base::out);

		for (int i = 0; i < Tasks.getSize(); i++) {
			myfile2 << cp.getStart(Tasks[i]) << "\t" << cp.getEnd(Tasks[i]) << "\t" << flush;
			//cout << cp.getStart(Tasks[i]) << "\t" << cp.getEnd(Tasks[i]) << "\t" << flush;
			for (int j = 0; j < fullModel->_vectorTasksPerAgents[i].size(); j++) {

				if (cp.isPresent(TasksV[i][j])) {
					usedAgents.insert(fullModel->_vectorTasksPerAgents[i][j]);
					myfile2 << fullModel->_vectorTasksPerAgents[i][j] << "\t" << fullModel->H[i + sigma] << endl;
					//cout << fullModel->_vectorTasksPerAgents[i][j] << "\t" << fullModel->H[i + sigma] << endl;
				}
			}
		}
		myfile2.close();


		//remove("c:\\Users\\bmc01\\Documents\\Visual Studio 2017\\Projects\\MR_MT_Solver\\MR_MT_Solver\\depot2.txt");
		//ofstream myfile3("depot2.txt", ios_base::app | ios_base::out);

		path2 = path + to_string(run) + "_depot2.txt";
		ofstream myfile3(path2, ios_base::app | ios_base::out);

		int deductable = 0;
		for (int i = 0; i < sigma; i++) {
			for (int j = 0; j < delta; j++) {
				if (cp.isPresent(DestDepot2[i][j])) {

					int dt = cp.getStart(DestDepot2[i][j]);
					auto it = find(usedAgents.begin(), usedAgents.end(), i);
					if (it == usedAgents.end()) {
						deductable += dt;
						dt = 0;
					}

					myfile3 << dt << endl;
					cout << "D[" << to_string(i) << "] " << dt << "\t" << flush;
				}
			}
		}
		myfile3 << cp.getStatus() << endl;
		myfile3 << cp.getObjGap() << endl;
		myfile3 << ga->pVirtual << endl;
		myfile3 << ga->pParallel << endl;

		myfile3.close();
		cout << endl; cout << endl;
		//cout << "Final obj. value: " << cp.getObjValue() - deductable << endl;
		cout << "Final obj. value: " << cp.getObjValue() << endl;
	}
	else { cp.out() << "No solution found." << std::endl; }
	cout << endl;
	env.end();
}

/***********/


/* Helper Methods */
void findParallel(vector<int> & _vector, vector<int> & _tasksPerAgents, vector<int> & _return, int sigma, int delta, int task) {

	for (int i = 0; i < _tasksPerAgents.size(); i++) {
		if (_vector[_tasksPerAgents[i]] == 0 && _tasksPerAgents[i] != task) {
			_return.emplace_back(_tasksPerAgents[i]);
		}
	}

	return;
}