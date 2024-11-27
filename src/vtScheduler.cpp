#include "include/vtScheduler.h"

vtScheduler::vtScheduler(shared_ptr<GA> Genetic, MODEL *_fullModel) {

	ga = Genetic; fullModel = _fullModel;

	/* Update the size of the chromosome for the full model */
	chromo = ga->population[ga->_bestIndividualsIndex];

	for (int i = 0; i < chromo.individual.size(); i++) {
		chromo.individual[i].resize(fullModel->V + fullModel->sigma);
	}

	// update destination depot for the full model
	for (int i = 0; i < fullModel->sigma; i++) {
		int _tmp = i;
		while (_tmp < chromo.individual[i].size()) {

			if (chromo.individual[i][_tmp] >= ga->population[ga->_bestIndividualsIndex].individual[i].size()) {
				chromo.individual[i][_tmp] = fullModel->V + fullModel->sigma + (chromo.individual[i][_tmp] - ga->population[ga->_bestIndividualsIndex].individual[i].size());
			}

			if (chromo.individual[i][_tmp] >= chromo.individual[i].size()) {
				break;
			}
			// Update tmp
			_tmp = chromo.individual[i][_tmp];
		}
	}

	// update the size of timeline for the full model 
	fullTimeline = ga->timeline;
	fullTimeline.resize(fullModel->V);

	//_depots.resize(fullModel->sigma);
}

void vtScheduler::solver()
{

	//vector<int> taskidx = { 23   , 28  ,  19  ,  22   , 30  ,  25  ,  26  ,  21  ,  18  ,  24  ,  20  ,  29  ,  17 , 27 };


	for (int ii = ga->model->V + fullModel->sigma; ii < fullModel->V + fullModel->sigma; ii++) {
		//int ii = taskidx[iii- (ga->model->V + fullModel->sigma)];

		if (fullTimeline[ii - fullModel->sigma].endTime == 0) {

			scheduleTask(ii);


		}
	}


	ga2 = shared_ptr<GA>(new GA(fullModel, 1));
	chromo.timeline = fullTimeline;
	ga2->evaluateXD(chromo);

	ga2->population.emplace_back(chromo);

	cout << "PC PENAL: " << ga2->population[0].penal.pc << endl;
	if (ga2->population[0].penal.pc != 0) {
		cout << "ovo govno opet jede govna " << endl;
	}
	ga2->_bestIndividualsIndex = 0;
	ga2->timeline = fullTimeline;
	ga2->depotStart = ga->depotStart;
	ganttLog(0);
	//cout << " basda " << endl;
}

bool vtScheduler::findTask(int _task, triplet &_tri, int ii)
{

	for (int i = 0; i < fullModel->sigma; i++) {
		int _tmp = i;
		while (_tmp < chromo.individual[i].size()) {

			if (chromo.individual[i][_tmp] >= chromo.individual[i].size()) { break; }
			//cout << chromo.individual[i][_tmp] << "  " << _task << endl;
			if (chromo.individual[i][_tmp] == _task) {
				if (checkColor(i, ii) == 0) {
					//&& checkPC(i, ii) == 0
					vector<pair<int, int>> startTimes(1);
					//TAKE PREC CONSTRAINTS INTO ACCOUNT
					_tri.agentID = i;
					bool flag = findLimits(_tri, startTimes, chromo.individual[i][_tmp], ii);
					if (!flag) {
						return false;
					}

					return true;
				}
				return false;
			}
			// Update tmp
			_tmp = chromo.individual[i][_tmp];
		}
	}

	return false;
	// the function should not reach this point
	cout << "Something went horribly wrong in vtScheduler findTask() function !" << endl;
}

triplet vtScheduler::findGap(vector<int> &a2t, int _task)
{
	vector<triplet> _pair(1);
	int diff = 0, minDiff = MAXINT, _tt = 0;;
	bool flag = false, flag2 = false;
	set<int> godHelpMe;

	auto _pcTime = findPClimits(_task);

	for (int i = 0; i < a2t.size(); i++) {
		startTime st = checkPC(a2t[i], _task - fullModel->sigma);

		for (int j = 0; j < transitionGaps[a2t[i]].size(); j++) {
			if (chromo.individual[a2t[i]][j] != 0) {
				diff = fullModel->taskDuration[_task] - transitionGaps[a2t[i]][j];
				cout << "i " << i << " j " << j << endl;


				if (diff < minDiff && checkColor(a2t[i], _task) == 0) {
					auto dur = fullModel->w[j][_task][a2t[i]];

					godHelpMe.insert(a2t[i]);
					if (j < fullModel->sigma) { _tt = 0; }
					else { _tt = fullTimeline[j - fullModel->sigma].endTime; }

					if (_tt >= _pcTime.first && (_tt + dur + fullModel->taskDuration[_task]) <= _pcTime.second) {


						int ft = 0;

						if (j >= fullModel->sigma) {
							ft = fullTimeline[j - fullModel->sigma].endTime;
						}


						if (ft >= st.earliest && (ft + fullModel->taskDuration[j]) <= st.latest) {
							minDiff = diff;
							flag2 = true;
							if (diff > 0) {
								_pair[0].agentID = a2t[i]; _pair[0].taskID = j; _pair[0].diff = diff;
								flag = false;
							}
							else {
								if (!flag) {
									_pair.clear();
									flag = true;
								}
								triplet tmp;
								tmp.agentID = a2t[i]; tmp.taskID = j; tmp.diff = diff;
								if (chromo.individual[tmp.agentID][tmp.taskID] == 0) {
									cout << "bla" << endl;
								}
								_pair.emplace_back(tmp);

							}
						}
					}
				}
			}
		}
	}
	if (!flag2) {
		triplet tmp;
		int rnd = getRandomIntegerInRange(0, (int)godHelpMe.size() - 1);

		tmp.agentID = *next(godHelpMe.begin(), rnd);
		tmp.diff = 0;

		int i = tmp.agentID;
		int _tmp = i;
		while (_tmp < chromo.individual[i].size()) {

			if (chromo.individual[i][_tmp] >= chromo.individual[i].size()) {
				tmp.taskID = _tmp;
				break;
			}

			// Update tmp
			_tmp = chromo.individual[i][_tmp];
		}

		fullTimeline[_task - fullModel->sigma].startTime = _pcTime.first;

		fullTimeline[_task - fullModel->sigma].endTime = fullTimeline[_task - fullModel->sigma].startTime + fullModel->taskDuration[_task];

		return tmp;
	}

	if (_pair.size() > 1) {
		int rnd = getRandomIntegerInRange(0, (int)_pair.size() - 1);
		return _pair[rnd];
	}
	else {
		return _pair[0];
	}

}
int vtScheduler::checkXD()
{
	for (int i = 0; i < fullTimeline.size(); i++) {

	}

	return 0;
}
//
//bool vtScheduler::tmpFindTask(int _task, triplet & _tri, int ii)
//{
//	return false;
//}
//
//void vtScheduler::tmpUpdateTimeline(int _task, triplet & _tri, int _tmpVal)
//{
//
//	// update timeline
//	tmpFullTimeline[_task - fullModel->sigma].agentIdx = _tri.agentID;
//	tmpFullTimeline[_task - fullModel->sigma].virtualTask = 1;
//	tmpFullTimeline[_task - fullModel->sigma].taskIdx = _task;
//
//	if (_tri.taskID < fullModel->sigma) {
//		tmpFullTimeline[_task - fullModel->sigma].startTime = 0;
//	}
//	else if (fullModel->R[_task][_tmpVal] == 1) {
//		tmpFullTimeline[_task - fullModel->sigma].startTime = tmpFullTimeline[_tmpVal - fullModel->sigma].startTime;
//	}
//	else {
//		tmpFullTimeline[_task - fullModel->sigma].startTime = tmpFullTimeline[_tri.taskID - fullModel->sigma].endTime;
//	}
//
//	tmpFullTimeline[_task - fullModel->sigma].endTime = tmpFullTimeline[_task - fullModel->sigma].startTime + fullModel->taskDuration[_task];
//
//	if (_tri.diff > 0) {
//
//		int mem_max_dur = 0;
//		vector<int> tasks2check;
//		int  old_tmp = _tri.taskID;
//
//
//		if (fullModel->H[old_tmp] == 1) {
//			int tmp = old_tmp;
//
//			while (fullModel->H[tmp] != 0) {
//				auto it = find(tmpChromo.individual[_tri.agentID].begin(), tmpChromo.individual[_tri.agentID].end(), tmp);
//				tmp = distance(tmpChromo.individual[_tri.agentID].begin(), it);
//				mem_task = tmp;
//
//			}
//			old_tmp = mem_task;
//		}
//
//		int _tmp = tmpChromo.individual[_tri.agentID][old_tmp];
//
//		while (_tmp < tmpChromo.individual[_tri.agentID].size()) {
//
//			int tmpEnd = 0;
//			int tmpStart = 0;
//			int memTmpEnd = 0;
//
//			if (old_tmp >= fullModel->sigma) {
//				tmpEnd = tmpFullTimeline[old_tmp - fullModel->sigma].endTime;
//				tmpStart = tmpFullTimeline[old_tmp - fullModel->sigma].startTime;
//
//			}
//
//			if (mem_task >= fullModel->sigma) {
//				memTmpEnd = tmpFullTimeline[mem_task - fullModel->sigma].endTime;
//			}
//
//			if (fullModel->H[old_tmp] == 0 && fullModel->H[_tmp] == 1) {
//				mem_task = old_tmp;
//
//				if (fullModel->R[old_tmp][_tmp] == 0) {
//					tmpFullTimeline[_tmp - fullModel->sigma].startTime = tmpEnd;
//					tmpFullTimeline[_tmp - fullModel->sigma].endTime = tmpFullTimeline[_tmp - fullModel->sigma].startTime + fullModel->taskDuration[_tmp];
//					mem_max_dur = fullModel->taskDuration[_tmp];
//				}
//				else {
//					tasks2check.emplace_back(old_tmp);
//					tmpFullTimeline[_tmp - fullModel->sigma].startTime = tmpStart;
//					tmpFullTimeline[_tmp - fullModel->sigma].endTime = tmpFullTimeline[_tmp - fullModel->sigma].startTime + fullModel->taskDuration[_tmp];
//					if (fullModel->taskDuration[_tmp] > mem_max_dur) {
//						mem_max_dur = fullModel->taskDuration[_tmp];
//					}
//				}
//			}
//			else if (fullModel->H[old_tmp] == 1 && fullModel->H[_tmp] == 1) {
//
//				if (fullModel->R[old_tmp][_tmp] == 0) {
//					tmpFullTimeline[_tmp - fullModel->sigma].startTime = tmpEnd;
//					tmpFullTimeline[_tmp - fullModel->sigma].endTime = tmpFullTimeline[_tmp - fullModel->sigma].startTime + fullModel->taskDuration[_tmp];
//					if (mem_max_dur < (tmpFullTimeline[_tmp - fullModel->sigma].endTime - memTmpEnd)) {
//						mem_max_dur = tmpFullTimeline[_tmp - fullModel->sigma].endTime - memTmpEnd;
//					}
//
//				}
//				else {
//					tasks2check.emplace_back(old_tmp);
//
//					int eTime = findEarliestSchedulingTime(tasks2check, _tmp);
//
//					tmpFullTimeline[_tmp - fullModel->sigma].startTime = eTime;
//					tmpFullTimeline[_tmp - fullModel->sigma].endTime = tmpFullTimeline[_tmp - fullModel->sigma].startTime + fullModel->taskDuration[_tmp];
//
//					if (fullModel->taskDuration[_tmp] > fullModel->taskDuration[old_tmp]) {
//						if (fullModel->taskDuration[_tmp] > mem_max_dur) {
//							mem_max_dur = fullModel->taskDuration[_tmp];
//						}
//					}
//					else {
//						if (fullModel->taskDuration[old_tmp] > mem_max_dur) {
//							mem_max_dur = fullModel->taskDuration[old_tmp];
//						}
//					}
//
//				}
//			}
//			else if (fullModel->H[old_tmp] == 1 && fullModel->H[_tmp] == 0) {
//				if (fullModel->R[old_tmp][_tmp] == 0) {
//					if (fullModel->w[mem_task][_tmp][_tri.agentID] > mem_max_dur) {
//						tmpFullTimeline[_tmp - fullModel->sigma].startTime = memTmpEnd + fullModel->w[mem_task][_tmp][_tri.agentID];
//					}
//					else {
//						tmpFullTimeline[_tmp - fullModel->sigma].startTime = memTmpEnd + mem_max_dur;
//					}
//
//					tmpFullTimeline[_tmp - fullModel->sigma].endTime = tmpFullTimeline[_tmp - fullModel->sigma].startTime + fullModel->taskDuration[_tmp];
//					mem_max_dur = 0;
//					tasks2check.clear();
//				}
//				else {
//					tasks2check.emplace_back(old_tmp);
//
//					int eTime = findEarliestSchedulingTime(tasks2check, _tmp);
//
//					tmpFullTimeline[_tmp - fullModel->sigma].startTime = eTime;
//					tmpFullTimeline[_tmp - fullModel->sigma].endTime = tmpFullTimeline[_tmp - fullModel->sigma].startTime + fullModel->taskDuration[_tmp];
//
//					if (fullModel->taskDuration[_tmp] > mem_max_dur) {
//						mem_max_dur += fullModel->taskDuration[_tmp];
//					}
//				}
//			}
//			else if (fullModel->H[old_tmp] == 0 && fullModel->H[_tmp] == 0) {
//				tasks2check.clear();
//				tmpFullTimeline[_tmp - fullModel->sigma].startTime = tmpEnd + fullModel->w[old_tmp][_tmp][_tri.agentID];
//				tmpFullTimeline[_tmp - fullModel->sigma].endTime = tmpFullTimeline[_tmp - fullModel->sigma].startTime + fullModel->taskDuration[_tmp];
//			}
//
//			if (chromo.individual[_tri.agentID][_tmp] >= chromo.individual[_tri.agentID].size()) {
//				if (fullModel->w[_tmp][ga->depots[_tri.agentID] + fullModel->sigma + fullModel->V][_tri.agentID] > mem_max_dur) {
//					ga->depotStart[_tri.agentID] = tmpFullTimeline[_tmp - fullModel->sigma].endTime + fullModel->w[_tmp][ga->depots[_tri.agentID] + fullModel->sigma + fullModel->V][_tri.agentID];
//				}
//				else {
//					ga->depotStart[_tri.agentID] = tmpFullTimeline[mem_task - fullModel->sigma].endTime + mem_max_dur;
//				}
//
//
//				break;
//			}
//
//			// Update tmp
//			old_tmp = _tmp;
//			_tmp = tmpChromo.individual[_tri.agentID][_tmp];
//		}
//	}
//}
//
//int vtScheduler::tmpInsertTask(int _task, triplet & _tri)
//{
//	// insert task
//	int _tmpVal = tmpChromo.individual[_tri.agentID][_tri.taskID];
//	if (_tmpVal == 0) {
//		cout << "bla" << endl;
//	}
//	tmpChromo.individual[_tri.agentID][_tri.taskID] = _task;
//	tmpChromo.individual[_tri.agentID][_task] = _tmpVal;
//
//	if (tmpChromo.individual[_tri.agentID][_task] == 0) {
//		cout << "bla" << endl;
//	}
//	return _tmpVal;
//}
//
//bool vtScheduler::tmpFindLimits(triplet & _tri, vector<pair<int, int>>& startTimes, int task, int insertTask)
//{
//
//	int min = tmpFullTimeline[task - fullModel->sigma].startTime;
//	int max = tmpFullTimeline[task - fullModel->sigma].endTime;
//
//	int tmpTime = 0, tmpTime2;
//
//	int index = -1;
//
//	startTimes[0].first = min;
//	startTimes[0].second = max;
//
//	for (int i = 0; i < tmpFullTimeline.size(); i++) {
//		if (i != (insertTask - fullModel->sigma) && i != (task - fullModel->sigma)) {
//
//			cout << i << endl;
//			if (tmpFullTimeline[i].startTime >= min && tmpFullTimeline[i].startTime < max && tmpFullTimeline[i].endTime > max) {
//				tmpTime = tmpFullTimeline[i].startTime - min;
//				if (tmpTime >= fullModel->taskDuration[insertTask] && fullModel->R[task][i + fullModel->sigma] == 0) {
//
//					/*updateTimes(startTimes, tmpFullTimeline[i].startTime, min);*/
//					updateTimes(startTimes, min, tmpFullTimeline[i].startTime);
//					index = i;
//				}
//				else if (tmpTime < fullModel->taskDuration[insertTask]) {
//					return false;
//				}
//			}
//			else if (tmpFullTimeline[i].endTime > min && tmpFullTimeline[i].endTime <= max && tmpFullTimeline[i].startTime < min) {
//				tmpTime = max - tmpFullTimeline[i].endTime;
//				if (tmpTime >= fullModel->taskDuration[insertTask] && fullModel->R[task][i + fullModel->sigma] == 0) {
//
//					updateTimes(startTimes, tmpFullTimeline[i].endTime, max);
//					index = i;
//				}
//				else if (tmpTime < fullModel->taskDuration[insertTask]) {
//					return false;
//				}
//			}
//			else if (tmpFullTimeline[i].startTime <= min && tmpFullTimeline[i].endTime > max) {
//				if (fullModel->R[task][i + fullModel->sigma] == 0) {
//					return false;
//				}
//			}
//			else if (tmpFullTimeline[i].startTime < min && tmpFullTimeline[i].endTime >= max) {
//				if (fullModel->R[task][i + fullModel->sigma] == 0) {
//					return false;
//				}
//			}
//			else if (tmpFullTimeline[i].startTime >= min && tmpFullTimeline[i].startTime <= max && tmpFullTimeline[i].endTime >= min && tmpFullTimeline[i].endTime <= max) {
//
//				if (fullModel->R[insertTask][i + fullModel->sigma] == 0) {
//					tmpTime = tmpFullTimeline[i].startTime - min;
//					tmpTime2 = max - tmpFullTimeline[i].endTime;
//
//					if (tmpTime > 0 && tmpTime2 > 0) {
//						if (tmpTime >= fullModel->taskDuration[insertTask]) {
//							updateTimes(startTimes, min, tmpFullTimeline[i].startTime);
//							index = i;
//						}
//						startTimes.emplace_back(make_pair(tmpFullTimeline[i].endTime, startTimes[startTimes.size() - 1].second));
//
//					}
//
//					else if (tmpTime > 0 && tmpTime >= fullModel->taskDuration[insertTask]) {
//						updateTimes(startTimes, min, tmpFullTimeline[i].startTime);
//						index = i;
//					}
//					else if (tmpTime2 > 0) {
//						updateTimes(startTimes, tmpFullTimeline[i].endTime, max);
//						index = i;
//					}
//				}
//			}
//
//		}
//		else {
//
//		}
//
//	}
//
//	if (index == -1) {
//		index = task - fullModel->sigma;
//	}
//
//	for (int i = 0; i < startTimes.size(); i++) {
//		if (startTimes[i].second == max || startTimes[i].second - startTimes[i].first >= fullModel->taskDuration[insertTask]) {
//			_tri.agentID = tmpFullTimeline[index].agentIdx;
//			_tri.taskID = index + fullModel->sigma;
//			_tri.diff = fullModel->taskDuration[insertTask] - (startTimes[i].second - startTimes[i].first);
//			return true;
//		}
//	}
//	return false;
//}

void vtScheduler::findParallel(vector<int> &_parallelVector) {

	_possibleParallelTasks.clear();
	for (int i = fullModel->sigma; i < _parallelVector.size(); i++) {
		if (_parallelVector[i] == 1 && fullTimeline[i - fullModel->sigma].endTime > 0) {
			_possibleParallelTasks.emplace_back(i);
		}
	}
	return;
}

void vtScheduler::updatePlan(int _task, triplet &_tri, int ii, vector<int> &a2t_ii)
{
	triplet  _tri2;
	createTransitionGaps();

	bool flag = findTask(_task, _tri, ii);
	if (!flag || _tri.diff > 0) {
		_tri2 = findGap(a2t_ii, ii);
		if (_tri.diff / 5 > _tri2.diff) { _tri = _tri2; }
	}


	cout << "insert task" << endl;
	int _tmpVal = insertTask(ii, _tri);
	cout << "update timeline " << endl;
	updateTimeline(ii, _tri, _tmpVal);

}

int vtScheduler::pickTask(int _taskD) {

	int _min = MAXINT, index = 0;

	for (int i = 0; i < _possibleParallelTasks.size(); i++) {

		int _currentMin = abs(fullModel->taskDuration[_possibleParallelTasks[i]] - fullModel->taskDuration[_taskD]);
		if (_currentMin < _min) {
			_min = _currentMin;
			index = i;
		}
	}

	return _possibleParallelTasks[index];
}

void vtScheduler::createTransitionGaps()
{
	/* Explanation */
	/* Task0->Task1, the gap between these two tasks is recorded in transitionGaps at index of Task0 */
	transitionGaps.clear();
	transitionGaps.resize(fullModel->sigma);

	for (int i = 0; i < fullModel->sigma; i++) {
		int _tmp = i, old_tmp = _tmp;
		int max_end_time = 0;
		transitionGaps[i].resize(fullModel->V + fullModel->sigma);

		while (_tmp < chromo.individual[i].size()) {

			//cout << _tmp << " " << chromo.individual[i][_tmp] << endl;

			if (chromo.individual[i][_tmp] >= chromo.individual[i].size()) {
				transitionGaps[i][_tmp] = ga->depotStart[i] - fullTimeline[_tmp - fullModel->sigma].endTime;
				break;
			}

			if (_tmp < fullModel->sigma) {
				transitionGaps[i][_tmp] = fullTimeline[chromo.individual[i][_tmp] - fullModel->sigma].startTime;
			}
			else {
				if (fullModel->R[_tmp][chromo.individual[i][_tmp]] == 0) {


					if (max_end_time < fullTimeline[_tmp - fullModel->sigma].endTime) {
						max_end_time = fullTimeline[_tmp - fullModel->sigma].endTime;

					}
					if (fullModel->H[_tmp] == 1 && fullModel->H[chromo.individual[i][_tmp]] == 0) {
						transitionGaps[i][_tmp] = fullTimeline[chromo.individual[i][_tmp] - fullModel->sigma].startTime - max_end_time;
						max_end_time = 0;
					}
					else {
						transitionGaps[i][_tmp] = fullTimeline[chromo.individual[i][_tmp] - fullModel->sigma].startTime - fullTimeline[_tmp - fullModel->sigma].endTime;
					}
				}
				else {
					transitionGaps[i][_tmp] = 0;
					if (max_end_time < fullTimeline[_tmp - fullModel->sigma].endTime) {
						max_end_time = fullTimeline[_tmp - fullModel->sigma].endTime;
					}
				}

			}

			// Update tmp
			old_tmp = _tmp;
			_tmp = chromo.individual[i][_tmp];
		}
	}
	return;
}

void vtScheduler::scheduleTask(int ii)
{
	cout << "------------ ii: " << ii - fullModel->sigma << " ------------" << endl;

	triplet _tri;
	pair<int, int> _pair;
	int _tmpVal; bool flag = false;

	createTransitionGaps();


	for (int i = 0; i < fullModel->vpc[ii - fullModel->sigma].after.size(); i++) {
		int __task = fullModel->vpc[ii - fullModel->sigma].after[i];
		//if (fullModel->vpc[__task].after.empty() && fullTimeline[__task].endTime == 0) {
		//	//_task = __task + fullModel->sigma;

		//	scheduleTask(__task + fullModel->sigma);


		//	break;
		//}
		//else {
		if (fullTimeline[__task].endTime == 0) {
			scheduleTask(__task + fullModel->sigma);
		}
		//}
	}
	/* Get agents that can perform task ii */
	auto a2t_ii = fullModel->_vectorTasksPerAgents[ii - fullModel->sigma];
	findParallel(fullModel->R[ii]);


	if (!_possibleParallelTasks.empty()) {
		originalChromo = chromo;
		originalFullTimeline = fullTimeline;
		int _task;
		while (true) {

			chromo = originalChromo;
			fullTimeline = originalFullTimeline;

			cout << "parallel to ii  " << flush;
			for (auto a : _possibleParallelTasks) { cout << a - fullModel->sigma << ", " << flush; }
			cout << endl;
			_task = pickTask(ii);

			cout << "------------ _task: " << _task - fullModel->sigma << " ------------" << endl;

			auto a2t_task = fullModel->_vectorTasksPerAgents[_task - fullModel->sigma];

			//_tri.diff = MAXINT;


			/*_tri = findGap(a2t_task, ii);

			_tmpVal = insertTask(ii, _tri);

			updateTimeline(ii, _tri, _tmpVal);*/
			_tri.diff = MAXINT;
			bool flag = findTask(_task, _tri, ii);

			auto it = find(_possibleParallelTasks.begin(), _possibleParallelTasks.end(), _task);
			std::iter_swap(it, _possibleParallelTasks.end() - 1);
			_possibleParallelTasks.pop_back();

			if (flag || _possibleParallelTasks.empty()) { break; }
			else { _task = pickTask(ii); }


		}

		updatePlan(_task, _tri, ii, a2t_ii);
	}
	else {

		_tri = findGap(a2t_ii, ii);

		_tmpVal = insertTask(ii, _tri);

		updateTimeline(ii, _tri, _tmpVal);

	}
	//if (fullTimeline[_task - fullModel->sigma].endTime == 0) {


	//	for (int i = 0; i < fullModel->vpc[_task - fullModel->sigma].after.size(); i++) {
	//		int __task = fullModel->vpc[_task - fullModel->sigma].after[i];
	//		if (fullModel->vpc[__task].after.empty() && fullTimeline[__task].endTime == 0) {
	//			//_task = __task + fullModel->sigma;

	//			scheduleTask(__task + fullModel->sigma);


	//			break;
	//		}
	//		else {
	//			if (fullTimeline[__task].endTime == 0) {
	//				scheduleTask(__task + fullModel->sigma);
	//			}
	//		}


	//	}

	//	_tri = findGap(a2t_task, _task);

	//	_tmpVal = insertTask(_task, _tri);

	//	updateTimeline(_task, _tri, _tmpVal);

	//	originalChromo = chromo;
	//	originalFullTimeline = fullTimeline;

	//	while (true) {

	//		chromo = originalChromo;
	//		fullTimeline = originalFullTimeline;

	//		//cout << " findGap " << endl;
	//		_tri = findGap(a2t_task, _task);

	//		_tmpVal = insertTask(_task, _tri);

	//		updateTimeline(_task, _tri, _tmpVal);

	//		_tri.diff = MAXINT;
	//		flag = findTask(_task, _tri, ii);

	//		auto it = find(_possibleParallelTasks.begin(), _possibleParallelTasks.end(), _task);
	//		std::iter_swap(it, _possibleParallelTasks.end() - 1);
	//		_possibleParallelTasks.pop_back();

	//		if (flag || !_possibleParallelTasks.empty()) { break; }
	//		else { _task = pickTask(ii); }

	//	}
	//
	//	updatePlan(_task, _tri, ii, a2t_ii);

	//}
	//else {
	//	cout << " findTask " << endl;

	//	_tri.diff = MAXINT;

	//	updatePlan(_task,_tri, ii, a2t_ii);
	//}

	for (int ag = 0; ag < fullModel->sigma; ag++) {
		printIndividual(chromo.individual[ag], ag); cout << endl;
		cout << "--------------------------------------------" << endl;
	}

	ga2 = shared_ptr<GA>(new GA(fullModel, 1));
	chromo.timeline = fullTimeline;
	ga2->evaluateXD(chromo);

	ga2->population.emplace_back(chromo);

	cout << "PC PENAL: " << ga2->population[0].penal.pc << endl;
	if (ga2->population[0].penal.pc != 0) {
		cout << "ovo govno opet jede govna " << endl;
	}
	ga2->_bestIndividualsIndex = 0;
	ga2->timeline = fullTimeline;
	ga2->depotStart = ga->depotStart;
	ganttLog(count);
	count++;

}

int vtScheduler::checkColor(int ag, int _task)
{
	auto *tmp = fullModel->getTasksPerAgents();

	if ((*tmp)[ag][_task - fullModel->sigma] == 1) {
		return 0;
	}

	return 1;
}

// Do this.
startTime vtScheduler::checkPC(int ag, int _task)
{
	/* Check Precedence Constraints */
	startTime st;
	st.earliest = 0; st.latest = MAXINT;

	/* Find Before */
	if (!fullModel->vpc[_task].before.empty()) {
		for (int j = 0; j < fullModel->vpc[_task].before.size(); j++) {
			if (fullTimeline[fullModel->vpc[_task].before[j]].endTime > st.latest) {
				st.latest = fullTimeline[fullModel->vpc[_task].before[j]].endTime;
			}
		}
	}

	/* Find After */
	if (!fullModel->vpc[_task].after.empty()) {
		for (int j = 0; j < fullModel->vpc[_task].after.size(); j++) {
			if (fullTimeline[fullModel->vpc[_task].after[j]].startTime > st.earliest) {
				st.earliest = fullTimeline[fullModel->vpc[_task].after[j]].startTime;
			}
		}
	}

	return st;
}

/* Helper Methods */

int vtScheduler::insertTask(int _task, triplet &_tri) {
	// insert task
	int _tmp = _tri.taskID;
	int i = _tri.agentID;
	
	while (_tmp < chromo.individual[i].size()) {

		if (chromo.individual[i][_tmp] >= chromo.individual[i].size()) { break; }
		//cout << chromo.individual[i][_tmp] << "  " << _task << endl;
		if (fullModel->R[_tmp][chromo.individual[i][_tmp]] == 1) {
			_tri.taskID = chromo.individual[i][_tmp];
			
		}
		else {
			break;
		}
		// Update tmp
		_tmp = chromo.individual[i][_tmp];
	}

	int _tmpVal = chromo.individual[_tri.agentID][_tri.taskID];
	if (_tmpVal == 0) {
		cout << "bla" << endl;
	}
	chromo.individual[_tri.agentID][_tri.taskID] = _task;
	chromo.individual[_tri.agentID][_task] = _tmpVal;

	if (chromo.individual[_tri.agentID][_task] == 0) {
		cout << "bla" << endl;
	}
	return _tmpVal;
}

pair<int, int> vtScheduler::findPClimits(int _task)
{
	cout << endl;
	cout << "Find pcLimits: ------------------- " << _task << endl;
	int pcStartTime = 0;
	int pcEndTime = MAXINT;

	if (!fullModel->_PV[_task].empty()) {
		for (int p = 0; p < fullModel->_PV[_task].size(); p++) {
			if (pcEndTime > fullTimeline[fullModel->_PV[_task][p] - fullModel->sigma].startTime && fullTimeline[fullModel->_PV[_task][p] - fullModel->sigma].endTime != 0) {
				pcEndTime = fullTimeline[fullModel->_PV[_task][p] - fullModel->sigma].startTime;
			}
		}
	}

	pair<int, int> new_pair;
	// this can be done more efficient, use fullModel->vpc
	for (int i = fullModel->sigma; i < fullModel->_PV.size(); i++) {
		auto it = find(fullModel->_PV[i].begin(), fullModel->_PV[i].end(), _task);
		if (it != fullModel->_PV[i].end()) {
			cout << i << flush;


			if (fullTimeline[i - fullModel->sigma].endTime > pcStartTime) {
				pcStartTime = fullTimeline[i - fullModel->sigma].endTime;
			}
			else {
				new_pair = findPClimits(i);
			}

			if ((new_pair.first + fullModel->taskDuration[i]) > pcStartTime) {
				pcStartTime = (new_pair.first + fullModel->taskDuration[i]);
			}

			cout << " startTime: " << pcStartTime << endl;

		}
	}



	return make_pair(pcStartTime, pcEndTime);
}

void vtScheduler::updateTimeline(int _task, triplet &_tri, int _tmpVal) {

	// update timeline
	fullTimeline[_task - fullModel->sigma].agentIdx = _tri.agentID;
	fullTimeline[_task - fullModel->sigma].virtualTask = 1;
	fullTimeline[_task - fullModel->sigma].taskIdx = _task;

	//int  old_tmp = _tri.taskID;

	//if (fullModel->H[_task] == 1) {
	//	int tmp = old_tmp;

	//	while (fullModel->H[tmp] != 0) {
	//		auto it = find(chromo.individual[_tri.agentID].begin(), chromo.individual[_tri.agentID].end(), tmp);
	//		tmp = distance(chromo.individual[_tri.agentID].begin(), it);
	//		mem_task = tmp;

	//	}
	//	old_tmp = mem_task;
	//}
	if (fullTimeline[_task - fullModel->sigma].endTime == 0) {
		if (_tri.taskID < fullModel->sigma) {
			fullTimeline[_task - fullModel->sigma].startTime = 0;
		}
		else if (fullModel->R[_task][_tri.taskID] == 1) {
			fullTimeline[_task - fullModel->sigma].startTime = fullTimeline[_tri.taskID - fullModel->sigma].startTime;
		}
		else {
			fullTimeline[_task - fullModel->sigma].startTime = fullTimeline[_tri.taskID - fullModel->sigma].endTime;
		}

		fullTimeline[_task - fullModel->sigma].endTime = fullTimeline[_task - fullModel->sigma].startTime + fullModel->taskDuration[_task];
	}

	if (_tri.diff > 0) {

		int mem_max_dur = 0;
		vector<int> tasks2check;
		int  old_tmp = _tri.taskID;


		if (fullModel->H[old_tmp] == 1) {
			int tmp = old_tmp;

			while (fullModel->H[tmp] != 0) {
				auto it = find(chromo.individual[_tri.agentID].begin(), chromo.individual[_tri.agentID].end(), tmp);
				tmp = distance(chromo.individual[_tri.agentID].begin(), it);
				mem_task = tmp;

			}
			old_tmp = mem_task;
		}

		int _tmp = chromo.individual[_tri.agentID][old_tmp];

		while (_tmp < chromo.individual[_tri.agentID].size()) {

			cout << _tmp << "\t" << chromo.individual[_tri.agentID][_tmp] << endl;

			int tmpEnd = 0;
			int tmpStart = 0;
			int memTmpEnd = 0;

			if (old_tmp >= fullModel->sigma) {
				tmpEnd = fullTimeline[old_tmp - fullModel->sigma].endTime;
				tmpStart = fullTimeline[old_tmp - fullModel->sigma].startTime;

			}

			if (mem_task >= fullModel->sigma) {
				memTmpEnd = fullTimeline[mem_task - fullModel->sigma].endTime;
			}

			if (fullModel->H[old_tmp] == 0 && fullModel->H[_tmp] == 1) {
				mem_task = old_tmp;

				if (fullModel->R[old_tmp][_tmp] == 0) {
					fullTimeline[_tmp - fullModel->sigma].startTime = tmpEnd;
					fullTimeline[_tmp - fullModel->sigma].endTime = fullTimeline[_tmp - fullModel->sigma].startTime + fullModel->taskDuration[_tmp];
					mem_max_dur = fullModel->taskDuration[_tmp];
				}
				else {
					tasks2check.emplace_back(old_tmp);
					fullTimeline[_tmp - fullModel->sigma].startTime = tmpStart;
					fullTimeline[_tmp - fullModel->sigma].endTime = fullTimeline[_tmp - fullModel->sigma].startTime + fullModel->taskDuration[_tmp];
					if (fullModel->taskDuration[_tmp] > mem_max_dur) {
						mem_max_dur = fullModel->taskDuration[_tmp];
					}
				}
			}
			else if (fullModel->H[old_tmp] == 1 && fullModel->H[_tmp] == 1) {

				if (fullModel->R[old_tmp][_tmp] == 0) {
					fullTimeline[_tmp - fullModel->sigma].startTime = tmpEnd;
					fullTimeline[_tmp - fullModel->sigma].endTime = fullTimeline[_tmp - fullModel->sigma].startTime + fullModel->taskDuration[_tmp];
					if (mem_max_dur < (fullTimeline[_tmp - fullModel->sigma].endTime - memTmpEnd)) {
						mem_max_dur = fullTimeline[_tmp - fullModel->sigma].endTime - memTmpEnd;
					}

				}
				else {
					tasks2check.emplace_back(old_tmp);

					int eTime = findEarliestSchedulingTime(tasks2check, _tmp);

					fullTimeline[_tmp - fullModel->sigma].startTime = eTime;
					fullTimeline[_tmp - fullModel->sigma].endTime = fullTimeline[_tmp - fullModel->sigma].startTime + fullModel->taskDuration[_tmp];

					if (fullModel->taskDuration[_tmp] > fullModel->taskDuration[old_tmp]) {
						if (fullModel->taskDuration[_tmp] > mem_max_dur) {
							mem_max_dur = fullModel->taskDuration[_tmp];
						}
					}
					else {
						if (fullModel->taskDuration[old_tmp] > mem_max_dur) {
							mem_max_dur = fullModel->taskDuration[old_tmp];
						}
					}

				}
			}
			else if (fullModel->H[old_tmp] == 1 && fullModel->H[_tmp] == 0) {
				if (fullModel->R[old_tmp][_tmp] == 0) {
					if (fullModel->w[mem_task][_tmp][_tri.agentID] > mem_max_dur) {
						fullTimeline[_tmp - fullModel->sigma].startTime = memTmpEnd + fullModel->w[mem_task][_tmp][_tri.agentID];
					}
					else {
						fullTimeline[_tmp - fullModel->sigma].startTime = memTmpEnd + mem_max_dur;
					}

					fullTimeline[_tmp - fullModel->sigma].endTime = fullTimeline[_tmp - fullModel->sigma].startTime + fullModel->taskDuration[_tmp];
					mem_max_dur = 0;
					tasks2check.clear();
				}
				else {
					tasks2check.emplace_back(old_tmp);

					int eTime = findEarliestSchedulingTime(tasks2check, _tmp);

					fullTimeline[_tmp - fullModel->sigma].startTime = eTime;
					fullTimeline[_tmp - fullModel->sigma].endTime = fullTimeline[_tmp - fullModel->sigma].startTime + fullModel->taskDuration[_tmp];

					if (fullModel->taskDuration[_tmp] > mem_max_dur) {
						mem_max_dur += fullModel->taskDuration[_tmp];
					}
				}
			}
			else if (fullModel->H[old_tmp] == 0 && fullModel->H[_tmp] == 0) {
				tasks2check.clear();
				fullTimeline[_tmp - fullModel->sigma].startTime = tmpEnd + fullModel->w[old_tmp][_tmp][_tri.agentID];
				fullTimeline[_tmp - fullModel->sigma].endTime = fullTimeline[_tmp - fullModel->sigma].startTime + fullModel->taskDuration[_tmp];
			}

			if (chromo.individual[_tri.agentID][_tmp] >= chromo.individual[_tri.agentID].size()) {
				if (fullModel->w[_tmp][ga->depots[_tri.agentID] + fullModel->sigma + fullModel->V][_tri.agentID] >= mem_max_dur) {
					ga->depotStart[_tri.agentID] = fullTimeline[_tmp - fullModel->sigma].endTime + fullModel->w[_tmp][ga->depots[_tri.agentID] + fullModel->sigma + fullModel->V][_tri.agentID];
				}
				else {
					ga->depotStart[_tri.agentID] = fullTimeline[mem_task - fullModel->sigma].endTime + mem_max_dur;
				}


				break;
			}

			// Update tmp
			old_tmp = _tmp;
			_tmp = chromo.individual[_tri.agentID][_tmp];
		}
	}
}

int vtScheduler::findEarliestSchedulingTime(vector<int> &tasks2check, int _tmp)
{
	int eTime = fullTimeline[tasks2check[0] - fullModel->sigma].startTime;
	for (int i = 0; i < tasks2check.size(); i++) {
		if (fullModel->R[tasks2check[i]][_tmp] == 0) {
			eTime = fullTimeline[tasks2check[i] - fullModel->sigma].endTime;
		}
	}

	return eTime;
}

bool vtScheduler::findLimits(triplet &_tri, vector<pair<int, int>> &startTimes, int task, int insertTask)
{
	int min = fullTimeline[task - fullModel->sigma].startTime;
	int max = fullTimeline[task - fullModel->sigma].endTime;

	int tmpTime = 0, tmpTime2;

	int index = -1;

	auto _pcTime = findPClimits(insertTask);

	if (_pcTime.first > min) {
		min = _pcTime.first;
	}
	if (_pcTime.second < max) {
		max = _pcTime.second;
	}

	startTimes[0].first = min;
	startTimes[0].second = max;

	for (int i = 0; i < fullTimeline.size(); i++) {
		if (_tri.agentID == fullTimeline[i].agentIdx && fullTimeline[i].endTime > 0) {
			if (i != (insertTask - fullModel->sigma) && i != (task - fullModel->sigma)) {

				cout << i << endl;

				if (fullTimeline[i].startTime >= min && fullTimeline[i].startTime < max && fullTimeline[i].endTime > max) {
					tmpTime = fullTimeline[i].startTime - min;
					if (tmpTime >= fullModel->taskDuration[insertTask] && fullModel->R[task][i + fullModel->sigma] == 0) {
						updateTimes(startTimes, min, fullTimeline[i].startTime);
						index = i;
					}
					else if (tmpTime < fullModel->taskDuration[insertTask]) { return false; }
				}
				else if (fullTimeline[i].endTime > min && fullTimeline[i].endTime <= max && fullTimeline[i].startTime < min) {
					tmpTime = max - fullTimeline[i].endTime;
					if (tmpTime >= fullModel->taskDuration[insertTask] && fullModel->R[task][i + fullModel->sigma] == 0) {
						updateTimes(startTimes, fullTimeline[i].endTime, max);
						index = i;
					}
					else if (tmpTime < fullModel->taskDuration[insertTask]) { return false; }
				}
				else if (fullTimeline[i].startTime <= min && fullTimeline[i].endTime > max) {
					if (fullModel->R[task][i + fullModel->sigma] == 0) { return false; }
				}
				else if (fullTimeline[i].startTime < min && fullTimeline[i].endTime >= max) {
					if (fullModel->R[task][i + fullModel->sigma] == 0) { return false; }
				}
				else if (fullTimeline[i].startTime >= min && fullTimeline[i].startTime <= max && fullTimeline[i].endTime >= min && fullTimeline[i].endTime <= max) {

					if (fullModel->R[insertTask][i + fullModel->sigma] == 0) {
						tmpTime = fullTimeline[i].startTime - min;
						tmpTime2 = max - fullTimeline[i].endTime;

						if (tmpTime > 0 && tmpTime2 > 0) {
							if (tmpTime >= fullModel->taskDuration[insertTask]) {
								updateTimes(startTimes, min, fullTimeline[i].startTime);
								index = i;
							}
							startTimes.emplace_back(make_pair(fullTimeline[i].endTime, startTimes[startTimes.size() - 1].second));
						}

						else if (tmpTime > 0 && tmpTime >= fullModel->taskDuration[insertTask]) {
							updateTimes(startTimes, min, fullTimeline[i].startTime);
							index = i;
						}
						else if (tmpTime2 > 0) {
							updateTimes(startTimes, fullTimeline[i].endTime, max);
							index = i;
						}
					}
				}
			}
		}
	}

	if (index == -1 && min < max) {
		
		_tri.taskID = task;
		_tri.diff = fullModel->taskDuration[insertTask] - fullModel->taskDuration[task];
			return true;
	}


	for (int i = 0; i < startTimes.size(); i++) {
		if (index > -1 && (startTimes[i].second == max || startTimes[i].second - startTimes[i].first >= fullModel->taskDuration[insertTask])) {
			_tri.agentID = fullTimeline[index].agentIdx;
			_tri.taskID = index + fullModel->sigma;
			_tri.diff = fullModel->taskDuration[insertTask] - (startTimes[i].second - startTimes[i].first);
			return true;
		} 
	}
	return false;
}

void vtScheduler::updateTimes(vector<pair<int, int>>& startTimes, int start, int end)
{
	for (int i = 0; i < startTimes.size(); i++) {
		if (startTimes[i].first < start) { startTimes[i].first = start; }
		if (startTimes[i].second > end) { startTimes[i].second = end; }
	}

	return;
}

void vtScheduler::printIndividual(vector<int> &vector, int _ag) {
	cout << endl;
	cout << endl;
	int _k = _ag;
	int count = 0;
	while (true) {

		cout << _k << "   " << flush;

		if (vector[_k] >= fullModel->sigma + fullModel->V) {
			cout << vector[_k] << "   " << flush;
			break;
		}

		_k = vector[_k];

		if (count > fullModel->Vtilde) {
			break;
		}
		count++;
	}
}

void vtScheduler::ganttLog(int count)
{
	//write tasks
	string path = ("c:\\Users\\bmc01\\Documents\\Visual Studio 2017\\Projects\\MR_MT_Solver\\MR_MT_Solver\\");
	string path2 = path + to_string(count) + "_drawGantt_withVT.txt";
	string path3 = path + to_string(count) + "_drawGantt_withVT_depot.txt";

	char* char_path = &path2[0];
	remove(char_path);

	ofstream myfile(path2, ios_base::app | ios_base::out);

	for (int i = 0; i < ga2->timeline.size(); i++) {
		myfile << ga2->timeline[i].startTime << "\t" << ga2->timeline[i].endTime << "\t" << flush;
		myfile << ga2->timeline[i].agentIdx << "\t" << ga2->timeline[i].virtualTask << "\t" << i << endl;
	}

	// write destination depots

	char_path = &path3[0];
	remove(char_path);

	ofstream myfile2(path3, ios_base::app | ios_base::out);

	for (int i = 0; i < ga2->depotStart.size(); i++) {
		myfile2 << ga2->depotStart[i] << endl;
	}

	myfile2.close();
}