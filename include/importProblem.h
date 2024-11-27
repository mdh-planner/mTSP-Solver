#pragma once
#include "definitions.h"
#include <string>
#include <fstream>
#include <iostream>
#include <cstring>

using namespace std;

struct Point
{
	double x;
	double y;
};

class importProblem
{

public:
	importProblem() {}

	vector<Problem> importAnders()
	{

		int numberOfInstances = 30;
		vector<int> numberOfNodes = {41, 81, 121, 161, 201, 401};
		vector<int> numberOfRobots = {2, 3, 4, 5, 6};

		vector<Problem> _return;

		string path = "AndersPaper\\input data\\";
		cout << "Reading Instances " << flush;
		for (int nr = 0; nr < numberOfRobots.size(); nr++)
		{
			cout << "." << flush;
			for (int nn = 0; nn < numberOfNodes.size(); nn++)
			{

				int nRoutes, Vtilde;
				vector<int> startNodes, deliveryNodes, actionCost;
				vector<vector<double>> routingCost;

				Problem object;
				// cout << "nr: " << nr << " | nn: " << nn << endl;

				// create path and filename
				string fname = "UC_" + to_string(numberOfRobots[nr]) + "_" + to_string(numberOfNodes[nn]) + "_8.dat";
				string fileName = path + fname;

				// import file into a string
				std::ifstream in(fileName);
				std::string fileContents((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

				string delimiter = "\n";
				size_t pos = 0;

				std::string token;

				while ((pos = fileContents.find(delimiter)) != std::string::npos)
				{
					token = fileContents.substr(0, pos);

					if (token.compare("NO_OF_ROUTES:") == 0)
					{
						fileContents.erase(0, pos + delimiter.length());
						pos = fileContents.find(delimiter);
						token = fileContents.substr(0, pos);
						fileContents.erase(0, pos + delimiter.length());

						nRoutes = stoi(token);
						continue;
					}

					if (token.compare("NO_OF_NODES:") == 0)
					{
						fileContents.erase(0, pos + delimiter.length());
						pos = fileContents.find(delimiter);
						token = fileContents.substr(0, pos);
						fileContents.erase(0, pos + delimiter.length());

						Vtilde = stoi(token);
						continue;
					}

					if (token.compare("START_NODES:") == 0)
					{
						fileContents.erase(0, pos + delimiter.length());
						pos = fileContents.find(delimiter);
						token = fileContents.substr(0, pos);
						fileContents.erase(0, pos + delimiter.length());

						for (int i = 0; i < token.length(); i++)
						{
							if (token[i] != ' ')
							{

								startNodes.push_back(token[i] - '0');
							}
						}
						continue;
					}

					if (token.compare("DELIVERY_NODES:") == 0)
					{
						fileContents.erase(0, pos + delimiter.length());
						pos = fileContents.find(delimiter);
						token = fileContents.substr(0, pos);
						fileContents.erase(0, pos + delimiter.length());

						/*	for (int i = 0; i < token.length(); i++) {
								if (token[i] != ' ') {

									deliveryNodes.push_back(token[i] - '0');
								}
							}*/
						string delimiter2 = " ";
						int pos2;
						string token2;
						while ((pos2 = token.find(delimiter2)) != std::string::npos)
						{

							token2 = token.substr(0, pos);

							int a = stoi(token2.c_str());
							deliveryNodes.push_back(a);
							token.erase(0, pos2 + delimiter2.length());
						}
						int a = stoi(token.c_str());
						deliveryNodes.push_back(a);

						continue;
					}

					if (token.compare("ACTION_COST:") == 0)
					{
						fileContents.erase(0, pos + delimiter.length());
						pos = fileContents.find(delimiter);
						token = fileContents.substr(0, pos);
						fileContents.erase(0, pos + delimiter.length());

						for (int i = 0; i < token.length(); i++)
						{
							if (token[i] != ' ')
							{

								actionCost.push_back(token[i] - '0');
							}
						}
						continue;
					}

					if (token.compare("ROUTING_COST:") == 0)
					{
						routingCost.resize(Vtilde);
						fileContents.erase(0, pos + delimiter.length());
						for (int r = 0; r < Vtilde; r++)
						{

							pos = fileContents.find(delimiter);
							token = fileContents.substr(0, pos);
							fileContents.erase(0, pos + delimiter.length());
							string delimiter2 = " ";
							int pos2;
							string token2;
							while ((pos2 = token.find(delimiter2)) != std::string::npos)
							{

								token2 = token.substr(0, pos);

								double a = atof(token2.c_str());
								routingCost[r].push_back(a);
								token.erase(0, pos2 + delimiter2.length());
							}
							double a = atof(token2.c_str());
							routingCost[r].push_back(a);
						}
						continue;
					}

					fileContents.erase(0, pos + delimiter.length());
				}

				// create task object
				int sigma = startNodes.size();
				int delta = deliveryNodes.size();
				for (int i = sigma; i < Vtilde - delta; i++)
				{

					Task T;

					/*Not in the text files */
					T.color = 0;
					T.prec = -1;
					T.index = 0;
					T.virt = 0;
					T.reqA = 1;
					T.x = -1;
					/*********/
					T.duration = actionCost[i];

					object.T.push_back(T);
				}

				for (int i = 0; i < sigma; i++)
				{

					Agent A;
					Source S;

					/*Not in the text files */
					A.color.push_back(0);
					A.v = 1;

					object.A.push_back(A);
					object.src.push_back(S);
				}

				for (int i = 0; i < delta; i++)
				{
					Destination D;
					object.dest.push_back(D);
				}
				object.transitionMatrix = routingCost;
				object.numberOfUsedAgents = nRoutes;
				_return.push_back(object);
			}
		}

		cout << endl;
		return _return;
	}

	static vector<Problem> importInstanceUVA()
	{

		int _numberOfInstances = 1;
		Problem object;
		vector<Problem> _return(_numberOfInstances);
		std::ifstream inputFile("../instance_3.txt");

		if (!inputFile.is_open())
		{
			std::cerr << "Failed to open the file." << std::endl;
			
		}

		std::vector<Point> A_points, T_points, D_points;
		std::vector<Point> *currentContainer = &A_points;

		double x, y;
		while (inputFile >> x >> y)
		{
			if (x == -1 && y == -1)
			{
				if (currentContainer == &A_points)
				{
					currentContainer = &T_points;
				}
				else if (currentContainer == &T_points)
				{
					currentContainer = &D_points;
				}
				continue;
			}

			currentContainer->push_back({x, y});
		}

		inputFile.close();

		// Assign the points to the corresponding objects
		Task T;
		Agent A;
		Source S;
		Destination D;

		// Print the points in each container
		std::cout << "Points in A:" << std::endl;
		for (const auto &point : A_points)
		{
			A.x = point.x;
			A.y = point.y;
			A.color = {0};
			A.v = 2;
			A.commRadius = 225;
			object.A.emplace_back(A);

			S.x = point.x;
			S.y = point.y;
			object.src.emplace_back(S);
			std::cout << point.x << "\t" << point.y << std::endl;
		}

		std::cout << "\nPoints in T:" << std::endl;
		for (const auto &point : T_points)
		{
			T.x = point.x;
			T.y = point.y;
			T.color = 0;
			T.duration = 0;
			T.prec = -1;
			T.index = 0;
			T.virt = 0;
			T.reqA = 1;
			object.T.emplace_back(T);
			std::cout << point.x << "\t" << point.y << std::endl;
		}

		std::cout << "\nPoints in D:" << std::endl;
		for (const auto &point : D_points)
		{
			D.x = point.x;
			D.y = point.y;
			object.dest.emplace_back(D);

			std::cout << point.x << "\t" << point.y << std::endl;
		}
		std::cout << std::endl;
		object.pParallel = 0;
		object.pVirtual = 0;
		_return[0] = object;
		return _return;
	}

	// vector<Problem> importInstanceECTSP()
	// {

	// 	int _numberOfInstances = 10;
	// 	Problem object;
	// 	vector<Problem> _return(_numberOfInstances);

	// 	string path = "Test Instances\\";

	// 	for (int i = 0; i < _numberOfInstances; i++)
	// 	{
	// 		string fileName = path + "Instance " + to_string(i) + "\\Cities_" + to_string(i) + ".txt";
	// 		ifstream infile;
	// 		infile.exceptions(ifstream::failbit | ifstream::badbit);
	// 		infile.open(fileName.c_str());

	// 		string str;
	// 		char *pch;
	// 		char *line;
	// 		Task T;
	// 		Agent A;
	// 		Source S;
	// 		Destination D;

	// 		while (infile.peek() != EOF)
	// 		{

	// 			getline(infile, str);
	// 			line = _strdup(str.c_str());
	// 			pch = strtok(line, " ");

	// 			if (strcmp(pch, "City") == 0)
	// 			{
	// 				getline(infile, str);
	// 				line = _strdup(str.c_str());
	// 				pch = strtok(line, " ");
	// 			}

	// 			pch = strtok(NULL, " ");
	// 			T.x = atoi(pch);

	// 			pch = strtok(NULL, " ");
	// 			T.y = atoi(pch);

	// 			pch = strtok(NULL, " ");
	// 			T.duration = atoi(pch);

	// 			pch = strtok(NULL, " ");
	// 			T.color = atoi(pch);

	// 			pch = strtok(NULL, " ");
	// 			if (atoi(pch) == -1)
	// 			{
	// 				T.prec = -1;
	// 			}
	// 			else
	// 			{
	// 				T.prec = atoi(pch);
	// 			}

	// 			/*Not in the text files, but added so I can use the model with CP planner */
	// 			T.index = 0;
	// 			T.virt = 0;
	// 			T.reqA = 1;
	// 			/*********/

	// 			object.T.emplace_back(T);
	// 		}

	// 		object.pParallel = 0;
	// 		object.pVirtual = 0;

	// 		fileName = path + "Instance " + to_string(i) + "\\Depots_" + to_string(i) + ".txt";
	// 		ifstream infile2;
	// 		infile2.exceptions(ifstream::failbit | ifstream::badbit);
	// 		infile2.open(fileName.c_str());

	// 		while (infile2.peek() != EOF)
	// 		{

	// 			getline(infile2, str);
	// 			line = _strdup(str.c_str());
	// 			pch = strtok(line, " ");

	// 			if (strcmp(pch, "destinationDepot") == 0)
	// 			{
	// 				getline(infile2, str);
	// 				line = _strdup(str.c_str());
	// 				pch = strtok(line, " ");
	// 			}

	// 			pch = strtok(NULL, " ");
	// 			D.x = atoi(pch);

	// 			pch = strtok(NULL, " ");
	// 			D.y = atoi(pch);

	// 			object.dest.emplace_back(D);
	// 		}

	// 		fileName = path + "Instance " + to_string(i) + "\\Salespersons_" + to_string(i) + ".txt";
	// 		ifstream infile3;
	// 		infile3.exceptions(ifstream::failbit | ifstream::badbit);
	// 		infile3.open(fileName.c_str());

	// 		while (infile3.peek() != EOF)
	// 		{

	// 			getline(infile3, str);
	// 			line = _strdup(str.c_str());
	// 			pch = strtok(line, " ");

	// 			if (strcmp(pch, "Salesperson") == 0)
	// 			{
	// 				getline(infile3, str);
	// 				line = _strdup(str.c_str());
	// 				pch = strtok(line, " ");
	// 			}

	// 			pch = strtok(NULL, " ");
	// 			A.x = atoi(pch);
	// 			S.x = A.x;

	// 			pch = strtok(NULL, " ");
	// 			A.y = atoi(pch);
	// 			S.y = A.y;

	// 			while (true)
	// 			{
	// 				pch = strtok(NULL, " ");
	// 				if (atoi(pch) < 10)
	// 				{
	// 					A.color.emplace_back(atoi(pch));
	// 				}
	// 				else
	// 				{
	// 					A.v = atoi(pch);
	// 					break;
	// 				}
	// 			}

	// 			object.src.emplace_back(S);
	// 			object.A.emplace_back(A);
	// 			A.color.clear();
	// 		}
	// 		_return[i] = object;

	// 		object = {};
	// 	}

	// 	return _return;
	// }
};
