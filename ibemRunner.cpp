/* ----------------------------------------------------------------------
   iBEM - inclusion-based Boundary Element Method
   This software is distributed under the GNU General Public License
   version 3 or any later version.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
    This file is part of iBEM.

    iBEM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    iBEM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with iBEM.  If not, see <https://www.gnu.org/licenses/>.

------------------------------------------------------------------------- */

#include "ibemRunner.h"
#include "inputGenerator.h"
#include<iostream>
#include<iterator>
#include<string>
#include<vector>
#include<map>
#include<queue>
#include "serviceWrapper.h"

using namespace std;

void ibemRunner::run(string mainService) {


	// Initiate variables
	vector<string> ALLServices{ "preprocessService", "iBEMSystemService", "systemSolverService", "postprocessService" };  // all the services involdes
	vector<string> preprocessorService;            // define the dependency for each service in order to form a DAG
	vector<string> iBEMSystemService;
	vector<string> systemSloverService;
	vector<string> postprocessService;
	std::map<string, int> indegree;                     // a mapping describes the number of prerequisites for each service
	std::map<string, vector<string>> dependency;        // dependency map
	std::map<string, vector<string>> graph;             // a mapping describes the dependents for each service


	// read the caltulation type from input fiels
	int ID = 0;

	ifstream myfile;
	myfile.open("calculationType.txt", ios::in);

	string T;
	getline(myfile, T);
	ID = stoi(T);

	// get the corresponding config based on ID
	shared_ptr<Config> config = inputGenerator(ID);


	// the following codes apply topological sorting algorithm to find the order of service to be execute.
	iBEMSystemService.push_back("preprocessService");
	systemSloverService.push_back("iBEMSystemService");
	postprocessService.push_back("systemSolverService");

	dependency.insert({ "preprocessService" , preprocessorService });
	dependency.insert({ "iBEMSystemService" , iBEMSystemService });
	dependency.insert({ "systemSolverService" , systemSloverService });
	dependency.insert({ "postprocessService" , postprocessService });

	for (size_t j = 0; j < ALLServices.size(); j++) {

		graph.insert({ ALLServices[j], vector<string> {} });
		indegree.insert({ ALLServices[j], 0 });

	}

	for (auto i = dependency.begin(); i != dependency.end(); i++)
	{
		auto k = i->first;
		auto v = i->second;

		for (size_t j = 0; j < v.size(); j++) {

			graph[v[j]].push_back(k);
			indegree[k] = indegree[k] + 1;

		}

	}

	vector<string> orders{};
	queue<string> frontqueue;


	for (auto i = indegree.begin(); i != indegree.end(); i++) {
		if (i->second == 0) {
			orders.push_back(i->first);
			frontqueue.push(i->first);
		}
	}


	while (!frontqueue.empty()) {
		string node = frontqueue.front();
		frontqueue.pop();
		auto dependList = graph[node];
		for (size_t i = 0; i < dependList.size(); i++) {
			indegree[dependList[i]] -= 1;
			if (indegree[dependList[i]] == 0) {
				orders.push_back(dependList[i]);
				frontqueue.push(dependList[i]);
			}
		}

	}

	// given the orders, the service excutes one by one (.run())
	cout << "-------------------\n";
	for (size_t i = 0; i < orders.size(); i++) {
		serviceWrapper serviceWrapper_{ orders[i],  *config };
		serviceWrapper_.run();
		if (orders[i] == mainService) {
			break;
		}
	}


}
