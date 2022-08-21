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

#pragma once
#include<string>
#include<iostream>
#include"preprocessService.h"
#include"postprocessService.h"
#include"iBEMSystemService.h"
#include "systemSolverService.h"

using namespace std;

template < const char* N> class service {
public:
	virtual void run() {};
};

char preprocessServiceName[] = "preprocessService";

template<> class service<preprocessServiceName> {
public:
	preprocessService* underlyingservice;
	service(Config& config_) {
		underlyingservice = new preprocessService(config_);
	};
	void run() {
		cout << "preprocessService\n";
		underlyingservice->run();
	};
	void test() {
		cout << "preprocessService\n";
	}
};

char iBEMSystemServiceName[] = "iBEMSystemService";

template<> class service<iBEMSystemServiceName> {
public:
	iBEMSystemService* underlyingservice;
	service(Config& config_, BEMBuilder& bembuilder_, integrator& integrator_) {
		underlyingservice = new iBEMSystemService(config_, bembuilder_, integrator_);
	};
	void run() {
		cout << "ibemSystemService\n";
		underlyingservice->run();
	};
	void test() {
		cout << "ibemSystemService\n";
	}
};

char systemSolverServiceName[] = "systemSolverService";

template<> class service<systemSolverServiceName> {
public:
	systemSolverService* underlyingservice;
	service(Config& config_, BEMBuilder& bembuilder_, systemSolver& systemsolver_) {
		underlyingservice = new systemSolverService(config_, bembuilder_, systemsolver_);
	};
	void run() {
		cout << "systemSolverService\n";
		underlyingservice->run();
	};
	void test() {
		cout << "systemSolverService\n";
	}
};

char postprocessServiceName[] = "postprocessService";

template<> class service<postprocessServiceName> {
public:
	postprocessService* underlyingservice;
	service(Config& config_, postProcessor& postProcesspr_) {
		underlyingservice = new postprocessService(config_, postProcesspr_);
	};
	void run() {
		cout << "postprocessService\n";
		underlyingservice->run();
	};
	void test() {
		cout << "postprocessService\n";
	}
};