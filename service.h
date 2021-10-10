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