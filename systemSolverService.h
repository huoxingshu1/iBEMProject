#pragma once
#include "Config.h"
#include "BEMBuilder.h"
#include "systemSolver.h"

class systemSolverService {
public:
	Config& config;
	BEMBuilder& bembuilder;
	systemSolver& systemsolver;

	systemSolverService(Config& config_, BEMBuilder& bembuilder_, systemSolver& systemSolver_) : config{ config_ }, bembuilder{ bembuilder_ }, systemsolver{ systemSolver_} {};
	void run();

};

void systemSolverService::run() {
	bembuilder.ApplyBC(config);            // apply boundary conditions
	systemsolver.solve(config);            // solve the iBEM system of equations
	bembuilder.getUTVector(config);
}
