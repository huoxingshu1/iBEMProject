#pragma once
#include "BEMBuilderElastic3D.h"
#include "integratorElastic3D.h"
#include "postProcessorElastic3D.h"
#include "systemSolverElastic3D.h"


shared_ptr<integrator> getIntegrator(string type) {

	if (type == "integratorElastic3D") {
		//return make_shared <integratorElastic3D> (new integratorElastic3D());
		return make_shared <integratorElastic3D>();
	}

	return NULL;
}


shared_ptr<BEMBuilder> getBEMBuilder(string type) {

	if (type == "BEMBuilderElastic3D") {
		return make_shared <BEMBuilderElastic3D>();
	}

	return NULL;
}

shared_ptr<postProcessor> getPostProcessor(string type) {
	if (type == "postProcessorElastic3D") {
		return make_shared <postProcessorElastic3D>();
	}
	return NULL;
}

shared_ptr<systemSolver> getSystemSolver(string type) {
	if (type == "systemSolverElastic3D") {
		return make_shared <systemSolverElastic3D>();
	}
	return NULL;
}