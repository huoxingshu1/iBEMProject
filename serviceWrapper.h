#pragma once
# include <string>
# include "service.h"
# include "Config.h"
# include "classFinder.h"

class serviceWrapper {
public:

	string serviceType;
	Config& config;
	shared_ptr<BEMBuilder> bembuilder;
	shared_ptr<integrator> Integrator;
	shared_ptr<postProcessor> Postprocessor;
	shared_ptr<systemSolver> SystemSolver;

	serviceWrapper(string type_, Config& config_) : serviceType{ type_ }, config{ config_ }
	{
		bembuilder = getBEMBuilder(config.getbembuilder());
		Integrator = getIntegrator(config.getintegrator());
		Postprocessor = getPostProcessor(config.getpostprocessor());
		SystemSolver = getSystemSolver(config.getsystemsolver());
	};

	void run();
};


void serviceWrapper::run() {

	if (serviceType == "preprocessService") {
		service<preprocessServiceName> s{ config };
		s.run();
	}

	if (serviceType == "iBEMSystemService") {
		service<iBEMSystemServiceName> s{ config , *bembuilder , *Integrator };
		s.run();
	}


	if (serviceType == "systemSolverService") {
		service<systemSolverServiceName> s{ config , *bembuilder , *SystemSolver };
		s.run();
	}

	if (serviceType == "postprocessService") {
		service<postprocessServiceName> s{ config, *Postprocessor };
		s.run();
	};


};

