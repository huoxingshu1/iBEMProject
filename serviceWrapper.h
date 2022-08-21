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

