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