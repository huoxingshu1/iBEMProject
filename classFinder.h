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
// 3D library
#include "BEMBuilderElastic3D.h"
#include "integratorElastic3D.h"
#include "postProcessorElastic3D.h"
#include "systemSolverElastic3D.h"

// 2D library
#include "integratorElastic2D.h"
#include "BEMBuilderElastic2D.h"
#include "postProcessorElastic2D.h"
#include "systemSolverElastic2D.h"

// 2D arbitrary library
#include "integratorElastic2D_arbitrary.h"
#include "BEMBuilderElastic2D_arbitrary.h"
#include"postProcessorElastic2D_arbitrary.h"
#include "systemSolverElastic2D_arbitrary.h"

/*
This header file allows iBEM to locate which specific problems to be run

The string names are defined in inputGenerator and config. h

*/

shared_ptr<integrator> getIntegrator(string type) {

	if (type == "integratorElastic3D") {
		//return make_shared <integratorElastic3D> (new integratorElastic3D());
		return make_shared <integratorElastic3D>();
	}

	else if (type == "integratorElastic2D") {
		return make_shared <integratorElastic2D>();
	}

	else if (type == "integratorElastic2D_arbitrary") {
		return make_shared <integratorElastic2D_arbitrary>();
	}

	return NULL;
}


shared_ptr<BEMBuilder> getBEMBuilder(string type) {

	if (type == "BEMBuilderElastic3D") {
		return make_shared <BEMBuilderElastic3D>();
	}

	else if (type == "BEMBuilderElastic2D") {
		return make_shared <BEMBuilderElastic2D>();
	}

	else if (type == "BEMBuilderElastic2D_arbitrary") {
		return make_shared <BEMBuilderElastic2D_arbitrary>();
	}

	return NULL;
}

shared_ptr<postProcessor> getPostProcessor(string type) {
	if (type == "postProcessorElastic3D") {
		return make_shared <postProcessorElastic3D>();
	}
	else if (type == "postProcessorElastic2D") {
		return make_shared <postProcessorElastic2D>();
	}
	else if (type == "postProcessorElastic2D_arbitrary") {
		return make_shared <postProcessorElastic2D_arbitrary>();
	}
	return NULL;
}

shared_ptr<systemSolver> getSystemSolver(string type) {
	if (type == "systemSolverElastic3D") {
		return make_shared <systemSolverElastic3D>();
	}
	else if (type == "systemSolverElastic2D") {
		return make_shared <systemSolverElastic2D>();
	}
	else if (type == "systemSolverElastic2D_arbitrary") {
		return make_shared <systemSolverElastic2D_arbitrary>();
	}

	return NULL;
}