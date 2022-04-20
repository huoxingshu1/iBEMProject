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
#include "Config.h"
#include "BEMBuilder.h"
#include "integrator.h"

class iBEMSystemService {
public:
	Config& config;
	BEMBuilder& bembuilder;
	integrator& Integrator;
	iBEMSystemService(Config& config_, BEMBuilder& bembuilder_, integrator& integrator_) : config{ config_ }, bembuilder{ bembuilder_ }, Integrator{ integrator_ }{};
	void run();
	void tradition_BEM();             // build up traditioanl boundary integral equations
	void addFieldToBEM();        // build up effects from inclusion to boundary integral equations
	void addFluxEquivalentBEM(); // build up constant term of flux equivalent equations from boundary integral
	void addFluxEquivalentFirstOrderBEM();  // build up first order term of flux equivalent equations from boundary integral
	void addFluxEquivalentSecondOrderBEM();  // build up second order term of flux equivalent equations from boundary integral
	void addFluxEquivalentInclusion();  // build up  stress quivalent equations for includsions (equivalent includsion method)


};

/*
	In ellipsoidal cases, the interchange of function:
	(a) addFluxEquivalentInclusion (b) addFluxEquivalentBEM
	will not change the results; 

	However, for arbitrary cases, we need to calculate the Gaussian weights
	for subdomain collocation method. Hence, it is necessary to run inclusion functions
	firtstly!

*/

void iBEMSystemService::run() {
	this->tradition_BEM();
	this->addFieldToBEM();
	this->addFluxEquivalentInclusion();
	this->addFluxEquivalentBEM();
	this->addFluxEquivalentFirstOrderBEM();
	this->addFluxEquivalentSecondOrderBEM();
	

}


void iBEMSystemService::tradition_BEM() {
	bembuilder.addBEM(config);
}

void iBEMSystemService::addFieldToBEM() {
	Integrator.addFieldToBEM(config);

}

void iBEMSystemService::addFluxEquivalentBEM() {
	bembuilder.addFluxEquivalentBEM(config);

}

void iBEMSystemService::addFluxEquivalentFirstOrderBEM() {
	bembuilder.addFluxEquivalentFirstOrderBEM(config);

}

void iBEMSystemService::addFluxEquivalentSecondOrderBEM() {
	bembuilder.addFluxEquivalentSecondOrderBEM(config);

}

void iBEMSystemService::addFluxEquivalentInclusion() {
	Integrator.addFluxEquivalentInclusion(config);
	Integrator.addFluxEquivalentFirstOrderInclusion(config);
	Integrator.addFluxEquivalentSecondOrderInclusion(config);

}