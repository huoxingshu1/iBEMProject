#pragma once
#include "Config.h"
class BEMBuilder {
public:
	virtual void addBEM(Config& config_) {};  // build up traditioanl boundary integral equations
	virtual void addFluxEquivalentBEM(Config& config_) {};  // build up constant term of flux equivalent equations from boundary integral
	virtual void addFluxEquivalentFirstOrderBEM(Config& config_) {};  // build up first order term of flux equivalent equations from boundary integral
	virtual void addFluxEquivalentSecondOrderBEM(Config& config_) {};  // build up second order term of flux equivalent equations from boundary integral
	virtual void ApplyBC(Config& config_) {};      // apply boundary conditions
	virtual void getUTVector(Config& config_) {};  // get the result vectors after solving system of linear equations
};