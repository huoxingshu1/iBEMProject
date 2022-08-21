#pragma once
# include "Config.h"

class integrator {
public:
	virtual void addFieldToBEM(Config& config_) {};                    // builds up the effect of inclusion into traditional BEM equations.
	virtual void addFluxEquivalentInclusion(Config& config_) {};       // builds up inclusion impact on flux equivalent equations in EIM.
	virtual void addFluxEquivalentFirstOrderInclusion(Config& config_) {}; // builds up inclusion impact on first order of flux equivalent equations in EIM.
	virtual void addFluxEquivalentSecondOrderInclusion(Config& config_) {}; // builds up inclusion impact on second order of flux equivalent equations in EIM.
};
