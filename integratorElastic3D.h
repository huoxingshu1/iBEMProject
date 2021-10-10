#pragma once
#include "integrator.h"

class integratorElastic3D : public integrator {
public:
	void addFieldToBEM(Config& config_);
	void addFluxEquivalentInclusion(Config& config_);
	void addFluxEquivalentFirstOrderInclusion(Config& config_);
	void addFluxEquivalentSecondOrderInclusion(Config& config_);
};

//void integratorElastic3D::addFieldToBEM(Config& config_) {};
//void integratorElastic3D::addFluxEquivalentInclusion(Config& config_) {};
//void integratorElastic3D::addFluxEquivalentFirstOrderInclusion(Config& config_) {};
//void integratorElastic3D::addFluxEquivalentSecondOrderInclusion(Config& config_) {};