
/*
 This file contains the integration of Green's function of
 Elliptical particles
 
 */

#pragma once
#include "integrator.h"

class integratorElastic2D : public integrator {
public:
    void addFieldToBEM(Config& config_);
    void addFluxEquivalentInclusion(Config& config_);
    void addFluxEquivalentFirstOrderInclusion(Config& config_);
    void addFluxEquivalentSecondOrderInclusion(Config& config_);
};
