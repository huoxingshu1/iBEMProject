
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

private:
    void get_addu_tensor(double* x, double a, double** DD_u3, int** index);
    void get_addu1_tensor(double* x, double a, double** DD_u4, int*** indexp);
    void get_addu2_tensor(double* x, double a, double** DD_u5, int**** indexpp);
};
