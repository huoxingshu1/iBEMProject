#pragma once
#include "BEMBuilder.h"
#include "Eigen/Dense"
using namespace Eigen;

class BEMBuilderElastic2D : public BEMBuilder {
public:
	double** record;
    void Integ3(double E, double nu, int* sctrb, double** Elcor, double** dUe, double** dTe, Ref<MatrixXd> xP, int Ncol);
    
	void addBEM(Config& config_);
    
	void addFluxEquivalentBEM(Config& config_);
    
	void addFluxEquivalentFirstOrderBEM(Config& config_);
    
	void addFluxEquivalentSecondOrderBEM(Config& config_);
    
	void ApplyBC(Config& config_);
    
	void getUTVector(Config& config_);
};

