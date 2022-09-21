#pragma once
#include "BEMBuilder.h"
#include "Eigen/Dense"
using namespace Eigen;

class BEMBuilderPotential3D : public BEMBuilder {
public:
    double** record;
    void Integ3(int* sctrb, double** Elcor, double** dUe, double** dTe, Ref<MatrixXd> xP, int Ncol);
    void addBEM(Config& config_);
    void addFluxEquivalentBEM(Config& config_);
    void addFluxEquivalentFirstOrderBEM(Config& config_);
    void addFluxEquivalentSecondOrderBEM(Config& config_);
    void ApplyBC(Config& config_);
    void getUTVector(Config& config_);

    // those functions have same name, to avoid link errors, they are set as private members
private:
    int ffloord(int i, int j);
    int Find(int h, int e, int i, double** record, Ref<MatrixXi> NConnect, int NE);
};
