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
#include "BEMBuilder.h"
#include "Eigen/Dense"
using namespace Eigen;

class BEMBuilderElastic3D_arbitrary : public BEMBuilder {
public:
    double** record;
    void addBEM(Config& config_);
    void addFluxEquivalentBEM(Config& config_);
    void ApplyBC(Config& config_);
    void getUTVector(Config& config_);

    // two empty function in arbitrary case
    void addFluxEquivalentFirstOrderBEM(Config& config_);
    void addFluxEquivalentSecondOrderBEM(Config& config_);

    // those functions have same name, to avoid link errors, they are set as private members
private:
    int ffloord(int i, int j);
    int Find(int h, int e, int i, double** record, Ref<MatrixXi> NConnect, int NE);

    void get_for_stress_BEM(double* Vnorm, double* point, double* GCcor, double** Dst, double** Dsu);
    void Integ3(double E, double nu, int* sctrb, double** Elcor, double** dUe, double** dTe, Ref<MatrixXd> xP, int Ncol);
    void Integ3_stress(double** Elcor, double** DDsu, double** DDst, double* point);
};
