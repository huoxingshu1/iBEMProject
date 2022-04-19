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
#include "BEMElasticHelper2D.h"
#include "Eigen/Dense"
using namespace Eigen;

class BEMBuilderElastic2D : public BEMBuilder {
public:
    double** record; BEM_2D B2D;    
    
	void addBEM(Config& config_);
    
	void addFluxEquivalentBEM(Config& config_);
    
	void addFluxEquivalentFirstOrderBEM(Config& config_);
    
	void addFluxEquivalentSecondOrderBEM(Config& config_);
    
	void ApplyBC(Config& config_);
    
	void getUTVector(Config& config_);

private:
    void Integ3(double E, double nu, int* sctrb, double** Elcor, double** dUe, double** dTe, Ref<MatrixXd> xP, int Ncol);
    int Find(int h, int e, int i, double** record, Ref<MatrixXi> NConnect, int NE);
    int  ffloord(int i, int j);

    // for functions addFluxEquivalentBEM to second order
    void sr_ijk_inclusion(double* Vnorm, double* point, double* GCcor, double** Dst, double** Dsu);
    void sr_ijk_inclusion_linear(double* Vnorm, double* point, double* GCcor, double** Dstw, double** Dsuw);
    void sr_ijk_inclusion_quadratic(double* Vnorm, double* point, double* GCcor, double** Dstwz, double** Dsuwz);

    // stress and to second order
    void Integ2e_stress(double** Elcor, double** DDsu, double** DDst, double* point);
    void Integ2e_stress_linear(double** Elcor, double** DDsuw, double** DDstw, double* point);
    void Integ2e_stress_quadratic(double** Elcor, double** DDsuwz, double** DDstwz, double* point);
};

