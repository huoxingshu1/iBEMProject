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
#include "integrator.h"
#include "integratorhelper.h"
#include "Uniform_Polyhedral_integral.h"
#include "Quadratic_Polyhedral_integral.h"
#include "Linear_Polyhedral_integral.h"

class integratorElastic3D_arbitrary : public integrator {
public:
    Uniform_Polyhedral U_P; Linear_Polyhedral L_P; Quadratic_Polyhedral Q_P;
    integratorhelper3D_arbitrary I3A;

    int number_in_array = 500;

    void addFieldToBEM(Config& config_);
    void addFluxEquivalentInclusion(Config& config_);
    
    // two empty functions
    void addFluxEquivalentFirstOrderInclusion(Config& config_);
    void addFluxEquivalentSecondOrderInclusion(Config& config_);
    void Polyhedral_Eshelby_strain_all_10(int i, int j, int k, int l, double** coord, double* x, double* result);
    void Polyhedral_Eshelby_displacement_all_10(int i, int j, int k, double** coord, double* x, double* result);
private:
 
    void displacement_collect_normal(double** coord, double* x, int** index, double*** ret_u3);

    void Eshelby_collect_normal_quadratic(int i, int j, double** coord, double* x, int** index_E_ij, double*** ret_u4);
    void Eshleby_strain_polyhedra(int num_NN, int num_NE, double** A, Ref<MatrixXd> eigen_point, Ref<MatrixXi> eigen_connect, \
        int** index_E_ij, double*** changed_coord, double** number_weight, int* count_number);

};