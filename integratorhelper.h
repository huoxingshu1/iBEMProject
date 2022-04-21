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
#include "Eigen/Dense"

using namespace Eigen;

class integratorhelper3D_arbitrary {
public:
    int number_in_array = 500;
    double Volume_tet(double** coord);
    void tet_shape(double** coord, int n, double* SH);
    double tet_shape_10(double** coord, int n, double* x);
    void surface_normal(double* Svnorm, double* x1, double* x2, double* x3);
    void edge_direction_normal(double* Svnorm, double* Vnorm, double* Vdir, double* x1, double* x2);
    void tet_shape_10_coeff_A(double* SH_0, double* SH_1, double* SH_2, double* SH_3, double** coeff);
    void gauss_10_point_shared(int s, int num_NE, Ref<MatrixXd> eigen_point, Ref<MatrixXi> eigen_connect, \
        int& num_array, double* x, double*** returned_coord, double**** record_weight, double*** changed_coord, double** number_weight, int* count_number);
    double func_surface(double* ko, double* Sv, double* x);
    int judge_tet_point(double* x, double** coord);
};