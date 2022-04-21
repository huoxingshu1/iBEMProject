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

#include <fstream>
#include "Eigen/Dense"
#include <iostream>
#include <sstream>
#include "Config.h"

using namespace std;
using namespace Eigen;

extern double mu_0, mu_1, nu0, nu1;

class configElastic3D_arbitrary : public Config {
public:
    int NN, NE, nump, num_NN, num_NE, ** index_E_ij;
    double XNU, XNU1, EE, EE1;

    // create three arrays to store Gaussian point weight and their changed coordinates
    double*** changed_coord, ** number_weight; int number_in_array, * count_number;

    int Dimension = 3; double** s_pt;

    string particleInput = "position_3D_arbitrary.txt";
    string bemInput = "BEM_3D_arbitrary.txt";
    string postProcessInput = "postprocess_3D_arbitrary.txt";

    MatrixXd eigen_point, NODES, BC, t_bc, HMAT, GMAT, A;
    MatrixXi eigen_connect, NConnect;
    VectorXd U, T, b;

    void Readparticle_num();
    void Readpartical_info();
    void ReadBEM_num();
    void ReadBEM_info();
    void ReadPostProcess_info();

    void Initialization_record();

    string getbembuilder();
    string getintegrator();
    string getsystemsolver();
    string getpostprocessor();


};
