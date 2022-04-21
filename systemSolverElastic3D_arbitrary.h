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
# include "systemSolver.h"
# include "Eigen/Dense"
# include "configElastic3D_arbitrary.h"

using namespace Eigen;

class systemSolverElastic3D_arbitrary : public systemSolver {

    void solve(Config& config_);  // solve the syetem of equations
    template <typename T> void _solve(Config& config_);
};

void systemSolverElastic3D_arbitrary::solve(Config& config_) {

    _solve<configElastic3D_arbitrary>(config_);

}

template <typename T>
void systemSolverElastic3D_arbitrary::_solve(Config& config_) {

    T& config = dynamic_cast<T&> (config_);
    MatrixXd& AA = config.A;
    VectorXd& BB = config.b;

    config.X = AA.fullPivLu().solve(BB);

    //    cout << config.X;
    double r = (AA * config.X - BB).norm();
    cout << "residual of solver" << " " << r << endl;
}
