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
#include "postProcessor.h"
#include "configElastic2D_arbitrary.h"
#include "Config.h"
#include "integratorElastic2D_arbitrary.h"
#include "BEMElasticHelper2D_arbitrary.h"
#include "Eigen/Dense"

using namespace Eigen;
class postProcessorElastic2D_arbitrary : public postProcessor {
public:
    void getFlux(Config& config_);
    void getField(Config& config_);
    BEM_2D_arbitrary B2D_A; integratorElastic2D_arbitrary INT_A;

private:
    void post_processor(int NN, int NE, Ref<MatrixXi> NConnect, Ref<MatrixXd> NODES, Ref<VectorXd> U, Ref<VectorXd> T,\
        int nump, int num, int num_NE, double** s_p, double** s_pb, double** e_p, int** index_E_ij, \
        Ref<MatrixXd> x_o, Ref<MatrixXi> eigen_connect, Ref<MatrixXd> x_p, double* Eigen);

    void post_bem_stress(int nump, int NE, int NN, Ref<MatrixXd> x_p, Ref<MatrixXi> NConnect, Ref<MatrixXd> NODES, double** s_pb, Ref<VectorXd> T, Ref<VectorXd> U);

    void Integ2_stress_post(double** Elcor, double** DDsup, double** DDstp, double* point);

    void get_for_stress_BEM_post(double* Vnorm, double* point, double* GCcor, double** Dst, double** Dsu);

    int pnpoly(int nvert, double* vertx, double* verty, double testx, double testy);
};
