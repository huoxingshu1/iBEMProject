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
#include "GlobalConstant.h"
#include "Eigen/Dense"


using namespace Eigen;

//////////////////////////////////////////////////////////

namespace Elastic3DSpherical {

    inline double I_n0(double lambda, double a, int n);

    inline double I_n1(double lambda, double a, double* x, int n, int p);

    inline double I_n2(double lambda, double a, double* x, int n, int p, int q);

    inline double I_n3(double lambda, double a, double* x, int n, int p, int q, int r);

    inline double I_n4(double lambda, double a, double* x, int n, int p, int q, int r, int t);

    inline double I_n5(double lambda, double a, double* x, int n, int p, int q, int r, int t, int h);

    ////////////////////////////////////////////////////////////
    inline double V_00(double* x, double a);

    inline double V_10(double* x, double a, int i);

    inline double V_20(double* x, double a, int i, int j);

    inline double V_n0(double* x, double a, int n);


    inline double V_n1(double* x, double a, int n, int p);

    inline double V_n2(double* x, double a, int n, int p, int q);

    inline double V_n3(double* x, double a, int n, int p, int q, int r);

    inline double V_n4(double* x, double a, int n, int p, int q, int r, int t);

    inline double V_n5(double* x, double a, int n, int p, int q, int r, int t, int s);

    inline double V_n6(double* x, double a, int n, int p, int q, int r, int t, int s, int h);

    ////////////////////////////////////////////////////////////
    inline double Phi_00(double* x, double a);

    inline double Phi_10(double* x, double a, int n);

    inline double Phi_20(double* x, double a, int m, int n);

    inline double Phi_01(double* x, double a, int i);

    inline double Phi_02(double* x, double a, int i, int j);

    inline double Phi_03(double* x, double a, int i, int j, int k);

    inline double Phi_04(double* x, double a, int i, int j, int k, int l);

    inline double Phi_05(double* x, double a, int i, int j, int k, int l, int p);

    inline double Phi_06(double* x, double a, int i, int j, int k, int l, int p, int q);

    inline double Phi_11(double* x, double a, int n, int i);

    inline double Phi_12(double* x, double a, int n, int i, int j);

    inline double Phi_13(double* x, double a, int n, int i, int j, int k);

    inline double Phi_14(double* x, double a, int n, int i, int j, int k, int l);

    inline double Phi_15(double* x, double a, int n, int i, int j, int k, int l, int p);

    inline double Phi_16(double* x, double a, int n, int i, int j, int k, int l, int p, int q);

    inline double Phi_21(double* x, double a, int m, int n, int i);

    inline double Phi_22(double* x, double a, int m, int n, int i, int j);

    inline double Phi_23(double* x, double a, int m, int n, int i, int j, int k);

    inline double Phi_24(double* x, double a, int m, int n, int i, int j, int k, int l);

    inline double Phi_25(double* x, double a, int m, int n, int i, int j, int k, int l, int p);

    inline double Phi_26(double* x, double a, int m, int n, int i, int j, int k, int l, int p, int q);

    inline double Psi_03(double* x, double a, int i, int j, int k);

    inline double Psi_04(double* x, double a, int i, int j, int k, int l);

    inline double Psi_05(double* x, double a, int i, int j, int k, int l, int s);

    inline double Psi_06(double* x, double a, int i, int j, int k, int l, int s, int t);

    inline double Psi_07(double* x, double a, int i, int j, int k, int l, int s, int t, int p);

    inline double Psi_13(double* x, double a, int n, int i, int j, int k);

    inline double Psi_14(double* x, double a, int n, int i, int j, int k, int l);

    inline double Psi_15(double* x, double a, int n, int i, int j, int k, int l, int s);

    inline double Psi_16(double* x, double a, int n, int i, int j, int k, int l, int s, int q);

    inline double Psi_17(double* x, double a, int n, int i, int j, int k, int l, int s, int q, int p);

    inline double Psi_23(double* x, double a, int m, int n, int i, int j, int k);

    inline double Psi_24(double* x, double a, int m, int n, int i, int j, int k, int l);

    inline double Psi_25(double* x, double a, int m, int n, int i, int j, int k, int l, int s);

    inline double Psi_26(double* x, double a, int m, int n, int i, int j, int k, int l, int s, int t);

    inline double Psi_27(double* x, double a, int m, int n, int i, int j, int k, int l, int s, int t, int h);

    //////////////////////////////////////////////////////////////////////////////////////

    double R(double* x);

    double R_01(double* x, int i);

    double R_02(double* x, int i, int j);

    double R_03(double* x, int i, int j, int k);

    double R_04(double* x, int i, int j, int k, int l);

    double R_05(double* x, int i, int j, int k, int l, int s);

    double R_06(double* x, int i, int j, int k, int l, int s, int t);

    double R_07(double* x, int i, int j, int k, int l, int s, int t, int p);

    double fR(double* x, int sign);

    double fR_01(double* x, int i, int sign);

    double fR_02(double* x, int i, int j, int sign);

    double fR_03(double* x, int i, int j, int k, int sign);

    double fR_04(double* x, int i, int j, int k, int l, int sign);

    double fR_05(double* x, int i, int j, int k, int l, int s, int sign);

    double fR_06(double* x, int i, int j, int k, int l, int s, int t, int sign);

    double fR_07(double* x, int i, int j, int k, int l, int s, int t, int p, int sign);

    double lnR_01(double* x, int i, int sign);

    double lnR_02(double* x, int i, int j, int sign);

    double lnR_03(double* x, int i, int j, int k, int sign);

    double lnR_04(double* x, int i, int j, int k, int l, int sign);

    double lnR_05(double* x, int i, int j, int k, int l, int s, int sign);

    double lnR_06(double* x, int i, int j, int k, int l, int s, int t, int sign);

    double lnR_07(double* x, int i, int j, int k, int l, int s, int t, int p, int sign);


    /////////////////////////Constitutive//////////////////////////////////
    double Constitutive(int i, int j, int k, int l);
    ////////////////////////////////////////////////////////////////

    ////////////////////////Constitutive//////////////////////////////////
    double Constitutive1(int i, int j, int k, int l);
    ////////////////////////////////


    ///////////////////////////////////////////////////////////////////////////////////////////
    double D_40(int i, int j, int k, int l, int s, int h, Ref<VectorXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p);


    double D_50(int i, int j, int k, int l, int p, int s, int h, Ref<VectorXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p);


    double D_60(int i, int j, int k, int l, int p, int q, int s, int h, Ref<VectorXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p);


    /////////////////////////

    double D_41(int i, int j, int k, int l, int r, int s, int h, Ref<VectorXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p);
    //////////////////////////////////


    double D_51(int i, int j, int k, int l, int p, int r, int s, int h, Ref<VectorXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p);

    //////////////////////////////////

    double D_61(int i, int j, int k, int l, int p, int q, int r, int s, int h, Ref<VectorXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p);

    //////////////////////////////////////////////////////
    double D_42(int i, int j, int k, int l, int r, int z, int s, int h, Ref<VectorXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p);

    //////////////////////////////////////////////////////
    double D_52(int i, int j, int k, int l, int p, int r, int z, int s, int h, Ref<VectorXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p);

    //////////////////////////////////////////////////////
    double D_62(int i, int j, int k, int l, int p, int q, int r, int z, int s, int h, Ref<VectorXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p);


    double Green(double* x, int i, int j);

    double Green_01(double* x, int i, int j, int k);

    double Green_02(double* x, int i, int j, int k, int l);

    double Green_03(double* x, int i, int j, int k, int l, int s);

    double Green_04(double* x, int i, int j, int k, int l, int s, int t);

    double Green_05(double* x, int i, int j, int k, int l, int s, int t, int r);

    double IGreen_01(double* x, double a, int i, int j, int l);

    double IGreen_11(double* x, double a, int i, int j, int p, int l);

    double IGreen_21(double* x, double a, int i, int j, int p, int q, int l);

}

