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
#include "Uniform_Polyhedral_integral.h"

class Linear_Polyhedral {
public:

    Uniform_Polyhedral U_PL;
    // Psi_I2
    double Psi_I2(double a, double b, double le);

    // first derivative
    double Psi_I2_a(double a, double b, double le);
    double Psi_I2_b(double a, double b, double le);
    double Psi_I2_le(double a, double b, double le);


    double Psi_I2_a_a(double a, double b, double le);
    double Psi_I2_a_b(double a, double b, double le);
    double Psi_I2_a_le(double a, double b, double le);
    double Psi_I2_b_b(double a, double b, double le);
    double Psi_I2_b_le(double a, double b, double le);
    double Psi_I2_le_le(double a, double b, double le);

    double Psi_I2_a_a_a(double a, double b, double le);
    double Psi_I2_a_a_b(double a, double b, double le);
    double Psi_I2_a_a_le(double a, double b, double le);
    double Psi_I2_a_b_b(double a, double b, double le);
    double Psi_I2_a_b_le(double a, double b, double le);
    double Psi_I2_a_le_le(double a, double b, double le);
    double Psi_I2_b_b_b(double a, double b, double le);
    double Psi_I2_b_b_le(double a, double b, double le);
    double Psi_I2_b_le_le(double a, double b, double le);
    double Psi_I2_le_le_le(double a, double b, double le);

    // Psi_I3
    double Psi_I3(double a, double b, double le);

    double Psi_I3_a(double a, double b, double le);
    double Psi_I3_b(double a, double b, double le);
    double Psi_I3_le(double a, double b, double le);

    double Psi_I3_a_a(double a, double b, double le);
    double Psi_I3_a_b(double a, double b, double le);
    double Psi_I3_a_le(double a, double b, double le);
    double Psi_I3_b_b(double a, double b, double le);
    double Psi_I3_b_le(double a, double b, double le);
    double Psi_I3_le_le(double a, double b, double le);

    double Psi_I3_a_a_a(double a, double b, double le);
    double Psi_I3_a_a_b(double a, double b, double le);
    double Psi_I3_a_a_le(double a, double b, double le);
    double Psi_I3_a_b_b(double a, double b, double le);
    double Psi_I3_a_b_le(double a, double b, double le);
    double Psi_I3_a_le_le(double a, double b, double le);
    double Psi_I3_b_b_b(double a, double b, double le);
    double Psi_I3_b_b_le(double a, double b, double le);
    double Psi_I3_b_le_le(double a, double b, double le);
    double Psi_I3_le_le_le(double a, double b, double le);

    // the non-derivative form of Phi_p potential
    double Phi_Lin_H(double a, double b, double le);
    double Phi_Lin_pcos(double a, double b, double le);
    double Phi_Lin_psin(double a, double b, double le);


    // Phi_I2, I3

    double Phi_I2_III(double a, double b, double le);
    double Phi_I2_III_a(double a, double b, double le);
    double Phi_I2_III_b(double a, double b, double le);
    double Phi_I2_III_le(double a, double b, double le);
    double Phi_I3_III(double a, double b, double le);
    double Phi_I3_III_a(double a, double b, double le);
    double Phi_I3_III_b(double a, double b, double le);
    double Phi_I3_III_le(double a, double b, double le);


    // define the combined potentials
    double I2(double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double I3(double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double I2_j(int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double I3_j(int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);

    double I5_j(int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double I6_j(int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double I5_jk(int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double I6_jk(int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double I5_jkl(int i, int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double I6_jkl(int i, int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);


    // the complete set of potentials
    double PHI_p(int p, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x);
    double PHI_p_i(int p, int i, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x);
    double PHI_p_ij(int p, int i, int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x);

    double PSI_p_ij(int p, int i, int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x);
    double PSI_p_ijk(int p, int i, int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x);
    double PSI_p_ijkl(int p, int i, int j, int k, int l, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x);

    void PSI_p_ijk_combine(int i, int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x, double* result, \
        double& ij, double& ik, double& jk, double& kj, double& ijk, double& jk_5, double& jk_6);
    void PSI_p_ijkl_combine(int i, int j, int k, int l, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x, double* result, double& ijk, double& ikl, double& ijl, \
        double& jkl, double& kjl, double& ljk, double& ijkl, double& jkl_5, double& jkl_6);

};