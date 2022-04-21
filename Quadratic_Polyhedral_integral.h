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
#include "Linear_Polyhedral_integral.h"

class Quadratic_Polyhedral {
public:
    Uniform_Polyhedral U_PL; Linear_Polyhedral L_PL;

    // I4_second derivatives
    double Psi_I4_a_a(double a, double b, double le);
    double Psi_I4_a_b(double a, double b, double le);
    double Psi_I4_a_le(double a, double b, double le);
    double Psi_I4_b_b(double a, double b, double le);
    double Psi_I4_b_le(double a, double b, double le);
    double Psi_I4_le_le(double a, double b, double le);


    // I4_ third derivatives
    double Psi_I4_a_a_a(double a, double b, double le);
    double Psi_I4_a_a_b(double a, double b, double le);
    double Psi_I4_a_a_le(double a, double b, double le);
    double Psi_I4_a_b_b(double a, double b, double le);
    double Psi_I4_a_b_le(double a, double b, double le);
    double Psi_I4_a_le_le(double a, double b, double le);
    double Psi_I4_b_b_b(double a, double b, double le);
    double Psi_I4_b_b_le(double a, double b, double le);
    double Psi_I4_b_le_le(double a, double b, double le);
    double Psi_I4_le_le_le(double a, double b, double le);

    // I4_second derivatives
    double Psi_I5_a_a(double a, double b, double le);
    double Psi_I5_a_b(double a, double b, double le);
    double Psi_I5_a_le(double a, double b, double le);
    double Psi_I5_b_b(double a, double b, double le);
    double Psi_I5_b_le(double a, double b, double le);
    double Psi_I5_le_le(double a, double b, double le);

    // I5_ third derivatives
    double Psi_I5_a_a_a(double a, double b, double le);
    double Psi_I5_a_a_b(double a, double b, double le);
    double Psi_I5_a_a_le(double a, double b, double le);
    double Psi_I5_a_b_b(double a, double b, double le);
    double Psi_I5_a_b_le(double a, double b, double le);
    double Psi_I5_a_le_le(double a, double b, double le);
    double Psi_I5_b_b_b(double a, double b, double le);
    double Psi_I5_b_b_le(double a, double b, double le);
    double Psi_I5_b_le_le(double a, double b, double le);
    double Psi_I5_le_le_le(double a, double b, double le);

    // I6_second derivatives
    double Psi_I6_a_a(double a, double b, double le);
    double Psi_I6_a_b(double a, double b, double le);
    double Psi_I6_a_le(double a, double b, double le);
    double Psi_I6_b_b(double a, double b, double le);
    double Psi_I6_b_le(double a, double b, double le);
    double Psi_I6_le_le(double a, double b, double le);

    // I6_ third derivatives
    double Psi_I6_a_a_a(double a, double b, double le);
    double Psi_I6_a_a_b(double a, double b, double le);
    double Psi_I6_a_a_le(double a, double b, double le);
    double Psi_I6_a_b_b(double a, double b, double le);
    double Psi_I6_a_b_le(double a, double b, double le);
    double Psi_I6_a_le_le(double a, double b, double le);
    double Psi_I6_b_b_b(double a, double b, double le);
    double Psi_I6_b_b_le(double a, double b, double le);
    double Psi_I6_b_le_le(double a, double b, double le);
    double Psi_I6_le_le_le(double a, double b, double le);

    // derivatives of the Phi
    double Phi_I4_III(double a, double b, double le);
    double Phi_I5_III(double a, double b, double le);
    double Phi_I6_III(double a, double b, double le);

    double Phi_I4_III_a(double a, double b, double le);
    double Phi_I4_III_b(double a, double b, double le);
    double Phi_I4_III_le(double a, double b, double le);
    double Phi_I5_III_a(double a, double b, double le);
    double Phi_I5_III_b(double a, double b, double le);
    double Phi_I5_III_le(double a, double b, double le);
    double Phi_I6_III_a(double a, double b, double le);
    double Phi_I6_III_b(double a, double b, double le);
    double Phi_I6_III_le(double a, double b, double le);

    // add the combined potentials

    double I7(double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double I8(double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double I9(double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);

    double I7_j(int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double I8_j(int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double I9_j(int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);

    double I10_jk(int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double I11_jk(int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double I12_jk(int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);

    double I10_jkl(int i, int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double I11_jkl(int i, int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double I12_jkl(int i, int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);

    // complete set of potentials
    double PHI_pq_i(int p, int q, int i, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x);
    double PHI_pq_ij(int p, int q, int i, int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x);

    double PSI_pq_ijk(int p, int q, int i, int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x);
    double PSI_pq_ijkl(int p, int q, int i, int j, int k, int l, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x);

    void PSI_pq_ijk_combine(int i, int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, \
        double* x, double* linear_result, double* result, double ij, double ik, double jk, double kj, double ijk, double jk_5, double jk_6);

    void PSI_pq_ijkl_combine(int i, int j, int k, int l, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x, double* linear_result, double* result, double ijk, double ikl, double ijl, \
        double jkl, double kjl, double ljk, double ijkl, double jkl_5, double jkl_6);
    void PSI_pq_ijkl_combine_II(int i, int j, int k, int l, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x, double* linear_result, double* result);
};