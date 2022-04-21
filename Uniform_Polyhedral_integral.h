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

class Uniform_Polyhedral {
public:
    // Psi I1
    double Psi_I1(double a, double b, double le);

    double Psi_I1_a(double a, double b, double le);
    double Psi_I1_b(double a, double b, double le);
    double Psi_I1_le(double a, double b, double le);

    double Psi_I1_a_a(double a, double b, double le);
    double Psi_I1_a_b(double a, double b, double le);
    double Psi_I1_a_le(double a, double b, double le);
    double Psi_I1_b_b(double a, double b, double le);
    double Psi_I1_b_le(double a, double b, double le);
    double Psi_I1_le_le(double a, double b, double le);


    double Psi_I1_a_a_a(double a, double b, double le);
    double Psi_I1_a_a_b(double a, double b, double le);
    double Psi_I1_a_a_le(double a, double b, double le);
    double Psi_I1_a_b_b(double a, double b, double le);
    double Psi_I1_a_b_le(double a, double b, double le);
    double Psi_I1_a_le_le(double a, double b, double le);
    double Psi_I1_b_b_b(double a, double b, double le);
    double Psi_I1_b_b_le(double a, double b, double le);
    double Psi_I1_b_le_le(double a, double b, double le);
    double Psi_I1_le_le_le(double a, double b, double le);

    //Phi I1
    double Phi_u(double a, double b, double le);
    double Phi_u_II(double a, double b, double le);
    double Phi_u_II_a(double a, double b, double le);
    double Phi_u_II_b(double a, double b, double le);
    double Phi_u_II_le(double a, double b, double le);

    // define the combined coeffcients
    double I1(double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double I1_j(int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);

    double I4(double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double I4_j(int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double I4_jk(int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double I4_jkl(int j, int k, int l, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);

    // define the complete set of potentials
    double PHI(double a, double b, double lm, double lp);
    double PHI_i(int i, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double PHI_ij(int i, int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);

    double PSI_i(int i, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double PSI_ij(int i, int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double PSI_ijk(int i, int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
    double PSI_ijkl(int i, int j, int k, int l, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir);
};
