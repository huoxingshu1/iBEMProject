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
#include "Uniform_polygon_integral.h"
class Linear_Polygon {
public:

	// 
	Uniform_Polygon U_P;

	double PHI_p(double b, double lm, double lp);
	double Phi_p_b(double b, double lm, double lp);
	double Phi_p_lp(double b, double lp);
	double Phi_p_lm(double b, double lm);
	double PHI_p_j(int j, double b, double lm, double lp, double* Vnorm, double* Vdir);

	// biharmonic potentials
	double Psi_p_b(double b, double lm, double lp);
	double Psi_p_lm(double b, double lm);
	double Psi_p_lp(double b, double lp);
	double PSI_p_j(int j, double b, double lm, double lp, double* Vnorm, double* Vdir);

	double Psi_p_2b(double b, double lm, double lp);
	double Psi_p_lp_b(double b, double lp);
	double Psi_p_lm_b(double b, double lm);
	double Psi_p_2lp(double b, double lp);
	double Psi_p_2lm(double b, double lm);
	double PSI_p_jk(int j, int k, double b, double lm, double lp, double* Vnorm, double* Vdir);

	double Psi_p_3b(double b, double lm, double lp);
	double Psi_p_lp_2b(double b, double lp);
	double Psi_p_lm_2b(double b, double lm);
	double Psi_p_2lp_b(double b, double lp);
	double Psi_p_2lm_b(double b, double lm);
	double Psi_p_3lp(double b, double lp);
	double Psi_p_3lm(double b, double lm);
	double PSI_p_jkl(int j, int k, int l, double b, double lm, double lp, double* Vnorm, double* Vdir);

	// combined derivatives:
	double PSI_p_2(int p, int i, int j, double b, double lm, double lp, double* Vnorm, double* Vdir, double* x);
	double PSI_p_3(int p, int i, int j, int k, double b, double lm, double lp, double* Vnorm, double* Vdir, double* x);
	double PSI_p_4(int p, int i, int j, int k, int l, double b, double lm, double lp, double* Vnorm, double* Vdir, double* x);

	double PHI_p_o(int p, double b, double lm, double lp, double* Vnorm, double* Vdir, double* x);
	double PHI_p_1(int p, int i, double b, double lm, double lp, double* Vnorm, double* Vdir, double* x);
	double PHI_p_2(int p, int i, int j, double b, double lm, double lp, double* Vnorm, double* Vdir, double* x);


};
