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

class Uniform_Polygon {
public:
	double PHI(double b, double lm, double lp);
	double Phi_b(double b, double lm, double lp);
	double Phi_lm(double b, double lm);
	double Phi_lp(double b, double lp);

	double PHI_j(int j, double b, double lm, double lp, double* Vnorm, double* Vdir);

	// biharmonic potentials
	double Psi(double b, double lm, double lp);

	double Psi_b(double b, double lm, double lp);
	double Psi_lp(double b, double lp);
	double Psi_lm(double b, double lm);

	double PSI_j(int j, double b, double lm, double lp, double* Vnorm, double* Vdir);

	double Psi_2b(double b, double lm, double lp);
	double Psi_lp_b(double b, double lp);
	double Psi_lm_b(double b, double lm);
	double Psi_2lp(double b, double lp);
	double Psi_2lm(double b, double lm);
	double PSI_jk(int j, int k, double b, double lm, double lp, double* Vnorm, double* Vdir);

	// third potentials;
	double Psi_3b(double b, double lm, double lp);
	double Psi_3lp(double b, double lp);
	double Psi_3lm(double b, double lm);
	double Psi_2lp_b(double b, double lp);
	double Psi_2lm_b(double b, double lm);
	double Psi_lp_2b(double b, double lp);
	double Psi_lm_2b(double b, double lm);
	double PSI_jkl(int j, int k, int l, double b, double lm, double lp, double* Vnorm, double* Vdir);

	double PHI_o(double b, double lm, double lp);
	double PHI_1(int i, double b, double lm, double lp, double* Vnorm, double* Vdir);
	double PHI_2(int i, int j, double b, double lm, double lp, double* Vnorm, double* Vdir);

	double PSI_1(int i, double b, double lm, double lp, double* Vnorm, double* Vdir);
	double PSI_2(int i, int j, double b, double lm, double lp, double* Vnorm, double* Vdir);
	double PSI_3(int i, int j, int k, double b, double lm, double lp, double* Vnorm, double* Vdir);
	double PSI_4(int i, int j, int k, int l, double b, double lm, double lp, double* Vnorm, double* Vdir);

};
