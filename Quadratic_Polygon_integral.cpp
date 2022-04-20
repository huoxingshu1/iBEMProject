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

#include <cmath>
#include "Uniform_polygon_integral.h"
#include "Linear_Polygon_integral.h"
#include "Quadratic_Polygon_integral.h"

// second derivative of \Psi
double Quadratic_Polygon::Psi_pq_2b(double b, double lm, double lp)
{
	double result = (1.0 / 9.0) * (
		2.0 * pow(lm, 3.0) - 2.0 * pow(lp, 3.0) + 24.0 * b * b * (-lm + lp)
		+ 24.0 * pow(b, 3.0) * (atan(lm / b) - atan(lp / b)) - 3.0 * pow(lm, 3.0) * log(b * b + lm * lm)
		+ 3.0 * pow(lp, 3.0) * log(b * b + lp * lp)
		);
	if (b == 0) {
		result = (1.0 / 9.0) * (2.0 * pow(lm, 3.0) - 2.0 * pow(lp, 3.0) - 3.0 * pow(lm, 3.0) * log(lm * lm) + 3.0 * pow(lp, 3.0) * log(lp * lp));
	}
	return result;
}

double Quadratic_Polygon::Psi_pq_lp_b(double b, double lp)
{
	double result;
	if (b != 0) {
		result = b * lp * lp * log(b * b + lp * lp);
	}
	else {
		result = 0.0;
	}

	return result;
}

double Quadratic_Polygon::Psi_pq_lm_b(double b, double lm)
{
	double result;
	if (b != 0) {
		result = -b * lm * lm * log(b * b + lm * lm);
	}
	else {
		result = 0.0;
	}

	return result;
}

double Quadratic_Polygon::Psi_pq_2lp(double b, double lp)
{
	double result;

	if (b != 0) {
		result = -lp * (b * b + lp * lp - (b * b + 2.0 * lp * lp) * log(b * b + lp * lp));
	}
	else {
		result = -lp * lp * lp * (1.0 - 2.0 * log(lp * lp));
	}

	return result;
}

double Quadratic_Polygon::Psi_pq_2lm(double b, double lm)
{
	double result;

	if (b != 0) {
		result = lm * (b * b + lm * lm - (b * b + 2.0 * lm * lm) * log(b * b + lm * lm));
	}
	else {
		result = lm * lm * lm * (1.0 - 2.0 * log(lm * lm));
	}

	return result;
}




// third derivative of \Psi
double Quadratic_Polygon::Psi_pq_3b(double b, double lm, double lp)
{
	double para = (2.0 * b) / ((pow(b, 2.0) + pow(lm, 2.0)) * (pow(b, 2.0) + pow(lp, 2.0)));
	double result = -4.0 * pow(b, 4.0) * lm - 3.0 * pow(b, 2.0) * pow(lm, 3.0) + 4.0 * pow(b, 4.0) * lp + 4.0 * pow(b, 2.0) * pow(lm, 2.0) * lp
		- 4.0 * pow(b, 2.0) * pow(lp, 2.0) * lm - 3.0 * pow(lm, 3.0) * pow(lp, 2.0) + 3.0 * pow(b, 2.0) * pow(lp, 3.0)
		+ 3.0 * pow(lm, 2.0) * pow(lp, 3.0) + 4.0 * b * (pow(b, 2.0) + pow(lm, 2.0)) * (pow(b, 2.0) + pow(lp, 2.0)) * (atan(lm / b) - atan(lp / b));

	if (b == 0) result = -4.0 * pow(b, 4.0) * lm - 3.0 * pow(b, 2.0) * pow(lm, 3.0) + 4.0 * pow(b, 4.0) * lp + 4.0 * pow(b, 2.0) * pow(lm, 2.0) * lp
		- 4.0 * pow(b, 2.0) * pow(lp, 2.0) * lm - 3.0 * pow(lm, 3.0) * pow(lp, 2.0) + 3.0 * pow(b, 2.0) * pow(lp, 3.0)
		+ 3.0 * pow(lm, 2.0) * pow(lp, 3.0);
	result *= para;
	return result;
}

double Quadratic_Polygon::Psi_pq_2lp_b(double b, double lp)
{
	double result = 2.0 * b * lp * (log(b * b + lp * lp) + (lp * lp) / (b * b + lp * lp));
	return result;
}

double Quadratic_Polygon::Psi_pq_2lm_b(double b, double lm)
{
	double result = -2.0 * b * lm * (log(b * b + lm * lm) + (lm * lm) / (b * b + lm * lm));
	return result;
}

double Quadratic_Polygon::Psi_pq_lp_2b(double b, double lp)
{
	double result = pow(lp, 2.0) * (log(b * b + lp * lp) + (2.0 * b * b) / (b * b + lp * lp));
	return result;
}

double Quadratic_Polygon::Psi_pq_lm_2b(double b, double lm)
{
	double result = -pow(lm, 2.0) * (log(b * b + lm * lm) + (2.0 * b * b) / (b * b + lm * lm));
	return result;
}

double Quadratic_Polygon::Psi_pq_3lp(double b, double lp)
{
	double result = -pow(b, 4.0) - 2.0 * pow(b, 2.0) * pow(lp, 2.0) + pow(lp, 4.0) + (pow(b, 4.0) + 7.0 * pow(b * lp, 2.0) + 6.0 * pow(lp, 4.0)) * log(b * b + lp * lp);
	result *= 1.0 / (b * b + lp * lp);
	return result;
}

double Quadratic_Polygon::Psi_pq_3lm(double b, double lm)
{
	double result = pow(b, 4.0) + 2.0 * pow(b, 2.0) * pow(lm, 2.0) - pow(lm, 4.0) - (pow(b, 4.0) + 7.0 * pow(b * lm, 2.0) + 6.0 * pow(lm, 4.0)) * log(b * b + lm * lm);
	result *= 1.0 / (b * b + lm * lm);
	return result;
}

/*----------------------------------------
	Harmonic potentials
------------------------------------------*/
double Quadratic_Polygon::Phi_pq(double b, double lm, double lp)
{
	double result = 0.0;
	if (b != 0) {
		result = 1.0 / 9.0 * (
			-6.0 * b * b * lm + 2.0 * lm * lm * lm + 6.0 * b * b * lp - 2.0 * lp * lp * lp
			+ 6.0 * b * b * b * atan(lm / b) - 6.0 * b * b * b * atan(lp / b)
			- 3.0 * lm * lm * lm * log(b * b + lm * lm) + 3.0 * lp * lp * lp * log(b * b + lp * lp)
			);
	}
	else {
		result = 1.0 / 9.0 * (
			2.0 * lm * lm * lm - 2.0 * lp * lp * lp - 3.0 * lm * lm * lm * log(lm * lm) + 3.0 * lp * lp * lp * log(lp * lp)
			);
	}
	return result;
}


double Quadratic_Polygon::Phi_pq_b(double b, double lm, double lp)
{
	double result = 2.0 * b * (-lm + lp + b * atan(lm / b) - b * atan(lp / b));
	if (b == 0) result = 0;
	return result;
}

double Quadratic_Polygon::Phi_pq_lp(double b, double lp)
{
	double result = pow(lp, 2.0) * log(b * b + lp * lp);
	return result;
}

double Quadratic_Polygon::Phi_pq_lm(double b, double lm)
{
	double result = -pow(lm, 2.0) * log(b * b + lm * lm);
	return result;
}

/*--------------------------------------------
		Combined potentials of functions
		for derivatives
---------------------------------------------*/

double Quadratic_Polygon::PSI_pq_jk(int j, int k, double b, double lm, double lp, double* Vnorm, double* Vdir)
{
	double result = 0.0;

	result = Vnorm[j] * Vnorm[k] * Psi_pq_2b(b, lm, lp)
		+ (Psi_pq_lm_b(b, lm) + Psi_pq_lp_b(b, lp)) * (Vnorm[j] * Vdir[k] + Vnorm[k] * Vdir[j])
		+ (Psi_pq_2lm(b, lm) + Psi_pq_2lp(b, lp)) * Vdir[j] * Vdir[k];

	return result;
}

double Quadratic_Polygon::PSI_pq_jkl(int j, int k, int l, double b, double lm, double lp, double* Vnorm, double* Vdir)
{
	double result = 0.0;

	result = -Psi_pq_3b(b, lm, lp) * Vnorm[j] * Vnorm[k] * Vnorm[l] - (Psi_pq_lm_2b(b, lm) + Psi_pq_lp_2b(b, lp)) * (Vnorm[j] * Vnorm[k] * Vdir[l] + Vnorm[j] * Vdir[k] * Vnorm[l] + Vdir[j] * Vnorm[k] * Vnorm[l]);
	result += -(Psi_pq_2lp_b(b, lp) + Psi_pq_2lm_b(b, lm)) * (Vnorm[j] * Vdir[k] * Vdir[l] + Vdir[j] * Vnorm[k] * Vdir[l] + Vdir[j] * Vdir[k] * Vnorm[l]);
	result += -(Psi_pq_3lp(b, lp) + Psi_pq_3lm(b, lm)) * Vdir[j] * Vdir[k] * Vdir[l];

	return result;
}

double Quadratic_Polygon::PHI_pq_j(int j, double b, double lm, double lp, double* Vnorm, double* Vdir)
{
	double result = 0.0;

	result = -Vnorm[j] * Phi_pq_b(b, lm, lp) - Vdir[j] * (Phi_pq_lp(b, lp) + Phi_pq_lm(b, lm));

	return result;
}

double Quadratic_Polygon::PSI_pq_3(int p, int q, int i, int j, int k, double b, double lm, double lp, double* Vnorm, double* Vdir, double* x)
{
	double result = 0.0; double d[2][2] = { 1, 0, 0, 1 };

	// Uniform part
	result = Vdir[p] * Vdir[q] * Vnorm[i] * PSI_pq_jk(j, k, b, lm, lp, Vnorm, Vdir);
	result += Vnorm[i] * Vnorm[p] * Vnorm[q] * (2.0 * Vnorm[j] * Vnorm[k] * U_P.Psi(b, lm, lp) - 2.0 * b * Vnorm[j] * U_P.PSI_j(k, b, lm, lp, Vnorm, Vdir) - 2.0 * b * Vnorm[k] * U_P.PSI_j(j, b, lm, lp, Vnorm, Vdir) + b * b * U_P.PSI_jk(j, k, b, lm, lp, Vnorm, Vdir));
	result += (Vnorm[i] * Vnorm[p] * Vdir[q] + Vnorm[i] * Vnorm[q] * Vdir[p]) * (-Vnorm[j] * L_P.PSI_p_j(k, b, lm, lp, Vnorm, Vdir) - Vnorm[k] * L_P.PSI_p_j(j, b, lm, lp, Vnorm, Vdir) + b * L_P.PSI_p_jk(j, k, b, lm, lp, Vnorm, Vdir));

	result += (d[p][i]) * L_P.PSI_p_2(q, j, k, b, lm, lp, Vnorm, Vdir, x);
	result += (d[p][j]) * L_P.PSI_p_2(q, i, k, b, lm, lp, Vnorm, Vdir, x);
	result += (d[p][k]) * L_P.PSI_p_2(q, i, j, b, lm, lp, Vnorm, Vdir, x);
	result += (x[p]) * L_P.PSI_p_3(q, i, j, k, b, lm, lp, Vnorm, Vdir, x);

	result += (d[q][i]) * L_P.PSI_p_2(p, j, k, b, lm, lp, Vnorm, Vdir, x);
	result += (d[q][j]) * L_P.PSI_p_2(p, i, k, b, lm, lp, Vnorm, Vdir, x);
	result += (d[q][k]) * L_P.PSI_p_2(p, i, j, b, lm, lp, Vnorm, Vdir, x);
	result += (x[q]) * L_P.PSI_p_3(p, i, j, k, b, lm, lp, Vnorm, Vdir, x);

	result += -(d[p][i] * d[q][j] + d[q][i] * d[p][j]) * Vnorm[k] * U_P.Psi(b, lm, lp);
	result += -(d[p][i] * d[q][k] + d[q][i] * d[p][k]) * Vnorm[j] * U_P.Psi(b, lm, lp);
	result += -(d[p][i] * x[q] + d[q][i] * x[p]) * Vnorm[j] * U_P.PSI_j(k, b, lm, lp, Vnorm, Vdir);
	result += -(d[p][j] * d[q][k] + d[q][j] * d[p][k]) * Vnorm[i] * U_P.Psi(b, lm, lp);
	result += -(d[p][j] * x[q] + d[q][j] * x[p]) * Vnorm[i] * U_P.PSI_j(k, b, lm, lp, Vnorm, Vdir);
	result += -(d[p][k] * x[q] + d[q][k] * x[p]) * Vnorm[i] * U_P.PSI_j(j, b, lm, lp, Vnorm, Vdir);
	result += -x[p] * x[q] * Vnorm[i] * U_P.PSI_jk(j, k, b, lm, lp, Vnorm, Vdir);

	return result;
}

double Quadratic_Polygon::PSI_pq_4(int p, int q, int i, int j, int k, int l, double b, double lm, double lp, double* Vnorm, double* Vdir, double* x)
{
	double result = 0.0;
	double d[2][2] = { 1,0,0,1 };

	//Uniform_polygon Uni; Linear_polygon Lin;

	result = Vdir[p] * Vdir[q] * Vnorm[i] * PSI_pq_jkl(j, k, l, b, lm, lp, Vnorm, Vdir);
	result += Vnorm[i] * Vnorm[p] * Vnorm[q] * (2.0 * Vnorm[j] * Vnorm[k] * U_P.PSI_j(l, b, lm, lp, Vnorm, Vdir) + 2.0 * Vnorm[j] * Vnorm[l] * U_P.PSI_j(k, b, lm, lp, Vnorm, Vdir) - 2.0 * Vnorm[j] * b * U_P.PSI_jk(k, l, b, lm, lp, Vnorm, Vdir));
	result += Vnorm[i] * Vnorm[p] * Vnorm[q] * (2.0 * Vnorm[k] * Vnorm[l] * U_P.PSI_j(j, b, lm, lp, Vnorm, Vdir) - 2.0 * Vnorm[k] * b * U_P.PSI_jk(j, l, b, lm, lp, Vnorm, Vdir) - 2.0 * b * Vnorm[l] * U_P.PSI_jk(j, k, b, lm, lp, Vnorm, Vdir) + b * b * U_P.PSI_jkl(j, k, l, b, lm, lp, Vnorm, Vdir));
	result += (Vnorm[i] * Vnorm[p] * Vdir[q] + Vnorm[i] * Vnorm[q] * Vdir[p]) * (-Vnorm[j] * L_P.PSI_p_jk(k, l, b, lm, lp, Vnorm, Vdir) - Vnorm[k] * L_P.PSI_p_jk(j, l, b, lm, lp, Vnorm, Vdir) - Vnorm[l] * L_P.PSI_p_jk(j, k, b, lm, lp, Vnorm, Vdir) + b * L_P.PSI_p_jkl(j, k, l, b, lm, lp, Vnorm, Vdir));

	///
	result += (d[p][i]) * L_P.PSI_p_3(q, j, k, l, b, lm, lp, Vnorm, Vdir, x);
	result += (d[p][j]) * L_P.PSI_p_3(q, i, k, l, b, lm, lp, Vnorm, Vdir, x);
	result += (d[p][k]) * L_P.PSI_p_3(q, i, j, l, b, lm, lp, Vnorm, Vdir, x);
	result += (d[p][l]) * L_P.PSI_p_3(q, i, j, k, b, lm, lp, Vnorm, Vdir, x);
	result += (x[p]) * L_P.PSI_p_4(q, i, j, k, l, b, lm, lp, Vnorm, Vdir, x);

	////

	result += (d[q][i]) * L_P.PSI_p_3(p, j, k, l, b, lm, lp, Vnorm, Vdir, x);
	result += (d[q][j]) * L_P.PSI_p_3(p, i, k, l, b, lm, lp, Vnorm, Vdir, x);
	result += (d[q][k]) * L_P.PSI_p_3(p, i, j, l, b, lm, lp, Vnorm, Vdir, x);
	result += (d[q][l]) * L_P.PSI_p_3(p, i, j, k, b, lm, lp, Vnorm, Vdir, x);
	result += (x[q]) * L_P.PSI_p_4(p, i, j, k, l, b, lm, lp, Vnorm, Vdir, x);

	////

	result += -(d[p][i] * d[q][j] + d[q][i] * d[p][j]) * Vnorm[k] * U_P.PSI_j(l, b, lm, lp, Vnorm, Vdir);
	result += -(d[p][i] * d[q][k] + d[q][i] * d[p][k]) * Vnorm[j] * U_P.PSI_j(l, b, lm, lp, Vnorm, Vdir);
	result += -(d[p][i] * d[q][l] + d[q][i] * d[p][l]) * Vnorm[j] * U_P.PSI_j(k, b, lm, lp, Vnorm, Vdir);
	result += -(d[p][i] * x[q] + x[p] * d[q][i]) * Vnorm[j] * U_P.PSI_jk(k, l, b, lm, lp, Vnorm, Vdir);
	result += -(d[p][j] * d[q][k] + d[q][j] * d[p][k]) * Vnorm[i] * U_P.PSI_j(l, b, lm, lp, Vnorm, Vdir);
	result += -(d[p][j] * d[q][l] + d[q][j] * d[p][l]) * Vnorm[i] * U_P.PSI_j(k, b, lm, lp, Vnorm, Vdir);
	result += -(d[p][j] * x[q] + x[p] * d[q][j]) * Vnorm[i] * U_P.PSI_jk(k, l, b, lm, lp, Vnorm, Vdir);
	result += -(d[p][k] * d[q][l] + d[q][k] * d[p][l]) * Vnorm[i] * U_P.PSI_j(j, b, lm, lp, Vnorm, Vdir);
	result += -(d[p][k] * x[q] + x[p] * d[q][k]) * Vnorm[i] * U_P.PSI_jk(j, l, b, lm, lp, Vnorm, Vdir);
	result += -(d[p][l] * x[q] + x[p] * d[q][l]) * Vnorm[i] * U_P.PSI_jk(j, k, b, lm, lp, Vnorm, Vdir);
	result += -(x[p] * x[q]) * Vnorm[i] * U_P.PSI_jkl(j, k, l, b, lm, lp, Vnorm, Vdir);

	return result;


}

double Quadratic_Polygon::PHI_pq_1(int p, int q, int i, double b, double lm, double lp, double* Vnorm, double* Vdir, double* x)
{
	double result = 0.0;
	double d[2][2] = { 1,0,0,1 };

	result = Vdir[p] * Vdir[q] * Vnorm[i] * Phi_pq(b, lm, lp);
	result += Vnorm[i] * Vnorm[p] * Vnorm[q] * b * b * U_P.PHI(b, lm, lp);
	result += Vnorm[i] * (Vnorm[p] * Vdir[q] + Vnorm[q] * Vdir[p]) * b * L_P.PHI_p(b, lm, lp);

	///
	result += (d[p][i]) * L_P.PHI_p_o(q, b, lm, lp, Vnorm, Vdir, x);
	result += x[p] * L_P.PHI_p_1(q, i, b, lm, lp, Vnorm, Vdir, x);

	result += (d[q][i]) * L_P.PHI_p_o(p, b, lm, lp, Vnorm, Vdir, x);
	result += (x[q]) * L_P.PHI_p_1(p, i, b, lm, lp, Vnorm, Vdir, x);

	result = result - (d[p][i] * x[q] + d[q][i] * x[p]) * U_P.PHI_o(b, lm, lp) - x[p] * x[q] * U_P.PHI_1(i, b, lm, lp, Vnorm, Vdir);

	return result;
}

double Quadratic_Polygon::PHI_pq_2(int p, int q, int i, int j, double b, double lm, double lp, double* Vnorm, double* Vdir, double* x)
{
	//Uniform_polygon Uni; Linear_polygon Lin;

	double result = 0.0;
	double d[2][2] = { 1,0,0,1 };

	result = Vdir[p] * Vdir[q] * Vnorm[i] * PHI_pq_j(j, b, lm, lp, Vnorm, Vdir);
	result += Vnorm[i] * Vnorm[p] * Vnorm[q] * (-2.0 * b * Vnorm[j] * U_P.PHI(b, lm, lp) + b * b * U_P.PHI_j(j, b, lm, lp, Vnorm, Vdir));
	result += (Vnorm[i] * Vnorm[p] * Vdir[q] + Vnorm[i] * Vnorm[q] * Vdir[p]) * (-Vnorm[j] * L_P.PHI_p(b, lm, lp) + b * L_P.PHI_p_j(j, b, lm, lp, Vnorm, Vdir));

	///
	result += (d[p][i]) * L_P.PHI_p_1(q, j, b, lm, lp, Vnorm, Vdir, x);
	result += (d[p][j]) * L_P.PHI_p_1(q, i, b, lm, lp, Vnorm, Vdir, x);
	result += (x[p]) * L_P.PHI_p_2(q, i, j, b, lm, lp, Vnorm, Vdir, x);

	///

	result += (d[q][i]) * L_P.PHI_p_1(p, j, b, lm, lp, Vnorm, Vdir, x);
	result += (d[q][j]) * L_P.PHI_p_1(p, i, b, lm, lp, Vnorm, Vdir, x);
	result += (x[q]) * L_P.PHI_p_2(p, i, j, b, lm, lp, Vnorm, Vdir, x);

	////

	result += -(d[p][i] * d[q][j] + d[p][j] * d[q][i]) * U_P.PHI_o(b, lm, lp);
	result += -(d[p][i] * x[q] + x[p] * d[q][i]) * Vnorm[j] * U_P.PHI(b, lm, lp);
	result += -(d[p][j] * x[q] + x[p] * d[q][j]) * Vnorm[i] * U_P.PHI(b, lm, lp);
	result += -(x[p] * x[q]) * Vnorm[i] * U_P.PHI_j(j, b, lm, lp, Vnorm, Vdir);

	return result;



}