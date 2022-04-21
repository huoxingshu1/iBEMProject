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

#include "integratorElastic3DHelper.h"
#include "configElastic3D.h"
#include <math.h>
#include"Eigen/Dense"

using namespace Eigen;
////////////////////////////////////////////////////////////new subroutine//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//double mu_0 = 1.0;
//double nu0 = 1.0;
//double mu_1 = 3.0;
//double nu1 = 3.0;

////////////////////////////////////////////////////////////new subroutine//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline double I_n0(double lambda, double a, int n)
{
	double I_n = 4 * Pi * a * a * a / ((2 * n + 1) * pow((a * a + lambda), (n + 0.5)));
	return I_n;
}

inline double I_n1(double lambda, double a, double* x, int n, int p)
{
	double delta = sqrt((a * a + lambda) * (a * a + lambda) * (a * a + lambda));
	double I_n1;
	if (lambda == 0)
		I_n1 = 0.0;
	else
		I_n1 = -2 * Pi * a * a * a / (pow((a * a + lambda), n) * delta) * 2 * x[p];

	return I_n1;
}

inline double I_n2(double lambda, double a, double* x, int n, int p, int q)
{
	double delta = sqrt((a * a + lambda) * (a * a + lambda) * (a * a + lambda));
	double Z1 = (double)(n + 1.5) / (a * a + lambda);

	double I_n2;
	if (lambda == 0)
		I_n2 = 0.0;
	else
		I_n2 = -2 * Pi * a * a * a / (pow((a * a + lambda), n) * delta) * (2 * d[p][q] - 2 * x[p] * 2 * x[q] * Z1);

	return I_n2;

}

inline double I_n3(double lambda, double a, double* x, int n, int p, int q, int r)
{
	double delta = sqrt((a * a + lambda) * (a * a + lambda) * (a * a + lambda));
	double Z1 = (double)(n + 1.5) / (a * a + lambda);
	double Z2 = (double)(n + 1.5) / ((a * a + lambda) * (a * a + lambda));

	double I_n3;
	if (lambda == 0)
		I_n3 = 0.0;
	else
		I_n3 = -2 * Pi * a * a * a / (pow((a * a + lambda), n) * delta) * (-(2 * d[p][q] * 2 * x[r] + 2 * d[p][r] * 2 * x[q] + 2 * d[q][r] * 2 * x[p]) * Z1 + 2 * x[p] * 2 * x[q] * 2 * x[r] * (Z2 + Z1 * Z1));
	return I_n3;

}

inline double I_n4(double lambda, double a, double* x, int n, int p, int q, int r, int t)
{
	double delta = sqrt((a * a + lambda) * (a * a + lambda) * (a * a + lambda));
	double Z1 = (double)(n + 1.5) / (a * a + lambda);
	double Z2 = (double)(n + 1.5) / ((a * a + lambda) * (a * a + lambda));
	double Z3 = (double)(n + 1.5) / ((a * a + lambda) * (a * a + lambda) * (a * a + lambda));

	double I_n4;
	if (lambda == 0)
		I_n4 = 0.0;
	else
		I_n4 = -2 * Pi * a * a * a / (pow((a * a + lambda), n) * delta) * (0 - (0 + 0 + 0 + 0 + 2 * d[p][q] * 2 * d[r][t] + 2 * d[p][r] * 2 * d[q][t] + 2 * d[q][r] * 2 * d[p][t]) * Z1
			+ (2 * d[p][q] * 2 * x[r] * 2 * x[t] + 2 * d[p][r] * 2 * x[q] * 2 * x[t] + 2 * d[q][r] * 2 * x[p] * 2 * x[t] + 2 * d[p][t] * 2 * x[q] * 2 * x[r] + 2 * d[q][t] * 2 * x[p] * 2 * x[r] + 2 * d[r][t] * 2 * x[p] * 2 * x[q]) * (Z2 + Z1 * Z1)
			- 2 * x[p] * 2 * x[q] * 2 * x[r] * 2 * x[t] * (2 * Z3 + 3 * Z1 * Z2 + Z1 * Z1 * Z1));
	return I_n4;
}

inline double I_n5(double lambda, double a, double* x, int n, int p, int q, int r, int t, int h)
{
	double delta = sqrt((a * a + lambda) * (a * a + lambda) * (a * a + lambda));
	double Z1 = (double)(n + 1.5) / (a * a + lambda);
	double Z2 = (double)(n + 1.5) / ((a * a + lambda) * (a * a + lambda));
	double Z3 = (double)(n + 1.5) / ((a * a + lambda) * (a * a + lambda) * (a * a + lambda));
	double Z4 = (double)(n + 1.5) / ((a * a + lambda) * (a * a + lambda) * (a * a + lambda) * (a * a + lambda));

	double I_n5;
	if (lambda == 0)
		I_n5 = 0.0;
	else
		I_n5 = -2 * x[h] * Z1 * I_n4(lambda, a, x, n, p, q, r, t) - 2 * Pi * a * a * a / (pow((a * a + lambda), n) * delta) * (

			0 + (0 + 0 + 0 + 0 + 2 * d[p][q] * 2 * d[r][t] + 2 * d[p][r] * 2 * d[q][t] + 2 * d[q][r] * 2 * d[p][t]) * 2 * x[h] * Z2

			+ (2 * d[p][q] * (2 * d[r][h] * 2 * x[t] + 2 * d[t][h] * 2 * x[r]) + 2 * d[p][r] * (2 * d[q][h] * 2 * x[t] + 2 * d[t][h] * 2 * x[q]) + 2 * d[q][r] * (2 * d[p][h] * 2 * x[t] + 2 * d[t][h] * 2 * x[p])
				+ 2 * d[p][t] * (2 * d[q][h] * 2 * x[r] + 2 * d[r][h] * 2 * x[q]) + 2 * d[q][t] * (2 * d[p][h] * 2 * x[r] + 2 * d[r][h] * 2 * x[p]) + 2 * d[r][t] * (2 * d[p][h] * 2 * x[q] + 2 * d[q][h] * 2 * x[p])) * (Z2 + Z1 * Z1)

			+ (2 * d[p][q] * 2 * x[r] * 2 * x[t] + 2 * d[p][r] * 2 * x[q] * 2 * x[t] + 2 * d[q][r] * 2 * x[p] * 2 * x[t] + 2 * d[p][t] * 2 * x[q] * 2 * x[r] + 2 * d[q][t] * 2 * x[p] * 2 * x[r] + 2 * d[r][t] * 2 * x[p] * 2 * x[q]) * (-2 * 2 * x[h] * Z3 - 2 * Z1 * 2 * x[h] * Z2)

			- (2 * d[p][h] * 2 * x[q] * 2 * x[r] * 2 * x[t] + 2 * d[q][h] * 2 * x[p] * 2 * x[r] * 2 * x[t] + 2 * d[r][h] * 2 * x[p] * 2 * x[q] * 2 * x[t] + 2 * d[t][h] * 2 * x[p] * 2 * x[q] * 2 * x[r]) * (2 * Z3 + 3 * Z1 * Z2 + Z1 * Z1 * Z1)
			- 2 * x[p] * 2 * x[q] * 2 * x[r] * 2 * x[t] * (0 - 6 * 2 * x[h] * Z4 + 3 * (0 - 1 * 2 * x[h] * Z2) * Z2 + 3 * Z1 * (0 - 2 * 2 * x[h] * Z3) + 3 * Z1 * Z1 * (0 - 1 * 2 * x[h] * Z2))

			);
	return I_n5;
}
////////////////////////////////////////////////////////////
inline double V_00(double* x, double a)
{
	double lambda;
	if ((x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - a * a) > 0)
		lambda = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - a * a;
	else
		lambda = 0;
	double V = 0.5 * (I_n0(lambda, a, 0) - x[0] * x[0] * I_n0(lambda, a, 1) - x[1] * x[1] * I_n0(lambda, a, 1) - x[2] * x[2] * I_n0(lambda, a, 1));
	return V;
}

inline double V_10(double* x, double a, int i)
{
	double lambda;
	if ((x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - a * a) > 0)
		lambda = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - a * a;
	else
		lambda = 0;

	double V_i = 0.5 * (I_n0(lambda, a, 1) - x[0] * x[0] * I_n0(lambda, a, 2) - x[1] * x[1] * I_n0(lambda, a, 2) - x[2] * x[2] * I_n0(lambda, a, 2));
	return V_i;
}

inline double V_20(double* x, double a, int i, int j)
{
	double lambda;
	if ((x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - a * a) > 0)
		lambda = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - a * a;
	else
		lambda = 0;

	double V_ij = 0.5 * (I_n0(lambda, a, 2) - x[0] * x[0] * I_n0(lambda, a, 3) - x[1] * x[1] * I_n0(lambda, a, 3) - x[2] * x[2] * I_n0(lambda, a, 3));
	return V_ij;
}

inline double V_n0(double* x, double a, int n)
{
	double lambda;
	if ((x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - a * a) > 0)
		lambda = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - a * a;
	else
		lambda = 0;

	double V_ijk = 0.5 * (I_n0(lambda, a, n) - x[0] * x[0] * I_n0(lambda, a, n + 1) - x[1] * x[1] * I_n0(lambda, a, n + 1) - x[2] * x[2] * I_n0(lambda, a, n + 1));
	return V_ijk;
}


inline double V_n1(double* x, double a, int n, int p)
{
	double lambda;
	if ((x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - a * a) > 0)
		lambda = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - a * a;
	else
		lambda = 0;

	double V_i_p = -x[p] * I_n0(lambda, a, n + 1);
	return V_i_p;

}

inline double V_n2(double* x, double a, int n, int p, int q)
{
	double lambda;
	if ((x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - a * a) > 0)
		lambda = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - a * a;
	else
		lambda = 0;

	double V_i_pq = -(d[p][q] * I_n0(lambda, a, n + 1) + x[p] * I_n1(lambda, a, x, n + 1, q));
	return V_i_pq;
}

inline double V_n3(double* x, double a, int n, int p, int q, int r)
{
	double lambda;
	if ((x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - a * a) > 0)
		lambda = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - a * a;
	else
		lambda = 0;

	double V_i_pqr = -(d[p][q] * I_n1(lambda, a, x, n + 1, r) + d[p][r] * I_n1(lambda, a, x, n + 1, q) + x[p] * I_n2(lambda, a, x, n + 1, q, r));
	return V_i_pqr;

}

inline double V_n4(double* x, double a, int n, int p, int q, int r, int t)
{
	double lambda;
	if ((x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - a * a) > 0)
		lambda = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - a * a;
	else
		lambda = 0;

	double V_i_pqrt = -(d[p][q] * I_n2(lambda, a, x, n + 1, r, t) + d[p][r] * I_n2(lambda, a, x, n + 1, q, t) + d[p][t] * I_n2(lambda, a, x, n + 1, q, r) + x[p] * I_n3(lambda, a, x, n + 1, q, r, t));
	return V_i_pqrt;

}

inline double V_n5(double* x, double a, int n, int p, int q, int r, int t, int s)
{
	double lambda;
	if ((x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - a * a) > 0)
		lambda = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - a * a;
	else
		lambda = 0;

	double V_i_pqrts = -(d[p][q] * I_n3(lambda, a, x, n + 1, r, t, s) + d[p][r] * I_n3(lambda, a, x, n + 1, q, t, s) + d[p][t] * I_n3(lambda, a, x, n + 1, q, r, s) + d[p][s] * I_n3(lambda, a, x, n + 1, q, r, t) + x[p] * I_n4(lambda, a, x, n + 1, q, r, t, s));
	return V_i_pqrts;

}

inline double V_n6(double* x, double a, int n, int p, int q, int r, int t, int s, int h)
{
	double lambda;
	if ((x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - a * a) > 0)
		lambda = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - a * a;
	else
		lambda = 0;

	double V_i_pqrtsh = -(d[p][q] * I_n4(lambda, a, x, n + 1, r, t, s, h) + d[p][r] * I_n4(lambda, a, x, n + 1, q, t, s, h) + d[p][t] * I_n4(lambda, a, x, n + 1, q, r, s, h) + d[p][s] * I_n4(lambda, a, x, n + 1, q, r, t, h) + d[p][h] * I_n4(lambda, a, x, n + 1, q, r, t, s) + x[p] * I_n5(lambda, a, x, n + 1, q, r, t, s, h));
	return V_i_pqrtsh;
}

////////////////////////////////////////////////////////////

inline double Phi_00(double* x, double a)
{
	double Phi = V_00(x, a);
	return Phi;
}

inline double Phi_10(double* x, double a, int n)
{
	double Phi_n = a * a * x[n] * V_10(x, a, n);
	return Phi_n;
}

inline double Phi_20(double* x, double a, int m, int n)
{
	double Phi_mn = a * a * (x[m] * x[n] * a * a * V_20(x, a, m, n) + 0.25 * d[m][n] * (V_00(x, a) - x[1] * x[1] * V_10(x, a, 1) - x[2] * x[2] * V_10(x, a, 2) - x[3] * x[3] * V_10(x, a, 3)
		- a * a * (V_10(x, a, m) - x[1] * x[1] * V_20(x, a, 1, m) - x[2] * x[2] * V_20(x, a, 2, m) - x[3] * x[3] * V_20(x, a, 3, m))));
	return Phi_mn;
}

inline double Phi_01(double* x, double a, int i)
{
	double Phi_i = V_n1(x, a, 0, i);
	return Phi_i;
}

inline double Phi_02(double* x, double a, int i, int j)
{
	double Phi_ij = V_n2(x, a, 0, i, j);
	return Phi_ij;
}

inline double Phi_03(double* x, double a, int i, int j, int k)
{
	double Phi_ijk = V_n3(x, a, 0, i, j, k);
	return Phi_ijk;
}

inline double Phi_04(double* x, double a, int i, int j, int k, int l)
{
	double Phi_ijkl = V_n4(x, a, 0, i, j, k, l);
	return Phi_ijkl;
}

inline double Phi_05(double* x, double a, int i, int j, int k, int l, int p)
{
	double Phi_ijklp = V_n5(x, a, 0, i, j, k, l, p);
	return Phi_ijklp;
}

inline double Phi_06(double* x, double a, int i, int j, int k, int l, int p, int q)
{
	double Phi_ijklp = V_n6(x, a, 0, i, j, k, l, p, q);
	return Phi_ijklp;

}

inline double Phi_11(double* x, double a, int n, int i)
{
	double Phi_n_i = a * a * (d[i][n] * V_10(x, a, n) + x[n] * V_n1(x, a, 1, i));
	return Phi_n_i;
}

inline double Phi_12(double* x, double a, int n, int i, int j)
{
	double Phi_n_ij = a * a * (d[i][n] * V_n1(x, a, 1, j) + d[j][n] * V_n1(x, a, 1, i) + x[n] * V_n2(x, a, 1, i, j));
	return Phi_n_ij;
}

inline double Phi_13(double* x, double a, int n, int i, int j, int k)
{
	double Phi_n_ijk = a * a * (d[i][n] * V_n2(x, a, 1, j, k) + d[j][n] * V_n2(x, a, 1, i, k) + d[n][k] * V_n2(x, a, 1, i, j) + x[n] * V_n3(x, a, 1, i, j, k));
	return Phi_n_ijk;
}

inline double Phi_14(double* x, double a, int n, int i, int j, int k, int l)
{
	double Phi_n_ijkl = a * a * (d[i][n] * V_n3(x, a, 1, j, k, l) + d[j][n] * V_n3(x, a, 1, i, k, l) + d[n][k] * V_n3(x, a, 1, i, j, l) + d[n][l] * V_n3(x, a, 1, i, j, k) + x[n] * V_n4(x, a, 1, i, j, k, l));
	return Phi_n_ijkl;
}

inline double Phi_15(double* x, double a, int n, int i, int j, int k, int l, int p)
{
	double Phi_n_ijklp = a * a * (d[i][n] * V_n4(x, a, 1, j, k, l, p) + d[j][n] * V_n4(x, a, 1, i, k, l, p) + d[n][k] * V_n4(x, a, 1, i, j, l, p) + d[n][l] * V_n4(x, a, 1, i, j, k, p) + d[n][p] * V_n4(x, a, 1, i, j, k, l) + x[n] * V_n5(x, a, 1, i, j, k, l, p));
	return Phi_n_ijklp;

}

inline double Phi_16(double* x, double a, int n, int i, int j, int k, int l, int p, int q)
{
	double Phi_n_ijklpq = a * a * (d[i][n] * V_n5(x, a, 1, j, k, l, p, q) + d[j][n] * V_n5(x, a, 1, i, k, l, p, q) + d[n][k] * V_n5(x, a, 1, i, j, l, p, q)
		+ d[n][l] * V_n5(x, a, 1, i, j, k, p, q) + d[n][p] * V_n5(x, a, 1, i, j, k, l, q) + d[n][q] * V_n5(x, a, 1, i, j, k, l, p)
		+ x[n] * V_n6(x, a, 1, i, j, k, l, p, q));
	return Phi_n_ijklpq;
}

inline double Phi_21(double* x, double a, int m, int n, int i)
{
	double Phi_mn_i = a * a * (a * a * (d[m][i] * x[n] + d[n][i] * x[m]) * V_20(x, a, m, n) + a * a * x[m] * x[n] * V_n1(x, a, 2, i) - d[m][n] * x[i] * (V_10(x, a, i) - a * a * V_20(x, a, i, m)));
	return Phi_mn_i;
}

inline double Phi_22(double* x, double a, int m, int n, int i, int j)
{
	double Phi_mn_ij = a * a * (a * a * (d[m][i] * d[n][j] + d[n][i] * d[m][j]) * V_20(x, a, m, n) + a * a * (d[m][i] * x[n] + d[n][i] * x[m]) * V_n1(x, a, 2, j)
		+ a * a * (d[m][j] * x[n] + d[n][j] * x[m]) * V_n1(x, a, 2, i) + a * a * x[m] * x[n] * V_n2(x, a, 2, i, j)
		- d[m][n] * (d[i][j] * (V_10(x, a, i) - a * a * V_20(x, a, i, m)) + x[i] * (V_n1(x, a, 1, j) - a * a * V_n1(x, a, 2, j))));
	return Phi_mn_ij;
}

inline double Phi_23(double* x, double a, int m, int n, int i, int j, int k)
{
	double Phi_mn_ijk = a * a * a * a * ((d[m][i] * d[n][j] + d[n][i] * d[m][j]) * V_n1(x, a, 2, k) + (d[m][i] * d[n][k] + d[n][i] * d[m][k]) * V_n1(x, a, 2, j)
		+ (d[m][j] * d[n][k] + d[n][j] * d[m][k]) * V_n1(x, a, 2, i) + (d[m][i] * x[n] + d[n][i] * x[m]) * V_n2(x, a, 2, j, k)
		+ (d[m][j] * x[n] + d[n][j] * x[m]) * V_n2(x, a, 2, i, k) + (d[m][k] * x[n] + d[n][k] * x[m]) * V_n2(x, a, 2, i, j) + x[m] * x[n] * V_n3(x, a, 2, i, j, k))
		- a * a * d[m][n] * (d[i][j] * (V_n1(x, a, 1, k) - a * a * V_n1(x, a, 2, k)) + d[i][k] * (V_n1(x, a, 1, j) - a * a * V_n1(x, a, 2, j)) + x[i] * (V_n2(x, a, 1, j, k) - a * a * V_n2(x, a, 2, j, k)));
	return Phi_mn_ijk;
}

inline double Phi_24(double* x, double a, int m, int n, int i, int j, int k, int l)
{
	double Phi_mn_ijkl = a * a * a * a * ((d[m][i] * d[n][j] + d[n][i] * d[m][j]) * V_n2(x, a, 2, k, l) + (d[m][i] * d[n][k] + d[n][i] * d[m][k]) * V_n2(x, a, 2, j, l)
		+ (d[m][j] * d[n][k] + d[n][j] * d[m][k]) * V_n2(x, a, 2, i, l) + (d[m][i] * d[n][l] + d[n][i] * d[m][l]) * V_n2(x, a, 2, j, k)
		+ (d[m][j] * d[n][l] + d[n][j] * d[m][l]) * V_n2(x, a, 2, i, k) + (d[m][k] * d[n][l] + d[n][k] * d[m][l]) * V_n2(x, a, 2, i, j)
		+ (d[m][i] * x[n] + d[n][i] * x[m]) * V_n3(x, a, 2, j, k, l) + (d[m][j] * x[n] + d[n][j] * x[m]) * V_n3(x, a, 2, i, k, l)
		+ (d[m][k] * x[n] + d[n][k] * x[m]) * V_n3(x, a, 2, i, j, l) + (d[m][l] * x[n] + d[n][l] * x[m]) * V_n3(x, a, 2, i, j, k) + x[m] * x[n] * V_n4(x, a, 2, i, j, k, l))
		- a * a * d[m][n] * (d[i][j] * (V_n2(x, a, 1, k, l) - a * a * V_n2(x, a, 2, k, l)) + d[i][k] * (V_n2(x, a, 1, j, l) - a * a * V_n2(x, a, 2, j, l))
			+ d[i][l] * (V_n2(x, a, 1, j, k) - a * a * V_n2(x, a, 2, j, k)) + x[i] * (V_n3(x, a, 1, j, k, l) - a * a * V_n3(x, a, 2, j, k, l)));
	return Phi_mn_ijkl;
}

inline double Phi_25(double* x, double a, int m, int n, int i, int j, int k, int l, int p)
{
	double Phi_mn_ijklp = a * a * a * a * ((d[m][i] * d[n][j] + d[n][i] * d[m][j]) * V_n3(x, a, 2, k, l, p) + (d[m][i] * d[n][k] + d[n][i] * d[m][k]) * V_n3(x, a, 2, j, l, p)
		+ (d[m][j] * d[n][k] + d[n][j] * d[m][k]) * V_n3(x, a, 2, i, l, p) + (d[m][i] * d[n][l] + d[n][i] * d[m][l]) * V_n3(x, a, 2, j, k, p)
		+ (d[m][j] * d[n][l] + d[n][j] * d[m][l]) * V_n3(x, a, 2, i, k, p) + (d[m][k] * d[n][l] + d[n][k] * d[m][l]) * V_n3(x, a, 2, i, j, p)
		+ (d[m][i] * x[n] + d[n][i] * x[m]) * V_n4(x, a, 2, j, k, l, p) + (d[m][j] * x[n] + d[n][j] * x[m]) * V_n4(x, a, 2, i, k, l, p)
		+ (d[m][i] * d[n][p] + d[n][i] * d[m][p]) * V_n3(x, a, 2, j, k, l) + (d[m][j] * d[n][p] + d[n][j] * d[m][p]) * V_n3(x, a, 2, i, k, l)
		+ (d[m][k] * d[n][p] + d[n][k] * d[m][p]) * V_n3(x, a, 2, i, j, l) + (d[m][l] * d[n][p] + d[n][l] * d[m][p]) * V_n3(x, a, 2, i, j, k) + (d[m][p] * x[n] + d[n][p] * x[m]) * V_n4(x, a, 2, i, j, k, l)
		+ (d[m][k] * x[n] + d[n][k] * x[m]) * V_n4(x, a, 2, i, j, l, p) + (d[m][l] * x[n] + d[n][l] * x[m]) * V_n4(x, a, 2, i, j, k, p) + x[m] * x[n] * V_n5(x, a, 2, i, j, k, l, p))

		- a * a * d[m][n] * (d[i][j] * (V_n3(x, a, 1, k, l, p) - a * a * V_n3(x, a, 2, k, l, p)) + d[i][k] * (V_n3(x, a, 1, j, l, p) - a * a * V_n3(x, a, 2, j, l, p))
			+ d[i][l] * (V_n3(x, a, 1, j, k, p) - a * a * V_n3(x, a, 2, j, k, p)) + d[i][p] * (V_n3(x, a, 1, j, k, l) - a * a * V_n3(x, a, 2, j, k, l)) + x[i] * (V_n4(x, a, 1, j, k, l, p) - a * a * V_n4(x, a, 2, j, k, l, p)));
	return Phi_mn_ijklp;
}

inline double Phi_26(double* x, double a, int m, int n, int i, int j, int k, int l, int p, int q)
{
	double Phi_mn_ijklpq = a * a * a * a * ((d[m][i] * d[n][j] + d[n][i] * d[m][j]) * V_n4(x, a, 2, k, l, p, q) + (d[m][i] * d[n][k] + d[n][i] * d[m][k]) * V_n4(x, a, 2, j, l, p, q)
		+ (d[m][j] * d[n][k] + d[n][j] * d[m][k]) * V_n4(x, a, 2, i, l, p, q) + (d[m][i] * d[n][l] + d[n][i] * d[m][l]) * V_n4(x, a, 2, j, k, p, q)
		+ (d[m][j] * d[n][l] + d[n][j] * d[m][l]) * V_n4(x, a, 2, i, k, p, q) + (d[m][k] * d[n][l] + d[n][k] * d[m][l]) * V_n4(x, a, 2, i, j, p, q)
		+ (d[m][i] * d[n][q] + d[n][i] * d[m][q]) * V_n4(x, a, 2, j, k, l, p) + (d[m][i] * x[n] + d[n][i] * x[m]) * V_n5(x, a, 2, j, k, l, p, q)
		+ (d[m][j] * d[n][q] + d[n][j] * d[m][q]) * V_n4(x, a, 2, i, k, l, p) + (d[m][j] * x[n] + d[n][j] * x[m]) * V_n5(x, a, 2, i, k, l, p, q)

		+ (d[m][i] * d[n][p] + d[n][i] * d[m][p]) * V_n4(x, a, 2, j, k, l, q) + (d[m][j] * d[n][p] + d[n][j] * d[m][p]) * V_n4(x, a, 2, i, k, l, q)
		+ (d[m][k] * d[n][p] + d[n][k] * d[m][p]) * V_n4(x, a, 2, i, j, l, q) + (d[m][l] * d[n][p] + d[n][l] * d[m][p]) * V_n4(x, a, 2, i, j, k, q)
		+ (d[m][p] * d[n][q] + d[n][p] * d[m][q]) * V_n4(x, a, 2, i, j, k, l) + (d[m][p] * x[n] + d[n][p] * x[m]) * V_n5(x, a, 2, i, j, k, l, q)
		+ (d[m][k] * d[n][q] + d[n][k] * d[m][q]) * V_n4(x, a, 2, i, j, l, p) + (d[m][k] * x[n] + d[n][k] * x[m]) * V_n5(x, a, 2, i, j, l, p, q)
		+ (d[m][l] * d[n][q] + d[n][l] * d[m][q]) * V_n4(x, a, 2, i, j, k, p) + (d[m][l] * x[n] + d[n][l] * x[m]) * V_n5(x, a, 2, i, j, k, p, q)
		+ (d[m][q] * x[n] + d[n][q] * x[m]) * V_n5(x, a, 2, i, j, k, l, p) + x[m] * x[n] * V_n6(x, a, 2, i, j, k, l, p, q))

		- a * a * d[m][n] * (d[i][j] * (V_n4(x, a, 1, k, l, p, q) - a * a * V_n4(x, a, 2, k, l, p, q)) + d[i][k] * (V_n4(x, a, 1, j, l, p, q) - a * a * V_n4(x, a, 2, j, l, p, q))
			+ d[i][l] * (V_n4(x, a, 1, j, k, p, q) - a * a * V_n4(x, a, 2, j, k, p, q)) + d[i][p] * (V_n4(x, a, 1, j, k, l, q) - a * a * V_n4(x, a, 2, j, k, l, q))
			+ d[i][q] * (V_n4(x, a, 1, j, k, l, p) - a * a * V_n4(x, a, 2, j, k, l, p)) + x[i] * (V_n5(x, a, 1, j, k, l, p, q) - a * a * V_n5(x, a, 2, j, k, l, p, q)));
	return Phi_mn_ijklpq;

}

////////////////////////////////////////////////////////////

inline double Psi_03(double* x, double a, int i, int j, int k)
{
	double Psi_ijk = d[i][j] * (V_n1(x, a, 0, k) - a * a * V_n1(x, a, 1, k)) + d[i][k] * (V_n1(x, a, 0, j) - a * a * V_n1(x, a, 1, j)) + x[i] * (V_n2(x, a, 0, j, k) - a * a * V_n2(x, a, 1, j, k));
	return Psi_ijk;
}

inline double Psi_04(double* x, double a, int i, int j, int k, int l)
{
	double Psi_ijkl = d[i][j] * (V_n2(x, a, 0, k, l) - a * a * V_n2(x, a, 1, k, l)) + d[i][k] * (V_n2(x, a, 0, j, l) - a * a * V_n2(x, a, 1, j, l)) + d[i][l] * (V_n2(x, a, 0, j, k) - a * a * V_n2(x, a, 1, j, k))
		+ x[i] * (V_n3(x, a, 0, j, k, l) - a * a * V_n3(x, a, 1, j, k, l));
	return Psi_ijkl;
}

inline double Psi_05(double* x, double a, int i, int j, int k, int l, int s)
{
	double Psi_ijkls = d[i][j] * (V_n3(x, a, 0, k, l, s) - a * a * V_n3(x, a, 1, k, l, s)) + d[i][k] * (V_n3(x, a, 0, j, l, s) - a * a * V_n3(x, a, 1, j, l, s)) + d[i][l] * (V_n3(x, a, 0, j, k, s) - a * a * V_n3(x, a, 1, j, k, s))
		+ d[i][s] * (V_n3(x, a, 0, j, k, l) - a * a * V_n3(x, a, 1, j, k, l)) + x[i] * (V_n4(x, a, 0, j, k, l, s) - a * a * V_n4(x, a, 1, j, k, l, s));
	return Psi_ijkls;
}

inline double Psi_06(double* x, double a, int i, int j, int k, int l, int s, int t)
{
	double Psi_ijklst = d[i][j] * (V_n4(x, a, 0, k, l, s, t) - a * a * V_n4(x, a, 1, k, l, s, t)) + d[i][k] * (V_n4(x, a, 0, j, l, s, t) - a * a * V_n4(x, a, 1, j, l, s, t)) + d[i][l] * (V_n4(x, a, 0, j, k, s, t) - a * a * V_n4(x, a, 1, j, k, s, t))
		+ d[i][s] * (V_n4(x, a, 0, j, k, l, t) - a * a * V_n4(x, a, 1, j, k, l, t)) + d[i][t] * (V_n4(x, a, 0, j, k, l, s) - a * a * V_n4(x, a, 1, j, k, l, s)) + x[i] * (V_n5(x, a, 0, j, k, l, s, t) - a * a * V_n5(x, a, 1, j, k, l, s, t));
	return Psi_ijklst;
}

inline double Psi_07(double* x, double a, int i, int j, int k, int l, int s, int t, int p)
{
	double Psi_ijklstp = d[i][j] * (V_n5(x, a, 0, k, l, s, t, p) - a * a * V_n5(x, a, 1, k, l, s, t, p)) + d[i][k] * (V_n5(x, a, 0, j, l, s, t, p) - a * a * V_n5(x, a, 1, j, l, s, t, p)) + d[i][l] * (V_n5(x, a, 0, j, k, s, t, p) - a * a * V_n5(x, a, 1, j, k, s, t, p))
		+ d[i][s] * (V_n5(x, a, 0, j, k, l, t, p) - a * a * V_n5(x, a, 1, j, k, l, t, p)) + d[i][t] * (V_n5(x, a, 0, j, k, l, s, p) - a * a * V_n5(x, a, 1, j, k, l, s, p)) + d[i][p] * (V_n5(x, a, 0, j, k, l, s, t) - a * a * V_n5(x, a, 1, j, k, l, s, t)) + x[i] * (V_n6(x, a, 0, j, k, l, s, t, p) - a * a * V_n6(x, a, 1, j, k, l, s, t, p));
	return Psi_ijklstp;

}

inline double Psi_13(double* x, double a, int n, int i, int j, int k)
{
	double Psi_n_ijk = a * a * (d[i][n] * d[j][k] * (V_n0(x, a, 1) - a * a * V_n0(x, a, 2)) + d[i][n] * x[j] * (V_n1(x, a, 1, k) - a * a * V_n1(x, a, 2, k))
		+ (d[n][j] * d[i][k] + d[i][j] * d[n][k]) * (V_n0(x, a, 1) - a * a * V_n0(x, a, 2)) + (d[n][j] * x[i] + d[i][j] * x[n]) * (V_n1(x, a, 1, k) - a * a * V_n1(x, a, 2, k))
		+ (d[n][k] * x[i] + d[i][k] * x[n]) * (V_n1(x, a, 1, j) - a * a * V_n1(x, a, 2, j)) + x[n] * x[i] * (V_n2(x, a, 1, j, k) - a * a * V_n2(x, a, 2, j, k)));

	return Psi_n_ijk;
}

inline double Psi_14(double* x, double a, int n, int i, int j, int k, int l)
{
	double Psi_n_ijkl = a * a * (d[i][n] * d[j][k] * (V_n1(x, a, 1, l) - a * a * V_n1(x, a, 2, l)) + d[i][n] * d[j][l] * (V_n1(x, a, 1, k) - a * a * V_n1(x, a, 2, k))
		+ d[i][n] * x[j] * (V_n2(x, a, 1, k, l) - a * a * V_n2(x, a, 2, k, l)) + (d[n][j] * d[i][k] + d[i][j] * d[n][k]) * (V_n1(x, a, 1, l) - a * a * V_n1(x, a, 2, l))
		+ (d[n][j] * d[i][l] + d[i][j] * d[n][l]) * (V_n1(x, a, 1, k) - a * a * V_n1(x, a, 2, k)) + (d[n][j] * x[i] + d[i][j] * x[n]) * (V_n2(x, a, 1, k, l) - a * a * V_n2(x, a, 2, k, l))
		+ (d[n][k] * d[i][l] + d[i][k] * d[n][l]) * (V_n1(x, a, 1, j) - a * a * V_n1(x, a, 2, j))
		+ (d[n][k] * x[i] + d[i][k] * x[n]) * (V_n2(x, a, 1, j, l) - a * a * V_n2(x, a, 2, j, l))
		+ (d[n][l] * x[i] + d[i][l] * x[n]) * (V_n2(x, a, 1, j, k) - a * a * V_n2(x, a, 2, j, k)) + x[n] * x[i] * (V_n3(x, a, 1, j, k, l) - a * a * V_n3(x, a, 2, j, k, l)));
	return Psi_n_ijkl;
}

inline double Psi_15(double* x, double a, int n, int i, int j, int k, int l, int s)
{
	double Psi_n_ijkls = a * a * (d[i][n] * d[j][k] * (V_n2(x, a, 1, l, s) - a * a * V_n2(x, a, 2, l, s))
		+ d[i][n] * d[j][l] * (V_n2(x, a, 1, k, s) - a * a * V_n2(x, a, 2, k, s)) + d[i][n] * d[j][s] * (V_n2(x, a, 1, k, l) - a * a * V_n2(x, a, 2, k, l)) + d[i][n] * x[j] * (V_n3(x, a, 1, k, l, s) - a * a * V_n3(x, a, 2, k, l, s))
		+ (d[n][j] * d[i][k] + d[i][j] * d[n][k]) * (V_n2(x, a, 1, l, s) - a * a * V_n2(x, a, 2, l, s)) + (d[n][j] * d[i][l] + d[i][j] * d[n][l]) * (V_n2(x, a, 1, k, s) - a * a * V_n2(x, a, 2, k, s))
		+ (d[n][j] * d[i][s] + d[i][j] * d[n][s]) * (V_n2(x, a, 1, k, l) - a * a * V_n2(x, a, 2, k, l)) + (d[n][j] * x[i] + d[i][j] * x[n]) * (V_n3(x, a, 1, k, l, s) - a * a * V_n3(x, a, 2, k, l, s))
		+ (d[n][k] * d[i][l] + d[i][k] * d[n][l]) * (V_n2(x, a, 1, j, s) - a * a * V_n2(x, a, 2, j, s)) + (d[n][k] * d[i][s] + d[i][k] * d[n][s]) * (V_n2(x, a, 1, j, l) - a * a * V_n2(x, a, 2, j, l))
		+ (d[n][k] * x[i] + d[i][k] * x[n]) * (V_n3(x, a, 1, j, l, s) - a * a * V_n3(x, a, 2, j, l, s)) + (d[n][l] * d[i][s] + d[i][l] * d[n][s]) * (V_n2(x, a, 1, j, k) - a * a * V_n2(x, a, 2, j, k))
		+ (d[n][l] * x[i] + d[i][l] * x[n]) * (V_n3(x, a, 1, j, k, s) - a * a * V_n3(x, a, 2, j, k, s)) + (d[n][s] * x[i] + d[i][s] * x[n]) * (V_n3(x, a, 1, j, k, l) - a * a * V_n3(x, a, 2, j, k, l))
		+ x[n] * x[i] * (V_n4(x, a, 1, j, k, l, s) - a * a * V_n4(x, a, 2, j, k, l, s)));
	return Psi_n_ijkls;
}

inline double Psi_16(double* x, double a, int n, int i, int j, int k, int l, int s, int q)
{
	double Psi_n_ijklsq = a * a * (d[i][n] * d[j][k] * (V_n3(x, a, 1, l, s, q) - a * a * V_n3(x, a, 2, l, s, q)) + d[i][n] * d[j][l] * (V_n3(x, a, 1, k, s, q) - a * a * V_n3(x, a, 2, k, s, q))
		+ d[i][n] * d[j][s] * (V_n3(x, a, 1, k, l, q) - a * a * V_n3(x, a, 2, k, l, q)) + d[i][n] * d[j][q] * (V_n3(x, a, 1, k, l, s) - a * a * V_n3(x, a, 2, k, l, s)) + d[i][n] * x[j] * (V_n4(x, a, 1, k, l, s, q) - a * a * V_n4(x, a, 2, k, l, s, q))
		+ (d[n][j] * d[i][k] + d[i][j] * d[n][k]) * (V_n3(x, a, 1, l, s, q) - a * a * V_n3(x, a, 2, l, s, q)) + (d[n][j] * d[i][l] + d[i][j] * d[n][l]) * (V_n3(x, a, 1, k, s, q) - a * a * V_n3(x, a, 2, k, s, q))
		+ (d[n][j] * d[i][s] + d[i][j] * d[n][s]) * (V_n3(x, a, 1, k, l, q) - a * a * V_n3(x, a, 2, k, l, q)) + (d[n][j] * d[i][q] + d[i][j] * d[n][q]) * (V_n3(x, a, 1, k, l, s) - a * a * V_n3(x, a, 2, k, l, s))
		+ (d[n][j] * x[i] + d[i][j] * x[n]) * (V_n4(x, a, 1, k, l, s, q) - a * a * V_n4(x, a, 2, k, l, s, q)) + (d[n][k] * d[i][l] + d[i][k] * d[n][l]) * (V_n3(x, a, 1, j, s, q) - a * a * V_n3(x, a, 2, j, s, q))
		+ (d[n][k] * d[i][s] + d[i][k] * d[n][s]) * (V_n3(x, a, 1, j, l, q) - a * a * V_n3(x, a, 2, j, l, q)) + (d[n][k] * d[i][q] + d[i][k] * d[n][q]) * (V_n3(x, a, 1, j, l, s) - a * a * V_n3(x, a, 2, j, l, s))
		+ (d[n][k] * x[i] + d[i][k] * x[n]) * (V_n4(x, a, 1, j, l, s, q) - a * a * V_n4(x, a, 2, j, l, s, q)) + (d[n][l] * d[i][s] + d[i][l] * d[n][s]) * (V_n3(x, a, 1, j, k, q) - a * a * V_n3(x, a, 2, j, k, q))
		+ (d[n][l] * d[i][q] + d[i][l] * d[n][q]) * (V_n3(x, a, 1, j, k, s) - a * a * V_n3(x, a, 2, j, k, s)) + (d[n][l] * x[i] + d[i][l] * x[n]) * (V_n4(x, a, 1, j, k, s, q) - a * a * V_n4(x, a, 2, j, k, s, q))
		+ (d[n][s] * d[i][q] + d[i][s] * d[n][q]) * (V_n3(x, a, 1, j, k, l) - a * a * V_n3(x, a, 2, j, k, l)) + (d[n][s] * x[i] + d[i][s] * x[n]) * (V_n4(x, a, 1, j, k, l, q) - a * a * V_n4(x, a, 2, j, k, l, q))
		+ (d[n][q] * x[i] + d[i][q] * x[n]) * (V_n4(x, a, 1, j, k, l, s) - a * a * V_n4(x, a, 2, j, k, l, s)) + x[n] * x[i] * (V_n5(x, a, 1, j, k, l, s, q) - a * a * V_n5(x, a, 2, j, k, l, s, q)));
	return Psi_n_ijklsq;

}

inline double Psi_17(double* x, double a, int n, int i, int j, int k, int l, int s, int q, int p)
{
	double Psi_n_ijklsqp = a * a * (d[i][n] * d[j][k] * (V_n4(x, a, 1, l, s, q, p) - a * a * V_n4(x, a, 2, l, s, q, p)) + d[i][n] * d[j][l] * (V_n4(x, a, 1, k, s, q, p) - a * a * V_n4(x, a, 2, k, s, q, p))
		+ d[i][n] * d[j][s] * (V_n4(x, a, 1, k, l, q, p) - a * a * V_n4(x, a, 2, k, l, q, p)) + d[i][n] * d[j][q] * (V_n4(x, a, 1, k, l, s, p) - a * a * V_n4(x, a, 2, k, l, s, p)) + d[i][n] * d[j][p] * (V_n4(x, a, 1, k, l, s, q) - a * a * V_n4(x, a, 2, k, l, s, q)) + d[i][n] * x[j] * (V_n5(x, a, 1, k, l, s, q, p) - a * a * V_n5(x, a, 2, k, l, s, q, p))
		+ (d[n][j] * d[i][k] + d[i][j] * d[n][k]) * (V_n4(x, a, 1, l, s, q, p) - a * a * V_n4(x, a, 2, l, s, q, p)) + (d[n][j] * d[i][l] + d[i][j] * d[n][l]) * (V_n4(x, a, 1, k, s, q, p) - a * a * V_n4(x, a, 2, k, s, q, p))
		+ (d[n][j] * d[i][s] + d[i][j] * d[n][s]) * (V_n4(x, a, 1, k, l, q, p) - a * a * V_n4(x, a, 2, k, l, q, p)) + (d[n][j] * d[i][q] + d[i][j] * d[n][q]) * (V_n4(x, a, 1, k, l, s, p) - a * a * V_n4(x, a, 2, k, l, s, p))
		+ (d[n][j] * x[i] + d[i][j] * x[n]) * (V_n5(x, a, 1, k, l, s, q, p) - a * a * V_n5(x, a, 2, k, l, s, q, p)) + (d[n][j] * d[i][p] + d[i][j] * d[n][p]) * (V_n4(x, a, 1, k, l, s, q) - a * a * V_n4(x, a, 2, k, l, s, q)) + (d[n][k] * d[i][l] + d[i][k] * d[n][l]) * (V_n4(x, a, 1, j, s, q, p) - a * a * V_n4(x, a, 2, j, s, q, p))
		+ (d[n][k] * d[i][s] + d[i][k] * d[n][s]) * (V_n4(x, a, 1, j, l, q, p) - a * a * V_n4(x, a, 2, j, l, q, p)) + (d[n][k] * d[i][q] + d[i][k] * d[n][q]) * (V_n4(x, a, 1, j, l, s, p) - a * a * V_n4(x, a, 2, j, l, s, p))
		+ (d[n][k] * x[i] + d[i][k] * x[n]) * (V_n5(x, a, 1, j, l, s, q, p) - a * a * V_n5(x, a, 2, j, l, s, q, p)) + (d[n][k] * d[i][p] + d[i][k] * d[n][p]) * (V_n4(x, a, 1, j, l, s, q) - a * a * V_n4(x, a, 2, j, l, s, q)) + (d[n][l] * d[i][s] + d[i][l] * d[n][s]) * (V_n4(x, a, 1, j, k, q, p) - a * a * V_n4(x, a, 2, j, k, q, p))
		+ (d[n][l] * d[i][q] + d[i][l] * d[n][q]) * (V_n4(x, a, 1, j, k, s, p) - a * a * V_n4(x, a, 2, j, k, s, p)) + (d[n][l] * x[i] + d[i][l] * x[n]) * (V_n5(x, a, 1, j, k, s, q, p) - a * a * V_n5(x, a, 2, j, k, s, q, p)) + (d[n][l] * d[i][p] + d[i][l] * d[n][p]) * (V_n4(x, a, 1, j, k, s, q) - a * a * V_n4(x, a, 2, j, k, s, q))
		+ (d[n][s] * d[i][q] + d[i][s] * d[n][q]) * (V_n4(x, a, 1, j, k, l, p) - a * a * V_n4(x, a, 2, j, k, l, p)) + (d[n][s] * x[i] + d[i][s] * x[n]) * (V_n5(x, a, 1, j, k, l, q, p) - a * a * V_n5(x, a, 2, j, k, l, q, p)) + (d[n][s] * d[i][p] + d[i][s] * d[n][p]) * (V_n4(x, a, 1, j, k, l, q) - a * a * V_n4(x, a, 2, j, k, l, q))
		+ (d[n][q] * x[i] + d[i][q] * x[n]) * (V_n5(x, a, 1, j, k, l, s, p) - a * a * V_n5(x, a, 2, j, k, l, s, p)) + (d[n][q] * d[i][p] + d[i][q] * d[n][p]) * (V_n4(x, a, 1, j, k, l, s) - a * a * V_n4(x, a, 2, j, k, l, s))
		+ x[n] * x[i] * (V_n6(x, a, 1, j, k, l, s, q, p) - a * a * V_n6(x, a, 2, j, k, l, s, q, p)) + (d[n][p] * x[i] + d[i][p] * x[n]) * (V_n5(x, a, 1, j, k, l, s, q) - a * a * V_n5(x, a, 2, j, k, l, s, q)));
	return Psi_n_ijklsqp;
}

inline double Psi_23(double* x, double a, int m, int n, int i, int j, int k)
{
	double Psi_mn_ijk = a * a * a * a * ((d[m][i] * d[n][j] + d[n][i] * d[m][j]) * x[k] * (V_n0(x, a, 2) - a * a * V_n0(x, a, 3))
		+ ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][k]) * (V_n0(x, a, 2) - a * a * V_n0(x, a, 3))
		+ (d[m][i] * x[n] + d[n][i] * x[m]) * x[j] * (V_n1(x, a, 2, k) - a * a * V_n1(x, a, 3, k))
		+ (d[m][j] * (d[n][k] * x[i] + d[i][k] * x[n]) + d[n][j] * (d[m][k] * x[i] + d[i][k] * x[m]) + d[i][j] * (d[m][k] * x[n] + d[n][k] * x[m])) * (V_n0(x, a, 2) - a * a * V_n0(x, a, 3))
		+ (d[m][j] * x[n] * x[i] + d[n][j] * x[m] * x[i] + d[i][j] * x[m] * x[n]) * (V_n1(x, a, 2, k) - a * a * V_n1(x, a, 3, k))
		+ (d[m][k] * x[n] * x[i] + d[n][k] * x[m] * x[i] + d[i][k] * x[m] * x[n]) * (V_n1(x, a, 2, j) - a * a * V_n1(x, a, 3, j))
		+ x[m] * x[n] * x[i] * (V_n2(x, a, 2, j, k) - a * a * V_n2(x, a, 3, j, k)))

		+ a * a * d[m][n] * (0 - d[i][j] * x[k] * (V_n0(x, a, 1) - 2 * a * a * V_n0(x, a, 2) + a * a * a * a * V_n0(x, a, 3))
			- (d[i][k] * x[j] + d[j][k] * x[i]) * (V_n0(x, a, 1) - 2 * a * a * V_n0(x, a, 2) + a * a * a * a * V_n0(x, a, 3))
			- x[i] * x[j] * (V_n1(x, a, 1, k) - 2 * a * a * V_n1(x, a, 2, k) + a * a * a * a * V_n1(x, a, 3, k)));
	return Psi_mn_ijk;
}

inline double Psi_24(double* x, double a, int m, int n, int i, int j, int k, int l)
{
	double Psi_mn_ijkl = a * a * a * a * ((d[m][i] * d[n][j] + d[n][i] * d[m][j]) * d[k][l] * (V_20(x, a, m, k) - a * a * V_n0(x, a, 3)) + (d[m][i] * d[n][j] + d[n][i] * d[m][j]) * x[k] * (V_n1(x, a, 2, l) - a * a * V_n1(x, a, 3, l))
		+ ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * d[j][l] + (d[m][i] * d[n][l] + d[n][i] * d[m][l]) * d[j][k]) * (V_20(x, a, m, j) - a * a * V_n0(x, a, 3))
		+ ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][k]) * (V_n1(x, a, 2, l) - a * a * V_n1(x, a, 3, l))
		+ ((d[m][i] * d[n][l] + d[n][i] * d[m][l]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][l]) * (V_n1(x, a, 2, k) - a * a * V_n1(x, a, 3, k))
		+ (d[m][i] * x[n] + d[n][i] * x[m]) * x[j] * (V_n2(x, a, 2, k, l) - a * a * V_n2(x, a, 3, k, l))
		+ (d[m][j] * (d[n][k] * d[i][l] + d[i][k] * d[n][l]) + d[n][j] * (d[m][k] * d[i][l] + d[i][k] * d[m][l]) + d[i][j] * (d[m][k] * d[n][l] + d[n][k] * d[m][l])) * (V_n0(x, a, 2) - a * a * V_n0(x, a, 3))
		+ (d[m][j] * (d[n][k] * x[i] + d[i][k] * x[n]) + d[n][j] * (d[m][k] * x[i] + d[i][k] * x[m]) + d[i][j] * (d[m][k] * x[n] + d[n][k] * x[m])) * (V_n1(x, a, 2, l) - a * a * V_n1(x, a, 3, l))
		+ (d[m][j] * (d[n][l] * x[i] + d[i][l] * x[n]) + d[n][j] * (d[m][l] * x[i] + d[i][l] * x[m]) + d[i][j] * (d[m][l] * x[n] + d[n][l] * x[m])) * (V_n1(x, a, 2, k) - a * a * V_n1(x, a, 3, k))
		+ (d[m][j] * x[n] * x[i] + d[n][j] * x[m] * x[i] + d[i][j] * x[m] * x[n]) * (V_n2(x, a, 2, k, l) - a * a * V_n2(x, a, 3, k, l))
		+ (d[m][k] * (d[n][l] * x[i] + d[i][l] * x[n]) + d[n][k] * (d[m][l] * x[i] + d[i][l] * x[m]) + d[i][k] * (d[m][l] * x[n] + d[n][l] * x[m])) * (V_n1(x, a, 2, j) - a * a * V_n1(x, a, 3, j))
		+ (d[m][k] * x[n] * x[i] + d[n][k] * x[m] * x[i] + d[i][k] * x[m] * x[n]) * (V_n2(x, a, 2, j, l) - a * a * V_n2(x, a, 3, j, l))
		+ (d[m][l] * x[n] * x[i] + d[n][l] * x[m] * x[i] + d[i][l] * x[m] * x[n]) * (V_n2(x, a, 2, j, k) - a * a * V_n2(x, a, 3, j, k))
		+ x[m] * x[n] * x[i] * (V_n3(x, a, 2, j, k, l) - a * a * V_n3(x, a, 3, j, k, l)))
		+ a * a * d[m][n] * (-d[i][j] * d[k][l] * (V_10(x, a, k) - (a * a + a * a) * V_20(x, a, m, k) + a * a * a * a * V_n0(x, a, 3))
			- d[i][j] * x[k] * (V_n1(x, a, 1, l) - (a * a + a * a) * V_n1(x, a, 2, l) + a * a * a * a * V_n1(x, a, 3, l))
			- (d[i][k] * d[j][l] + d[j][k] * d[i][l]) * (V_10(x, a, j) - (a * a + a * a) * V_20(x, a, m, j) + a * a * a * a * V_n0(x, a, 3))
			- (d[i][k] * x[j] + d[j][k] * x[i]) * (V_n1(x, a, 1, l) - (a * a + a * a) * V_n1(x, a, 2, l) + a * a * a * a * V_n1(x, a, 3, l))
			- (d[i][l] * x[j] + d[j][l] * x[i]) * (V_n1(x, a, 1, k) - (a * a + a * a) * V_n1(x, a, 2, k) + a * a * a * a * V_n1(x, a, 3, k))
			- x[i] * x[j] * (V_n2(x, a, 1, k, l) - (a * a + a * a) * V_n2(x, a, 2, k, l) + a * a * a * a * V_n2(x, a, 3, k, l)));
	return Psi_mn_ijkl;
}

inline double Psi_25(double* x, double a, int m, int n, int i, int j, int k, int l, int s)
{
	double Psi_mn_ijkls = a * a * a * a * ((d[m][i] * d[n][j] + d[n][i] * d[m][j]) * (d[k][l] * (V_n1(x, a, 2, s) - a * a * V_n1(x, a, 3, s)) + d[k][s] * (V_n1(x, a, 2, l) - a * a * V_n1(x, a, 3, l)) + x[k] * (V_n2(x, a, 2, l, s) - a * a * V_n2(x, a, 3, l, s)))
		+ ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * d[j][l] + (d[m][i] * d[n][l] + d[n][i] * d[m][l]) * d[j][k]) * (V_n1(x, a, 2, s) - a * a * V_n1(x, a, 3, s))
		+ ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * d[j][s] + (d[m][i] * d[n][s] + d[n][i] * d[m][s]) * d[j][k]) * (V_n1(x, a, 2, l) - a * a * V_n1(x, a, 3, l))
		+ ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][k]) * (V_n2(x, a, 2, l, s) - a * a * V_n2(x, a, 3, l, s))
		+ ((d[m][i] * d[n][l] + d[n][i] * d[m][l]) * d[j][s] + (d[m][i] * d[n][s] + d[n][i] * d[m][s]) * d[j][l]) * (V_n1(x, a, 2, k) - a * a * V_n1(x, a, 3, k))
		+ ((d[m][i] * d[n][l] + d[n][i] * d[m][l]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][l]) * (V_n2(x, a, 2, k, s) - a * a * V_n2(x, a, 3, k, s))
		+ ((d[m][i] * d[n][s] + d[n][i] * d[m][s]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][s]) * (V_n2(x, a, 2, k, l) - a * a * V_n2(x, a, 3, k, l))
		+ (d[m][i] * x[n] + d[n][i] * x[m]) * x[j] * (V_n3(x, a, 2, k, l, s) - a * a * V_n3(x, a, 3, k, l, s))
		+ (d[m][s] * x[n] * x[i] + d[n][s] * x[m] * x[i] + d[i][s] * x[m] * x[n]) * (V_n3(x, a, 2, j, k, l) - a * a * V_n3(x, a, 3, j, k, l)) + x[m] * x[n] * x[i] * (V_n4(x, a, 2, j, k, l, s) - a * a * V_n4(x, a, 3, j, k, l, s))
		+ (d[m][l] * (d[n][s] * x[i] + d[i][s] * x[n]) + d[n][l] * (d[m][s] * x[i] + d[i][s] * x[m]) + d[i][l] * (d[n][s] * x[m] + d[m][s] * x[n])) * (V_n2(x, a, 2, j, k) - a * a * V_n2(x, a, 3, j, k))
		+ (d[m][l] * x[n] * x[i] + d[n][l] * x[m] * x[i] + d[i][l] * x[n] * x[m]) * (V_n3(x, a, 2, j, k, s) - a * a * V_n3(x, a, 3, j, k, s))
		+ (d[m][k] * (d[n][l] * d[i][s] + d[i][l] * d[n][s]) + d[n][k] * (d[m][l] * d[i][s] + d[i][l] * d[m][s]) + d[i][k] * (d[m][l] * d[n][s] + d[n][l] * d[m][s])) * (V_n1(x, a, 2, j) - a * a * V_n1(x, a, 3, j))
		+ (d[m][k] * (d[n][l] * x[i] + d[i][l] * x[n]) + d[n][k] * (d[m][l] * x[i] + d[i][l] * x[m]) + d[i][k] * (d[m][l] * x[n] + d[n][l] * x[m])) * (V_n2(x, a, 2, j, s) - a * a * V_n2(x, a, 3, j, s))
		+ (d[m][k] * (d[n][s] * x[i] + d[i][s] * x[n]) + d[n][k] * (d[m][s] * x[i] + d[i][s] * x[m]) + d[i][k] * (d[m][s] * x[n] + d[n][s] * x[m])) * (V_n2(x, a, 2, j, l) - a * a * V_n2(x, a, 3, j, l))
		+ (d[m][k] * x[n] * x[i] + d[n][k] * x[m] * x[i] + d[i][k] * x[m] * x[n]) * (V_n3(x, a, 2, j, l, s) - a * a * V_n3(x, a, 3, j, l, s))
		+ (d[m][j] * (d[n][k] * d[i][l] + d[i][k] * d[n][l]) + d[n][j] * (d[m][k] * d[i][l] + d[i][k] * d[m][l]) + d[i][j] * (d[m][k] * d[n][l] + d[n][k] * d[m][l])) * (V_n1(x, a, 2, s) - a * a * V_n1(x, a, 3, s))
		+ (d[m][j] * (d[n][k] * d[i][s] + d[i][k] * d[n][s]) + d[n][j] * (d[m][k] * d[i][s] + d[i][k] * d[m][s]) + d[i][j] * (d[m][k] * d[n][s] + d[n][k] * d[m][s])) * (V_n1(x, a, 2, l) - a * a * V_n1(x, a, 3, l))
		+ (d[m][j] * (d[n][l] * d[i][s] + d[i][l] * d[n][s]) + d[n][j] * (d[m][l] * d[i][s] + d[i][l] * d[m][s]) + d[i][j] * (d[m][l] * d[n][s] + d[n][l] * d[m][s])) * (V_n1(x, a, 2, k) - a * a * V_n1(x, a, 3, k))
		+ (d[m][j] * (d[n][k] * x[i] + d[i][k] * x[n]) + d[n][j] * (d[m][k] * x[i] + d[i][k] * x[m]) + d[i][j] * (d[m][k] * x[n] + d[n][k] * x[m])) * (V_n2(x, a, 2, l, s) - a * a * V_n2(x, a, 3, l, s))
		+ (d[m][j] * (d[n][l] * x[i] + d[i][l] * x[n]) + d[n][j] * (d[m][l] * x[i] + d[i][l] * x[m]) + d[i][j] * (d[m][l] * x[n] + d[n][l] * x[m])) * (V_n2(x, a, 2, k, s) - a * a * V_n2(x, a, 3, k, s))
		+ (d[m][j] * (d[n][s] * x[i] + d[i][s] * x[n]) + d[n][j] * (d[m][s] * x[i] + d[i][s] * x[m]) + d[i][j] * (d[m][s] * x[n] + d[n][s] * x[m])) * (V_n2(x, a, 2, k, l) - a * a * V_n2(x, a, 3, k, l))
		+ (d[m][j] * x[n] * x[i] + d[n][j] * x[m] * x[i] + d[i][j] * x[m] * x[n]) * (V_n3(x, a, 2, k, l, s) - a * a * V_n3(x, a, 3, k, l, s)))
		+ a * a * d[m][n] * (-d[i][j] * d[k][l] * (V_n1(x, a, 1, s) - (a * a + a * a) * V_n1(x, a, 2, s) + a * a * a * a * V_n1(x, a, 3, s))
			- d[i][j] * d[k][s] * (V_n1(x, a, 1, l) - (a * a + a * a) * V_n1(x, a, 2, l) + a * a * a * a * V_n1(x, a, 3, l))
			- d[i][j] * x[k] * (V_n2(x, a, 1, l, s) - (a * a + a * a) * V_n2(x, a, 2, l, s) + a * a * a * a * V_n2(x, a, 3, l, s))
			- (d[i][k] * d[j][l] + d[j][k] * d[i][l]) * (V_n1(x, a, 1, s) - (a * a + a * a) * V_n1(x, a, 2, s) + a * a * a * a * V_n1(x, a, 3, s))
			- (d[i][k] * d[j][s] + d[j][k] * d[i][s]) * (V_n1(x, a, 1, l) - (a * a + a * a) * V_n1(x, a, 2, l) + a * a * a * a * V_n1(x, a, 3, l))
			- (d[i][k] * x[j] + d[j][k] * x[i]) * (V_n2(x, a, 1, l, s) - (a * a + a * a) * V_n2(x, a, 2, l, s) + a * a * a * a * V_n2(x, a, 3, l, s))
			- (d[i][l] * d[j][s] + d[j][l] * d[i][s]) * (V_n1(x, a, 1, k) - (a * a + a * a) * V_n1(x, a, 2, k) + a * a * a * a * V_n1(x, a, 3, k))
			- (d[i][l] * x[j] + d[j][l] * x[i]) * (V_n2(x, a, 1, k, s) - (a * a + a * a) * V_n2(x, a, 2, k, s) + a * a * a * a * V_n2(x, a, 3, k, s))
			- (d[i][s] * x[j] + d[j][s] * x[i]) * (V_n2(x, a, 1, k, l) - (a * a + a * a) * V_n2(x, a, 2, k, l) + a * a * a * a * V_n2(x, a, 3, k, l))
			- x[i] * x[j] * (V_n3(x, a, 1, k, l, s) - (a * a + a * a) * V_n3(x, a, 2, k, l, s) + a * a * a * a * V_n3(x, a, 3, k, l, s)));
	return Psi_mn_ijkls;
}

inline double Psi_26(double* x, double a, int m, int n, int i, int j, int k, int l, int s, int t)
{
	double Psi_mn_ijklst = a * a * a * a * ((d[m][i] * d[n][j] + d[n][i] * d[m][j]) * (d[k][l] * (V_n2(x, a, 2, s, t) - a * a * V_n2(x, a, 3, s, t)) + d[k][s] * (V_n2(x, a, 2, l, t) - a * a * V_n2(x, a, 3, l, t)) + x[k] * (V_n3(x, a, 2, l, s, t) - a * a * V_n3(x, a, 3, l, s, t)) + d[k][t] * (V_n2(x, a, 2, l, s) - a * a * V_n2(x, a, 3, l, s)))
		+ ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * d[j][l] + (d[m][i] * d[n][l] + d[n][i] * d[m][l]) * d[j][k]) * (V_n2(x, a, 2, s, t) - a * a * V_n2(x, a, 3, s, t))
		+ ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * d[j][s] + (d[m][i] * d[n][s] + d[n][i] * d[m][s]) * d[j][k]) * (V_n2(x, a, 2, l, t) - a * a * V_n2(x, a, 3, l, t))
		+ ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * d[j][t] + (d[m][i] * d[n][t] + d[n][i] * d[m][t]) * d[j][k]) * (V_n2(x, a, 2, l, s) - a * a * V_n2(x, a, 3, l, s))
		+ ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][k]) * (V_n3(x, a, 2, l, s, t) - a * a * V_n3(x, a, 3, l, s, t))
		+ ((d[m][i] * d[n][l] + d[n][i] * d[m][l]) * d[j][s] + (d[m][i] * d[n][s] + d[n][i] * d[m][s]) * d[j][l]) * (V_n2(x, a, 2, k, t) - a * a * V_n2(x, a, 3, k, t))
		+ ((d[m][i] * d[n][l] + d[n][i] * d[m][l]) * d[j][t] + (d[m][i] * d[n][t] + d[n][i] * d[m][t]) * d[j][l]) * (V_n2(x, a, 2, k, s) - a * a * V_n2(x, a, 3, k, s))
		+ ((d[m][i] * d[n][l] + d[n][i] * d[m][l]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][l]) * (V_n3(x, a, 2, k, s, t) - a * a * V_n3(x, a, 3, k, s, t))
		+ ((d[m][i] * d[n][s] + d[n][i] * d[m][s]) * d[j][t] + (d[m][i] * d[n][t] + d[n][i] * d[m][t]) * d[j][s]) * (V_n2(x, a, 2, k, l) - a * a * V_n2(x, a, 3, k, l))
		+ ((d[m][i] * d[n][s] + d[n][i] * d[m][s]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][s]) * (V_n3(x, a, 2, k, l, t) - a * a * V_n3(x, a, 3, k, l, t))
		+ ((d[m][i] * d[n][t] + d[n][i] * d[m][t]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][t]) * (V_n3(x, a, 2, k, l, s) - a * a * V_n3(x, a, 3, k, l, s))
		+ (d[m][i] * x[n] + d[n][i] * x[m]) * x[j] * (V_n4(x, a, 2, k, l, s, t) - a * a * V_n4(x, a, 3, k, l, s, t))
		+ (d[m][s] * (d[n][t] * x[i] + d[i][t] * x[n]) + d[n][s] * (d[m][t] * x[i] + d[i][t] * x[m]) + d[i][s] * (d[m][t] * x[n] + d[n][t] * x[m])) * (V_n3(x, a, 2, j, k, l) - a * a * V_n3(x, a, 3, j, k, l))
		+ (d[m][s] * x[n] * x[i] + d[n][s] * x[m] * x[i] + d[i][s] * x[m] * x[n]) * (V_n4(x, a, 2, j, k, l, t) - a * a * V_n4(x, a, 3, j, k, l, t))
		+ x[m] * x[n] * x[i] * (V_n5(x, a, 2, j, k, l, s, t) - a * a * V_n5(x, a, 3, j, k, l, s, t))
		+ (d[m][t] * x[n] * x[i] + d[n][t] * x[m] * x[i] + d[i][t] * x[m] * x[n]) * (V_n4(x, a, 2, j, k, l, s) - a * a * V_n4(x, a, 3, j, k, l, s))
		+ (d[m][l] * (d[n][s] * d[i][t] + d[i][s] * d[n][t]) + d[n][l] * (d[m][s] * d[i][t] + d[i][s] * d[m][t]) + d[i][l] * (d[n][s] * d[m][t] + d[m][s] * d[n][t])) * (V_n2(x, a, 2, j, k) - a * a * V_n2(x, a, 3, j, k))
		+ (d[m][l] * (d[n][s] * x[i] + d[i][s] * x[n]) + d[n][l] * (d[m][s] * x[i] + d[i][s] * x[m]) + d[i][l] * (d[n][s] * x[m] + d[m][s] * x[n])) * (V_n3(x, a, 2, j, k, t) - a * a * V_n3(x, a, 3, j, k, t))
		+ (d[m][l] * (d[n][t] * x[i] + d[i][t] * x[n]) + d[n][l] * (d[m][t] * x[i] + d[i][t] * x[m]) + d[i][l] * (d[n][t] * x[m] + d[m][t] * x[n])) * (V_n3(x, a, 2, j, k, s) - a * a * V_n3(x, a, 3, j, k, s))
		+ (d[m][l] * x[n] * x[i] + d[n][l] * x[m] * x[i] + d[i][l] * x[n] * x[m]) * (V_n4(x, a, 2, j, k, s, t) - a * a * V_n4(x, a, 3, j, k, s, t))
		+ (d[m][k] * (d[n][l] * d[i][s] + d[i][l] * d[n][s]) + d[n][k] * (d[m][l] * d[i][s] + d[i][l] * d[m][s]) + d[i][k] * (d[m][l] * d[n][s] + d[n][l] * d[m][s])) * (V_n2(x, a, 2, j, t) - a * a * V_n2(x, a, 3, j, t))
		+ (d[m][k] * (d[n][l] * d[i][t] + d[i][l] * d[n][t]) + d[n][k] * (d[m][l] * d[i][t] + d[i][l] * d[m][t]) + d[i][k] * (d[m][l] * d[n][t] + d[n][l] * d[m][t])) * (V_n2(x, a, 2, j, s) - a * a * V_n2(x, a, 3, j, s))
		+ (d[m][k] * (d[n][l] * x[i] + d[i][l] * x[n]) + d[n][k] * (d[m][l] * x[i] + d[i][l] * x[m]) + d[i][k] * (d[m][l] * x[n] + d[n][l] * x[m])) * (V_n3(x, a, 2, j, s, t) - a * a * V_n3(x, a, 3, j, s, t))
		+ (d[m][k] * (d[n][s] * d[i][t] + d[i][s] * d[n][t]) + d[n][k] * (d[m][s] * d[i][t] + d[i][s] * d[m][t]) + d[i][k] * (d[m][s] * d[n][t] + d[n][s] * d[m][t])) * (V_n2(x, a, 2, j, l) - a * a * V_n2(x, a, 3, j, l))
		+ (d[m][k] * (d[n][s] * x[i] + d[i][s] * x[n]) + d[n][k] * (d[m][s] * x[i] + d[i][s] * x[m]) + d[i][k] * (d[m][s] * x[n] + d[n][s] * x[m])) * (V_n3(x, a, 2, j, l, t) - a * a * V_n3(x, a, 3, j, l, t))
		+ (d[m][k] * (d[n][t] * x[i] + d[i][t] * x[n]) + d[n][k] * (d[m][t] * x[i] + d[i][t] * x[m]) + d[i][k] * (d[m][t] * x[n] + d[n][t] * x[m])) * (V_n3(x, a, 2, j, l, s) - a * a * V_n3(x, a, 3, j, l, s))
		+ (d[m][k] * x[n] * x[i] + d[n][k] * x[m] * x[i] + d[i][k] * x[m] * x[n]) * (V_n4(x, a, 2, j, l, s, t) - a * a * V_n4(x, a, 3, j, l, s, t))
		+ (d[m][j] * (d[n][k] * d[i][l] + d[i][k] * d[n][l]) + d[n][j] * (d[m][k] * d[i][l] + d[i][k] * d[m][l]) + d[i][j] * (d[m][k] * d[n][l] + d[n][k] * d[m][l])) * (V_n2(x, a, 2, s, t) - a * a * V_n2(x, a, 3, s, t))
		+ (d[m][j] * (d[n][k] * d[i][s] + d[i][k] * d[n][s]) + d[n][j] * (d[m][k] * d[i][s] + d[i][k] * d[m][s]) + d[i][j] * (d[m][k] * d[n][s] + d[n][k] * d[m][s])) * (V_n2(x, a, 2, l, t) - a * a * V_n2(x, a, 3, l, t))
		+ (d[m][j] * (d[n][l] * d[i][s] + d[i][l] * d[n][s]) + d[n][j] * (d[m][l] * d[i][s] + d[i][l] * d[m][s]) + d[i][j] * (d[m][l] * d[n][s] + d[n][l] * d[m][s])) * (V_n2(x, a, 2, k, t) - a * a * V_n2(x, a, 3, k, t))
		+ (d[m][j] * (d[n][k] * d[i][t] + d[i][k] * d[n][t]) + d[n][j] * (d[m][k] * d[i][t] + d[i][k] * d[m][t]) + d[i][j] * (d[m][k] * d[n][t] + d[n][k] * d[m][t])) * (V_n2(x, a, 2, l, s) - a * a * V_n2(x, a, 3, l, s))
		+ (d[m][j] * (d[n][k] * x[i] + d[i][k] * x[n]) + d[n][j] * (d[m][k] * x[i] + d[i][k] * x[m]) + d[i][j] * (d[m][k] * x[n] + d[n][k] * x[m])) * (V_n3(x, a, 2, l, s, t) - a * a * V_n3(x, a, 3, l, s, t))
		+ (d[m][j] * (d[n][l] * d[i][t] + d[i][l] * d[n][t]) + d[n][j] * (d[m][l] * d[i][t] + d[i][l] * d[m][t]) + d[i][j] * (d[m][l] * d[n][t] + d[n][l] * d[m][t])) * (V_n2(x, a, 2, k, s) - a * a * V_n2(x, a, 3, k, s))
		+ (d[m][j] * (d[n][l] * x[i] + d[i][l] * x[n]) + d[n][j] * (d[m][l] * x[i] + d[i][l] * x[m]) + d[i][j] * (d[m][l] * x[n] + d[n][l] * x[m])) * (V_n3(x, a, 2, k, s, t) - a * a * V_n3(x, a, 3, k, s, t))
		+ (d[m][j] * (d[n][s] * d[i][t] + d[i][s] * d[n][t]) + d[n][j] * (d[m][s] * d[i][t] + d[i][s] * d[m][t]) + d[i][j] * (d[m][s] * d[n][t] + d[n][s] * d[m][t])) * (V_n2(x, a, 2, k, l) - a * a * V_n2(x, a, 3, k, l))
		+ (d[m][j] * (d[n][s] * x[i] + d[i][s] * x[n]) + d[n][j] * (d[m][s] * x[i] + d[i][s] * x[m]) + d[i][j] * (d[m][s] * x[n] + d[n][s] * x[m])) * (V_n3(x, a, 2, k, l, t) - a * a * V_n3(x, a, 3, k, l, t))
		+ (d[m][j] * (d[n][t] * x[i] + d[i][t] * x[n]) + d[n][j] * (d[m][t] * x[i] + d[i][t] * x[m]) + d[i][j] * (d[m][t] * x[n] + d[n][t] * x[m])) * (V_n3(x, a, 2, k, l, s) - a * a * V_n3(x, a, 3, k, l, s))
		+ (d[m][j] * x[n] * x[i] + d[n][j] * x[m] * x[i] + d[i][j] * x[m] * x[n]) * (V_n4(x, a, 2, k, l, s, t) - a * a * V_n4(x, a, 3, k, l, s, t)))
		+ a * a * d[m][n] * (-d[i][j] * d[k][l] * (V_n2(x, a, 1, s, t) - (a * a + a * a) * V_n2(x, a, 2, s, t) + a * a * a * a * V_n2(x, a, 3, s, t))
			- d[i][j] * d[k][s] * (V_n2(x, a, 1, l, t) - (a * a + a * a) * V_n2(x, a, 2, l, t) + a * a * a * a * V_n2(x, a, 3, l, t))
			- d[i][j] * d[k][t] * (V_n2(x, a, 1, l, s) - (a * a + a * a) * V_n2(x, a, 2, l, s) + a * a * a * a * V_n2(x, a, 3, l, s))
			- d[i][j] * x[k] * (V_n3(x, a, 1, l, s, t) - (a * a + a * a) * V_n3(x, a, 2, l, s, t) + a * a * a * a * V_n3(x, a, 3, l, s, t))
			- (d[i][k] * d[j][l] + d[j][k] * d[i][l]) * (V_n2(x, a, 1, s, t) - (a * a + a * a) * V_n2(x, a, 2, s, t) + a * a * a * a * V_n2(x, a, 3, s, t))
			- (d[i][k] * d[j][s] + d[j][k] * d[i][s]) * (V_n2(x, a, 1, l, t) - (a * a + a * a) * V_n2(x, a, 2, l, t) + a * a * a * a * V_n2(x, a, 3, l, t))
			- (d[i][k] * d[j][t] + d[j][k] * d[i][t]) * (V_n2(x, a, 1, l, s) - (a * a + a * a) * V_n2(x, a, 2, l, s) + a * a * a * a * V_n2(x, a, 3, l, s))
			- (d[i][k] * x[j] + d[j][k] * x[i]) * (V_n3(x, a, 1, l, s, t) - (a * a + a * a) * V_n3(x, a, 2, l, s, t) + a * a * a * a * V_n3(x, a, 3, l, s, t))
			- (d[i][l] * d[j][s] + d[j][l] * d[i][s]) * (V_n2(x, a, 1, k, t) - (a * a + a * a) * V_n2(x, a, 2, k, t) + a * a * a * a * V_n2(x, a, 3, k, t))
			- (d[i][l] * d[j][t] + d[j][l] * d[i][t]) * (V_n2(x, a, 1, k, s) - (a * a + a * a) * V_n2(x, a, 2, k, s) + a * a * a * a * V_n2(x, a, 3, k, s))
			- (d[i][l] * x[j] + d[j][l] * x[i]) * (V_n3(x, a, 1, k, s, t) - (a * a + a * a) * V_n3(x, a, 2, k, s, t) + a * a * a * a * V_n3(x, a, 3, k, s, t))
			- (d[i][s] * d[j][t] + d[j][s] * d[i][t]) * (V_n2(x, a, 1, k, l) - (a * a + a * a) * V_n2(x, a, 2, k, l) + a * a * a * a * V_n2(x, a, 3, k, l))
			- (d[i][s] * x[j] + d[j][s] * x[i]) * (V_n3(x, a, 1, k, l, t) - (a * a + a * a) * V_n3(x, a, 2, k, l, t) + a * a * a * a * V_n3(x, a, 3, k, l, t))
			- (d[i][t] * x[j] + d[j][t] * x[i]) * (V_n3(x, a, 1, k, l, s) - (a * a + a * a) * V_n3(x, a, 2, k, l, s) + a * a * a * a * V_n3(x, a, 3, k, l, s))
			- x[i] * x[j] * (V_n4(x, a, 1, k, l, s, t) - (a * a + a * a) * V_n4(x, a, 2, k, l, s, t) + a * a * a * a * V_n4(x, a, 3, k, l, s, t)));
	return Psi_mn_ijklst;

}

inline double Psi_27(double* x, double a, int m, int n, int i, int j, int k, int l, int s, int t, int h)
{
	double Psi_mn_ijklsth = a * a * a * a * ((d[m][i] * d[n][j] + d[n][i] * d[m][j]) * (d[k][l] * (V_n3(x, a, 2, s, t, h) - a * a * V_n3(x, a, 3, s, t, h)) + d[k][s] * (V_n3(x, a, 2, l, t, h) - a * a * V_n3(x, a, 3, l, t, h))
		+ d[k][h] * (V_n3(x, a, 2, l, s, t) - a * a * V_n3(x, a, 3, l, s, t)) + x[k] * (V_n4(x, a, 2, l, s, t, h) - a * a * V_n4(x, a, 3, l, s, t, h)) + d[k][t] * (V_n3(x, a, 2, l, s, h) - a * a * V_n3(x, a, 3, l, s, h)))
		+ ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * d[j][l] + (d[m][i] * d[n][l] + d[n][i] * d[m][l]) * d[j][k]) * (V_n3(x, a, 2, s, t, h) - a * a * V_n3(x, a, 3, s, t, h))
		+ ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * d[j][s] + (d[m][i] * d[n][s] + d[n][i] * d[m][s]) * d[j][k]) * (V_n3(x, a, 2, l, t, h) - a * a * V_n3(x, a, 3, l, t, h))
		+ ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * d[j][t] + (d[m][i] * d[n][t] + d[n][i] * d[m][t]) * d[j][k]) * (V_n3(x, a, 2, l, s, h) - a * a * V_n3(x, a, 3, l, s, h))
		+ ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * d[j][h] + (d[m][i] * d[n][h] + d[n][i] * d[m][h]) * d[j][k]) * (V_n3(x, a, 2, l, s, t) - a * a * V_n3(x, a, 3, l, s, t))
		+ ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][k]) * (V_n4(x, a, 2, l, s, t, h) - a * a * V_n4(x, a, 3, l, s, t, h))
		+ ((d[m][i] * d[n][l] + d[n][i] * d[m][l]) * d[j][s] + (d[m][i] * d[n][s] + d[n][i] * d[m][s]) * d[j][l]) * (V_n3(x, a, 2, k, t, h) - a * a * V_n3(x, a, 3, k, t, h))
		+ ((d[m][i] * d[n][l] + d[n][i] * d[m][l]) * d[j][t] + (d[m][i] * d[n][t] + d[n][i] * d[m][t]) * d[j][l]) * (V_n3(x, a, 2, k, s, h) - a * a * V_n3(x, a, 3, k, s, h))
		+ ((d[m][i] * d[n][l] + d[n][i] * d[m][l]) * d[j][h] + (d[m][i] * d[n][h] + d[n][i] * d[m][h]) * d[j][l]) * (V_n3(x, a, 2, k, s, t) - a * a * V_n3(x, a, 3, k, s, t))
		+ ((d[m][i] * d[n][l] + d[n][i] * d[m][l]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][l]) * (V_n4(x, a, 2, k, s, t, h) - a * a * V_n4(x, a, 3, k, s, t, h))
		+ ((d[m][i] * d[n][s] + d[n][i] * d[m][s]) * d[j][t] + (d[m][i] * d[n][t] + d[n][i] * d[m][t]) * d[j][s]) * (V_n3(x, a, 2, k, l, h) - a * a * V_n3(x, a, 3, k, l, h))
		+ ((d[m][i] * d[n][s] + d[n][i] * d[m][s]) * d[j][h] + (d[m][i] * d[n][h] + d[n][i] * d[m][h]) * d[j][s]) * (V_n3(x, a, 2, k, l, t) - a * a * V_n3(x, a, 3, k, l, t))
		+ ((d[m][i] * d[n][s] + d[n][i] * d[m][s]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][s]) * (V_n4(x, a, 2, k, l, t, h) - a * a * V_n4(x, a, 3, k, l, t, h))
		+ ((d[m][i] * d[n][t] + d[n][i] * d[m][t]) * d[j][h] + (d[m][i] * d[n][h] + d[n][i] * d[m][h]) * d[j][t]) * (V_n3(x, a, 2, k, l, s) - a * a * V_n3(x, a, 3, k, l, s))
		+ ((d[m][i] * d[n][t] + d[n][i] * d[m][t]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][t]) * (V_n4(x, a, 2, k, l, s, h) - a * a * V_n4(x, a, 3, k, l, s, h))
		+ ((d[m][i] * x[n] + d[n][i] * x[m]) * d[j][h] + (d[m][i] * d[n][h] + d[m][h] * d[n][i]) * x[j]) * (V_n4(x, a, 2, k, l, s, t) - a * a * V_n4(x, a, 3, k, l, s, t))
		+ (d[m][i] * x[n] + d[n][i] * x[m]) * x[j] * (V_n5(x, a, 2, k, l, s, t, h) - a * a * V_n5(x, a, 3, k, l, s, t, h))
		+ (d[m][s] * (d[n][t] * d[i][h] + d[i][t] * d[n][h]) + d[n][s] * (d[m][t] * d[i][h] + d[i][t] * d[m][h]) + d[i][s] * (d[m][t] * d[n][h] + d[n][t] * d[m][h])) * (V_n3(x, a, 2, j, k, l) - a * a * V_n3(x, a, 3, j, k, l))
		+ (d[m][s] * (d[n][t] * x[i] + d[i][t] * x[n]) + d[n][s] * (d[m][t] * x[i] + d[i][t] * x[m]) + d[i][s] * (d[m][t] * x[n] + d[n][t] * x[m])) * (V_n4(x, a, 2, j, k, l, h) - a * a * V_n4(x, a, 3, j, k, l, h))
		+ (d[m][s] * (d[n][h] * x[i] + d[i][h] * x[n]) + d[n][s] * (d[m][h] * x[i] + d[i][h] * x[m]) + d[i][s] * (d[m][h] * x[n] + d[n][h] * x[m])) * (V_n4(x, a, 2, j, k, l, t) - a * a * V_n4(x, a, 3, j, k, l, t))
		+ (d[m][s] * x[n] * x[i] + d[n][s] * x[m] * x[i] + d[i][s] * x[m] * x[n]) * (V_n5(x, a, 2, j, k, l, t, h) - a * a * V_n5(x, a, 3, j, k, l, t, h))
		+ (d[m][h] * x[n] * x[i] + d[n][h] * x[m] * x[i] + d[i][h] * x[m] * x[n]) * (V_n5(x, a, 2, j, k, l, s, t) - a * a * V_n5(x, a, 3, j, k, l, s, t))
		+ x[m] * x[n] * x[i] * (V_n6(x, a, 2, j, k, l, s, t, h) - a * a * V_n6(x, a, 3, j, k, l, s, t, h))
		+ (d[m][t] * (d[n][h] * x[i] + d[i][h] * x[n]) + d[n][t] * (d[m][h] * x[i] + d[i][h] * x[m]) + d[i][t] * (d[m][h] * x[n] + d[n][h] * x[m])) * (V_n4(x, a, 2, j, k, l, s) - a * a * V_n4(x, a, 3, j, k, l, s))
		+ (d[m][t] * x[n] * x[i] + d[n][t] * x[m] * x[i] + d[i][t] * x[m] * x[n]) * (V_n5(x, a, 2, j, k, l, s, h) - a * a * V_n5(x, a, 3, j, k, l, s, h))
		+ (d[m][l] * (d[n][s] * d[i][t] + d[i][s] * d[n][t]) + d[n][l] * (d[m][s] * d[i][t] + d[i][s] * d[m][t]) + d[i][l] * (d[n][s] * d[m][t] + d[m][s] * d[n][t])) * (V_n3(x, a, 2, j, k, h) - a * a * V_n3(x, a, 3, j, k, h))
		+ (d[m][l] * (d[n][s] * d[i][h] + d[i][s] * d[n][h]) + d[n][l] * (d[m][s] * d[i][h] + d[i][s] * d[m][h]) + d[i][l] * (d[n][s] * d[m][h] + d[m][s] * d[n][h])) * (V_n3(x, a, 2, j, k, t) - a * a * V_n3(x, a, 3, j, k, t))
		+ (d[m][l] * (d[n][s] * x[i] + d[i][s] * x[n]) + d[n][l] * (d[m][s] * x[i] + d[i][s] * x[m]) + d[i][l] * (d[n][s] * x[m] + d[m][s] * x[n])) * (V_n4(x, a, 2, j, k, t, h) - a * a * V_n4(x, a, 3, j, k, t, h))
		+ (d[m][l] * (d[n][t] * d[i][h] + d[i][t] * d[n][h]) + d[n][l] * (d[m][t] * d[i][h] + d[i][t] * d[m][h]) + d[i][l] * (d[n][t] * d[m][h] + d[m][t] * d[n][h])) * (V_n3(x, a, 2, j, k, s) - a * a * V_n3(x, a, 3, j, k, s))
		+ (d[m][l] * (d[n][t] * x[i] + d[i][t] * x[n]) + d[n][l] * (d[m][t] * x[i] + d[i][t] * x[m]) + d[i][l] * (d[n][t] * x[m] + d[m][t] * x[n])) * (V_n4(x, a, 2, j, k, s, h) - a * a * V_n4(x, a, 3, j, k, s, h))
		+ (d[m][l] * (d[n][h] * x[i] + d[i][h] * x[n]) + d[n][l] * (d[m][h] * x[i] + d[i][h] * x[m]) + d[i][l] * (d[n][h] * x[m] + d[m][h] * x[n])) * (V_n4(x, a, 2, j, k, s, t) - a * a * V_n4(x, a, 3, j, k, s, t))
		+ (d[m][l] * x[n] * x[i] + d[n][l] * x[m] * x[i] + d[i][l] * x[n] * x[m]) * (V_n5(x, a, 2, j, k, s, t, h) - a * a * V_n5(x, a, 3, j, k, s, t, h))
		+ (d[m][k] * (d[n][l] * d[i][s] + d[i][l] * d[n][s]) + d[n][k] * (d[m][l] * d[i][s] + d[i][l] * d[m][s]) + d[i][k] * (d[m][l] * d[n][s] + d[n][l] * d[m][s])) * (V_n3(x, a, 2, j, t, h) - a * a * V_n3(x, a, 3, j, t, h))
		+ (d[m][k] * (d[n][l] * d[i][t] + d[i][l] * d[n][t]) + d[n][k] * (d[m][l] * d[i][t] + d[i][l] * d[m][t]) + d[i][k] * (d[m][l] * d[n][t] + d[n][l] * d[m][t])) * (V_n3(x, a, 2, j, s, h) - a * a * V_n3(x, a, 3, j, s, h))
		+ (d[m][k] * (d[n][l] * d[i][h] + d[i][l] * d[n][h]) + d[n][k] * (d[m][l] * d[i][h] + d[i][l] * d[m][h]) + d[i][k] * (d[m][l] * d[n][h] + d[n][l] * d[m][h])) * (V_n3(x, a, 2, j, s, t) - a * a * V_n3(x, a, 3, j, s, t))
		+ (d[m][k] * (d[n][l] * x[i] + d[i][l] * x[n]) + d[n][k] * (d[m][l] * x[i] + d[i][l] * x[m]) + d[i][k] * (d[m][l] * x[n] + d[n][l] * x[m])) * (V_n4(x, a, 2, j, s, t, h) - a * a * V_n4(x, a, 3, j, s, t, h))
		+ (d[m][k] * (d[n][s] * d[i][t] + d[i][s] * d[n][t]) + d[n][k] * (d[m][s] * d[i][t] + d[i][s] * d[m][t]) + d[i][k] * (d[m][s] * d[n][t] + d[n][s] * d[m][t])) * (V_n3(x, a, 2, j, l, h) - a * a * V_n3(x, a, 3, j, l, h))
		+ (d[m][k] * (d[n][s] * d[i][h] + d[i][s] * d[n][h]) + d[n][k] * (d[m][s] * d[i][h] + d[i][s] * d[m][h]) + d[i][k] * (d[m][s] * d[n][h] + d[n][s] * d[m][h])) * (V_n3(x, a, 2, j, l, t) - a * a * V_n3(x, a, 3, j, l, t))
		+ (d[m][k] * (d[n][s] * x[i] + d[i][s] * x[n]) + d[n][k] * (d[m][s] * x[i] + d[i][s] * x[m]) + d[i][k] * (d[m][s] * x[n] + d[n][s] * x[m])) * (V_n4(x, a, 2, j, l, t, h) - a * a * V_n4(x, a, 3, j, l, t, h))
		+ (d[m][k] * (d[n][t] * d[i][h] + d[i][t] * d[n][h]) + d[n][k] * (d[m][t] * d[i][h] + d[i][t] * d[m][h]) + d[i][k] * (d[m][t] * d[n][h] + d[n][t] * d[m][h])) * (V_n3(x, a, 2, j, l, s) - a * a * V_n3(x, a, 3, j, l, s))
		+ (d[m][k] * (d[n][t] * x[i] + d[i][t] * x[n]) + d[n][k] * (d[m][t] * x[i] + d[i][t] * x[m]) + d[i][k] * (d[m][t] * x[n] + d[n][t] * x[m])) * (V_n4(x, a, 2, j, l, s, h) - a * a * V_n4(x, a, 3, j, l, s, h))
		+ (d[m][k] * (d[n][h] * x[i] + d[i][h] * x[n]) + d[n][k] * (d[m][h] * x[i] + d[i][h] * x[m]) + d[i][k] * (d[m][h] * x[n] + d[n][h] * x[m])) * (V_n4(x, a, 2, j, l, s, t) - a * a * V_n4(x, a, 3, j, l, s, t))
		+ (d[m][k] * x[n] * x[i] + d[n][k] * x[m] * x[i] + d[i][k] * x[m] * x[n]) * (V_n5(x, a, 2, j, l, s, t, h) - a * a * V_n5(x, a, 3, j, l, s, t, h))
		+ (d[m][j] * (d[n][k] * d[i][l] + d[i][k] * d[n][l]) + d[n][j] * (d[m][k] * d[i][l] + d[i][k] * d[m][l]) + d[i][j] * (d[m][k] * d[n][l] + d[n][k] * d[m][l])) * (V_n3(x, a, 2, s, t, h) - a * a * V_n3(x, a, 3, s, t, h))
		+ (d[m][j] * (d[n][k] * d[i][s] + d[i][k] * d[n][s]) + d[n][j] * (d[m][k] * d[i][s] + d[i][k] * d[m][s]) + d[i][j] * (d[m][k] * d[n][s] + d[n][k] * d[m][s])) * (V_n3(x, a, 2, l, t, h) - a * a * V_n3(x, a, 3, l, t, h))
		+ (d[m][j] * (d[n][l] * d[i][s] + d[i][l] * d[n][s]) + d[n][j] * (d[m][l] * d[i][s] + d[i][l] * d[m][s]) + d[i][j] * (d[m][l] * d[n][s] + d[n][l] * d[m][s])) * (V_n3(x, a, 2, k, t, h) - a * a * V_n3(x, a, 3, k, t, h))
		+ (d[m][j] * (d[n][k] * d[i][t] + d[i][k] * d[n][t]) + d[n][j] * (d[m][k] * d[i][t] + d[i][k] * d[m][t]) + d[i][j] * (d[m][k] * d[n][t] + d[n][k] * d[m][t])) * (V_n3(x, a, 2, l, s, h) - a * a * V_n3(x, a, 3, l, s, h))
		+ (d[m][j] * (d[n][k] * d[i][h] + d[i][k] * d[n][h]) + d[n][j] * (d[m][k] * d[i][h] + d[i][k] * d[m][h]) + d[i][j] * (d[m][k] * d[n][h] + d[n][k] * d[m][h])) * (V_n3(x, a, 2, l, s, t) - a * a * V_n3(x, a, 3, l, s, t))
		+ (d[m][j] * (d[n][k] * x[i] + d[i][k] * x[n]) + d[n][j] * (d[m][k] * x[i] + d[i][k] * x[m]) + d[i][j] * (d[m][k] * x[n] + d[n][k] * x[m])) * (V_n4(x, a, 2, l, s, t, h) - a * a * V_n4(x, a, 3, l, s, t, h))
		+ (d[m][j] * (d[n][l] * d[i][t] + d[i][l] * d[n][t]) + d[n][j] * (d[m][l] * d[i][t] + d[i][l] * d[m][t]) + d[i][j] * (d[m][l] * d[n][t] + d[n][l] * d[m][t])) * (V_n3(x, a, 2, k, s, h) - a * a * V_n3(x, a, 3, k, s, h))
		+ (d[m][j] * (d[n][l] * d[i][h] + d[i][l] * d[n][h]) + d[n][j] * (d[m][l] * d[i][h] + d[i][l] * d[m][h]) + d[i][j] * (d[m][l] * d[n][h] + d[n][l] * d[m][h])) * (V_n3(x, a, 2, k, s, t) - a * a * V_n3(x, a, 3, k, s, t))
		+ (d[m][j] * (d[n][l] * x[i] + d[i][l] * x[n]) + d[n][j] * (d[m][l] * x[i] + d[i][l] * x[m]) + d[i][j] * (d[m][l] * x[n] + d[n][l] * x[m])) * (V_n4(x, a, 2, k, s, t, h) - a * a * V_n4(x, a, 3, k, s, t, h))
		+ (d[m][j] * (d[n][s] * d[i][t] + d[i][s] * d[n][t]) + d[n][j] * (d[m][s] * d[i][t] + d[i][s] * d[m][t]) + d[i][j] * (d[m][s] * d[n][t] + d[n][s] * d[m][t])) * (V_n3(x, a, 2, k, l, h) - a * a * V_n3(x, a, 3, k, l, h))
		+ (d[m][j] * (d[n][s] * d[i][h] + d[i][s] * d[n][h]) + d[n][j] * (d[m][s] * d[i][h] + d[i][s] * d[m][h]) + d[i][j] * (d[m][s] * d[n][h] + d[n][s] * d[m][h])) * (V_n3(x, a, 2, k, l, t) - a * a * V_n3(x, a, 3, k, l, t))
		+ (d[m][j] * (d[n][s] * x[i] + d[i][s] * x[n]) + d[n][j] * (d[m][s] * x[i] + d[i][s] * x[m]) + d[i][j] * (d[m][s] * x[n] + d[n][s] * x[m])) * (V_n4(x, a, 2, k, l, t, h) - a * a * V_n4(x, a, 3, k, l, t, h))
		+ (d[m][j] * (d[n][t] * d[i][h] + d[i][t] * d[n][h]) + d[n][j] * (d[m][t] * d[i][h] + d[i][t] * d[m][h]) + d[i][j] * (d[m][t] * d[n][h] + d[n][t] * d[m][h])) * (V_n3(x, a, 2, k, l, s) - a * a * V_n3(x, a, 3, k, l, s))
		+ (d[m][j] * (d[n][t] * x[i] + d[i][t] * x[n]) + d[n][j] * (d[m][t] * x[i] + d[i][t] * x[m]) + d[i][j] * (d[m][t] * x[n] + d[n][t] * x[m])) * (V_n4(x, a, 2, k, l, s, h) - a * a * V_n4(x, a, 3, k, l, s, h))
		+ (d[m][j] * (d[n][h] * x[i] + d[i][h] * x[n]) + d[n][j] * (d[m][h] * x[i] + d[i][h] * x[m]) + d[i][j] * (d[m][h] * x[n] + d[n][h] * x[m])) * (V_n4(x, a, 2, k, l, s, t) - a * a * V_n4(x, a, 3, k, l, s, t))
		+ (d[m][j] * x[n] * x[i] + d[n][j] * x[m] * x[i] + d[i][j] * x[m] * x[n]) * (V_n5(x, a, 2, k, l, s, t, h) - a * a * V_n5(x, a, 3, k, l, s, t, h)))
		+ a * a * d[m][n] * (-d[i][j] * d[k][l] * (V_n3(x, a, 1, s, t, h) - (a * a + a * a) * V_n3(x, a, 2, s, t, h) + a * a * a * a * V_n3(x, a, 3, s, t, h))
			- d[i][j] * d[k][s] * (V_n3(x, a, 1, l, t, h) - (a * a + a * a) * V_n3(x, a, 2, l, t, h) + a * a * a * a * V_n3(x, a, 3, l, t, h))
			- d[i][j] * d[k][t] * (V_n3(x, a, 1, l, s, h) - (a * a + a * a) * V_n3(x, a, 2, l, s, h) + a * a * a * a * V_n3(x, a, 3, l, s, h))
			- d[i][j] * d[k][h] * (V_n3(x, a, 1, l, s, t) - (a * a + a * a) * V_n3(x, a, 2, l, s, t) + a * a * a * a * V_n3(x, a, 3, l, s, t))
			- d[i][j] * x[k] * (V_n4(x, a, 1, l, s, t, h) - (a * a + a * a) * V_n4(x, a, 2, l, s, t, h) + a * a * a * a * V_n4(x, a, 3, l, s, t, h))
			- (d[i][k] * d[j][l] + d[j][k] * d[i][l]) * (V_n3(x, a, 1, s, t, h) - (a * a + a * a) * V_n3(x, a, 2, s, t, h) + a * a * a * a * V_n3(x, a, 3, s, t, h))
			- (d[i][k] * d[j][s] + d[j][k] * d[i][s]) * (V_n3(x, a, 1, l, t, h) - (a * a + a * a) * V_n3(x, a, 2, l, t, h) + a * a * a * a * V_n3(x, a, 3, l, t, h))
			- (d[i][k] * d[j][t] + d[j][k] * d[i][t]) * (V_n3(x, a, 1, l, s, h) - (a * a + a * a) * V_n3(x, a, 2, l, s, h) + a * a * a * a * V_n3(x, a, 3, l, s, h))
			- (d[i][k] * d[j][h] + d[j][k] * d[i][h]) * (V_n3(x, a, 1, l, s, t) - (a * a + a * a) * V_n3(x, a, 2, l, s, t) + a * a * a * a * V_n3(x, a, 3, l, s, t))
			- (d[i][k] * x[j] + d[j][k] * x[i]) * (V_n4(x, a, 1, l, s, t, h) - (a * a + a * a) * V_n4(x, a, 2, l, s, t, h) + a * a * a * a * V_n4(x, a, 3, l, s, t, h))
			- (d[i][l] * d[j][s] + d[j][l] * d[i][s]) * (V_n3(x, a, 1, k, t, h) - (a * a + a * a) * V_n3(x, a, 2, k, t, h) + a * a * a * a * V_n3(x, a, 3, k, t, h))
			- (d[i][l] * d[j][t] + d[j][l] * d[i][t]) * (V_n3(x, a, 1, k, s, h) - (a * a + a * a) * V_n3(x, a, 2, k, s, h) + a * a * a * a * V_n3(x, a, 3, k, s, h))
			- (d[i][l] * d[j][h] + d[j][l] * d[i][h]) * (V_n3(x, a, 1, k, s, t) - (a * a + a * a) * V_n3(x, a, 2, k, s, t) + a * a * a * a * V_n3(x, a, 3, k, s, t))
			- (d[i][l] * x[j] + d[j][l] * x[i]) * (V_n4(x, a, 1, k, s, t, h) - (a * a + a * a) * V_n4(x, a, 2, k, s, t, h) + a * a * a * a * V_n4(x, a, 3, k, s, t, h))
			- (d[i][s] * d[j][t] + d[j][s] * d[i][t]) * (V_n3(x, a, 1, k, l, h) - (a * a + a * a) * V_n3(x, a, 2, k, l, h) + a * a * a * a * V_n3(x, a, 3, k, l, h))
			- (d[i][s] * d[j][h] + d[j][s] * d[i][h]) * (V_n3(x, a, 1, k, l, t) - (a * a + a * a) * V_n3(x, a, 2, k, l, t) + a * a * a * a * V_n3(x, a, 3, k, l, t))
			- (d[i][s] * x[j] + d[j][s] * x[i]) * (V_n4(x, a, 1, k, l, t, h) - (a * a + a * a) * V_n4(x, a, 2, k, l, t, h) + a * a * a * a * V_n4(x, a, 3, k, l, t, h))
			- (d[i][t] * d[j][h] + d[j][t] * d[i][h]) * (V_n3(x, a, 1, k, l, s) - (a * a + a * a) * V_n3(x, a, 2, k, l, s) + a * a * a * a * V_n3(x, a, 3, k, l, s))
			- (d[i][t] * x[j] + d[j][t] * x[i]) * (V_n4(x, a, 1, k, l, s, h) - (a * a + a * a) * V_n4(x, a, 2, k, l, s, h) + a * a * a * a * V_n4(x, a, 3, k, l, s, h))
			- (d[i][h] * x[j] + d[j][h] * x[i]) * (V_n4(x, a, 1, k, l, s, t) - (a * a + a * a) * V_n4(x, a, 2, k, l, s, t) + a * a * a * a * V_n4(x, a, 3, k, l, s, t))
			- x[i] * x[j] * (V_n5(x, a, 1, k, l, s, t, h) - (a * a + a * a) * V_n5(x, a, 2, k, l, s, t, h) + a * a * a * a * V_n5(x, a, 3, k, l, s, t, h)));
	return Psi_mn_ijklsth;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double R(double* x)
{
	double R = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];

	R = sqrt(R);
	return R;
}

double R_01(double* x, int i)
{
	double R_i = x[i] / R(x);
	return R_i;
}

double R_02(double* x, int i, int j)
{
	double R_ij = (d[i][j] - (x[i] / R(x)) * (x[j] / R(x))) / R(x);
	return R_ij;
}

double R_03(double* x, int i, int j, int k)
{
	double R_ijk = (3 * (x[i] / R(x)) * (x[j] / R(x)) * (x[k] / R(x)) - d[i][j] * x[k] / R(x) - d[i][k] * x[j] / R(x) - d[j][k] * x[i] / R(x)) / (R(x) * R(x));
	return R_ijk;
}

double R_04(double* x, int i, int j, int k, int l)
{
	double R_ijkl = 1 / (R(x) * R(x) * R(x)) * (-(d[i][j] * d[k][l] + d[i][k] * d[j][l] + d[i][l] * d[j][k])
		+ 3 * (d[i][j] * (x[k] / R(x)) * (x[l] / R(x)) + d[i][k] * (x[j] / R(x)) * (x[l] / R(x)) + d[i][l] * (x[j] / R(x)) * (x[k] / R(x))
			+ d[j][k] * (x[i] / R(x)) * (x[l] / R(x)) + d[j][l] * (x[i] / R(x)) * (x[k] / R(x)) + d[k][l] * (x[i] / R(x)) * (x[j] / R(x)))
		- 15 * (x[i] / R(x)) * (x[j] / R(x)) * (x[k] / R(x)) * (x[l] / R(x)));
	return R_ijkl;
}

double R_05(double* x, int i, int j, int k, int l, int s)
{
	double R_ijkls = 1 / (R(x) * R(x) * R(x) * R(x)) * (105 * R_01(x, i) * R_01(x, j) * R_01(x, k) * R_01(x, l) * R_01(x, s)
		- 15 * (d[i][j] * R_01(x, k) * R_01(x, l) * R_01(x, s) + d[i][k] * R_01(x, j) * R_01(x, l) * R_01(x, s) + d[i][l] * R_01(x, j) * R_01(x, k) * R_01(x, s)
			+ d[i][s] * R_01(x, j) * R_01(x, k) * R_01(x, l) + d[j][k] * R_01(x, i) * R_01(x, l) * R_01(x, s) + d[j][l] * R_01(x, i) * R_01(x, k) * R_01(x, s)
			+ d[j][s] * R_01(x, i) * R_01(x, k) * R_01(x, l) + d[k][l] * R_01(x, i) * R_01(x, j) * R_01(x, s) + d[k][s] * R_01(x, i) * R_01(x, j) * R_01(x, l)
			+ d[l][s] * R_01(x, i) * R_01(x, j) * R_01(x, k))
		+ 3 * (d[i][j] * d[k][l] * R_01(x, s) + d[i][k] * d[j][l] * R_01(x, s) + d[i][l] * d[j][k] * R_01(x, s)
			+ d[i][j] * d[k][s] * R_01(x, l) + d[i][k] * d[j][s] * R_01(x, l) + d[i][s] * d[k][j] * R_01(x, l)
			+ d[i][j] * d[l][s] * R_01(x, k) + d[i][l] * d[j][s] * R_01(x, k) + d[i][s] * d[j][l] * R_01(x, k)
			+ d[i][k] * d[l][s] * R_01(x, j) + d[i][l] * d[k][s] * R_01(x, j) + d[i][s] * d[k][l] * R_01(x, j)
			+ d[j][k] * d[l][s] * R_01(x, i) + d[j][l] * d[k][s] * R_01(x, i) + d[j][s] * d[k][l] * R_01(x, i)));
	return R_ijkls;

}

double R_06(double* x, int i, int j, int k, int l, int s, int t)
{
	double R_ijklst = 1 / (R(x) * R(x) * R(x) * R(x) * R(x)) * (-945 * R_01(x, i) * R_01(x, j) * R_01(x, k) * R_01(x, l) * R_01(x, s) * R_01(x, t)
		+ 105 * (d[i][j] * R_01(x, k) * R_01(x, l) * R_01(x, s) * R_01(x, t) + d[i][k] * R_01(x, j) * R_01(x, l) * R_01(x, s) * R_01(x, t)
			+ d[i][l] * R_01(x, k) * R_01(x, j) * R_01(x, s) * R_01(x, t) + d[i][s] * R_01(x, k) * R_01(x, l) * R_01(x, j) * R_01(x, t)
			+ d[i][t] * R_01(x, k) * R_01(x, l) * R_01(x, s) * R_01(x, j) + d[j][k] * R_01(x, i) * R_01(x, l) * R_01(x, s) * R_01(x, t)
			+ d[l][j] * R_01(x, k) * R_01(x, i) * R_01(x, s) * R_01(x, t) + d[s][j] * R_01(x, k) * R_01(x, l) * R_01(x, i) * R_01(x, t)
			+ d[t][j] * R_01(x, k) * R_01(x, l) * R_01(x, s) * R_01(x, i) + d[k][l] * R_01(x, i) * R_01(x, j) * R_01(x, s) * R_01(x, t)
			+ d[k][s] * R_01(x, i) * R_01(x, j) * R_01(x, l) * R_01(x, t) + d[k][t] * R_01(x, i) * R_01(x, j) * R_01(x, s) * R_01(x, l)
			+ d[l][s] * R_01(x, i) * R_01(x, j) * R_01(x, k) * R_01(x, t) + d[l][t] * R_01(x, i) * R_01(x, j) * R_01(x, k) * R_01(x, s)
			+ d[s][t] * R_01(x, i) * R_01(x, j) * R_01(x, k) * R_01(x, l))
		- 15 * (d[s][t] * d[k][l] * R_01(x, i) * R_01(x, j) + d[s][k] * d[t][l] * R_01(x, i) * R_01(x, j) + d[s][l] * d[t][k] * R_01(x, i) * R_01(x, j)
			+ d[l][j] * d[s][t] * R_01(x, i) * R_01(x, k) + d[l][s] * d[j][t] * R_01(x, i) * R_01(x, k) + d[l][t] * d[j][s] * R_01(x, i) * R_01(x, k)
			+ d[j][k] * d[s][t] * R_01(x, i) * R_01(x, l) + d[j][s] * d[k][t] * R_01(x, i) * R_01(x, l) + d[j][t] * d[k][s] * R_01(x, i) * R_01(x, l)
			+ d[j][k] * d[t][l] * R_01(x, i) * R_01(x, s) + d[j][l] * d[k][t] * R_01(x, i) * R_01(x, s) + d[j][t] * d[k][l] * R_01(x, i) * R_01(x, s)
			+ d[j][k] * d[l][s] * R_01(x, i) * R_01(x, t) + d[j][l] * d[k][s] * R_01(x, i) * R_01(x, t) + d[j][s] * d[k][l] * R_01(x, i) * R_01(x, t)
			+ d[i][l] * d[s][t] * R_01(x, j) * R_01(x, k) + d[i][s] * d[t][l] * R_01(x, j) * R_01(x, k) + d[i][t] * d[l][s] * R_01(x, j) * R_01(x, k)
			+ d[i][k] * d[s][t] * R_01(x, j) * R_01(x, l) + d[i][s] * d[k][t] * R_01(x, j) * R_01(x, l) + d[i][t] * d[k][s] * R_01(x, j) * R_01(x, l)
			+ d[i][k] * d[l][t] * R_01(x, j) * R_01(x, s) + d[i][l] * d[k][t] * R_01(x, j) * R_01(x, s) + d[i][t] * d[k][l] * R_01(x, j) * R_01(x, s)
			+ d[i][k] * d[l][s] * R_01(x, j) * R_01(x, t) + d[i][l] * d[k][s] * R_01(x, j) * R_01(x, t) + d[i][s] * d[k][l] * R_01(x, j) * R_01(x, t)
			+ d[i][j] * d[s][t] * R_01(x, k) * R_01(x, l) + d[i][s] * d[j][t] * R_01(x, k) * R_01(x, l) + d[i][t] * d[j][s] * R_01(x, k) * R_01(x, l)
			+ d[i][j] * d[l][t] * R_01(x, k) * R_01(x, s) + d[i][l] * d[j][t] * R_01(x, k) * R_01(x, s) + d[i][t] * d[j][l] * R_01(x, k) * R_01(x, s)
			+ d[i][j] * d[l][s] * R_01(x, k) * R_01(x, t) + d[i][l] * d[j][s] * R_01(x, k) * R_01(x, t) + d[i][s] * d[j][l] * R_01(x, k) * R_01(x, t)
			+ d[i][j] * d[k][t] * R_01(x, l) * R_01(x, s) + d[i][k] * d[j][t] * R_01(x, l) * R_01(x, s) + d[i][t] * d[k][j] * R_01(x, l) * R_01(x, s)
			+ d[i][j] * d[k][s] * R_01(x, l) * R_01(x, t) + d[i][k] * d[j][s] * R_01(x, l) * R_01(x, t) + d[i][s] * d[k][j] * R_01(x, l) * R_01(x, t)
			+ d[i][j] * d[k][l] * R_01(x, s) * R_01(x, t) + d[i][k] * d[j][l] * R_01(x, s) * R_01(x, t) + d[i][l] * d[k][j] * R_01(x, s) * R_01(x, t))
		+ 3 * (d[i][j] * (d[k][l] * d[s][t] + d[k][s] * d[l][t] + d[k][t] * d[l][s]) + d[i][k] * (d[j][l] * d[s][t] + d[j][s] * d[l][t] + d[j][t] * d[s][l])
			+ d[i][l] * (d[j][k] * d[s][t] + d[j][s] * d[k][t] + d[j][t] * d[k][s]) + d[i][s] * (d[j][k] * d[l][t] + d[j][l] * d[k][t] + d[j][t] * d[l][k])
			+ d[i][t] * (d[j][k] * d[l][s] + d[j][l] * d[k][s] + d[j][s] * d[k][l])));
	return R_ijklst;


}

double R_07(double* x, int i, int j, int k, int l, int s, int t, int p)
{
	double R_ijklstp = 1 / (R(x) * R(x) * R(x) * R(x) * R(x) * R(x)) * (10395 * R_01(x, i) * R_01(x, j) * R_01(x, k) * R_01(x, l) * R_01(x, s) * R_01(x, t) * R_01(x, p)

		- 945 * (d[i][j] * R_01(x, k) * R_01(x, l) * R_01(x, s) * R_01(x, t) * R_01(x, p) + d[i][k] * R_01(x, j) * R_01(x, l) * R_01(x, s) * R_01(x, t) * R_01(x, p)
			+ d[i][l] * R_01(x, k) * R_01(x, j) * R_01(x, s) * R_01(x, t) * R_01(x, p) + d[i][s] * R_01(x, k) * R_01(x, l) * R_01(x, j) * R_01(x, t) * R_01(x, p)
			+ d[i][t] * R_01(x, k) * R_01(x, l) * R_01(x, s) * R_01(x, j) * R_01(x, p) + d[i][p] * R_01(x, k) * R_01(x, l) * R_01(x, s) * R_01(x, j) * R_01(x, t)

			+ d[j][k] * R_01(x, i) * R_01(x, l) * R_01(x, s) * R_01(x, t) * R_01(x, p) + d[l][j] * R_01(x, k) * R_01(x, i) * R_01(x, s) * R_01(x, t) * R_01(x, p)
			+ d[s][j] * R_01(x, k) * R_01(x, l) * R_01(x, i) * R_01(x, t) * R_01(x, p) + d[t][j] * R_01(x, k) * R_01(x, l) * R_01(x, s) * R_01(x, i) * R_01(x, p)
			+ d[p][j] * R_01(x, k) * R_01(x, l) * R_01(x, i) * R_01(x, t) * R_01(x, s)

			+ d[k][l] * R_01(x, i) * R_01(x, j) * R_01(x, s) * R_01(x, t) * R_01(x, p) + d[k][s] * R_01(x, i) * R_01(x, j) * R_01(x, l) * R_01(x, t) * R_01(x, p)
			+ d[k][t] * R_01(x, i) * R_01(x, j) * R_01(x, s) * R_01(x, l) * R_01(x, p) + d[k][p] * R_01(x, i) * R_01(x, j) * R_01(x, s) * R_01(x, l) * R_01(x, t)

			+ d[l][s] * R_01(x, i) * R_01(x, j) * R_01(x, k) * R_01(x, t) * R_01(x, p) + d[l][t] * R_01(x, i) * R_01(x, j) * R_01(x, k) * R_01(x, s) * R_01(x, p)
			+ d[l][p] * R_01(x, i) * R_01(x, j) * R_01(x, k) * R_01(x, t) * R_01(x, s)

			+ d[s][t] * R_01(x, i) * R_01(x, j) * R_01(x, k) * R_01(x, l) * R_01(x, p) + d[s][p] * R_01(x, i) * R_01(x, j) * R_01(x, k) * R_01(x, l) * R_01(x, t)

			+ d[t][p] * R_01(x, i) * R_01(x, j) * R_01(x, k) * R_01(x, l) * R_01(x, s))


		+ 105 * (d[s][t] * d[p][l] * R_01(x, i) * R_01(x, j) * R_01(x, k) + d[s][p] * d[t][l] * R_01(x, i) * R_01(x, j) * R_01(x, k) + d[s][l] * d[t][p] * R_01(x, i) * R_01(x, j) * R_01(x, k)
			+ d[s][t] * d[k][p] * R_01(x, i) * R_01(x, j) * R_01(x, l) + d[s][k] * d[t][p] * R_01(x, i) * R_01(x, j) * R_01(x, l) + d[s][p] * d[t][k] * R_01(x, i) * R_01(x, j) * R_01(x, l)
			+ d[p][t] * d[k][l] * R_01(x, i) * R_01(x, j) * R_01(x, s) + d[p][k] * d[t][l] * R_01(x, i) * R_01(x, j) * R_01(x, s) + d[p][l] * d[t][k] * R_01(x, i) * R_01(x, j) * R_01(x, s)
			+ d[s][p] * d[k][l] * R_01(x, i) * R_01(x, j) * R_01(x, t) + d[s][k] * d[p][l] * R_01(x, i) * R_01(x, j) * R_01(x, t) + d[s][l] * d[p][k] * R_01(x, i) * R_01(x, j) * R_01(x, t)
			+ d[s][t] * d[k][l] * R_01(x, i) * R_01(x, j) * R_01(x, p) + d[s][k] * d[t][l] * R_01(x, i) * R_01(x, j) * R_01(x, p) + d[s][l] * d[t][k] * R_01(x, i) * R_01(x, j) * R_01(x, p)

			+ d[s][t] * d[j][p] * R_01(x, i) * R_01(x, k) * R_01(x, l) + d[s][j] * d[t][p] * R_01(x, i) * R_01(x, k) * R_01(x, l) + d[s][p] * d[t][j] * R_01(x, i) * R_01(x, k) * R_01(x, l)
			+ d[p][t] * d[j][l] * R_01(x, i) * R_01(x, k) * R_01(x, s) + d[p][j] * d[t][l] * R_01(x, i) * R_01(x, k) * R_01(x, s) + d[p][l] * d[t][j] * R_01(x, i) * R_01(x, k) * R_01(x, s)
			+ d[s][p] * d[j][l] * R_01(x, i) * R_01(x, k) * R_01(x, t) + d[s][j] * d[p][l] * R_01(x, i) * R_01(x, k) * R_01(x, t) + d[s][l] * d[p][j] * R_01(x, i) * R_01(x, k) * R_01(x, t)
			+ d[s][t] * d[j][l] * R_01(x, i) * R_01(x, k) * R_01(x, p) + d[s][j] * d[t][l] * R_01(x, i) * R_01(x, k) * R_01(x, p) + d[s][l] * d[t][j] * R_01(x, i) * R_01(x, k) * R_01(x, p)

			+ d[p][t] * d[j][k] * R_01(x, i) * R_01(x, l) * R_01(x, s) + d[p][j] * d[t][k] * R_01(x, i) * R_01(x, l) * R_01(x, s) + d[p][k] * d[t][j] * R_01(x, i) * R_01(x, l) * R_01(x, s)
			+ d[s][p] * d[j][k] * R_01(x, i) * R_01(x, l) * R_01(x, t) + d[s][j] * d[p][k] * R_01(x, i) * R_01(x, l) * R_01(x, t) + d[s][k] * d[p][j] * R_01(x, i) * R_01(x, l) * R_01(x, t)
			+ d[s][t] * d[j][k] * R_01(x, i) * R_01(x, l) * R_01(x, p) + d[s][j] * d[t][k] * R_01(x, i) * R_01(x, l) * R_01(x, p) + d[s][k] * d[t][j] * R_01(x, i) * R_01(x, l) * R_01(x, p)

			+ d[l][p] * d[j][k] * R_01(x, i) * R_01(x, s) * R_01(x, t) + d[l][j] * d[p][k] * R_01(x, i) * R_01(x, s) * R_01(x, t) + d[l][k] * d[p][j] * R_01(x, i) * R_01(x, s) * R_01(x, t)
			+ d[l][t] * d[j][k] * R_01(x, i) * R_01(x, s) * R_01(x, p) + d[l][j] * d[t][k] * R_01(x, i) * R_01(x, s) * R_01(x, p) + d[l][k] * d[t][j] * R_01(x, i) * R_01(x, s) * R_01(x, p)

			+ d[l][s] * d[j][k] * R_01(x, i) * R_01(x, t) * R_01(x, p) + d[l][j] * d[s][k] * R_01(x, i) * R_01(x, t) * R_01(x, p) + d[l][k] * d[s][j] * R_01(x, i) * R_01(x, t) * R_01(x, p)


			+ d[i][p] * d[s][t] * R_01(x, j) * R_01(x, k) * R_01(x, l) + d[i][s] * d[t][p] * R_01(x, j) * R_01(x, k) * R_01(x, l) + d[i][t] * d[p][s] * R_01(x, j) * R_01(x, k) * R_01(x, l)
			+ d[i][p] * d[l][t] * R_01(x, j) * R_01(x, k) * R_01(x, s) + d[i][l] * d[t][p] * R_01(x, j) * R_01(x, k) * R_01(x, s) + d[i][t] * d[p][l] * R_01(x, j) * R_01(x, k) * R_01(x, s)
			+ d[i][p] * d[l][s] * R_01(x, j) * R_01(x, k) * R_01(x, t) + d[i][l] * d[s][p] * R_01(x, j) * R_01(x, k) * R_01(x, t) + d[i][s] * d[p][l] * R_01(x, j) * R_01(x, k) * R_01(x, t)
			+ d[i][l] * d[s][t] * R_01(x, j) * R_01(x, k) * R_01(x, p) + d[i][s] * d[t][l] * R_01(x, j) * R_01(x, k) * R_01(x, p) + d[i][t] * d[l][s] * R_01(x, j) * R_01(x, k) * R_01(x, p)

			+ d[i][p] * d[k][t] * R_01(x, j) * R_01(x, l) * R_01(x, s) + d[i][k] * d[t][p] * R_01(x, j) * R_01(x, l) * R_01(x, s) + d[i][t] * d[p][k] * R_01(x, j) * R_01(x, l) * R_01(x, s)
			+ d[i][p] * d[k][s] * R_01(x, j) * R_01(x, l) * R_01(x, t) + d[i][k] * d[s][p] * R_01(x, j) * R_01(x, l) * R_01(x, t) + d[i][s] * d[p][k] * R_01(x, j) * R_01(x, l) * R_01(x, t)
			+ d[i][k] * d[s][t] * R_01(x, j) * R_01(x, l) * R_01(x, p) + d[i][s] * d[t][k] * R_01(x, j) * R_01(x, l) * R_01(x, p) + d[i][t] * d[k][s] * R_01(x, j) * R_01(x, l) * R_01(x, p)

			+ d[i][p] * d[k][l] * R_01(x, j) * R_01(x, s) * R_01(x, t) + d[i][k] * d[l][p] * R_01(x, j) * R_01(x, s) * R_01(x, t) + d[i][l] * d[p][k] * R_01(x, j) * R_01(x, s) * R_01(x, t)
			+ d[i][k] * d[l][t] * R_01(x, j) * R_01(x, s) * R_01(x, p) + d[i][l] * d[t][k] * R_01(x, j) * R_01(x, s) * R_01(x, p) + d[i][t] * d[k][l] * R_01(x, j) * R_01(x, s) * R_01(x, p)

			+ d[i][k] * d[l][s] * R_01(x, j) * R_01(x, t) * R_01(x, p) + d[i][l] * d[s][k] * R_01(x, j) * R_01(x, t) * R_01(x, p) + d[i][s] * d[k][l] * R_01(x, j) * R_01(x, t) * R_01(x, p)


			+ d[i][j] * d[p][t] * R_01(x, k) * R_01(x, l) * R_01(x, s) + d[i][p] * d[j][t] * R_01(x, k) * R_01(x, l) * R_01(x, s) + d[i][t] * d[j][p] * R_01(x, k) * R_01(x, l) * R_01(x, s)
			+ d[i][j] * d[p][s] * R_01(x, k) * R_01(x, l) * R_01(x, t) + d[i][p] * d[j][s] * R_01(x, k) * R_01(x, l) * R_01(x, t) + d[i][s] * d[j][p] * R_01(x, k) * R_01(x, l) * R_01(x, t)
			+ d[i][j] * d[s][t] * R_01(x, k) * R_01(x, l) * R_01(x, p) + d[i][s] * d[j][t] * R_01(x, k) * R_01(x, l) * R_01(x, p) + d[i][t] * d[j][s] * R_01(x, k) * R_01(x, l) * R_01(x, p)

			+ d[i][j] * d[p][l] * R_01(x, k) * R_01(x, s) * R_01(x, t) + d[i][p] * d[j][l] * R_01(x, k) * R_01(x, s) * R_01(x, t) + d[i][l] * d[j][p] * R_01(x, k) * R_01(x, s) * R_01(x, t)
			+ d[i][j] * d[l][t] * R_01(x, k) * R_01(x, s) * R_01(x, p) + d[i][l] * d[j][t] * R_01(x, k) * R_01(x, s) * R_01(x, p) + d[i][t] * d[j][l] * R_01(x, k) * R_01(x, s) * R_01(x, p)

			+ d[i][j] * d[l][s] * R_01(x, k) * R_01(x, t) * R_01(x, p) + d[i][l] * d[j][s] * R_01(x, k) * R_01(x, t) * R_01(x, p) + d[i][s] * d[j][l] * R_01(x, k) * R_01(x, t) * R_01(x, p)


			+ d[i][j] * d[k][p] * R_01(x, l) * R_01(x, s) * R_01(x, t) + d[i][k] * d[j][p] * R_01(x, l) * R_01(x, s) * R_01(x, t) + d[i][p] * d[k][j] * R_01(x, l) * R_01(x, s) * R_01(x, t)
			+ d[i][j] * d[k][t] * R_01(x, l) * R_01(x, s) * R_01(x, p) + d[i][k] * d[j][t] * R_01(x, l) * R_01(x, s) * R_01(x, p) + d[i][t] * d[k][j] * R_01(x, l) * R_01(x, s) * R_01(x, p)

			+ d[i][j] * d[k][s] * R_01(x, l) * R_01(x, t) * R_01(x, p) + d[i][k] * d[j][s] * R_01(x, l) * R_01(x, t) * R_01(x, p) + d[i][s] * d[k][j] * R_01(x, l) * R_01(x, t) * R_01(x, p)


			+ d[i][j] * d[k][l] * R_01(x, s) * R_01(x, t) * R_01(x, p) + d[i][k] * d[j][l] * R_01(x, s) * R_01(x, t) * R_01(x, p) + d[i][l] * d[k][j] * R_01(x, s) * R_01(x, t) * R_01(x, p))


		- 15 * ((d[i][j] * (d[k][l] * d[s][t] + d[k][s] * d[l][t] + d[k][t] * d[l][s]) + d[i][k] * (d[j][l] * d[s][t] + d[j][s] * d[l][t] + d[j][t] * d[s][l])
			+ d[i][l] * (d[j][k] * d[s][t] + d[j][s] * d[k][t] + d[j][t] * d[k][s]) + d[i][s] * (d[j][k] * d[l][t] + d[j][l] * d[k][t] + d[j][t] * d[l][k])
			+ d[i][t] * (d[j][k] * d[l][s] + d[j][l] * d[k][s] + d[j][s] * d[k][l])) * R_01(x, p)

			+ (d[p][j] * (d[k][l] * d[s][t] + d[k][s] * d[l][t] + d[k][t] * d[l][s]) + d[p][k] * (d[j][l] * d[s][t] + d[j][s] * d[l][t] + d[j][t] * d[s][l])
				+ d[p][l] * (d[j][k] * d[s][t] + d[j][s] * d[k][t] + d[j][t] * d[k][s]) + d[p][s] * (d[j][k] * d[l][t] + d[j][l] * d[k][t] + d[j][t] * d[l][k])
				+ d[p][t] * (d[j][k] * d[l][s] + d[j][l] * d[k][s] + d[j][s] * d[k][l])) * R_01(x, i)

			+ (d[p][i] * (d[k][l] * d[s][t] + d[k][s] * d[l][t] + d[k][t] * d[l][s]) + d[p][k] * (d[i][l] * d[s][t] + d[i][s] * d[l][t] + d[i][t] * d[s][l])
				+ d[p][l] * (d[i][k] * d[s][t] + d[i][s] * d[k][t] + d[i][t] * d[k][s]) + d[p][s] * (d[i][k] * d[l][t] + d[i][l] * d[k][t] + d[i][t] * d[l][k])
				+ d[p][t] * (d[i][k] * d[l][s] + d[i][l] * d[k][s] + d[i][s] * d[k][l])) * R_01(x, j)

			+ (d[p][i] * (d[j][l] * d[s][t] + d[j][s] * d[l][t] + d[j][t] * d[l][s]) + d[p][j] * (d[i][l] * d[s][t] + d[i][s] * d[l][t] + d[i][t] * d[s][l])
				+ d[p][l] * (d[i][j] * d[s][t] + d[i][s] * d[j][t] + d[i][t] * d[j][s]) + d[p][s] * (d[i][j] * d[l][t] + d[i][l] * d[j][t] + d[i][t] * d[l][j])
				+ d[p][t] * (d[i][j] * d[l][s] + d[i][l] * d[j][s] + d[i][s] * d[j][l])) * R_01(x, k)

			+ (d[p][i] * (d[j][k] * d[s][t] + d[j][s] * d[k][t] + d[j][t] * d[k][s]) + d[p][j] * (d[i][k] * d[s][t] + d[i][s] * d[k][t] + d[i][t] * d[s][k])
				+ d[p][k] * (d[i][j] * d[s][t] + d[i][s] * d[j][t] + d[i][t] * d[j][s]) + d[p][s] * (d[i][j] * d[k][t] + d[i][k] * d[j][t] + d[i][t] * d[k][j])
				+ d[p][t] * (d[i][j] * d[k][s] + d[i][k] * d[j][s] + d[i][s] * d[j][k])) * R_01(x, l)

			+ (d[p][i] * (d[j][k] * d[l][t] + d[j][l] * d[k][t] + d[j][t] * d[k][l]) + d[p][j] * (d[i][k] * d[l][t] + d[i][l] * d[k][t] + d[i][t] * d[l][k])
				+ d[p][k] * (d[i][j] * d[l][t] + d[i][l] * d[j][t] + d[i][t] * d[j][l]) + d[p][l] * (d[i][j] * d[k][t] + d[i][k] * d[j][t] + d[i][t] * d[k][j])
				+ d[p][t] * (d[i][j] * d[k][l] + d[i][k] * d[j][l] + d[i][l] * d[j][k])) * R_01(x, s)

			+ (d[p][i] * (d[j][k] * d[l][s] + d[j][l] * d[k][s] + d[j][s] * d[k][l]) + d[p][j] * (d[i][k] * d[l][s] + d[i][l] * d[k][s] + d[i][s] * d[l][k])
				+ d[p][k] * (d[i][j] * d[l][s] + d[i][l] * d[j][s] + d[i][s] * d[j][l]) + d[p][l] * (d[i][j] * d[k][s] + d[i][k] * d[j][s] + d[i][s] * d[k][j])
				+ d[p][s] * (d[i][j] * d[k][l] + d[i][k] * d[j][l] + d[i][l] * d[j][k])) * R_01(x, t))


		);
	return R_ijklstp;


}

double fR(double* x, int sign)
{
	double R1 = sign * x[2] + R(x);
	double fR = 1 / R1;
	return fR;
}

double fR_01(double* x, int i, int sign)
{
	double R1 = sign * x[2] + R(x);
	double fR_i = -(R_01(x, i) + sign * d[2][i]) / (R1 * R1);
	return fR_i;
}

double fR_02(double* x, int i, int j, int sign)
{
	double R1 = R(x) + sign * x[2];
	double fR_ij = -1.00 / (R1 * R1 * R1) * (R1 * R_02(x, i, j) - 2 * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]));
	return fR_ij;
}

double fR_03(double* x, int i, int j, int k, int sign)
{
	double R1 = R(x) + sign * x[2];
	double fR_ijk = -1.00 / (R1 * R1 * R1 * R1) * (R1 * R1 * R_03(x, i, j, k) - 2 * R1 * ((R_01(x, i) + sign * d[2][i]) * R_02(x, j, k) + (R_01(x, j) + sign * d[2][j]) * R_02(x, i, k) + (R_01(x, k) + sign * d[2][k]) * R_02(x, i, j))
		+ 6 * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]));
	return fR_ijk;
}

double fR_04(double* x, int i, int j, int k, int l, int sign)
{
	double R1 = R(x) + sign * x[2];
	double fR_ijkl = -1.00 / (R1 * R1 * R1 * R1 * R1) * (R1 * R1 * R1 * R_04(x, i, j, k, l) - 2 * R1 * R1 * ((R_01(x, i) + sign * d[2][i]) * R_03(x, j, k, l) + (R_01(x, j) + sign * d[2][j]) * R_03(x, i, k, l)
		+ (R_01(x, k) + sign * d[2][k]) * R_03(x, i, j, l) + (R_01(x, l) + sign * d[2][l]) * R_03(x, i, j, k)

		+ R_02(x, i, j) * R_02(x, k, l) + R_02(x, i, k) * R_02(x, j, l) + R_02(x, i, l) * R_02(x, j, k)
		)
		+ 6 * R1 * ((R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * R_02(x, k, l) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * R_02(x, j, l)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * R_02(x, k, j) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * R_02(x, i, l)
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * R_02(x, k, i) + (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * R_02(x, i, j))
		- 24 * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]));
	return fR_ijkl;

}

double fR_05(double* x, int i, int j, int k, int l, int s, int sign)
{
	double R1 = R(x) + sign * x[2];
	double fR_ijkls = -1.00 / (R1 * R1 * R1 * R1 * R1 * R1) * (R1 * R1 * R1 * R1 * R_05(x, i, j, k, l, s)

		- 2 * R1 * R1 * R1 * ((R_01(x, i) + sign * d[2][i]) * R_04(x, j, k, l, s) + (R_01(x, j) + sign * d[2][j]) * R_04(x, i, k, l, s)
			+ (R_01(x, k) + sign * d[2][k]) * R_04(x, i, j, l, s) + (R_01(x, l) + sign * d[2][l]) * R_04(x, i, j, k, s)
			+ (R_01(x, s) + sign * d[2][s]) * R_04(x, i, j, k, l)

			+ R_02(x, i, j) * R_03(x, k, l, s) + R_02(x, i, k) * R_03(x, j, l, s) + R_02(x, i, l) * R_03(x, k, j, s)
			+ R_02(x, i, s) * R_03(x, k, l, j)

			+ R_02(x, j, k) * R_03(x, i, l, s) + R_02(x, j, l) * R_03(x, k, i, s) + R_02(x, j, s) * R_03(x, k, i, l)

			+ R_02(x, k, l) * R_03(x, i, j, s) + R_02(x, k, s) * R_03(x, i, j, l)

			+ R_02(x, l, s) * R_03(x, k, i, j)

			)
		+ 6 * R1 * R1 * ((R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * R_03(x, k, l, s) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * R_03(x, j, l, s)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * R_03(x, k, j, s) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s]) * R_03(x, k, j, l)
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * R_03(x, i, l, s) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * R_03(x, k, i, s)
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * R_03(x, i, l, k) + (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * R_03(x, i, j, s)
			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * R_03(x, i, l, j) + (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * R_03(x, i, j, k)

			+ (R_01(x, i) + sign * d[2][i]) * (R_02(x, k, j) * R_02(x, l, s) + R_02(x, j, l) * R_02(x, k, s) + R_02(x, s, j) * R_02(x, l, k))
			+ (R_01(x, j) + sign * d[2][j]) * (R_02(x, k, i) * R_02(x, l, s) + R_02(x, i, l) * R_02(x, k, s) + R_02(x, s, i) * R_02(x, l, k))
			+ (R_01(x, k) + sign * d[2][k]) * (R_02(x, j, i) * R_02(x, l, s) + R_02(x, i, l) * R_02(x, j, s) + R_02(x, s, i) * R_02(x, l, j))
			+ (R_01(x, l) + sign * d[2][l]) * (R_02(x, j, i) * R_02(x, k, s) + R_02(x, i, k) * R_02(x, j, s) + R_02(x, s, i) * R_02(x, k, j))
			+ (R_01(x, s) + sign * d[2][s]) * (R_02(x, j, i) * R_02(x, k, l) + R_02(x, i, k) * R_02(x, j, l) + R_02(x, l, i) * R_02(x, k, j))


			)
		- 24 * R1 * (R_02(x, i, j) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l])
			+ R_02(x, i, k) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l])
			+ R_02(x, i, l) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, j) + sign * d[2][j])
			+ R_02(x, i, s) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l])
			+ R_02(x, k, j) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l])
			+ R_02(x, j, l) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k])
			+ R_02(x, s, j) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l])
			+ R_02(x, l, k) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, i) + sign * d[2][i])
			+ R_02(x, k, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l])
			+ R_02(x, l, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]))
		+ 120 * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]));
	return fR_ijkls;
}

double fR_06(double* x, int i, int j, int k, int l, int s, int t, int sign)
{
	double R1 = R(x) + sign * x[2];
	double fR_ijklst = -1.00 / (R1 * R1 * R1 * R1 * R1 * R1 * R1) * (R1 * R1 * R1 * R1 * R1 * R_06(x, i, j, k, l, s, t)

		- 2 * R1 * R1 * R1 * R1 * ((R_01(x, t) + sign * d[2][t]) * R_05(x, i, j, k, l, s) + (R_01(x, i) + sign * d[2][i]) * R_05(x, j, k, l, s, t) + (R_01(x, j) + sign * d[2][j]) * R_05(x, i, k, l, s, t)
			+ (R_01(x, k) + sign * d[2][k]) * R_05(x, i, j, l, s, t) + (R_01(x, l) + sign * d[2][l]) * R_05(x, i, j, k, s, t) + (R_01(x, s) + sign * d[2][s]) * R_05(x, i, j, k, l, t)


			+ R_02(x, i, j) * R_04(x, k, l, s, t) + R_02(x, i, k) * R_04(x, j, l, s, t) + R_02(x, i, l) * R_04(x, k, j, s, t) + R_02(x, i, s) * R_04(x, k, l, j, t) + R_02(x, i, t) * R_04(x, k, l, s, j)
			+ R_02(x, j, k) * R_04(x, i, l, s, t) + R_02(x, j, l) * R_04(x, i, k, s, t) + R_02(x, j, s) * R_04(x, i, l, k, t) + +R_02(x, j, t) * R_04(x, i, l, s, k)
			+ R_02(x, k, l) * R_04(x, i, j, s, t) + R_02(x, k, s) * R_04(x, i, j, l, t) + R_02(x, k, t) * R_04(x, i, j, s, l)
			+ R_02(x, l, s) * R_04(x, i, j, k, t) + R_02(x, l, t) * R_04(x, i, j, k, s) + R_02(x, s, t) * R_04(x, i, j, k, l)


			+ R_03(x, i, j, k) * R_03(x, l, s, t) + R_03(x, i, j, l) * R_03(x, k, s, t) + R_03(x, i, j, s) * R_03(x, l, k, t) + R_03(x, i, j, t) * R_03(x, l, s, k)
			+ R_03(x, i, k, l) * R_03(x, j, s, t) + R_03(x, i, k, s) * R_03(x, l, j, t) + R_03(x, i, k, t) * R_03(x, j, s, l)
			+ R_03(x, i, l, s) * R_03(x, j, k, t) + R_03(x, i, l, t) * R_03(x, j, k, s) + R_03(x, i, s, t) * R_03(x, j, k, l)

			)

		+ 6 * R1 * R1 * R1 * ((R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * R_04(x, k, l, s, t) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * R_04(x, j, l, s, t) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * R_04(x, j, k, s, t)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s]) * R_04(x, j, k, l, t) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, j, k, l, s) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * R_04(x, i, l, s, t)
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * R_04(x, i, k, s, t) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * R_04(x, i, k, l, t) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, i, k, l, s)
			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * R_04(x, i, j, s, t) + (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * R_04(x, i, j, l, t) + (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, i, j, l, s)
			+ (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * R_04(x, i, j, k, t) + (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, i, j, k, s) + (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, i, j, k, l)

			+ (R_01(x, i) + sign * d[2][i]) * (R_02(x, j, k) * R_03(x, l, s, t) + R_02(x, j, l) * R_03(x, k, s, t) + R_02(x, j, s) * R_03(x, l, k, t) + R_02(x, j, t) * R_03(x, l, s, k)
				+ R_02(x, k, l) * R_03(x, j, s, t) + R_02(x, k, s) * R_03(x, l, j, t) + R_02(x, k, t) * R_03(x, l, j, s)
				+ R_02(x, s, l) * R_03(x, j, k, t) + R_02(x, l, t) * R_03(x, k, j, s) + R_02(x, s, t) * R_03(x, l, j, k)
				)
			+ (R_01(x, j) + sign * d[2][j]) * (R_02(x, i, k) * R_03(x, l, s, t) + R_02(x, i, l) * R_03(x, k, s, t) + R_02(x, i, s) * R_03(x, l, k, t) + R_02(x, i, t) * R_03(x, l, s, k)
				+ R_02(x, k, l) * R_03(x, i, s, t) + R_02(x, k, s) * R_03(x, l, i, t) + R_02(x, k, t) * R_03(x, l, i, s)
				+ R_02(x, s, l) * R_03(x, i, k, t) + R_02(x, l, t) * R_03(x, k, i, s) + R_02(x, s, t) * R_03(x, l, i, k)
				)
			+ (R_01(x, k) + sign * d[2][k]) * (R_02(x, i, j) * R_03(x, l, s, t) + R_02(x, i, l) * R_03(x, j, s, t) + R_02(x, i, s) * R_03(x, l, j, t) + R_02(x, i, t) * R_03(x, l, s, j)
				+ R_02(x, j, l) * R_03(x, i, s, t) + R_02(x, j, s) * R_03(x, l, i, t) + R_02(x, j, t) * R_03(x, l, i, s)
				+ R_02(x, s, l) * R_03(x, i, j, t) + R_02(x, l, t) * R_03(x, j, i, s) + R_02(x, s, t) * R_03(x, l, i, j)
				)
			+ (R_01(x, l) + sign * d[2][l]) * (R_02(x, i, j) * R_03(x, k, s, t) + R_02(x, i, k) * R_03(x, j, s, t) + R_02(x, i, s) * R_03(x, k, j, t) + R_02(x, i, t) * R_03(x, k, s, j)
				+ R_02(x, j, k) * R_03(x, i, s, t) + R_02(x, j, s) * R_03(x, k, i, t) + R_02(x, j, t) * R_03(x, k, i, s)
				+ R_02(x, s, k) * R_03(x, i, j, t) + R_02(x, k, t) * R_03(x, j, i, s) + R_02(x, s, t) * R_03(x, k, i, j)
				)
			+ (R_01(x, s) + sign * d[2][s]) * (R_02(x, i, j) * R_03(x, k, l, t) + R_02(x, i, k) * R_03(x, j, l, t) + R_02(x, i, l) * R_03(x, k, j, t) + R_02(x, i, t) * R_03(x, k, l, j)
				+ R_02(x, j, k) * R_03(x, i, l, t) + R_02(x, j, l) * R_03(x, k, i, t) + R_02(x, j, t) * R_03(x, k, i, l)
				+ R_02(x, l, k) * R_03(x, i, j, t) + R_02(x, k, t) * R_03(x, j, i, l) + R_02(x, l, t) * R_03(x, k, i, j)
				)
			+ (R_01(x, t) + sign * d[2][t]) * (R_02(x, i, j) * R_03(x, k, l, s) + R_02(x, i, k) * R_03(x, j, l, s) + R_02(x, i, l) * R_03(x, k, j, s) + R_02(x, i, s) * R_03(x, k, l, j)
				+ R_02(x, j, k) * R_03(x, i, l, s) + R_02(x, j, l) * R_03(x, k, i, s) + R_02(x, j, s) * R_03(x, k, i, l)
				+ R_02(x, l, k) * R_03(x, i, j, s) + R_02(x, k, s) * R_03(x, j, i, l) + R_02(x, l, s) * R_03(x, k, i, j)
				)

			+ R_02(x, i, j) * (R_02(x, k, l) * R_02(x, s, t) + R_02(x, k, s) * R_02(x, l, t) + R_02(x, k, t) * R_02(x, s, l))
			+ R_02(x, i, k) * (R_02(x, j, l) * R_02(x, s, t) + R_02(x, j, s) * R_02(x, l, t) + R_02(x, j, t) * R_02(x, s, l))
			+ R_02(x, i, l) * (R_02(x, k, j) * R_02(x, s, t) + R_02(x, k, s) * R_02(x, j, t) + R_02(x, k, t) * R_02(x, s, j))
			+ R_02(x, i, s) * (R_02(x, k, l) * R_02(x, j, t) + R_02(x, k, j) * R_02(x, l, t) + R_02(x, k, t) * R_02(x, j, l))
			+ R_02(x, i, t) * (R_02(x, k, l) * R_02(x, s, j) + R_02(x, k, s) * R_02(x, l, j) + R_02(x, k, j) * R_02(x, s, l))

			)


		- 24 * R1 * R1 * ((R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * R_03(x, l, s, t) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * R_03(x, k, s, t)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * R_03(x, k, l, t) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, t) + sign * d[2][t]) * R_03(x, k, l, s)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * R_03(x, j, s, t) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * R_03(x, j, l, t)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * R_03(x, j, l, s) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * R_03(x, k, j, t)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * R_03(x, j, k, s) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * R_03(x, j, l, k)

			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * R_03(x, i, s, t) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * R_03(x, i, l, t)
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * R_03(x, i, l, s) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * R_03(x, i, k, t)
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * R_03(x, i, k, s) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * R_03(x, i, l, k)

			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * R_03(x, i, j, t) + (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * R_03(x, i, j, s)
			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * R_03(x, i, j, l) + (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * R_03(x, i, j, k)

			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_02(x, k, l) * R_02(x, s, t) + R_02(x, k, s) * R_02(x, l, t) + R_02(x, k, t) * R_02(x, s, l))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_02(x, j, l) * R_02(x, s, t) + R_02(x, j, s) * R_02(x, l, t) + R_02(x, j, t) * R_02(x, s, l))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_02(x, j, k) * R_02(x, s, t) + R_02(x, j, s) * R_02(x, k, t) + R_02(x, j, t) * R_02(x, s, k))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s]) * (R_02(x, j, k) * R_02(x, l, t) + R_02(x, j, l) * R_02(x, k, t) + R_02(x, j, t) * R_02(x, l, k))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, j, k) * R_02(x, l, s) + R_02(x, j, l) * R_02(x, k, s) + R_02(x, j, s) * R_02(x, l, k))

			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_02(x, i, l) * R_02(x, s, t) + R_02(x, i, s) * R_02(x, l, t) + R_02(x, i, t) * R_02(x, s, l))
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_02(x, i, k) * R_02(x, s, t) + R_02(x, i, s) * R_02(x, k, t) + R_02(x, i, t) * R_02(x, s, k))
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_02(x, i, k) * R_02(x, l, t) + R_02(x, i, l) * R_02(x, k, t) + R_02(x, i, t) * R_02(x, l, k))
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, i, k) * R_02(x, l, s) + R_02(x, i, l) * R_02(x, k, s) + R_02(x, i, s) * R_02(x, l, k))

			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_02(x, i, j) * R_02(x, s, t) + R_02(x, i, s) * R_02(x, j, t) + R_02(x, i, t) * R_02(x, s, j))
			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_02(x, i, j) * R_02(x, l, t) + R_02(x, i, l) * R_02(x, j, t) + R_02(x, i, t) * R_02(x, l, j))
			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, i, j) * R_02(x, l, s) + R_02(x, i, l) * R_02(x, j, s) + R_02(x, i, s) * R_02(x, l, j))

			+ (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_02(x, i, j) * R_02(x, k, t) + R_02(x, i, k) * R_02(x, j, t) + R_02(x, i, t) * R_02(x, k, j))
			+ (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, i, j) * R_02(x, k, s) + R_02(x, i, k) * R_02(x, j, s) + R_02(x, i, s) * R_02(x, k, j))

			+ (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, i, j) * R_02(x, k, l) + R_02(x, i, k) * R_02(x, j, l) + R_02(x, i, l) * R_02(x, k, j))


			)


		+ 120 * R1 * (R_02(x, i, j) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) + R_02(x, i, k) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])
			+ R_02(x, i, l) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) + R_02(x, i, s) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, t) + sign * d[2][t])
			+ R_02(x, i, t) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, j, k) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])
			+ R_02(x, j, l) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) + R_02(x, j, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t])
			+ R_02(x, j, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) + R_02(x, k, l) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])
			+ R_02(x, k, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) + R_02(x, k, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s])
			+ R_02(x, l, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) + R_02(x, l, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s])
			+ R_02(x, t, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]))

		- 720 * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]));
	return fR_ijklst;

}

double fR_07(double* x, int i, int j, int k, int l, int s, int t, int p, int sign)
{
	double R1 = R(x) + sign * x[2];
	double fR_ijklstp = -1.00 / (R1 * R1 * R1 * R1 * R1 * R1 * R1 * R1) * (R1 * R1 * R1 * R1 * R1 * R1 * R_07(x, i, j, k, l, s, t, p)
		- 2 * R1 * R1 * R1 * R1 * R1 * ((R_01(x, i) + sign * d[2][i]) * R_06(x, j, k, l, s, t, p) + (R_01(x, j) + sign * d[2][j]) * R_06(x, i, k, l, s, t, p) + (R_01(x, k) + sign * d[2][k]) * R_06(x, i, j, l, s, t, p)
			+ (R_01(x, l) + sign * d[2][l]) * R_06(x, i, j, k, s, t, p) + (R_01(x, s) + sign * d[2][s]) * R_06(x, i, j, k, l, t, p) + (R_01(x, t) + sign * d[2][t]) * R_06(x, i, j, k, l, s, p)
			+ (R_01(x, p) + sign * d[2][p]) * R_06(x, i, j, k, l, s, t)


			+ R_02(x, i, j) * R_05(x, k, l, s, t, p) + R_02(x, i, k) * R_05(x, j, l, s, t, p) + R_02(x, i, l) * R_05(x, j, k, s, t, p) + R_02(x, i, s) * R_05(x, j, k, l, t, p)
			+ R_02(x, i, t) * R_05(x, j, k, l, s, p) + R_02(x, i, p) * R_05(x, j, k, l, s, t)

			+ R_02(x, j, k) * R_05(x, i, l, s, t, p) + R_02(x, j, l) * R_05(x, i, k, s, t, p) + R_02(x, j, s) * R_05(x, i, k, l, t, p) + R_02(x, j, t) * R_05(x, i, k, l, s, p)
			+ R_02(x, j, p) * R_05(x, i, k, l, s, t)

			+ R_02(x, k, l) * R_05(x, i, j, s, t, p) + R_02(x, k, s) * R_05(x, i, j, l, t, p) + R_02(x, k, t) * R_05(x, i, j, l, s, p) + R_02(x, k, p) * R_05(x, i, j, l, s, t)

			+ R_02(x, l, s) * R_05(x, i, j, k, t, p) + R_02(x, l, t) * R_05(x, i, j, k, s, p) + R_02(x, l, p) * R_05(x, i, j, k, s, t)

			+ R_02(x, s, t) * R_05(x, i, j, k, l, p) + R_02(x, s, p) * R_05(x, i, j, k, l, t)

			+ R_02(x, t, p) * R_05(x, i, j, k, l, s)


			+ R_03(x, i, j, k) * R_04(x, l, s, t, p) + R_03(x, i, j, l) * R_04(x, k, s, t, p) + R_03(x, i, j, s) * R_04(x, l, k, t, p) + R_03(x, i, j, t) * R_04(x, l, s, k, p)
			+ R_03(x, i, j, p) * R_04(x, l, s, t, k)
			+ R_03(x, i, k, l) * R_04(x, j, s, t, p) + R_03(x, i, k, s) * R_04(x, l, j, t, p) + R_03(x, i, k, t) * R_04(x, l, s, j, p) + R_03(x, i, k, p) * R_04(x, l, s, t, j)
			+ R_03(x, i, l, s) * R_04(x, k, j, t, p) + R_03(x, i, l, t) * R_04(x, k, s, j, p) + R_03(x, i, l, p) * R_04(x, k, s, t, j)
			+ R_03(x, i, s, t) * R_04(x, k, l, j, p) + R_03(x, i, s, p) * R_04(x, k, l, t, j)
			+ R_03(x, i, t, p) * R_04(x, k, l, s, j)

			+ R_03(x, j, k, l) * R_04(x, i, s, t, p) + R_03(x, j, k, s) * R_04(x, l, i, t, p) + R_03(x, j, k, t) * R_04(x, l, s, i, p) + R_03(x, j, k, p) * R_04(x, l, s, t, i)
			+ R_03(x, j, l, s) * R_04(x, k, i, t, p) + R_03(x, j, l, t) * R_04(x, k, s, i, p) + R_03(x, j, l, p) * R_04(x, k, s, t, i)
			+ R_03(x, j, s, t) * R_04(x, k, l, i, p) + R_03(x, j, s, p) * R_04(x, k, l, t, i)
			+ R_03(x, j, t, p) * R_04(x, k, l, s, i)

			+ R_03(x, k, l, s) * R_04(x, j, i, t, p) + R_03(x, k, l, t) * R_04(x, j, s, i, p) + R_03(x, k, l, p) * R_04(x, j, s, t, i)
			+ R_03(x, k, s, t) * R_04(x, j, l, i, p) + R_03(x, k, s, p) * R_04(x, j, l, t, i)
			+ R_03(x, k, t, p) * R_04(x, j, l, s, i)

			+ R_03(x, l, s, t) * R_04(x, j, k, i, p) + R_03(x, l, s, p) * R_04(x, j, k, t, i)
			+ R_03(x, l, t, p) * R_04(x, j, k, s, i)

			+ R_03(x, s, t, p) * R_04(x, j, k, l, i)


			)

		+ 6 * R1 * R1 * R1 * R1 * ((R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * R_05(x, k, l, s, t, p) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * R_05(x, j, l, s, t, p)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * R_05(x, k, j, s, t, p) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s]) * R_05(x, k, j, l, t, p)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, t) + sign * d[2][t]) * R_05(x, k, j, l, s, p) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) * R_05(x, k, j, l, s, t)

			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * R_05(x, i, l, s, t, p) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * R_05(x, k, i, s, t, p)
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * R_05(x, i, l, k, t, p) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, t) + sign * d[2][t]) * R_05(x, k, i, s, l, p)
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) * R_05(x, i, l, k, t, s)

			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * R_05(x, i, j, s, t, p) + (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * R_05(x, i, l, j, t, p)
			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * R_05(x, i, l, j, s, p) + (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) * R_05(x, i, l, j, s, t)

			+ (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * R_05(x, i, j, k, t, p) + (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * R_05(x, i, j, k, s, p)
			+ (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) * R_05(x, i, j, k, t, s)

			+ (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * R_05(x, i, j, k, l, p) + (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) * R_05(x, i, j, k, l, t)

			+ (R_01(x, p) + sign * d[2][p]) * (R_01(x, t) + sign * d[2][t]) * R_05(x, i, j, k, l, s)


			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, j, k) + (R_01(x, j) + sign * d[2][j]) * R_02(x, i, k) + (R_01(x, k) + sign * d[2][k]) * R_02(x, j, i)) * R_04(x, p, l, s, t)
			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, j, l) + (R_01(x, j) + sign * d[2][j]) * R_02(x, i, l) + (R_01(x, l) + sign * d[2][l]) * R_02(x, j, i)) * R_04(x, p, k, s, t)
			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, j, s) + (R_01(x, j) + sign * d[2][j]) * R_02(x, i, s) + (R_01(x, s) + sign * d[2][s]) * R_02(x, j, i)) * R_04(x, p, k, l, t)
			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, j, t) + (R_01(x, j) + sign * d[2][j]) * R_02(x, i, t) + (R_01(x, t) + sign * d[2][t]) * R_02(x, j, i)) * R_04(x, p, k, l, s)
			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, j, p) + (R_01(x, j) + sign * d[2][j]) * R_02(x, i, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, j, i)) * R_04(x, t, k, l, s)

			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, k, l) + (R_01(x, k) + sign * d[2][k]) * R_02(x, i, l) + (R_01(x, l) + sign * d[2][l]) * R_02(x, k, i)) * R_04(x, p, j, s, t)
			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, k, s) + (R_01(x, k) + sign * d[2][k]) * R_02(x, i, s) + (R_01(x, s) + sign * d[2][s]) * R_02(x, k, i)) * R_04(x, p, j, l, t)
			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, k, t) + (R_01(x, k) + sign * d[2][k]) * R_02(x, i, t) + (R_01(x, t) + sign * d[2][t]) * R_02(x, k, i)) * R_04(x, p, j, l, s)
			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, k, p) + (R_01(x, k) + sign * d[2][k]) * R_02(x, i, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, k, i)) * R_04(x, t, j, l, s)

			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, l, s) + (R_01(x, l) + sign * d[2][l]) * R_02(x, i, s) + (R_01(x, s) + sign * d[2][s]) * R_02(x, l, i)) * R_04(x, p, j, k, t)
			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, l, t) + (R_01(x, l) + sign * d[2][l]) * R_02(x, i, t) + (R_01(x, t) + sign * d[2][t]) * R_02(x, l, i)) * R_04(x, p, j, k, s)
			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, l, p) + (R_01(x, l) + sign * d[2][l]) * R_02(x, i, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, l, i)) * R_04(x, t, j, k, s)

			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, s, t) + (R_01(x, s) + sign * d[2][s]) * R_02(x, i, t) + (R_01(x, t) + sign * d[2][t]) * R_02(x, s, i)) * R_04(x, p, j, k, l)
			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, s, p) + (R_01(x, s) + sign * d[2][s]) * R_02(x, i, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, s, i)) * R_04(x, t, j, k, l)

			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, t, p) + (R_01(x, t) + sign * d[2][t]) * R_02(x, i, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, t, i)) * R_04(x, s, j, k, l)


			+ ((R_01(x, j) + sign * d[2][j]) * R_02(x, k, l) + (R_01(x, k) + sign * d[2][k]) * R_02(x, j, l) + (R_01(x, l) + sign * d[2][l]) * R_02(x, k, j)) * R_04(x, p, i, s, t)
			+ ((R_01(x, j) + sign * d[2][j]) * R_02(x, k, s) + (R_01(x, k) + sign * d[2][k]) * R_02(x, j, s) + (R_01(x, s) + sign * d[2][s]) * R_02(x, k, j)) * R_04(x, p, i, l, t)
			+ ((R_01(x, j) + sign * d[2][j]) * R_02(x, k, t) + (R_01(x, k) + sign * d[2][k]) * R_02(x, j, t) + (R_01(x, t) + sign * d[2][t]) * R_02(x, k, j)) * R_04(x, p, i, l, s)
			+ ((R_01(x, j) + sign * d[2][j]) * R_02(x, k, p) + (R_01(x, k) + sign * d[2][k]) * R_02(x, j, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, k, j)) * R_04(x, t, i, l, s)

			+ ((R_01(x, j) + sign * d[2][j]) * R_02(x, l, s) + (R_01(x, l) + sign * d[2][l]) * R_02(x, j, s) + (R_01(x, s) + sign * d[2][s]) * R_02(x, l, j)) * R_04(x, p, i, k, t)
			+ ((R_01(x, j) + sign * d[2][j]) * R_02(x, l, t) + (R_01(x, l) + sign * d[2][l]) * R_02(x, j, t) + (R_01(x, t) + sign * d[2][t]) * R_02(x, l, j)) * R_04(x, p, i, k, s)
			+ ((R_01(x, j) + sign * d[2][j]) * R_02(x, l, p) + (R_01(x, l) + sign * d[2][l]) * R_02(x, j, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, l, j)) * R_04(x, t, i, k, s)

			+ ((R_01(x, j) + sign * d[2][j]) * R_02(x, s, t) + (R_01(x, s) + sign * d[2][s]) * R_02(x, j, t) + (R_01(x, t) + sign * d[2][t]) * R_02(x, s, j)) * R_04(x, p, i, k, l)
			+ ((R_01(x, j) + sign * d[2][j]) * R_02(x, s, p) + (R_01(x, s) + sign * d[2][s]) * R_02(x, j, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, s, j)) * R_04(x, t, i, k, l)

			+ ((R_01(x, j) + sign * d[2][j]) * R_02(x, t, p) + (R_01(x, t) + sign * d[2][t]) * R_02(x, j, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, t, j)) * R_04(x, s, i, k, l)


			+ ((R_01(x, k) + sign * d[2][k]) * R_02(x, l, s) + (R_01(x, l) + sign * d[2][l]) * R_02(x, k, s) + (R_01(x, s) + sign * d[2][s]) * R_02(x, l, k)) * R_04(x, p, i, j, t)
			+ ((R_01(x, k) + sign * d[2][k]) * R_02(x, l, t) + (R_01(x, l) + sign * d[2][l]) * R_02(x, k, t) + (R_01(x, t) + sign * d[2][t]) * R_02(x, l, k)) * R_04(x, p, i, j, s)
			+ ((R_01(x, k) + sign * d[2][k]) * R_02(x, l, p) + (R_01(x, l) + sign * d[2][l]) * R_02(x, k, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, l, k)) * R_04(x, t, i, j, s)

			+ ((R_01(x, k) + sign * d[2][k]) * R_02(x, s, t) + (R_01(x, s) + sign * d[2][s]) * R_02(x, k, t) + (R_01(x, t) + sign * d[2][t]) * R_02(x, s, k)) * R_04(x, p, i, j, l)
			+ ((R_01(x, k) + sign * d[2][k]) * R_02(x, s, p) + (R_01(x, s) + sign * d[2][s]) * R_02(x, k, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, s, k)) * R_04(x, t, i, j, l)

			+ ((R_01(x, k) + sign * d[2][k]) * R_02(x, t, p) + (R_01(x, t) + sign * d[2][t]) * R_02(x, k, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, t, k)) * R_04(x, s, i, j, l)


			+ ((R_01(x, l) + sign * d[2][l]) * R_02(x, s, t) + (R_01(x, s) + sign * d[2][s]) * R_02(x, l, t) + (R_01(x, t) + sign * d[2][t]) * R_02(x, s, l)) * R_04(x, p, i, j, k)
			+ ((R_01(x, l) + sign * d[2][l]) * R_02(x, s, p) + (R_01(x, s) + sign * d[2][s]) * R_02(x, l, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, s, l)) * R_04(x, t, i, j, k)

			+ ((R_01(x, l) + sign * d[2][l]) * R_02(x, t, p) + (R_01(x, t) + sign * d[2][t]) * R_02(x, l, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, t, l)) * R_04(x, s, i, j, k)

			+ ((R_01(x, s) + sign * d[2][s]) * R_02(x, t, p) + (R_01(x, t) + sign * d[2][t]) * R_02(x, s, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, t, s)) * R_04(x, l, i, j, k)


			+ (R_01(x, i) + sign * d[2][i]) * (R_03(x, j, k, l) * R_03(x, s, t, p) + R_03(x, j, k, s) * R_03(x, l, t, p) + R_03(x, j, k, t) * R_03(x, s, l, p) + R_03(x, j, k, p) * R_03(x, s, t, l)
				+ R_03(x, j, l, s) * R_03(x, k, t, p) + R_03(x, j, l, t) * R_03(x, s, k, p) + R_03(x, j, l, p) * R_03(x, s, t, k)
				+ R_03(x, j, s, t) * R_03(x, l, k, p) + R_03(x, j, s, p) * R_03(x, l, t, k) + R_03(x, j, t, p) * R_03(x, l, s, k))
			+ (R_01(x, j) + sign * d[2][j]) * (R_03(x, i, k, l) * R_03(x, s, t, p) + R_03(x, i, k, s) * R_03(x, l, t, p) + R_03(x, i, k, t) * R_03(x, s, l, p) + R_03(x, i, k, p) * R_03(x, s, t, l)
				+ R_03(x, i, l, s) * R_03(x, k, t, p) + R_03(x, i, l, t) * R_03(x, s, k, p) + R_03(x, i, l, p) * R_03(x, s, t, k)
				+ R_03(x, i, s, t) * R_03(x, l, k, p) + R_03(x, i, s, p) * R_03(x, l, t, k) + R_03(x, i, t, p) * R_03(x, l, s, k))
			+ (R_01(x, k) + sign * d[2][k]) * (R_03(x, i, j, l) * R_03(x, s, t, p) + R_03(x, i, j, s) * R_03(x, l, t, p) + R_03(x, i, j, t) * R_03(x, s, l, p) + R_03(x, i, j, p) * R_03(x, s, t, l)
				+ R_03(x, i, l, s) * R_03(x, j, t, p) + R_03(x, i, l, t) * R_03(x, s, j, p) + R_03(x, i, l, p) * R_03(x, s, t, j)
				+ R_03(x, i, s, t) * R_03(x, l, j, p) + R_03(x, i, s, p) * R_03(x, l, t, j) + R_03(x, i, t, p) * R_03(x, l, s, j))
			+ (R_01(x, l) + sign * d[2][l]) * (R_03(x, i, j, k) * R_03(x, s, t, p) + R_03(x, i, j, s) * R_03(x, k, t, p) + R_03(x, i, j, t) * R_03(x, s, k, p) + R_03(x, i, j, p) * R_03(x, s, t, k)
				+ R_03(x, i, k, s) * R_03(x, j, t, p) + R_03(x, i, k, t) * R_03(x, s, j, p) + R_03(x, i, k, p) * R_03(x, s, t, j)
				+ R_03(x, i, s, t) * R_03(x, k, j, p) + R_03(x, i, s, p) * R_03(x, k, t, j) + R_03(x, i, t, p) * R_03(x, k, s, j))
			+ (R_01(x, s) + sign * d[2][s]) * (R_03(x, i, j, k) * R_03(x, l, t, p) + R_03(x, i, j, l) * R_03(x, k, t, p) + R_03(x, i, j, t) * R_03(x, l, k, p) + R_03(x, i, j, p) * R_03(x, l, t, k)
				+ R_03(x, i, k, l) * R_03(x, j, t, p) + R_03(x, i, k, t) * R_03(x, l, j, p) + R_03(x, i, k, p) * R_03(x, l, t, j)
				+ R_03(x, i, l, t) * R_03(x, k, j, p) + R_03(x, i, l, p) * R_03(x, k, t, j) + R_03(x, i, t, p) * R_03(x, k, l, j))
			+ (R_01(x, t) + sign * d[2][t]) * (R_03(x, i, j, k) * R_03(x, l, s, p) + R_03(x, i, j, l) * R_03(x, k, s, p) + R_03(x, i, j, s) * R_03(x, l, k, p) + R_03(x, i, j, p) * R_03(x, l, s, k)
				+ R_03(x, i, k, l) * R_03(x, j, s, p) + R_03(x, i, k, s) * R_03(x, l, j, p) + R_03(x, i, k, p) * R_03(x, l, s, j)
				+ R_03(x, i, l, s) * R_03(x, k, j, p) + R_03(x, i, l, p) * R_03(x, k, s, j) + R_03(x, i, s, p) * R_03(x, k, l, j))
			+ (R_01(x, p) + sign * d[2][p]) * (R_03(x, i, j, k) * R_03(x, l, s, t) + R_03(x, i, j, l) * R_03(x, k, s, t) + R_03(x, i, j, s) * R_03(x, l, k, t) + R_03(x, i, j, t) * R_03(x, l, s, k)
				+ R_03(x, i, k, l) * R_03(x, j, s, t) + R_03(x, i, k, s) * R_03(x, l, j, t) + R_03(x, i, k, t) * R_03(x, l, s, j)
				+ R_03(x, i, l, s) * R_03(x, k, j, t) + R_03(x, i, l, t) * R_03(x, k, s, j) + R_03(x, i, s, t) * R_03(x, k, l, j))


			+ R_03(x, i, j, k) * (R_02(x, l, s) * R_02(x, t, p) + R_02(x, l, t) * R_02(x, s, p) + R_02(x, l, p) * R_02(x, s, t))
			+ R_03(x, i, j, l) * (R_02(x, k, s) * R_02(x, t, p) + R_02(x, k, t) * R_02(x, s, p) + R_02(x, k, p) * R_02(x, s, t))
			+ R_03(x, i, j, s) * (R_02(x, k, l) * R_02(x, t, p) + R_02(x, k, t) * R_02(x, l, p) + R_02(x, k, p) * R_02(x, l, t))
			+ R_03(x, i, j, t) * (R_02(x, k, l) * R_02(x, s, p) + R_02(x, k, s) * R_02(x, l, p) + R_02(x, k, p) * R_02(x, l, s))
			+ R_03(x, i, j, p) * (R_02(x, k, l) * R_02(x, s, t) + R_02(x, k, s) * R_02(x, l, t) + R_02(x, k, t) * R_02(x, l, s))

			+ R_03(x, i, k, l) * (R_02(x, j, s) * R_02(x, t, p) + R_02(x, j, t) * R_02(x, s, p) + R_02(x, j, p) * R_02(x, s, t))
			+ R_03(x, i, k, s) * (R_02(x, j, l) * R_02(x, t, p) + R_02(x, j, t) * R_02(x, l, p) + R_02(x, j, p) * R_02(x, l, t))
			+ R_03(x, i, k, t) * (R_02(x, j, l) * R_02(x, s, p) + R_02(x, j, s) * R_02(x, l, p) + R_02(x, j, p) * R_02(x, l, s))
			+ R_03(x, i, k, p) * (R_02(x, j, l) * R_02(x, s, t) + R_02(x, j, s) * R_02(x, l, t) + R_02(x, j, t) * R_02(x, l, s))

			+ R_03(x, i, l, s) * (R_02(x, j, k) * R_02(x, t, p) + R_02(x, j, t) * R_02(x, k, p) + R_02(x, j, p) * R_02(x, k, t))
			+ R_03(x, i, l, t) * (R_02(x, j, k) * R_02(x, s, p) + R_02(x, j, s) * R_02(x, k, p) + R_02(x, j, p) * R_02(x, k, s))
			+ R_03(x, i, l, p) * (R_02(x, j, k) * R_02(x, s, t) + R_02(x, j, s) * R_02(x, k, t) + R_02(x, j, t) * R_02(x, k, s))

			+ R_03(x, i, s, t) * (R_02(x, j, k) * R_02(x, l, p) + R_02(x, j, l) * R_02(x, k, p) + R_02(x, j, p) * R_02(x, k, l))
			+ R_03(x, i, s, p) * (R_02(x, j, k) * R_02(x, l, t) + R_02(x, j, l) * R_02(x, k, t) + R_02(x, j, t) * R_02(x, k, l))

			+ R_03(x, i, t, p) * (R_02(x, j, k) * R_02(x, l, s) + R_02(x, j, l) * R_02(x, k, s) + R_02(x, j, s) * R_02(x, k, l))

			+ R_03(x, j, k, l) * (R_02(x, i, s) * R_02(x, t, p) + R_02(x, i, t) * R_02(x, s, p) + R_02(x, i, p) * R_02(x, s, t))
			+ R_03(x, j, k, s) * (R_02(x, i, l) * R_02(x, t, p) + R_02(x, i, t) * R_02(x, l, p) + R_02(x, i, p) * R_02(x, l, t))
			+ R_03(x, j, k, t) * (R_02(x, i, l) * R_02(x, s, p) + R_02(x, i, s) * R_02(x, l, p) + R_02(x, i, p) * R_02(x, l, s))
			+ R_03(x, j, k, p) * (R_02(x, i, l) * R_02(x, s, t) + R_02(x, i, s) * R_02(x, l, t) + R_02(x, i, t) * R_02(x, l, s))

			+ R_03(x, j, l, s) * (R_02(x, i, k) * R_02(x, t, p) + R_02(x, i, t) * R_02(x, k, p) + R_02(x, i, p) * R_02(x, k, t))
			+ R_03(x, j, l, t) * (R_02(x, i, k) * R_02(x, s, p) + R_02(x, i, s) * R_02(x, k, p) + R_02(x, i, p) * R_02(x, k, s))
			+ R_03(x, j, l, p) * (R_02(x, i, k) * R_02(x, s, t) + R_02(x, i, s) * R_02(x, k, t) + R_02(x, i, t) * R_02(x, k, s))

			+ R_03(x, j, s, t) * (R_02(x, i, k) * R_02(x, l, p) + R_02(x, i, l) * R_02(x, k, p) + R_02(x, i, p) * R_02(x, k, l))
			+ R_03(x, j, s, p) * (R_02(x, i, k) * R_02(x, l, t) + R_02(x, i, l) * R_02(x, k, t) + R_02(x, i, t) * R_02(x, k, l))

			+ R_03(x, j, t, p) * (R_02(x, i, k) * R_02(x, l, s) + R_02(x, i, l) * R_02(x, k, s) + R_02(x, i, s) * R_02(x, k, l))

			+ R_03(x, k, l, s) * (R_02(x, i, j) * R_02(x, t, p) + R_02(x, i, t) * R_02(x, j, p) + R_02(x, i, p) * R_02(x, j, t))
			+ R_03(x, k, l, t) * (R_02(x, i, j) * R_02(x, s, p) + R_02(x, i, s) * R_02(x, j, p) + R_02(x, i, p) * R_02(x, j, s))
			+ R_03(x, k, l, p) * (R_02(x, i, j) * R_02(x, s, t) + R_02(x, i, s) * R_02(x, j, t) + R_02(x, i, t) * R_02(x, j, s))

			+ R_03(x, k, s, t) * (R_02(x, i, j) * R_02(x, l, p) + R_02(x, i, l) * R_02(x, j, p) + R_02(x, i, p) * R_02(x, j, l))
			+ R_03(x, k, s, p) * (R_02(x, i, j) * R_02(x, l, t) + R_02(x, i, l) * R_02(x, j, t) + R_02(x, i, t) * R_02(x, j, l))

			+ R_03(x, k, t, p) * (R_02(x, i, j) * R_02(x, l, s) + R_02(x, i, l) * R_02(x, j, s) + R_02(x, i, s) * R_02(x, j, l))

			+ R_03(x, l, s, t) * (R_02(x, i, j) * R_02(x, k, p) + R_02(x, i, k) * R_02(x, j, p) + R_02(x, i, p) * R_02(x, j, k))
			+ R_03(x, l, s, p) * (R_02(x, i, j) * R_02(x, k, t) + R_02(x, i, k) * R_02(x, j, t) + R_02(x, i, t) * R_02(x, j, k))

			+ R_03(x, l, t, p) * (R_02(x, i, j) * R_02(x, k, s) + R_02(x, i, k) * R_02(x, j, s) + R_02(x, i, s) * R_02(x, j, k))

			+ R_03(x, s, t, p) * (R_02(x, i, j) * R_02(x, k, l) + R_02(x, i, k) * R_02(x, j, l) + R_02(x, i, l) * R_02(x, j, k))

			)

		- 24 * R1 * R1 * R1 * ((R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * R_04(x, l, s, t, p) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * R_04(x, k, s, t, p)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * R_04(x, k, l, t, p) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, k, l, s, p)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, k, l, s, t)

			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * R_04(x, j, s, t, p) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * R_04(x, j, l, t, p)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, j, l, s, p) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, j, l, s, t)

			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * R_04(x, j, k, t, p) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, j, k, s, p)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, j, k, s, t)

			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, j, k, l, p) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, j, k, l, t)

			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, j, k, l, s)


			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * R_04(x, i, s, t, p) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * R_04(x, i, l, t, p)
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, i, l, s, p) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, i, l, s, t)

			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * R_04(x, i, k, t, p) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, i, k, s, p)
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, i, k, s, t)

			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, i, k, l, p) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, i, k, l, t)

			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, i, k, l, s)


			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * R_04(x, i, j, t, p) + (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, i, j, s, p)
			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, i, j, s, t)

			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, i, j, l, p) + (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, i, j, l, t)

			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, i, j, l, s)


			+ (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, i, j, k, p) + (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, i, j, k, t)

			+ (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, i, j, k, s)


			+ (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, i, j, k, l)

			+ R_03(x, i, j, k) * (R_02(x, l, s) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, t) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, t, s) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, s, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, l) + sign * d[2][l]) + R_02(x, t, p) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]))
			+ R_03(x, i, j, l) * (R_02(x, k, s) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, t) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, t, s) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, s, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, t, p) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]))
			+ R_03(x, i, j, s) * (R_02(x, l, k) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, t) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, k) + sign * d[2][k])
				+ R_02(x, t, k) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, l) + sign * d[2][l]) + R_02(x, t, p) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, k) + sign * d[2][k]))
			+ R_03(x, i, j, t) * (R_02(x, l, s) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, k) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, p) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, k, s) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, s, p) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) + R_02(x, k, p) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]))
			+ R_03(x, i, j, p) * (R_02(x, l, s) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, l, t) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, l, k) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, t, s) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, s, k) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, l) + sign * d[2][l]) + R_02(x, t, k) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]))


			+ R_03(x, i, k, l) * (R_02(x, j, s) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, t) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, t, s) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, s, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, t, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]))
			+ R_03(x, i, k, s) * (R_02(x, l, j) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, t) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, j) + sign * d[2][j])
				+ R_02(x, t, j) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, l) + sign * d[2][l]) + R_02(x, t, p) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, j) + sign * d[2][j]))
			+ R_03(x, i, k, t) * (R_02(x, l, s) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, j) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, j, s) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, s, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) + R_02(x, j, p) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]))
			+ R_03(x, i, k, p) * (R_02(x, l, s) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, l, t) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, l, j) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, t, s) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, s, j) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, l) + sign * d[2][l]) + R_02(x, t, j) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]))


			+ R_03(x, i, l, s) * (R_02(x, k, j) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, t) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, j) + sign * d[2][j])
				+ R_02(x, t, j) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, t, p) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, j) + sign * d[2][j]))
			+ R_03(x, i, l, t) * (R_02(x, k, s) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, j) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, j, s) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, s, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, j, p) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]))
			+ R_03(x, i, l, p) * (R_02(x, k, s) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, k, t) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, k, j) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, t, s) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, s, j) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, t, j) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]))


			+ R_03(x, i, s, t) * (R_02(x, k, l) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, j) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l])
				+ R_02(x, j, l) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, j, p) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]))
			+ R_03(x, i, s, p) * (R_02(x, k, l) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, k, t) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, k, j) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, l) + sign * d[2][l])
				+ R_02(x, t, l) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, l, j) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, t, j) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]))


			+ R_03(x, i, t, p) * (R_02(x, k, l) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, k, s) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, k, j) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, l) + sign * d[2][l])
				+ R_02(x, s, l) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, l, j) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, s, j) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]))


			+ R_03(x, j, k, l) * (R_02(x, i, s) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, i, t) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, i, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, t, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, s, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, t, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s]))
			+ R_03(x, j, k, s) * (R_02(x, l, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, i) + sign * d[2][i])
				+ R_02(x, t, i) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, i, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, l) + sign * d[2][l]) + R_02(x, t, p) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, i) + sign * d[2][i]))
			+ R_03(x, j, k, t) * (R_02(x, l, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, i) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, i, s) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, s, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) + R_02(x, i, p) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]))
			+ R_03(x, j, k, p) * (R_02(x, l, s) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, l, t) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, l, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, t, s) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, s, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, l) + sign * d[2][l]) + R_02(x, t, i) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]))


			+ R_03(x, j, l, s) * (R_02(x, k, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, i) + sign * d[2][i])
				+ R_02(x, t, i) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, i, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, t, p) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, i) + sign * d[2][i]))
			+ R_03(x, j, l, t) * (R_02(x, k, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, i) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, i, s) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, s, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, i, p) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]))
			+ R_03(x, j, l, p) * (R_02(x, k, s) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, k, t) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, k, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, t, s) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, s, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, t, i) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]))


			+ R_03(x, j, s, t) * (R_02(x, k, l) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, i) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l])
				+ R_02(x, i, l) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, i, p) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]))
			+ R_03(x, j, s, p) * (R_02(x, k, l) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, k, t) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, k, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, l) + sign * d[2][l])
				+ R_02(x, t, l) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, l, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, t, i) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]))


			+ R_03(x, j, t, p) * (R_02(x, k, l) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, k, s) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, k, i) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, l) + sign * d[2][l])
				+ R_02(x, s, l) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, l, i) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, s, i) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]))


			+ R_03(x, k, l, s) * (R_02(x, j, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, i) + sign * d[2][i])
				+ R_02(x, t, i) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, i, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, t, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, i) + sign * d[2][i]))
			+ R_03(x, k, l, t) * (R_02(x, j, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, i) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, i, s) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, s, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, i, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]))
			+ R_03(x, k, l, p) * (R_02(x, j, s) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, t) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, t, s) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, s, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, t, i) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]))


			+ R_03(x, k, s, t) * (R_02(x, j, l) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, i) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l])
				+ R_02(x, i, l) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, i, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]))
			+ R_03(x, k, s, p) * (R_02(x, j, l) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, t) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, l) + sign * d[2][l])
				+ R_02(x, t, l) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, l, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, t, i) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]))


			+ R_03(x, k, t, p) * (R_02(x, j, l) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, s) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, i) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, l) + sign * d[2][l])
				+ R_02(x, s, l) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, l, i) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, s, i) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]))


			+ R_03(x, l, s, t) * (R_02(x, j, k) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, i) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k])
				+ R_02(x, i, k) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, i, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]))
			+ R_03(x, l, s, p) * (R_02(x, j, k) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, t) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, k) + sign * d[2][k])
				+ R_02(x, t, k) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, k, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, t, i) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]))


			+ R_03(x, l, t, p) * (R_02(x, j, k) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, s) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, i) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, k) + sign * d[2][k])
				+ R_02(x, s, k) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, k, i) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, s, i) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]))


			+ R_03(x, s, t, p) * (R_02(x, j, k) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, l) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, i) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, k) + sign * d[2][k])
				+ R_02(x, l, k) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, k, i) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, l, i) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]))

			+ (R_01(x, i) + sign * d[2][i]) * (R_02(x, j, k) * (R_02(x, s, t) * R_02(x, l, p) + R_02(x, s, l) * R_02(x, p, t) + R_02(x, s, p) * R_02(x, l, t))
				+ R_02(x, j, l) * (R_02(x, s, t) * R_02(x, k, p) + R_02(x, s, k) * R_02(x, p, t) + R_02(x, s, p) * R_02(x, k, t))
				+ R_02(x, j, s) * (R_02(x, k, l) * R_02(x, t, p) + R_02(x, k, t) * R_02(x, l, p) + R_02(x, k, p) * R_02(x, l, t))
				+ R_02(x, j, t) * (R_02(x, k, l) * R_02(x, s, p) + R_02(x, k, s) * R_02(x, l, p) + R_02(x, k, p) * R_02(x, l, s))
				+ R_02(x, j, p) * (R_02(x, k, l) * R_02(x, s, t) + R_02(x, k, s) * R_02(x, l, t) + R_02(x, k, t) * R_02(x, l, s)))
			+ (R_01(x, j) + sign * d[2][j]) * (R_02(x, i, k) * (R_02(x, s, t) * R_02(x, l, p) + R_02(x, s, l) * R_02(x, p, t) + R_02(x, s, p) * R_02(x, l, t))
				+ R_02(x, i, l) * (R_02(x, s, t) * R_02(x, k, p) + R_02(x, s, k) * R_02(x, p, t) + R_02(x, s, p) * R_02(x, k, t))
				+ R_02(x, i, s) * (R_02(x, k, l) * R_02(x, t, p) + R_02(x, k, t) * R_02(x, l, p) + R_02(x, k, p) * R_02(x, l, t))
				+ R_02(x, i, t) * (R_02(x, k, l) * R_02(x, s, p) + R_02(x, k, s) * R_02(x, l, p) + R_02(x, k, p) * R_02(x, l, s))
				+ R_02(x, i, p) * (R_02(x, k, l) * R_02(x, s, t) + R_02(x, k, s) * R_02(x, l, t) + R_02(x, k, t) * R_02(x, l, s)))
			+ (R_01(x, k) + sign * d[2][k]) * (R_02(x, i, j) * (R_02(x, s, t) * R_02(x, l, p) + R_02(x, s, l) * R_02(x, p, t) + R_02(x, s, p) * R_02(x, l, t))
				+ R_02(x, i, l) * (R_02(x, s, t) * R_02(x, j, p) + R_02(x, s, j) * R_02(x, p, t) + R_02(x, s, p) * R_02(x, j, t))
				+ R_02(x, i, s) * (R_02(x, j, l) * R_02(x, t, p) + R_02(x, j, t) * R_02(x, l, p) + R_02(x, j, p) * R_02(x, l, t))
				+ R_02(x, i, t) * (R_02(x, j, l) * R_02(x, s, p) + R_02(x, j, s) * R_02(x, l, p) + R_02(x, j, p) * R_02(x, l, s))
				+ R_02(x, i, p) * (R_02(x, j, l) * R_02(x, s, t) + R_02(x, j, s) * R_02(x, l, t) + R_02(x, j, t) * R_02(x, l, s)))
			+ (R_01(x, l) + sign * d[2][l]) * (R_02(x, i, j) * (R_02(x, s, t) * R_02(x, k, p) + R_02(x, s, k) * R_02(x, p, t) + R_02(x, s, p) * R_02(x, k, t))
				+ R_02(x, i, k) * (R_02(x, s, t) * R_02(x, j, p) + R_02(x, s, j) * R_02(x, p, t) + R_02(x, s, p) * R_02(x, j, t))
				+ R_02(x, i, s) * (R_02(x, j, k) * R_02(x, t, p) + R_02(x, j, t) * R_02(x, k, p) + R_02(x, j, p) * R_02(x, k, t))
				+ R_02(x, i, t) * (R_02(x, j, k) * R_02(x, s, p) + R_02(x, j, s) * R_02(x, k, p) + R_02(x, j, p) * R_02(x, k, s))
				+ R_02(x, i, p) * (R_02(x, j, k) * R_02(x, s, t) + R_02(x, j, s) * R_02(x, k, t) + R_02(x, j, t) * R_02(x, k, s)))
			+ (R_01(x, s) + sign * d[2][s]) * (R_02(x, i, j) * (R_02(x, l, t) * R_02(x, k, p) + R_02(x, l, k) * R_02(x, p, t) + R_02(x, l, p) * R_02(x, k, t))
				+ R_02(x, i, k) * (R_02(x, l, t) * R_02(x, j, p) + R_02(x, l, j) * R_02(x, p, t) + R_02(x, l, p) * R_02(x, j, t))
				+ R_02(x, i, l) * (R_02(x, j, k) * R_02(x, t, p) + R_02(x, j, t) * R_02(x, k, p) + R_02(x, j, p) * R_02(x, k, t))
				+ R_02(x, i, t) * (R_02(x, j, k) * R_02(x, l, p) + R_02(x, j, l) * R_02(x, k, p) + R_02(x, j, p) * R_02(x, k, l))
				+ R_02(x, i, p) * (R_02(x, j, k) * R_02(x, l, t) + R_02(x, j, l) * R_02(x, k, t) + R_02(x, j, t) * R_02(x, k, l)))
			+ (R_01(x, t) + sign * d[2][t]) * (R_02(x, i, j) * (R_02(x, l, s) * R_02(x, k, p) + R_02(x, l, k) * R_02(x, p, s) + R_02(x, l, p) * R_02(x, k, s))
				+ R_02(x, i, k) * (R_02(x, l, s) * R_02(x, j, p) + R_02(x, l, j) * R_02(x, p, s) + R_02(x, l, p) * R_02(x, j, s))
				+ R_02(x, i, l) * (R_02(x, j, k) * R_02(x, s, p) + R_02(x, j, s) * R_02(x, k, p) + R_02(x, j, p) * R_02(x, k, s))
				+ R_02(x, i, s) * (R_02(x, j, k) * R_02(x, l, p) + R_02(x, j, l) * R_02(x, k, p) + R_02(x, j, p) * R_02(x, k, l))
				+ R_02(x, i, p) * (R_02(x, j, k) * R_02(x, l, s) + R_02(x, j, l) * R_02(x, k, s) + R_02(x, j, s) * R_02(x, k, l)))
			+ (R_01(x, p) + sign * d[2][p]) * (R_02(x, i, j) * (R_02(x, l, s) * R_02(x, k, t) + R_02(x, l, k) * R_02(x, t, s) + R_02(x, l, t) * R_02(x, k, s))
				+ R_02(x, i, k) * (R_02(x, l, s) * R_02(x, j, t) + R_02(x, l, j) * R_02(x, t, s) + R_02(x, l, t) * R_02(x, j, s))
				+ R_02(x, i, l) * (R_02(x, j, k) * R_02(x, s, t) + R_02(x, j, s) * R_02(x, k, t) + R_02(x, j, t) * R_02(x, k, s))
				+ R_02(x, i, s) * (R_02(x, j, k) * R_02(x, l, t) + R_02(x, j, l) * R_02(x, k, t) + R_02(x, j, t) * R_02(x, k, l))
				+ R_02(x, i, t) * (R_02(x, j, k) * R_02(x, l, s) + R_02(x, j, l) * R_02(x, k, s) + R_02(x, j, s) * R_02(x, k, l)))


			)

		+ 120 * R1 * R1 * (R_03(x, i, j, k) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) + R_03(x, i, j, l) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])
			+ R_03(x, i, j, s) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) + R_03(x, i, j, t) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s])
			+ R_03(x, i, j, p) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])

			+ R_03(x, i, k, l) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) + R_03(x, i, k, s) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t])
			+ R_03(x, i, k, t) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) + R_03(x, i, k, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])

			+ R_03(x, i, l, s) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) + R_03(x, i, l, t) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s])
			+ R_03(x, i, l, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])

			+ R_03(x, i, s, t) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) + R_03(x, i, s, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t])

			+ R_03(x, i, t, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s])


			+ R_03(x, j, k, l) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) + R_03(x, j, k, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t])
			+ R_03(x, j, k, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) + R_03(x, j, k, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])

			+ R_03(x, j, l, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) + R_03(x, j, l, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s])
			+ R_03(x, j, l, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])

			+ R_03(x, j, s, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) + R_03(x, j, s, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t])

			+ R_03(x, j, t, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s])


			+ R_03(x, k, l, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, t) + sign * d[2][t]) + R_03(x, k, l, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, s) + sign * d[2][s])
			+ R_03(x, k, l, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])

			+ R_03(x, k, s, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, l) + sign * d[2][l]) + R_03(x, k, s, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t])

			+ R_03(x, k, t, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s])


			+ R_03(x, l, s, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) + R_03(x, l, s, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t])

			+ R_03(x, l, t, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s])


			+ R_03(x, t, s, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l])




			+(R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_02(x, l, s) * R_02(x, t, p) + R_02(x, l, t) * R_02(x, s, p) + R_02(x, l, p) * R_02(x, s, t))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_02(x, k, s) * R_02(x, t, p) + R_02(x, k, t) * R_02(x, s, p) + R_02(x, k, p) * R_02(x, s, t))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_02(x, l, k) * R_02(x, t, p) + R_02(x, l, t) * R_02(x, k, p) + R_02(x, l, p) * R_02(x, k, t))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, l, s) * R_02(x, k, p) + R_02(x, l, k) * R_02(x, s, p) + R_02(x, l, p) * R_02(x, s, k))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, l, s) * R_02(x, t, k) + R_02(x, l, t) * R_02(x, s, k) + R_02(x, l, k) * R_02(x, s, t))

			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_02(x, j, s) * R_02(x, t, p) + R_02(x, j, t) * R_02(x, s, p) + R_02(x, j, p) * R_02(x, s, t))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_02(x, l, j) * R_02(x, t, p) + R_02(x, l, t) * R_02(x, j, p) + R_02(x, l, p) * R_02(x, j, t))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, l, s) * R_02(x, j, p) + R_02(x, l, j) * R_02(x, s, p) + R_02(x, l, p) * R_02(x, s, j))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, l, s) * R_02(x, t, j) + R_02(x, l, t) * R_02(x, s, j) + R_02(x, l, j) * R_02(x, s, t))


			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_02(x, k, j) * R_02(x, t, p) + R_02(x, k, t) * R_02(x, j, p) + R_02(x, k, p) * R_02(x, j, t))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, k, s) * R_02(x, j, p) + R_02(x, k, j) * R_02(x, s, p) + R_02(x, k, p) * R_02(x, s, j))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, k, s) * R_02(x, t, j) + R_02(x, k, t) * R_02(x, s, j) + R_02(x, k, j) * R_02(x, s, t))



			+(R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, k, l) * R_02(x, j, p) + R_02(x, k, j) * R_02(x, l, p) + R_02(x, k, p) * R_02(x, l, j))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, k, l) * R_02(x, t, j) + R_02(x, k, t) * R_02(x, l, j) + R_02(x, k, j) * R_02(x, l, t))

			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, k, l) * R_02(x, s, j) + R_02(x, k, s) * R_02(x, l, j) + R_02(x, k, j) * R_02(x, l, s))


			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_02(x, i, s) * R_02(x, t, p) + R_02(x, i, t) * R_02(x, s, p) + R_02(x, i, p) * R_02(x, s, t))
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_02(x, l, i) * R_02(x, t, p) + R_02(x, l, t) * R_02(x, i, p) + R_02(x, l, p) * R_02(x, i, t))
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, l, s) * R_02(x, i, p) + R_02(x, l, i) * R_02(x, s, p) + R_02(x, l, p) * R_02(x, s, i))
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, l, s) * R_02(x, t, i) + R_02(x, l, t) * R_02(x, s, i) + R_02(x, l, i) * R_02(x, s, t))


			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_02(x, k, i) * R_02(x, t, p) + R_02(x, k, t) * R_02(x, i, p) + R_02(x, k, p) * R_02(x, i, t))
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, k, s) * R_02(x, i, p) + R_02(x, k, i) * R_02(x, s, p) + R_02(x, k, p) * R_02(x, s, i))
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, k, s) * R_02(x, t, i) + R_02(x, k, t) * R_02(x, s, i) + R_02(x, k, i) * R_02(x, s, t))



			+(R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, k, l) * R_02(x, i, p) + R_02(x, k, i) * R_02(x, l, p) + R_02(x, k, p) * R_02(x, l, i))
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, k, l) * R_02(x, t, i) + R_02(x, k, t) * R_02(x, l, i) + R_02(x, k, i) * R_02(x, l, t))

			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, k, l) * R_02(x, s, i) + R_02(x, k, s) * R_02(x, l, i) + R_02(x, k, i) * R_02(x, l, s))



			+(R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_02(x, j, i) * R_02(x, t, p) + R_02(x, j, t) * R_02(x, i, p) + R_02(x, j, p) * R_02(x, i, t))
			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, j, s) * R_02(x, i, p) + R_02(x, j, i) * R_02(x, s, p) + R_02(x, j, p) * R_02(x, s, i))
			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, j, s) * R_02(x, t, i) + R_02(x, j, t) * R_02(x, s, i) + R_02(x, j, i) * R_02(x, s, t))


			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, j, l) * R_02(x, i, p) + R_02(x, j, i) * R_02(x, l, p) + R_02(x, j, p) * R_02(x, l, i))
			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, j, l) * R_02(x, t, i) + R_02(x, j, t) * R_02(x, l, i) + R_02(x, j, i) * R_02(x, l, t))

			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, j, l) * R_02(x, s, i) + R_02(x, j, s) * R_02(x, l, i) + R_02(x, j, i) * R_02(x, l, s))


			+ (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, j, k) * R_02(x, i, p) + R_02(x, j, i) * R_02(x, k, p) + R_02(x, j, p) * R_02(x, k, i))
			+ (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, j, k) * R_02(x, t, i) + R_02(x, j, t) * R_02(x, k, i) + R_02(x, j, i) * R_02(x, k, t))

			+ (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, j, k) * R_02(x, s, i) + R_02(x, j, s) * R_02(x, k, i) + R_02(x, j, i) * R_02(x, k, s))

			+ (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, j, k) * R_02(x, l, i) + R_02(x, j, l) * R_02(x, k, i) + R_02(x, j, i) * R_02(x, k, l))

			)

		- 720 * R1 * (R_02(x, i, j) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])
			+ R_02(x, i, k) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])
			+ R_02(x, i, l) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])
			+ R_02(x, i, s) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, t) + sign * d[2][t])
			+ R_02(x, i, t) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p])
			+ R_02(x, i, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])

			+ R_02(x, j, k) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])
			+ R_02(x, j, l) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])
			+ R_02(x, j, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, t) + sign * d[2][t])
			+ R_02(x, j, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p])
			+ R_02(x, j, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])

			+ R_02(x, k, l) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])
			+ R_02(x, k, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, t) + sign * d[2][t])
			+ R_02(x, k, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p])
			+ R_02(x, k, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])


			+ R_02(x, l, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, t) + sign * d[2][t])
			+ R_02(x, l, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p])
			+ R_02(x, l, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])


			+ R_02(x, s, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p])
			+ R_02(x, s, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t])

			+ R_02(x, t, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]))

		+ 5040 * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]));
	return fR_ijklstp;
}

double lnR_01(double* x, int i, int sign)
{
	double R1 = R(x) + sign * x[2];
	double lnR_i = 1.0 / (R1) * (R_01(x, i) + d[2][i]);
	return lnR_i;
}
double lnR_02(double* x, int i, int j, int sign)
{
	double R1 = R(x) + sign * x[2];

	double lnR_ij = 1.0 / (R1 * R1) * (R1 * R_02(x, i, j) - (R_01(x, i) + sign * d[i][2]) * (R_01(x, j) + sign * d[j][2]));

	return lnR_ij;
}

double lnR_03(double* x, int i, int j, int k, int sign)
{
	double R1 = R(x) + sign * x[2];

	double lnR_ijk = 1 / (R1 * R1 * R1) * (R1 * R1 * R_03(x, i, j, k)

		- R1 * ((R_01(x, i) + sign * d[i][2]) * R_02(x, j, k) + (R_01(x, j) + sign * d[j][2]) * R_02(x, i, k) + (R_01(x, k) + sign * d[k][2]) * R_02(x, i, j))
		+ 2 * (R_01(x, i) + sign * d[i][2]) * (R_01(x, j) + sign * d[j][2]) * (R_01(x, k) + sign * d[k][2]));

	return lnR_ijk;
}

double lnR_04(double* x, int i, int j, int k, int l, int sign)
{
	double R1 = R(x) + sign * x[2];
	double lnR_ijkl = 1 / (R1 * R1 * R1 * R1) * (R1 * R1 * R1 * R_04(x, i, j, k, l)
		- R1 * R1 * ((R_01(x, i) + sign * d[2][i]) * R_03(x, j, k, l) + (R_01(x, j) + sign * d[2][j]) * R_03(x, i, k, l) + (R_01(x, k) + sign * d[2][k]) * R_03(x, i, j, l) + (R_01(x, l) + sign * d[2][l]) * R_03(x, i, j, k)

			+ R_02(x, i, j) * R_02(x, k, l) + R_02(x, i, k) * R_02(x, j, l) + R_02(x, i, l) * R_02(x, j, k))

		+ 2 * R1 * ((R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * R_02(x, k, l) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * R_02(x, j, l) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * R_02(x, j, k)
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * R_02(x, i, l) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * R_02(x, i, k) + (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * R_02(x, i, j))

		- 6 * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]));
	return lnR_ijkl;
}

double lnR_05(double* x, int i, int j, int k, int l, int s, int sign)
{
	double R1 = R(x) + sign * x[2];
	double lnR_ijkls = 1 / (R1 * R1 * R1 * R1 * R1) * (R1 * R1 * R1 * R1 * R_05(x, i, j, k, l, s)

		- R1 * R1 * R1 * ((R_01(x, i) + sign * d[2][i]) * R_04(x, j, k, l, s) + (R_01(x, j) + sign * d[2][j]) * R_04(x, i, k, l, s)
			+ (R_01(x, k) + sign * d[2][k]) * R_04(x, i, j, l, s) + (R_01(x, l) + sign * d[2][l]) * R_04(x, i, j, k, s)
			+ (R_01(x, s) + sign * d[2][s]) * R_04(x, i, j, k, l)

			+ R_02(x, i, j) * R_03(x, k, l, s) + R_02(x, i, k) * R_03(x, j, l, s) + R_02(x, i, l) * R_03(x, k, j, s)
			+ R_02(x, i, s) * R_03(x, k, l, j)

			+ R_02(x, j, k) * R_03(x, i, l, s) + R_02(x, j, l) * R_03(x, k, i, s) + R_02(x, j, s) * R_03(x, k, i, l)

			+ R_02(x, k, l) * R_03(x, i, j, s) + R_02(x, k, s) * R_03(x, i, j, l)

			+ R_02(x, l, s) * R_03(x, k, i, j)

			)
		+ 2 * R1 * R1 * ((R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * R_03(x, k, l, s) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * R_03(x, j, l, s)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * R_03(x, k, j, s) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s]) * R_03(x, k, j, l)
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * R_03(x, i, l, s) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * R_03(x, k, i, s)
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * R_03(x, i, l, k) + (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * R_03(x, i, j, s)
			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * R_03(x, i, l, j) + (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * R_03(x, i, j, k)

			+ (R_01(x, i) + sign * d[2][i]) * (R_02(x, k, j) * R_02(x, l, s) + R_02(x, j, l) * R_02(x, k, s) + R_02(x, s, j) * R_02(x, l, k))
			+ (R_01(x, j) + sign * d[2][j]) * (R_02(x, k, i) * R_02(x, l, s) + R_02(x, i, l) * R_02(x, k, s) + R_02(x, s, i) * R_02(x, l, k))
			+ (R_01(x, k) + sign * d[2][k]) * (R_02(x, j, i) * R_02(x, l, s) + R_02(x, i, l) * R_02(x, j, s) + R_02(x, s, i) * R_02(x, l, j))
			+ (R_01(x, l) + sign * d[2][l]) * (R_02(x, j, i) * R_02(x, k, s) + R_02(x, i, k) * R_02(x, j, s) + R_02(x, s, i) * R_02(x, k, j))
			+ (R_01(x, s) + sign * d[2][s]) * (R_02(x, j, i) * R_02(x, k, l) + R_02(x, i, k) * R_02(x, j, l) + R_02(x, l, i) * R_02(x, k, j))

			)

		- 6 * R1 * (R_02(x, i, j) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l])
			+ R_02(x, i, k) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l])
			+ R_02(x, i, l) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, j) + sign * d[2][j])
			+ R_02(x, i, s) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l])
			+ R_02(x, k, j) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l])
			+ R_02(x, j, l) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k])
			+ R_02(x, s, j) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l])
			+ R_02(x, l, k) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, i) + sign * d[2][i])
			+ R_02(x, k, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l])
			+ R_02(x, l, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]))
		+ 24 * (R_01(x, i) + sign * d[i][2]) * (R_01(x, j) + sign * d[j][2]) * (R_01(x, k) + sign * d[k][2]) * (R_01(x, l) + sign * d[l][2]) * (R_01(x, s) + sign * d[s][2]));
	return lnR_ijkls;
}

double lnR_06(double* x, int i, int j, int k, int l, int s, int t, int sign)
{
	double R1 = R(x) + sign * x[2];
	double lnR_ijklst = 1 / (R1 * R1 * R1 * R1 * R1 * R1) * (R1 * R1 * R1 * R1 * R1 * R_06(x, i, j, k, l, s, t)
		- R1 * R1 * R1 * R1 * ((R_01(x, t) + sign * d[2][t]) * R_05(x, i, j, k, l, s) + (R_01(x, i) + sign * d[2][i]) * R_05(x, j, k, l, s, t) + (R_01(x, j) + sign * d[2][j]) * R_05(x, i, k, l, s, t)
			+ (R_01(x, k) + sign * d[2][k]) * R_05(x, i, j, l, s, t) + (R_01(x, l) + sign * d[2][l]) * R_05(x, i, j, k, s, t) + (R_01(x, s) + sign * d[2][s]) * R_05(x, i, j, k, l, t)


			+ R_02(x, i, j) * R_04(x, k, l, s, t) + R_02(x, i, k) * R_04(x, j, l, s, t) + R_02(x, i, l) * R_04(x, k, j, s, t) + R_02(x, i, s) * R_04(x, k, l, j, t) + R_02(x, i, t) * R_04(x, k, l, s, j)
			+ R_02(x, j, k) * R_04(x, i, l, s, t) + R_02(x, j, l) * R_04(x, i, k, s, t) + R_02(x, j, s) * R_04(x, i, l, k, t) + +R_02(x, j, t) * R_04(x, i, l, s, k)
			+ R_02(x, k, l) * R_04(x, i, j, s, t) + R_02(x, k, s) * R_04(x, i, j, l, t) + R_02(x, k, t) * R_04(x, i, j, s, l)
			+ R_02(x, l, s) * R_04(x, i, j, k, t) + R_02(x, l, t) * R_04(x, i, j, k, s) + R_02(x, s, t) * R_04(x, i, j, k, l)


			+ R_03(x, i, j, k) * R_03(x, l, s, t) + R_03(x, i, j, l) * R_03(x, k, s, t) + R_03(x, i, j, s) * R_03(x, l, k, t) + R_03(x, i, j, t) * R_03(x, l, s, k)
			+ R_03(x, i, k, l) * R_03(x, j, s, t) + R_03(x, i, k, s) * R_03(x, l, j, t) + R_03(x, i, k, t) * R_03(x, j, s, l)
			+ R_03(x, i, l, s) * R_03(x, j, k, t) + R_03(x, i, l, t) * R_03(x, j, k, s) + R_03(x, i, s, t) * R_03(x, j, k, l)
			)
		+ 2 * R1 * R1 * R1 * ((R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * R_04(x, k, l, s, t) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * R_04(x, j, l, s, t) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * R_04(x, j, k, s, t)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s]) * R_04(x, j, k, l, t) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, j, k, l, s) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * R_04(x, i, l, s, t)
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * R_04(x, i, k, s, t) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * R_04(x, i, k, l, t) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, i, k, l, s)
			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * R_04(x, i, j, s, t) + (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * R_04(x, i, j, l, t) + (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, i, j, l, s)
			+ (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * R_04(x, i, j, k, t) + (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, i, j, k, s) + (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, i, j, k, l)

			+ (R_01(x, i) + sign * d[2][i]) * (R_02(x, j, k) * R_03(x, l, s, t) + R_02(x, j, l) * R_03(x, k, s, t) + R_02(x, j, s) * R_03(x, l, k, t) + R_02(x, j, t) * R_03(x, l, s, k)
				+ R_02(x, k, l) * R_03(x, j, s, t) + R_02(x, k, s) * R_03(x, l, j, t) + R_02(x, k, t) * R_03(x, l, j, s)
				+ R_02(x, s, l) * R_03(x, j, k, t) + R_02(x, l, t) * R_03(x, k, j, s) + R_02(x, s, t) * R_03(x, l, j, k)
				)
			+ (R_01(x, j) + sign * d[2][j]) * (R_02(x, i, k) * R_03(x, l, s, t) + R_02(x, i, l) * R_03(x, k, s, t) + R_02(x, i, s) * R_03(x, l, k, t) + R_02(x, i, t) * R_03(x, l, s, k)
				+ R_02(x, k, l) * R_03(x, i, s, t) + R_02(x, k, s) * R_03(x, l, i, t) + R_02(x, k, t) * R_03(x, l, i, s)
				+ R_02(x, s, l) * R_03(x, i, k, t) + R_02(x, l, t) * R_03(x, k, i, s) + R_02(x, s, t) * R_03(x, l, i, k)
				)
			+ (R_01(x, k) + sign * d[2][k]) * (R_02(x, i, j) * R_03(x, l, s, t) + R_02(x, i, l) * R_03(x, j, s, t) + R_02(x, i, s) * R_03(x, l, j, t) + R_02(x, i, t) * R_03(x, l, s, j)
				+ R_02(x, j, l) * R_03(x, i, s, t) + R_02(x, j, s) * R_03(x, l, i, t) + R_02(x, j, t) * R_03(x, l, i, s)
				+ R_02(x, s, l) * R_03(x, i, j, t) + R_02(x, l, t) * R_03(x, j, i, s) + R_02(x, s, t) * R_03(x, l, i, j)
				)
			+ (R_01(x, l) + sign * d[2][l]) * (R_02(x, i, j) * R_03(x, k, s, t) + R_02(x, i, k) * R_03(x, j, s, t) + R_02(x, i, s) * R_03(x, k, j, t) + R_02(x, i, t) * R_03(x, k, s, j)
				+ R_02(x, j, k) * R_03(x, i, s, t) + R_02(x, j, s) * R_03(x, k, i, t) + R_02(x, j, t) * R_03(x, k, i, s)
				+ R_02(x, s, k) * R_03(x, i, j, t) + R_02(x, k, t) * R_03(x, j, i, s) + R_02(x, s, t) * R_03(x, k, i, j)
				)
			+ (R_01(x, s) + sign * d[2][s]) * (R_02(x, i, j) * R_03(x, k, l, t) + R_02(x, i, k) * R_03(x, j, l, t) + R_02(x, i, l) * R_03(x, k, j, t) + R_02(x, i, t) * R_03(x, k, l, j)
				+ R_02(x, j, k) * R_03(x, i, l, t) + R_02(x, j, l) * R_03(x, k, i, t) + R_02(x, j, t) * R_03(x, k, i, l)
				+ R_02(x, l, k) * R_03(x, i, j, t) + R_02(x, k, t) * R_03(x, j, i, l) + R_02(x, l, t) * R_03(x, k, i, j)
				)
			+ (R_01(x, t) + sign * d[2][t]) * (R_02(x, i, j) * R_03(x, k, l, s) + R_02(x, i, k) * R_03(x, j, l, s) + R_02(x, i, l) * R_03(x, k, j, s) + R_02(x, i, s) * R_03(x, k, l, j)
				+ R_02(x, j, k) * R_03(x, i, l, s) + R_02(x, j, l) * R_03(x, k, i, s) + R_02(x, j, s) * R_03(x, k, i, l)
				+ R_02(x, l, k) * R_03(x, i, j, s) + R_02(x, k, s) * R_03(x, j, i, l) + R_02(x, l, s) * R_03(x, k, i, j)
				)

			+ R_02(x, i, j) * (R_02(x, k, l) * R_02(x, s, t) + R_02(x, k, s) * R_02(x, l, t) + R_02(x, k, t) * R_02(x, s, l))
			+ R_02(x, i, k) * (R_02(x, j, l) * R_02(x, s, t) + R_02(x, j, s) * R_02(x, l, t) + R_02(x, j, t) * R_02(x, s, l))
			+ R_02(x, i, l) * (R_02(x, k, j) * R_02(x, s, t) + R_02(x, k, s) * R_02(x, j, t) + R_02(x, k, t) * R_02(x, s, j))
			+ R_02(x, i, s) * (R_02(x, k, l) * R_02(x, j, t) + R_02(x, k, j) * R_02(x, l, t) + R_02(x, k, t) * R_02(x, j, l))
			+ R_02(x, i, t) * (R_02(x, k, l) * R_02(x, s, j) + R_02(x, k, s) * R_02(x, l, j) + R_02(x, k, j) * R_02(x, s, l))


			)
		- 6 * R1 * R1 * ((R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * R_03(x, l, s, t) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * R_03(x, k, s, t)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * R_03(x, k, l, t) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, t) + sign * d[2][t]) * R_03(x, k, l, s)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * R_03(x, j, s, t) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * R_03(x, j, l, t)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * R_03(x, j, l, s) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * R_03(x, k, j, t)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * R_03(x, j, k, s) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * R_03(x, j, l, k)

			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * R_03(x, i, s, t) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * R_03(x, i, l, t)
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * R_03(x, i, l, s) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * R_03(x, i, k, t)
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * R_03(x, i, k, s) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * R_03(x, i, l, k)

			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * R_03(x, i, j, t) + (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * R_03(x, i, j, s)
			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * R_03(x, i, j, l) + (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * R_03(x, i, j, k)

			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_02(x, k, l) * R_02(x, s, t) + R_02(x, k, s) * R_02(x, l, t) + R_02(x, k, t) * R_02(x, s, l))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_02(x, j, l) * R_02(x, s, t) + R_02(x, j, s) * R_02(x, l, t) + R_02(x, j, t) * R_02(x, s, l))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_02(x, j, k) * R_02(x, s, t) + R_02(x, j, s) * R_02(x, k, t) + R_02(x, j, t) * R_02(x, s, k))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s]) * (R_02(x, j, k) * R_02(x, l, t) + R_02(x, j, l) * R_02(x, k, t) + R_02(x, j, t) * R_02(x, l, k))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, j, k) * R_02(x, l, s) + R_02(x, j, l) * R_02(x, k, s) + R_02(x, j, s) * R_02(x, l, k))

			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_02(x, i, l) * R_02(x, s, t) + R_02(x, i, s) * R_02(x, l, t) + R_02(x, i, t) * R_02(x, s, l))
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_02(x, i, k) * R_02(x, s, t) + R_02(x, i, s) * R_02(x, k, t) + R_02(x, i, t) * R_02(x, s, k))
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_02(x, i, k) * R_02(x, l, t) + R_02(x, i, l) * R_02(x, k, t) + R_02(x, i, t) * R_02(x, l, k))
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, i, k) * R_02(x, l, s) + R_02(x, i, l) * R_02(x, k, s) + R_02(x, i, s) * R_02(x, l, k))

			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_02(x, i, j) * R_02(x, s, t) + R_02(x, i, s) * R_02(x, j, t) + R_02(x, i, t) * R_02(x, s, j))
			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_02(x, i, j) * R_02(x, l, t) + R_02(x, i, l) * R_02(x, j, t) + R_02(x, i, t) * R_02(x, l, j))
			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, i, j) * R_02(x, l, s) + R_02(x, i, l) * R_02(x, j, s) + R_02(x, i, s) * R_02(x, l, j))

			+ (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_02(x, i, j) * R_02(x, k, t) + R_02(x, i, k) * R_02(x, j, t) + R_02(x, i, t) * R_02(x, k, j))
			+ (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, i, j) * R_02(x, k, s) + R_02(x, i, k) * R_02(x, j, s) + R_02(x, i, s) * R_02(x, k, j))

			+ (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, i, j) * R_02(x, k, l) + R_02(x, i, k) * R_02(x, j, l) + R_02(x, i, l) * R_02(x, k, j))



			)
		+ 24 * R1 * (R_02(x, i, j) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) + R_02(x, i, k) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])
			+ R_02(x, i, l) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) + R_02(x, i, s) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, t) + sign * d[2][t])
			+ R_02(x, i, t) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, j, k) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])
			+ R_02(x, j, l) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) + R_02(x, j, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t])
			+ R_02(x, j, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) + R_02(x, k, l) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])
			+ R_02(x, k, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) + R_02(x, k, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s])
			+ R_02(x, l, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) + R_02(x, l, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s])
			+ R_02(x, t, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]))
		- 120 * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]));
	return lnR_ijklst;
}

double lnR_07(double* x, int i, int j, int k, int l, int s, int t, int p, int sign)
{
	double R1 = R(x) + sign * x[2];
	double lnR_ijklstp = 1 / (R1 * R1 * R1 * R1 * R1 * R1 * R1) * (R1 * R1 * R1 * R1 * R1 * R1 * R_07(x, i, j, k, l, s, t, p)
		- R1 * R1 * R1 * R1 * R1 * ((R_01(x, p) + sign * d[2][p]) * R_06(x, i, j, k, l, s, t) + (R_01(x, t) + sign * d[2][t]) * R_06(x, i, j, k, l, s, p) + (R_01(x, i) + sign * d[2][i]) * R_06(x, j, k, l, s, t, p) + (R_01(x, j) + sign * d[2][j]) * R_06(x, i, k, l, s, t, p)
			+ (R_01(x, k) + sign * d[2][k]) * R_06(x, i, j, l, s, t, p) + (R_01(x, l) + sign * d[2][l]) * R_06(x, i, j, k, s, t, p) + (R_01(x, s) + sign * d[2][s]) * R_06(x, i, j, k, l, t, p)



			+R_02(x, i, j) * R_05(x, k, l, s, t, p) + R_02(x, i, k) * R_05(x, j, l, s, t, p) + R_02(x, i, l) * R_05(x, j, k, s, t, p) + R_02(x, i, s) * R_05(x, j, k, l, t, p)
			+ R_02(x, i, t) * R_05(x, j, k, l, s, p) + R_02(x, i, p) * R_05(x, j, k, l, s, t)

			+ R_02(x, j, k) * R_05(x, i, l, s, t, p) + R_02(x, j, l) * R_05(x, i, k, s, t, p) + R_02(x, j, s) * R_05(x, i, k, l, t, p) + R_02(x, j, t) * R_05(x, i, k, l, s, p)
			+ R_02(x, j, p) * R_05(x, i, k, l, s, t)

			+ R_02(x, k, l) * R_05(x, i, j, s, t, p) + R_02(x, k, s) * R_05(x, i, j, l, t, p) + R_02(x, k, t) * R_05(x, i, j, l, s, p) + R_02(x, k, p) * R_05(x, i, j, l, s, t)

			+ R_02(x, l, s) * R_05(x, i, j, k, t, p) + R_02(x, l, t) * R_05(x, i, j, k, s, p) + R_02(x, l, p) * R_05(x, i, j, k, s, t)

			+ R_02(x, s, t) * R_05(x, i, j, k, l, p) + R_02(x, s, p) * R_05(x, i, j, k, l, t)

			+ R_02(x, t, p) * R_05(x, i, j, k, l, s)


			+ R_03(x, i, j, k) * R_04(x, l, s, t, p) + R_03(x, i, j, l) * R_04(x, k, s, t, p) + R_03(x, i, j, s) * R_04(x, l, k, t, p) + R_03(x, i, j, t) * R_04(x, l, s, k, p)
			+ R_03(x, i, j, p) * R_04(x, l, s, t, k)
			+ R_03(x, i, k, l) * R_04(x, j, s, t, p) + R_03(x, i, k, s) * R_04(x, l, j, t, p) + R_03(x, i, k, t) * R_04(x, l, s, j, p) + R_03(x, i, k, p) * R_04(x, l, s, t, j)
			+ R_03(x, i, l, s) * R_04(x, k, j, t, p) + R_03(x, i, l, t) * R_04(x, k, s, j, p) + R_03(x, i, l, p) * R_04(x, k, s, t, j)
			+ R_03(x, i, s, t) * R_04(x, k, l, j, p) + R_03(x, i, s, p) * R_04(x, k, l, t, j)
			+ R_03(x, i, t, p) * R_04(x, k, l, s, j)

			+ R_03(x, j, k, l) * R_04(x, i, s, t, p) + R_03(x, j, k, s) * R_04(x, l, i, t, p) + R_03(x, j, k, t) * R_04(x, l, s, i, p) + R_03(x, j, k, p) * R_04(x, l, s, t, i)
			+ R_03(x, j, l, s) * R_04(x, k, i, t, p) + R_03(x, j, l, t) * R_04(x, k, s, i, p) + R_03(x, j, l, p) * R_04(x, k, s, t, i)
			+ R_03(x, j, s, t) * R_04(x, k, l, i, p) + R_03(x, j, s, p) * R_04(x, k, l, t, i)
			+ R_03(x, j, t, p) * R_04(x, k, l, s, i)

			+ R_03(x, k, l, s) * R_04(x, j, i, t, p) + R_03(x, k, l, t) * R_04(x, j, s, i, p) + R_03(x, k, l, p) * R_04(x, j, s, t, i)
			+ R_03(x, k, s, t) * R_04(x, j, l, i, p) + R_03(x, k, s, p) * R_04(x, j, l, t, i)
			+ R_03(x, k, t, p) * R_04(x, j, l, s, i)

			+ R_03(x, l, s, t) * R_04(x, j, k, i, p) + R_03(x, l, s, p) * R_04(x, j, k, t, i)
			+ R_03(x, l, t, p) * R_04(x, j, k, s, i)

			+ R_03(x, s, t, p) * R_04(x, j, k, l, i)


			)
		+ 2 * R1 * R1 * R1 * R1 * ((R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * R_05(x, k, l, s, t, p) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * R_05(x, j, l, s, t, p) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * R_05(x, j, k, s, t, p)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s]) * R_05(x, j, k, l, t, p) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, t) + sign * d[2][t]) * R_05(x, j, k, l, s, p) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) * R_05(x, j, k, l, s, t)
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * R_05(x, i, l, s, t, p) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * R_05(x, i, k, s, t, p) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * R_05(x, i, k, l, t, p)
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, t) + sign * d[2][t]) * R_05(x, i, k, l, s, p) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) * R_05(x, i, k, l, s, t)
			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * R_05(x, i, j, s, t, p) + (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * R_05(x, i, j, l, t, p) + (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * R_05(x, i, j, l, s, p)
			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) * R_05(x, i, j, l, s, t)
			+ (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * R_05(x, i, j, k, t, p) + (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * R_05(x, i, j, k, s, p) + (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) * R_05(x, i, j, k, s, t)
			+ (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * R_05(x, i, j, k, l, p) + (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) * R_05(x, i, j, k, l, t) + (R_01(x, p) + sign * d[2][p]) * (R_01(x, t) + sign * d[2][t]) * R_05(x, i, j, k, l, s)

			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, j, k) + (R_01(x, j) + sign * d[2][j]) * R_02(x, i, k) + (R_01(x, k) + sign * d[2][k]) * R_02(x, j, i)) * R_04(x, p, l, s, t)
			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, j, l) + (R_01(x, j) + sign * d[2][j]) * R_02(x, i, l) + (R_01(x, l) + sign * d[2][l]) * R_02(x, j, i)) * R_04(x, p, k, s, t)
			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, j, s) + (R_01(x, j) + sign * d[2][j]) * R_02(x, i, s) + (R_01(x, s) + sign * d[2][s]) * R_02(x, j, i)) * R_04(x, p, k, l, t)
			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, j, t) + (R_01(x, j) + sign * d[2][j]) * R_02(x, i, t) + (R_01(x, t) + sign * d[2][t]) * R_02(x, j, i)) * R_04(x, p, k, l, s)
			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, j, p) + (R_01(x, j) + sign * d[2][j]) * R_02(x, i, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, j, i)) * R_04(x, t, k, l, s)

			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, k, l) + (R_01(x, k) + sign * d[2][k]) * R_02(x, i, l) + (R_01(x, l) + sign * d[2][l]) * R_02(x, k, i)) * R_04(x, p, j, s, t)
			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, k, s) + (R_01(x, k) + sign * d[2][k]) * R_02(x, i, s) + (R_01(x, s) + sign * d[2][s]) * R_02(x, k, i)) * R_04(x, p, j, l, t)
			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, k, t) + (R_01(x, k) + sign * d[2][k]) * R_02(x, i, t) + (R_01(x, t) + sign * d[2][t]) * R_02(x, k, i)) * R_04(x, p, j, l, s)
			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, k, p) + (R_01(x, k) + sign * d[2][k]) * R_02(x, i, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, k, i)) * R_04(x, t, j, l, s)

			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, l, s) + (R_01(x, l) + sign * d[2][l]) * R_02(x, i, s) + (R_01(x, s) + sign * d[2][s]) * R_02(x, l, i)) * R_04(x, p, j, k, t)
			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, l, t) + (R_01(x, l) + sign * d[2][l]) * R_02(x, i, t) + (R_01(x, t) + sign * d[2][t]) * R_02(x, l, i)) * R_04(x, p, j, k, s)
			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, l, p) + (R_01(x, l) + sign * d[2][l]) * R_02(x, i, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, l, i)) * R_04(x, t, j, k, s)

			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, s, t) + (R_01(x, s) + sign * d[2][s]) * R_02(x, i, t) + (R_01(x, t) + sign * d[2][t]) * R_02(x, s, i)) * R_04(x, p, j, k, l)
			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, s, p) + (R_01(x, s) + sign * d[2][s]) * R_02(x, i, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, s, i)) * R_04(x, t, j, k, l)

			+ ((R_01(x, i) + sign * d[2][i]) * R_02(x, t, p) + (R_01(x, t) + sign * d[2][t]) * R_02(x, i, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, t, i)) * R_04(x, s, j, k, l)


			+ ((R_01(x, j) + sign * d[2][j]) * R_02(x, k, l) + (R_01(x, k) + sign * d[2][k]) * R_02(x, j, l) + (R_01(x, l) + sign * d[2][l]) * R_02(x, k, j)) * R_04(x, p, i, s, t)
			+ ((R_01(x, j) + sign * d[2][j]) * R_02(x, k, s) + (R_01(x, k) + sign * d[2][k]) * R_02(x, j, s) + (R_01(x, s) + sign * d[2][s]) * R_02(x, k, j)) * R_04(x, p, i, l, t)
			+ ((R_01(x, j) + sign * d[2][j]) * R_02(x, k, t) + (R_01(x, k) + sign * d[2][k]) * R_02(x, j, t) + (R_01(x, t) + sign * d[2][t]) * R_02(x, k, j)) * R_04(x, p, i, l, s)
			+ ((R_01(x, j) + sign * d[2][j]) * R_02(x, k, p) + (R_01(x, k) + sign * d[2][k]) * R_02(x, j, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, k, j)) * R_04(x, t, i, l, s)

			+ ((R_01(x, j) + sign * d[2][j]) * R_02(x, l, s) + (R_01(x, l) + sign * d[2][l]) * R_02(x, j, s) + (R_01(x, s) + sign * d[2][s]) * R_02(x, l, j)) * R_04(x, p, i, k, t)
			+ ((R_01(x, j) + sign * d[2][j]) * R_02(x, l, t) + (R_01(x, l) + sign * d[2][l]) * R_02(x, j, t) + (R_01(x, t) + sign * d[2][t]) * R_02(x, l, j)) * R_04(x, p, i, k, s)
			+ ((R_01(x, j) + sign * d[2][j]) * R_02(x, l, p) + (R_01(x, l) + sign * d[2][l]) * R_02(x, j, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, l, j)) * R_04(x, t, i, k, s)

			+ ((R_01(x, j) + sign * d[2][j]) * R_02(x, s, t) + (R_01(x, s) + sign * d[2][s]) * R_02(x, j, t) + (R_01(x, t) + sign * d[2][t]) * R_02(x, s, j)) * R_04(x, p, i, k, l)
			+ ((R_01(x, j) + sign * d[2][j]) * R_02(x, s, p) + (R_01(x, s) + sign * d[2][s]) * R_02(x, j, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, s, j)) * R_04(x, t, i, k, l)

			+ ((R_01(x, j) + sign * d[2][j]) * R_02(x, t, p) + (R_01(x, t) + sign * d[2][t]) * R_02(x, j, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, t, j)) * R_04(x, s, i, k, l)


			+ ((R_01(x, k) + sign * d[2][k]) * R_02(x, l, s) + (R_01(x, l) + sign * d[2][l]) * R_02(x, k, s) + (R_01(x, s) + sign * d[2][s]) * R_02(x, l, k)) * R_04(x, p, i, j, t)
			+ ((R_01(x, k) + sign * d[2][k]) * R_02(x, l, t) + (R_01(x, l) + sign * d[2][l]) * R_02(x, k, t) + (R_01(x, t) + sign * d[2][t]) * R_02(x, l, k)) * R_04(x, p, i, j, s)
			+ ((R_01(x, k) + sign * d[2][k]) * R_02(x, l, p) + (R_01(x, l) + sign * d[2][l]) * R_02(x, k, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, l, k)) * R_04(x, t, i, j, s)

			+ ((R_01(x, k) + sign * d[2][k]) * R_02(x, s, t) + (R_01(x, s) + sign * d[2][s]) * R_02(x, k, t) + (R_01(x, t) + sign * d[2][t]) * R_02(x, s, k)) * R_04(x, p, i, j, l)
			+ ((R_01(x, k) + sign * d[2][k]) * R_02(x, s, p) + (R_01(x, s) + sign * d[2][s]) * R_02(x, k, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, s, k)) * R_04(x, t, i, j, l)

			+ ((R_01(x, k) + sign * d[2][k]) * R_02(x, t, p) + (R_01(x, t) + sign * d[2][t]) * R_02(x, k, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, t, k)) * R_04(x, s, i, j, l)


			+ ((R_01(x, l) + sign * d[2][l]) * R_02(x, s, t) + (R_01(x, s) + sign * d[2][s]) * R_02(x, l, t) + (R_01(x, t) + sign * d[2][t]) * R_02(x, s, l)) * R_04(x, p, i, j, k)
			+ ((R_01(x, l) + sign * d[2][l]) * R_02(x, s, p) + (R_01(x, s) + sign * d[2][s]) * R_02(x, l, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, s, l)) * R_04(x, t, i, j, k)

			+ ((R_01(x, l) + sign * d[2][l]) * R_02(x, t, p) + (R_01(x, t) + sign * d[2][t]) * R_02(x, l, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, t, l)) * R_04(x, s, i, j, k)

			+ ((R_01(x, s) + sign * d[2][s]) * R_02(x, t, p) + (R_01(x, t) + sign * d[2][t]) * R_02(x, s, p) + (R_01(x, p) + sign * d[2][p]) * R_02(x, t, s)) * R_04(x, l, i, j, k)


			+ (R_01(x, i) + sign * d[2][i]) * (R_03(x, j, k, l) * R_03(x, s, t, p) + R_03(x, j, k, s) * R_03(x, l, t, p) + R_03(x, j, k, t) * R_03(x, s, l, p) + R_03(x, j, k, p) * R_03(x, s, t, l)
				+ R_03(x, j, l, s) * R_03(x, k, t, p) + R_03(x, j, l, t) * R_03(x, s, k, p) + R_03(x, j, l, p) * R_03(x, s, t, k)
				+ R_03(x, j, s, t) * R_03(x, l, k, p) + R_03(x, j, s, p) * R_03(x, l, t, k) + R_03(x, j, t, p) * R_03(x, l, s, k))
			+ (R_01(x, j) + sign * d[2][j]) * (R_03(x, i, k, l) * R_03(x, s, t, p) + R_03(x, i, k, s) * R_03(x, l, t, p) + R_03(x, i, k, t) * R_03(x, s, l, p) + R_03(x, i, k, p) * R_03(x, s, t, l)
				+ R_03(x, i, l, s) * R_03(x, k, t, p) + R_03(x, i, l, t) * R_03(x, s, k, p) + R_03(x, i, l, p) * R_03(x, s, t, k)
				+ R_03(x, i, s, t) * R_03(x, l, k, p) + R_03(x, i, s, p) * R_03(x, l, t, k) + R_03(x, i, t, p) * R_03(x, l, s, k))
			+ (R_01(x, k) + sign * d[2][k]) * (R_03(x, i, j, l) * R_03(x, s, t, p) + R_03(x, i, j, s) * R_03(x, l, t, p) + R_03(x, i, j, t) * R_03(x, s, l, p) + R_03(x, i, j, p) * R_03(x, s, t, l)
				+ R_03(x, i, l, s) * R_03(x, j, t, p) + R_03(x, i, l, t) * R_03(x, s, j, p) + R_03(x, i, l, p) * R_03(x, s, t, j)
				+ R_03(x, i, s, t) * R_03(x, l, j, p) + R_03(x, i, s, p) * R_03(x, l, t, j) + R_03(x, i, t, p) * R_03(x, l, s, j))
			+ (R_01(x, l) + sign * d[2][l]) * (R_03(x, i, j, k) * R_03(x, s, t, p) + R_03(x, i, j, s) * R_03(x, k, t, p) + R_03(x, i, j, t) * R_03(x, s, k, p) + R_03(x, i, j, p) * R_03(x, s, t, k)
				+ R_03(x, i, k, s) * R_03(x, j, t, p) + R_03(x, i, k, t) * R_03(x, s, j, p) + R_03(x, i, k, p) * R_03(x, s, t, j)
				+ R_03(x, i, s, t) * R_03(x, k, j, p) + R_03(x, i, s, p) * R_03(x, k, t, j) + R_03(x, i, t, p) * R_03(x, k, s, j))
			+ (R_01(x, s) + sign * d[2][s]) * (R_03(x, i, j, k) * R_03(x, l, t, p) + R_03(x, i, j, l) * R_03(x, k, t, p) + R_03(x, i, j, t) * R_03(x, l, k, p) + R_03(x, i, j, p) * R_03(x, l, t, k)
				+ R_03(x, i, k, l) * R_03(x, j, t, p) + R_03(x, i, k, t) * R_03(x, l, j, p) + R_03(x, i, k, p) * R_03(x, l, t, j)
				+ R_03(x, i, l, t) * R_03(x, k, j, p) + R_03(x, i, l, p) * R_03(x, k, t, j) + R_03(x, i, t, p) * R_03(x, k, l, j))
			+ (R_01(x, t) + sign * d[2][t]) * (R_03(x, i, j, k) * R_03(x, l, s, p) + R_03(x, i, j, l) * R_03(x, k, s, p) + R_03(x, i, j, s) * R_03(x, l, k, p) + R_03(x, i, j, p) * R_03(x, l, s, k)
				+ R_03(x, i, k, l) * R_03(x, j, s, p) + R_03(x, i, k, s) * R_03(x, l, j, p) + R_03(x, i, k, p) * R_03(x, l, s, j)
				+ R_03(x, i, l, s) * R_03(x, k, j, p) + R_03(x, i, l, p) * R_03(x, k, s, j) + R_03(x, i, s, p) * R_03(x, k, l, j))
			+ (R_01(x, p) + sign * d[2][p]) * (R_03(x, i, j, k) * R_03(x, l, s, t) + R_03(x, i, j, l) * R_03(x, k, s, t) + R_03(x, i, j, s) * R_03(x, l, k, t) + R_03(x, i, j, t) * R_03(x, l, s, k)
				+ R_03(x, i, k, l) * R_03(x, j, s, t) + R_03(x, i, k, s) * R_03(x, l, j, t) + R_03(x, i, k, t) * R_03(x, l, s, j)
				+ R_03(x, i, l, s) * R_03(x, k, j, t) + R_03(x, i, l, t) * R_03(x, k, s, j) + R_03(x, i, s, t) * R_03(x, k, l, j))


			+ R_03(x, i, j, k) * (R_02(x, l, s) * R_02(x, t, p) + R_02(x, l, t) * R_02(x, s, p) + R_02(x, l, p) * R_02(x, s, t))
			+ R_03(x, i, j, l) * (R_02(x, k, s) * R_02(x, t, p) + R_02(x, k, t) * R_02(x, s, p) + R_02(x, k, p) * R_02(x, s, t))
			+ R_03(x, i, j, s) * (R_02(x, k, l) * R_02(x, t, p) + R_02(x, k, t) * R_02(x, l, p) + R_02(x, k, p) * R_02(x, l, t))
			+ R_03(x, i, j, t) * (R_02(x, k, l) * R_02(x, s, p) + R_02(x, k, s) * R_02(x, l, p) + R_02(x, k, p) * R_02(x, l, s))
			+ R_03(x, i, j, p) * (R_02(x, k, l) * R_02(x, s, t) + R_02(x, k, s) * R_02(x, l, t) + R_02(x, k, t) * R_02(x, l, s))

			+ R_03(x, i, k, l) * (R_02(x, j, s) * R_02(x, t, p) + R_02(x, j, t) * R_02(x, s, p) + R_02(x, j, p) * R_02(x, s, t))
			+ R_03(x, i, k, s) * (R_02(x, j, l) * R_02(x, t, p) + R_02(x, j, t) * R_02(x, l, p) + R_02(x, j, p) * R_02(x, l, t))
			+ R_03(x, i, k, t) * (R_02(x, j, l) * R_02(x, s, p) + R_02(x, j, s) * R_02(x, l, p) + R_02(x, j, p) * R_02(x, l, s))
			+ R_03(x, i, k, p) * (R_02(x, j, l) * R_02(x, s, t) + R_02(x, j, s) * R_02(x, l, t) + R_02(x, j, t) * R_02(x, l, s))

			+ R_03(x, i, l, s) * (R_02(x, j, k) * R_02(x, t, p) + R_02(x, j, t) * R_02(x, k, p) + R_02(x, j, p) * R_02(x, k, t))
			+ R_03(x, i, l, t) * (R_02(x, j, k) * R_02(x, s, p) + R_02(x, j, s) * R_02(x, k, p) + R_02(x, j, p) * R_02(x, k, s))
			+ R_03(x, i, l, p) * (R_02(x, j, k) * R_02(x, s, t) + R_02(x, j, s) * R_02(x, k, t) + R_02(x, j, t) * R_02(x, k, s))

			+ R_03(x, i, s, t) * (R_02(x, j, k) * R_02(x, l, p) + R_02(x, j, l) * R_02(x, k, p) + R_02(x, j, p) * R_02(x, k, l))
			+ R_03(x, i, s, p) * (R_02(x, j, k) * R_02(x, l, t) + R_02(x, j, l) * R_02(x, k, t) + R_02(x, j, t) * R_02(x, k, l))

			+ R_03(x, i, t, p) * (R_02(x, j, k) * R_02(x, l, s) + R_02(x, j, l) * R_02(x, k, s) + R_02(x, j, s) * R_02(x, k, l))

			+ R_03(x, j, k, l) * (R_02(x, i, s) * R_02(x, t, p) + R_02(x, i, t) * R_02(x, s, p) + R_02(x, i, p) * R_02(x, s, t))
			+ R_03(x, j, k, s) * (R_02(x, i, l) * R_02(x, t, p) + R_02(x, i, t) * R_02(x, l, p) + R_02(x, i, p) * R_02(x, l, t))
			+ R_03(x, j, k, t) * (R_02(x, i, l) * R_02(x, s, p) + R_02(x, i, s) * R_02(x, l, p) + R_02(x, i, p) * R_02(x, l, s))
			+ R_03(x, j, k, p) * (R_02(x, i, l) * R_02(x, s, t) + R_02(x, i, s) * R_02(x, l, t) + R_02(x, i, t) * R_02(x, l, s))

			+ R_03(x, j, l, s) * (R_02(x, i, k) * R_02(x, t, p) + R_02(x, i, t) * R_02(x, k, p) + R_02(x, i, p) * R_02(x, k, t))
			+ R_03(x, j, l, t) * (R_02(x, i, k) * R_02(x, s, p) + R_02(x, i, s) * R_02(x, k, p) + R_02(x, i, p) * R_02(x, k, s))
			+ R_03(x, j, l, p) * (R_02(x, i, k) * R_02(x, s, t) + R_02(x, i, s) * R_02(x, k, t) + R_02(x, i, t) * R_02(x, k, s))

			+ R_03(x, j, s, t) * (R_02(x, i, k) * R_02(x, l, p) + R_02(x, i, l) * R_02(x, k, p) + R_02(x, i, p) * R_02(x, k, l))
			+ R_03(x, j, s, p) * (R_02(x, i, k) * R_02(x, l, t) + R_02(x, i, l) * R_02(x, k, t) + R_02(x, i, t) * R_02(x, k, l))

			+ R_03(x, j, t, p) * (R_02(x, i, k) * R_02(x, l, s) + R_02(x, i, l) * R_02(x, k, s) + R_02(x, i, s) * R_02(x, k, l))

			+ R_03(x, k, l, s) * (R_02(x, i, j) * R_02(x, t, p) + R_02(x, i, t) * R_02(x, j, p) + R_02(x, i, p) * R_02(x, j, t))
			+ R_03(x, k, l, t) * (R_02(x, i, j) * R_02(x, s, p) + R_02(x, i, s) * R_02(x, j, p) + R_02(x, i, p) * R_02(x, j, s))
			+ R_03(x, k, l, p) * (R_02(x, i, j) * R_02(x, s, t) + R_02(x, i, s) * R_02(x, j, t) + R_02(x, i, t) * R_02(x, j, s))

			+ R_03(x, k, s, t) * (R_02(x, i, j) * R_02(x, l, p) + R_02(x, i, l) * R_02(x, j, p) + R_02(x, i, p) * R_02(x, j, l))
			+ R_03(x, k, s, p) * (R_02(x, i, j) * R_02(x, l, t) + R_02(x, i, l) * R_02(x, j, t) + R_02(x, i, t) * R_02(x, j, l))

			+ R_03(x, k, t, p) * (R_02(x, i, j) * R_02(x, l, s) + R_02(x, i, l) * R_02(x, j, s) + R_02(x, i, s) * R_02(x, j, l))

			+ R_03(x, l, s, t) * (R_02(x, i, j) * R_02(x, k, p) + R_02(x, i, k) * R_02(x, j, p) + R_02(x, i, p) * R_02(x, j, k))
			+ R_03(x, l, s, p) * (R_02(x, i, j) * R_02(x, k, t) + R_02(x, i, k) * R_02(x, j, t) + R_02(x, i, t) * R_02(x, j, k))

			+ R_03(x, l, t, p) * (R_02(x, i, j) * R_02(x, k, s) + R_02(x, i, k) * R_02(x, j, s) + R_02(x, i, s) * R_02(x, j, k))

			+ R_03(x, s, t, p) * (R_02(x, i, j) * R_02(x, k, l) + R_02(x, i, k) * R_02(x, j, l) + R_02(x, i, l) * R_02(x, j, k))



			)


		- 6 * R1 * R1 * R1 * ((R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * R_04(x, l, s, t, p) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * R_04(x, k, s, t, p)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * R_04(x, k, l, t, p) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, k, l, s, p)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, k, l, s, t)

			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * R_04(x, j, s, t, p) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * R_04(x, j, l, t, p)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, j, l, s, p) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, j, l, s, t)

			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * R_04(x, j, k, t, p) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, j, k, s, p)
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, j, k, s, t)

			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, j, k, l, p) + (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, j, k, l, t)

			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, j, k, l, s)


			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * R_04(x, i, s, t, p) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * R_04(x, i, l, t, p)
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, i, l, s, p) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, i, l, s, t)

			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * R_04(x, i, k, t, p) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, i, k, s, p)
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, i, k, s, t)

			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, i, k, l, p) + (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, i, k, l, t)

			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, i, k, l, s)


			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * R_04(x, i, j, t, p) + (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, i, j, s, p)
			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, i, j, s, t)

			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, i, j, l, p) + (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, i, j, l, t)

			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, i, j, l, s)


			+ (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * R_04(x, i, j, k, p) + (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, i, j, k, t)

			+ (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, i, j, k, s)


			+ (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) * R_04(x, i, j, k, l)

			+ R_03(x, i, j, k) * (R_02(x, l, s) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, t) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, t, s) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, s, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, l) + sign * d[2][l]) + R_02(x, t, p) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]))
			+ R_03(x, i, j, l) * (R_02(x, k, s) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, t) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, t, s) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, s, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, t, p) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]))
			+ R_03(x, i, j, s) * (R_02(x, l, k) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, t) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, k) + sign * d[2][k])
				+ R_02(x, t, k) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, l) + sign * d[2][l]) + R_02(x, t, p) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, k) + sign * d[2][k]))
			+ R_03(x, i, j, t) * (R_02(x, l, s) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, k) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, p) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, k, s) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, s, p) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) + R_02(x, k, p) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]))
			+ R_03(x, i, j, p) * (R_02(x, l, s) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, l, t) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, l, k) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, t, s) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, s, k) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, l) + sign * d[2][l]) + R_02(x, t, k) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]))


			+ R_03(x, i, k, l) * (R_02(x, j, s) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, t) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, t, s) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, s, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, t, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]))
			+ R_03(x, i, k, s) * (R_02(x, l, j) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, t) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, j) + sign * d[2][j])
				+ R_02(x, t, j) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, l) + sign * d[2][l]) + R_02(x, t, p) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, j) + sign * d[2][j]))
			+ R_03(x, i, k, t) * (R_02(x, l, s) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, j) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, j, s) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, s, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) + R_02(x, j, p) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]))
			+ R_03(x, i, k, p) * (R_02(x, l, s) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, l, t) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, l, j) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, t, s) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, s, j) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, l) + sign * d[2][l]) + R_02(x, t, j) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]))


			+ R_03(x, i, l, s) * (R_02(x, k, j) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, t) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, j) + sign * d[2][j])
				+ R_02(x, t, j) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, t, p) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, j) + sign * d[2][j]))
			+ R_03(x, i, l, t) * (R_02(x, k, s) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, j) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, j, s) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, s, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, j, p) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]))
			+ R_03(x, i, l, p) * (R_02(x, k, s) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, k, t) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, k, j) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, t, s) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, s, j) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, t, j) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]))


			+ R_03(x, i, s, t) * (R_02(x, k, l) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, j) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l])
				+ R_02(x, j, l) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, j, p) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]))
			+ R_03(x, i, s, p) * (R_02(x, k, l) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, k, t) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, k, j) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, l) + sign * d[2][l])
				+ R_02(x, t, l) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, l, j) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, t, j) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]))


			+ R_03(x, i, t, p) * (R_02(x, k, l) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, k, s) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, k, j) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, l) + sign * d[2][l])
				+ R_02(x, s, l) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, l, j) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, s, j) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]))


			+ R_03(x, j, k, l) * (R_02(x, i, s) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, i, t) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, i, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, t, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, s, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, t, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s]))
			+ R_03(x, j, k, s) * (R_02(x, l, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, i) + sign * d[2][i])
				+ R_02(x, t, i) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, i, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, l) + sign * d[2][l]) + R_02(x, t, p) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, i) + sign * d[2][i]))
			+ R_03(x, j, k, t) * (R_02(x, l, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, i) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, i, s) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, s, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) + R_02(x, i, p) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]))
			+ R_03(x, j, k, p) * (R_02(x, l, s) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, l, t) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, l, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, t, s) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, s, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, l) + sign * d[2][l]) + R_02(x, t, i) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]))


			+ R_03(x, j, l, s) * (R_02(x, k, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, i) + sign * d[2][i])
				+ R_02(x, t, i) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, i, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, t, p) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, i) + sign * d[2][i]))
			+ R_03(x, j, l, t) * (R_02(x, k, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, i) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, i, s) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, s, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, i, p) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]))
			+ R_03(x, j, l, p) * (R_02(x, k, s) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, k, t) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, k, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, t, s) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, s, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, t, i) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]))


			+ R_03(x, j, s, t) * (R_02(x, k, l) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, i) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l])
				+ R_02(x, i, l) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, i, p) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]))
			+ R_03(x, j, s, p) * (R_02(x, k, l) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, k, t) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, k, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, l) + sign * d[2][l])
				+ R_02(x, t, l) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, l, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, t, i) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]))


			+ R_03(x, j, t, p) * (R_02(x, k, l) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, k, s) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, k, i) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, l) + sign * d[2][l])
				+ R_02(x, s, l) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, l, i) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, k) + sign * d[2][k]) + R_02(x, s, i) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]))


			+ R_03(x, k, l, s) * (R_02(x, j, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, i) + sign * d[2][i])
				+ R_02(x, t, i) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, i, p) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, t, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, i) + sign * d[2][i]))
			+ R_03(x, k, l, t) * (R_02(x, j, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, i) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, i, s) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, s, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, i, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]))
			+ R_03(x, k, l, p) * (R_02(x, j, s) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, t) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, s) + sign * d[2][s])
				+ R_02(x, t, s) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, s, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, t, i) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]))


			+ R_03(x, k, s, t) * (R_02(x, j, l) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, i) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l])
				+ R_02(x, i, l) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, l, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, i, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]))
			+ R_03(x, k, s, p) * (R_02(x, j, l) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, t) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, l) + sign * d[2][l])
				+ R_02(x, t, l) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, l, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, t, i) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]))


			+ R_03(x, k, t, p) * (R_02(x, j, l) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, s) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, i) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, l) + sign * d[2][l])
				+ R_02(x, s, l) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, l, i) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, s, i) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]))


			+ R_03(x, l, s, t) * (R_02(x, j, k) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, i) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, j, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k])
				+ R_02(x, i, k) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) + R_02(x, k, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, i, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]))
			+ R_03(x, l, s, p) * (R_02(x, j, k) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, t) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, k) + sign * d[2][k])
				+ R_02(x, t, k) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, k, i) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, t, i) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]))


			+ R_03(x, l, t, p) * (R_02(x, j, k) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, s) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, i) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, k) + sign * d[2][k])
				+ R_02(x, s, k) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, k, i) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, s, i) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]))


			+ R_03(x, s, t, p) * (R_02(x, j, k) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, l) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, j, i) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, k) + sign * d[2][k])
				+ R_02(x, l, k) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, i) + sign * d[2][i]) + R_02(x, k, i) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, j) + sign * d[2][j]) + R_02(x, l, i) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]))

			+ (R_01(x, i) + sign * d[2][i]) * (R_02(x, j, k) * (R_02(x, s, t) * R_02(x, l, p) + R_02(x, s, l) * R_02(x, p, t) + R_02(x, s, p) * R_02(x, l, t))
				+ R_02(x, j, l) * (R_02(x, s, t) * R_02(x, k, p) + R_02(x, s, k) * R_02(x, p, t) + R_02(x, s, p) * R_02(x, k, t))
				+ R_02(x, j, s) * (R_02(x, k, l) * R_02(x, t, p) + R_02(x, k, t) * R_02(x, l, p) + R_02(x, k, p) * R_02(x, l, t))
				+ R_02(x, j, t) * (R_02(x, k, l) * R_02(x, s, p) + R_02(x, k, s) * R_02(x, l, p) + R_02(x, k, p) * R_02(x, l, s))
				+ R_02(x, j, p) * (R_02(x, k, l) * R_02(x, s, t) + R_02(x, k, s) * R_02(x, l, t) + R_02(x, k, t) * R_02(x, l, s)))
			+ (R_01(x, j) + sign * d[2][j]) * (R_02(x, i, k) * (R_02(x, s, t) * R_02(x, l, p) + R_02(x, s, l) * R_02(x, p, t) + R_02(x, s, p) * R_02(x, l, t))
				+ R_02(x, i, l) * (R_02(x, s, t) * R_02(x, k, p) + R_02(x, s, k) * R_02(x, p, t) + R_02(x, s, p) * R_02(x, k, t))
				+ R_02(x, i, s) * (R_02(x, k, l) * R_02(x, t, p) + R_02(x, k, t) * R_02(x, l, p) + R_02(x, k, p) * R_02(x, l, t))
				+ R_02(x, i, t) * (R_02(x, k, l) * R_02(x, s, p) + R_02(x, k, s) * R_02(x, l, p) + R_02(x, k, p) * R_02(x, l, s))
				+ R_02(x, i, p) * (R_02(x, k, l) * R_02(x, s, t) + R_02(x, k, s) * R_02(x, l, t) + R_02(x, k, t) * R_02(x, l, s)))
			+ (R_01(x, k) + sign * d[2][k]) * (R_02(x, i, j) * (R_02(x, s, t) * R_02(x, l, p) + R_02(x, s, l) * R_02(x, p, t) + R_02(x, s, p) * R_02(x, l, t))
				+ R_02(x, i, l) * (R_02(x, s, t) * R_02(x, j, p) + R_02(x, s, j) * R_02(x, p, t) + R_02(x, s, p) * R_02(x, j, t))
				+ R_02(x, i, s) * (R_02(x, j, l) * R_02(x, t, p) + R_02(x, j, t) * R_02(x, l, p) + R_02(x, j, p) * R_02(x, l, t))
				+ R_02(x, i, t) * (R_02(x, j, l) * R_02(x, s, p) + R_02(x, j, s) * R_02(x, l, p) + R_02(x, j, p) * R_02(x, l, s))
				+ R_02(x, i, p) * (R_02(x, j, l) * R_02(x, s, t) + R_02(x, j, s) * R_02(x, l, t) + R_02(x, j, t) * R_02(x, l, s)))
			+ (R_01(x, l) + sign * d[2][l]) * (R_02(x, i, j) * (R_02(x, s, t) * R_02(x, k, p) + R_02(x, s, k) * R_02(x, p, t) + R_02(x, s, p) * R_02(x, k, t))
				+ R_02(x, i, k) * (R_02(x, s, t) * R_02(x, j, p) + R_02(x, s, j) * R_02(x, p, t) + R_02(x, s, p) * R_02(x, j, t))
				+ R_02(x, i, s) * (R_02(x, j, k) * R_02(x, t, p) + R_02(x, j, t) * R_02(x, k, p) + R_02(x, j, p) * R_02(x, k, t))
				+ R_02(x, i, t) * (R_02(x, j, k) * R_02(x, s, p) + R_02(x, j, s) * R_02(x, k, p) + R_02(x, j, p) * R_02(x, k, s))
				+ R_02(x, i, p) * (R_02(x, j, k) * R_02(x, s, t) + R_02(x, j, s) * R_02(x, k, t) + R_02(x, j, t) * R_02(x, k, s)))
			+ (R_01(x, s) + sign * d[2][s]) * (R_02(x, i, j) * (R_02(x, l, t) * R_02(x, k, p) + R_02(x, l, k) * R_02(x, p, t) + R_02(x, l, p) * R_02(x, k, t))
				+ R_02(x, i, k) * (R_02(x, l, t) * R_02(x, j, p) + R_02(x, l, j) * R_02(x, p, t) + R_02(x, l, p) * R_02(x, j, t))
				+ R_02(x, i, l) * (R_02(x, j, k) * R_02(x, t, p) + R_02(x, j, t) * R_02(x, k, p) + R_02(x, j, p) * R_02(x, k, t))
				+ R_02(x, i, t) * (R_02(x, j, k) * R_02(x, l, p) + R_02(x, j, l) * R_02(x, k, p) + R_02(x, j, p) * R_02(x, k, l))
				+ R_02(x, i, p) * (R_02(x, j, k) * R_02(x, l, t) + R_02(x, j, l) * R_02(x, k, t) + R_02(x, j, t) * R_02(x, k, l)))
			+ (R_01(x, t) + sign * d[2][t]) * (R_02(x, i, j) * (R_02(x, l, s) * R_02(x, k, p) + R_02(x, l, k) * R_02(x, p, s) + R_02(x, l, p) * R_02(x, k, s))
				+ R_02(x, i, k) * (R_02(x, l, s) * R_02(x, j, p) + R_02(x, l, j) * R_02(x, p, s) + R_02(x, l, p) * R_02(x, j, s))
				+ R_02(x, i, l) * (R_02(x, j, k) * R_02(x, s, p) + R_02(x, j, s) * R_02(x, k, p) + R_02(x, j, p) * R_02(x, k, s))
				+ R_02(x, i, s) * (R_02(x, j, k) * R_02(x, l, p) + R_02(x, j, l) * R_02(x, k, p) + R_02(x, j, p) * R_02(x, k, l))
				+ R_02(x, i, p) * (R_02(x, j, k) * R_02(x, l, s) + R_02(x, j, l) * R_02(x, k, s) + R_02(x, j, s) * R_02(x, k, l)))
			+ (R_01(x, p) + sign * d[2][p]) * (R_02(x, i, j) * (R_02(x, l, s) * R_02(x, k, t) + R_02(x, l, k) * R_02(x, t, s) + R_02(x, l, t) * R_02(x, k, s))
				+ R_02(x, i, k) * (R_02(x, l, s) * R_02(x, j, t) + R_02(x, l, j) * R_02(x, t, s) + R_02(x, l, t) * R_02(x, j, s))
				+ R_02(x, i, l) * (R_02(x, j, k) * R_02(x, s, t) + R_02(x, j, s) * R_02(x, k, t) + R_02(x, j, t) * R_02(x, k, s))
				+ R_02(x, i, s) * (R_02(x, j, k) * R_02(x, l, t) + R_02(x, j, l) * R_02(x, k, t) + R_02(x, j, t) * R_02(x, k, l))
				+ R_02(x, i, t) * (R_02(x, j, k) * R_02(x, l, s) + R_02(x, j, l) * R_02(x, k, s) + R_02(x, j, s) * R_02(x, k, l)))




			)

		+ 24 * R1 * R1 * (R_03(x, i, j, k) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) + R_03(x, i, j, l) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p])
			+ R_03(x, i, j, s) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, t) + sign * d[2][t]) + R_03(x, i, j, t) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p])
			+ R_03(x, i, j, p) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])

			+ R_03(x, i, k, l) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) + R_03(x, i, k, s) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p])
			+ R_03(x, i, k, t) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, s) + sign * d[2][s]) + R_03(x, i, k, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])

			+ R_03(x, i, l, s) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) + R_03(x, i, l, t) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p])
			+ R_03(x, i, l, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, s) + sign * d[2][s])

			+ R_03(x, i, s, t) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) + R_03(x, i, s, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, k) + sign * d[2][k])

			+ R_03(x, i, t, p) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s])


			+ R_03(x, j, k, l) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) + R_03(x, j, k, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t])
			+ R_03(x, j, k, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, l) + sign * d[2][l]) + R_03(x, j, k, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t])

			+ R_03(x, j, l, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) + R_03(x, j, l, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s])
			+ R_03(x, j, l, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])

			+ R_03(x, j, s, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) + R_03(x, j, s, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t])

			+ R_03(x, j, p, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s])


			+ R_03(x, k, l, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, t) + sign * d[2][t]) + R_03(x, k, l, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, p) + sign * d[2][p]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s])
			+ R_03(x, k, l, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])

			+ R_03(x, k, s, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) + R_03(x, k, s, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t])

			+ R_03(x, k, t, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s])


			+ R_03(x, l, s, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) + R_03(x, l, s, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t])

			+ R_03(x, l, t, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s])


			+ R_03(x, s, t, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l])


			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_02(x, l, s) * R_02(x, t, p) + R_02(x, l, t) * R_02(x, s, p) + R_02(x, l, p) * R_02(x, s, t))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_02(x, k, s) * R_02(x, t, p) + R_02(x, k, t) * R_02(x, s, p) + R_02(x, k, p) * R_02(x, s, t))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_02(x, l, k) * R_02(x, t, p) + R_02(x, l, t) * R_02(x, k, p) + R_02(x, l, p) * R_02(x, k, t))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, l, s) * R_02(x, k, p) + R_02(x, l, k) * R_02(x, s, p) + R_02(x, l, p) * R_02(x, s, k))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, l, s) * R_02(x, t, k) + R_02(x, l, t) * R_02(x, s, k) + R_02(x, l, k) * R_02(x, s, t))

			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_02(x, j, s) * R_02(x, t, p) + R_02(x, j, t) * R_02(x, s, p) + R_02(x, j, p) * R_02(x, s, t))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_02(x, l, j) * R_02(x, t, p) + R_02(x, l, t) * R_02(x, j, p) + R_02(x, l, p) * R_02(x, j, t))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, l, s) * R_02(x, j, p) + R_02(x, l, j) * R_02(x, s, p) + R_02(x, l, p) * R_02(x, s, j))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, l, s) * R_02(x, t, j) + R_02(x, l, t) * R_02(x, s, j) + R_02(x, l, j) * R_02(x, s, t))


			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_02(x, k, j) * R_02(x, t, p) + R_02(x, k, t) * R_02(x, j, p) + R_02(x, k, p) * R_02(x, j, t))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, k, s) * R_02(x, j, p) + R_02(x, k, j) * R_02(x, s, p) + R_02(x, k, p) * R_02(x, s, j))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, k, s) * R_02(x, t, j) + R_02(x, k, t) * R_02(x, s, j) + R_02(x, k, j) * R_02(x, s, t))



			+(R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, k, l) * R_02(x, j, p) + R_02(x, k, j) * R_02(x, l, p) + R_02(x, k, p) * R_02(x, l, j))
			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, k, l) * R_02(x, t, j) + R_02(x, k, t) * R_02(x, l, j) + R_02(x, k, j) * R_02(x, l, t))

			+ (R_01(x, i) + sign * d[2][i]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, k, l) * R_02(x, s, j) + R_02(x, k, s) * R_02(x, l, j) + R_02(x, k, j) * R_02(x, l, s))


			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_02(x, i, s) * R_02(x, t, p) + R_02(x, i, t) * R_02(x, s, p) + R_02(x, i, p) * R_02(x, s, t))
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_02(x, l, i) * R_02(x, t, p) + R_02(x, l, t) * R_02(x, i, p) + R_02(x, l, p) * R_02(x, i, t))
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, l, s) * R_02(x, i, p) + R_02(x, l, i) * R_02(x, s, p) + R_02(x, l, p) * R_02(x, s, i))
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, l, s) * R_02(x, t, i) + R_02(x, l, t) * R_02(x, s, i) + R_02(x, l, i) * R_02(x, s, t))


			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_02(x, k, i) * R_02(x, t, p) + R_02(x, k, t) * R_02(x, i, p) + R_02(x, k, p) * R_02(x, i, t))
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, k, s) * R_02(x, i, p) + R_02(x, k, i) * R_02(x, s, p) + R_02(x, k, p) * R_02(x, s, i))
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, k, s) * R_02(x, t, i) + R_02(x, k, t) * R_02(x, s, i) + R_02(x, k, i) * R_02(x, s, t))



			+(R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, k, l) * R_02(x, i, p) + R_02(x, k, i) * R_02(x, l, p) + R_02(x, k, p) * R_02(x, l, i))
			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, k, l) * R_02(x, t, i) + R_02(x, k, t) * R_02(x, l, i) + R_02(x, k, i) * R_02(x, l, t))

			+ (R_01(x, j) + sign * d[2][j]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, k, l) * R_02(x, s, i) + R_02(x, k, s) * R_02(x, l, i) + R_02(x, k, i) * R_02(x, l, s))



			+(R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_02(x, j, i) * R_02(x, t, p) + R_02(x, j, t) * R_02(x, i, p) + R_02(x, j, p) * R_02(x, i, t))
			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, j, s) * R_02(x, i, p) + R_02(x, j, i) * R_02(x, s, p) + R_02(x, j, p) * R_02(x, s, i))
			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, j, s) * R_02(x, t, i) + R_02(x, j, t) * R_02(x, s, i) + R_02(x, j, i) * R_02(x, s, t))


			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, j, l) * R_02(x, i, p) + R_02(x, j, i) * R_02(x, l, p) + R_02(x, j, p) * R_02(x, l, i))
			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, j, l) * R_02(x, t, i) + R_02(x, j, t) * R_02(x, l, i) + R_02(x, j, i) * R_02(x, l, t))

			+ (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, j, l) * R_02(x, s, i) + R_02(x, j, s) * R_02(x, l, i) + R_02(x, j, i) * R_02(x, l, s))


			+ (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_02(x, j, k) * R_02(x, i, p) + R_02(x, j, i) * R_02(x, k, p) + R_02(x, j, p) * R_02(x, k, i))
			+ (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, j, k) * R_02(x, t, i) + R_02(x, j, t) * R_02(x, k, i) + R_02(x, j, i) * R_02(x, k, t))

			+ (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, j, k) * R_02(x, s, i) + R_02(x, j, s) * R_02(x, k, i) + R_02(x, j, i) * R_02(x, k, s))

			+ (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]) * (R_02(x, j, k) * R_02(x, l, i) + R_02(x, j, l) * R_02(x, k, i) + R_02(x, j, i) * R_02(x, k, l))

			)

		- 120 * R1 * (R_02(x, i, j) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p])
			+ R_02(x, i, k) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p])
			+ R_02(x, i, l) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p])
			+ R_02(x, i, s) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p])
			+ R_02(x, i, t) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, p) + sign * d[2][p])
			+ R_02(x, i, p) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, t) + sign * d[2][t])

			+ R_02(x, j, k) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p])
			+ R_02(x, j, l) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p])
			+ R_02(x, j, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p])
			+ R_02(x, j, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p])
			+ R_02(x, j, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, s) + sign * d[2][s])

			+ R_02(x, k, l) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p])
			+ R_02(x, k, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p])
			+ R_02(x, k, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p])
			+ R_02(x, k, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])

			+ R_02(x, l, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p])
			+ R_02(x, l, t) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, p) + sign * d[2][p])
			+ R_02(x, l, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t])

			+ R_02(x, t, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, p) + sign * d[2][p])
			+ R_02(x, p, s) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, t) + sign * d[2][t])

			+ R_02(x, t, p) * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]))

		+ 720 * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]) * (R_01(x, l) + sign * d[2][l]) * (R_01(x, s) + sign * d[2][s]) * (R_01(x, t) + sign * d[2][t]) * (R_01(x, p) + sign * d[2][p]));
	return lnR_ijklstp;
}





////////////////////////////////////////////////////////////////////////
/////////////////////////Constitutive//////////////////////////////////
double Constitutive(int i, int j, int k, int l)
{
	double c_ijkl = 2.0 * mu_0 * nu0 / (1 - 2 * nu0) * d[i][j] * d[k][l] + 1.0 * mu_0 * (d[i][k] * d[j][l] + d[i][l] * d[j][k]);
	return c_ijkl;
}
////////////////////////////////////////////////////////////////


////////////////////////Constitutive//////////////////////////////////
double Constitutive1(int i, int j, int k, int l)
{
	double c_ijkl = 2.0 * mu_1 * nu1 / (1 - 2 * nu1) * d[i][j] * d[k][l] + 1.0 * mu_1 * (d[i][k] * d[j][l] + d[i][l] * d[j][k]);
	return c_ijkl;
}
////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////
double D_40(int i, int j, int k, int l, int s, int h, Ref<VectorXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p)
{
	double x1[3];

	x1[0] = x_p(s, 0) - x_o(h, 0);
	x1[1] = x_p(s, 1) - x_o(h, 1);
	x1[2] = x_p(s, 2) - x_o(h, 2);

	double a = radius[h];
	double D_4 = 1 / (8 * Pi * (1 - nu0)) * (Psi_04(x1, a, k, l, i, j) - 2 * nu0 * d[k][l] * Phi_02(x1, a, i, j) - (1 - nu0) * (Phi_02(x1, a, k, j) * d[i][l] + Phi_02(x1, a, k, i) * d[j][l] + Phi_02(x1, a, l, j) * d[i][k] + Phi_02(x1, a, l, i) * d[j][k]));


	return D_4;
}


double D_50(int i, int j, int k, int l, int p, int s, int h, Ref<VectorXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p)
{
	double x1[3];

	x1[0] = x_p(s, 0) - x_o(h, 0);
	x1[1] = x_p(s, 1) - x_o(h, 1);
	x1[2] = x_p(s, 2) - x_o(h, 2);


	double a = radius[h];
	double D_5 = 1 / (8 * Pi * (1 - nu0)) * (Psi_14(x1, a, p, k, l, i, j) - 2 * nu0 * d[k][l] * Phi_12(x1, a, p, i, j) - (1 - nu0) * (Phi_12(x1, a, p, k, j) * d[i][l] + Phi_12(x1, a, p, k, i) * d[j][l] + Phi_12(x1, a, p, l, j) * d[i][k] + Phi_12(x1, a, p, l, i) * d[j][k]));

	return D_5;
}


double D_60(int i, int j, int k, int l, int p, int q, int s, int h, Ref<VectorXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p)
{

	double x1[3];

	x1[0] = x_p(s, 0) - x_o(h, 0);
	x1[1] = x_p(s, 1) - x_o(h, 1);
	x1[2] = x_p(s, 2) - x_o(h, 2);


	double a = radius[h];
	double D_6 = 1 / (8 * Pi * (1 - nu0)) * (Psi_24(x1, a, p, q, k, l, i, j) - 2 * nu0 * d[k][l] * Phi_22(x1, a, p, q, i, j) - (1 - nu0) * (Phi_22(x1, a, p, q, k, j) * d[i][l] + Phi_22(x1, a, p, q, k, i) * d[j][l] + Phi_22(x1, a, p, q, l, j) * d[i][k] + Phi_22(x1, a, p, q, l, i) * d[j][k]));

	return D_6;
}


/////////////////////////

double D_41(int i, int j, int k, int l, int r, int s, int h, Ref<VectorXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p)
{

	double x1[3];

	x1[0] = x_p(s, 0) - x_o(h, 0);
	x1[1] = x_p(s, 1) - x_o(h, 1);
	x1[2] = x_p(s, 2) - x_o(h, 2);




	double a = radius[h];
	double D_4_r = 1 / (8 * Pi * (1 - nu0)) * (Psi_05(x1, a, k, l, i, j, r) - 2 * nu0 * d[k][l] * Phi_03(x1, a, i, j, r) - (1 - nu0) * (Phi_03(x1, a, k, j, r) * d[i][l] + Phi_03(x1, a, k, i, r) * d[j][l] + Phi_03(x1, a, l, j, r) * d[i][k] + Phi_03(x1, a, l, i, r) * d[j][k]));

	return D_4_r;
}
//////////////////////////////////
//////////////////////////////////


double D_51(int i, int j, int k, int l, int p, int r, int s, int h, Ref<VectorXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p)
{
	double x1[3];

	x1[0] = x_p(s, 0) - x_o(h, 0);
	x1[1] = x_p(s, 1) - x_o(h, 1);
	x1[2] = x_p(s, 2) - x_o(h, 2);



	double a = radius[h];
	double D_5_r = 1 / (8 * Pi * (1 - nu0)) * (Psi_15(x1, a, p, k, l, i, j, r) - 2 * nu0 * d[k][l] * Phi_13(x1, a, p, i, j, r) - (1 - nu0) * (Phi_13(x1, a, p, k, j, r) * d[i][l] + Phi_13(x1, a, p, k, i, r) * d[j][l] + Phi_13(x1, a, p, l, j, r) * d[i][k] + Phi_13(x1, a, p, l, i, r) * d[j][k]));


	return D_5_r;
}

//////////////////////////////////


double D_61(int i, int j, int k, int l, int p, int q, int r, int s, int h, Ref<VectorXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p)
{

	double x1[3];

	x1[0] = x_p(s, 0) - x_o(h, 0);
	x1[1] = x_p(s, 1) - x_o(h, 1);
	x1[2] = x_p(s, 2) - x_o(h, 2);



	double a = radius[h];
	double D_6_r = 1 / (8 * Pi * (1 - nu0)) * (Psi_25(x1, a, p, q, k, l, i, j, r) - 2 * nu0 * d[k][l] * Phi_23(x1, a, p, q, i, j, r) - (1 - nu0) * (Phi_23(x1, a, p, q, k, j, r) * d[i][l] + Phi_23(x1, a, p, q, k, i, r) * d[j][l] + Phi_23(x1, a, p, q, l, j, r) * d[i][k] + Phi_23(x1, a, p, q, l, i, r) * d[j][k]));

	return D_6_r;
}


//////////////////////////////////////////////////////
double D_42(int i, int j, int k, int l, int r, int z, int s, int h, Ref<VectorXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p)
{

	double x1[3];

	x1[0] = x_p(s, 0) - x_o(h, 0);
	x1[1] = x_p(s, 1) - x_o(h, 1);
	x1[2] = x_p(s, 2) - x_o(h, 2);


	double a = radius[h];
	double D_4_rz = 1 / (8 * Pi * (1 - nu0)) * (Psi_06(x1, a, k, l, i, j, r, z) - 2 * nu0 * d[k][l] * Phi_04(x1, a, i, j, r, z) - (1 - nu0) * (Phi_04(x1, a, k, j, r, z) * d[i][l] + Phi_04(x1, a, k, i, r, z) * d[j][l] + Phi_04(x1, a, l, j, r, z) * d[i][k] + Phi_04(x1, a, l, i, r, z) * d[j][k]));

	return D_4_rz;
}
///////////////////////////////////////////////////

//////////////////////////////////////////////////////
double D_52(int i, int j, int k, int l, int p, int r, int z, int s, int h, Ref<VectorXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p)
{


	double x1[3];

	x1[0] = x_p(s, 0) - x_o(h, 0);
	x1[1] = x_p(s, 1) - x_o(h, 1);
	x1[2] = x_p(s, 2) - x_o(h, 2);

	double a = radius[h];
	double D_5_rz = 1 / (8 * Pi * (1 - nu0)) * (Psi_16(x1, a, p, k, l, i, j, r, z) - 2 * nu0 * d[k][l] * Phi_14(x1, a, p, i, j, r, z) - (1 - nu0) * (Phi_14(x1, a, p, k, j, r, z) * d[i][l] + Phi_14(x1, a, p, k, i, r, z) * d[j][l] + Phi_14(x1, a, p, l, j, r, z) * d[i][k] + Phi_14(x1, a, p, l, i, r, z) * d[j][k]));

	return D_5_rz;
}
/////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////
double D_62(int i, int j, int k, int l, int p, int q, int r, int z, int s, int h, Ref<VectorXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p)
{


	double x1[3];

	x1[0] = x_p(s, 0) - x_o(h, 0);
	x1[1] = x_p(s, 1) - x_o(h, 1);
	x1[2] = x_p(s, 2) - x_o(h, 2);

	double a = radius[h];
	double D_6_rz = 1 / (8 * Pi * (1 - nu0)) * (Psi_26(x1, a, p, q, k, l, i, j, r, z) - 2 * nu0 * d[k][l] * Phi_24(x1, a, p, q, i, j, r, z) - (1 - nu0) * (Phi_24(x1, a, p, q, k, j, r, z) * d[i][l] + Phi_24(x1, a, p, q, k, i, r, z) * d[j][l] + Phi_24(x1, a, p, q, l, j, r, z) * d[i][k] + Phi_24(x1, a, p, q, l, i, r, z) * d[j][k]));

	return D_6_rz;
}


///////////////////////////////////////////////////////////////////////////////////////

double Green(double* x, int i, int j) {

	double G_ij = d[i][j] / (4 * Pi * mu_0) * fR(x, 0) - 1 / (16 * Pi * mu_0 * (1 - nu0)) * R_02(x, i, j);

	return G_ij;
}

double Green_01(double* x, int i, int j, int k) {

	double G_ijk = d[i][j] / (4 * Pi * mu_0) * fR_01(x, k, 0) - 1 / (16 * Pi * mu_0 * (1 - nu0)) * R_03(x, i, j, k);

	return G_ijk;
}

double Green_02(double* x, int i, int j, int k, int l) {

	double G_ijkl = d[i][j] / (4 * Pi * mu_0) * fR_02(x, k, l, 0) - 1 / (16 * Pi * mu_0 * (1 - nu0)) * R_04(x, i, j, k, l);

	return G_ijkl;
}

double Green_03(double* x, int i, int j, int k, int l, int s) {

	double G_ijkls = d[i][j] / (4 * Pi * mu_0) * fR_03(x, k, l, s, 0) - 1 / (16 * Pi * mu_0 * (1 - nu0)) * R_05(x, i, j, k, l, s);

	return G_ijkls;
}


double Green_04(double* x, int i, int j, int k, int l, int s, int t) {

	double G_ijklst = d[i][j] / (4 * Pi * mu_0) * fR_04(x, k, l, s, t, 0) - 1 / (16 * Pi * mu_0 * (1 - nu0)) * R_06(x, i, j, k, l, s, t);

	return G_ijklst;
}

double Green_05(double* x, int i, int j, int k, int l, int s, int t, int r) {

	double G_ijklstr = d[i][j] / (4 * Pi * mu_0) * fR_05(x, k, l, s, t, r, 0) - 1 / (16 * Pi * mu_0 * (1 - nu0)) * R_07(x, i, j, k, l, s, t, r);

	return G_ijklstr;
}

///////////////////////////////////////////////////////

double IGreen_01(double* x, double a, int i, int j, int l) {

	double IGreen_ij_l = 1 / (4 * Pi * mu_0) * (0.5 * d[i][j] * Phi_01(x, a, l) + 0.5 * d[i][l] * Phi_01(x, a, j)) - 1 / (16 * Pi * mu_0 * (1 - nu0)) * Psi_03(x, a, i, j, l);

	return IGreen_ij_l;
}


double IGreen_11(double* x, double a, int i, int j, int p, int l) {

	double IGreen_ijp_l = 1 / (4 * Pi * mu_0) * (0.5 * d[i][j] * Phi_11(x, a, p, l) + 0.5 * d[i][l] * Phi_11(x, a, p, j)) - 1 / (16 * Pi * mu_0 * (1 - nu0)) * Psi_13(x, a, p, i, j, l);

	return IGreen_ijp_l;
}

double IGreen_21(double* x, double a, int i, int j, int p, int q, int l) {

	double IGreen_ijpq_l = 1 / (4 * Pi * mu_0) * (0.5 * d[i][j] * Phi_21(x, a, p, q, l) + 0.5 * d[i][l] * Phi_21(x, a, p, q, j)) - 1 / (16 * Pi * mu_0 * (1 - nu0)) * Psi_23(x, a, p, q, i, j, l);

	return IGreen_ijpq_l;
}

