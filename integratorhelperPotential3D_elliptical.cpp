
#include "integratorhelperPotential3D_elliptical.h"
#include "configPotential3D_elliptical.h"
#include <math.h>
#include"Eigen/Dense"
#include"GlobalConstant.h"

using namespace Eigen;
using namespace Potential3DElliptical;

double Potential3DElliptical::I(double lambda, Ref<VectorXd> a, int orient, int type) {

	double a_l = a.maxCoeff();
	double a_s = a.minCoeff();
	double a_u = a(orient);
	double a_e = a(dd[orient + 1]);

	double b_bar = sqrt((a_u * a_u + lambda) / (a_e * a_e + lambda));

	double value = type == 0 ? acosh(b_bar) : acos(b_bar);

	double I = 4 * Pi * a_u * a_e * a_e * value / sqrt((a_l * a_l - a_s * a_s));

	return	I;

}


double  Potential3DElliptical::I_i(int i, int location, double lambda, Ref<VectorXd> a, int orient, int type)
{

	double I_i;

	double a_l = a.maxCoeff();
	double a_s = a.minCoeff();
	double a_u = a(orient);
	double a_e = a(dd[orient + 1]);

	if (location == 0) {
		lambda = 0;
	}

	double  b_bar = sqrt((a_u * a_u + lambda) / (a_e * a_e + lambda));
	double	d_bar = sqrt((a_l * a_l - a_s * a_s) / (a_e * a_e + lambda));
	double value = type == 0 ? acosh(b_bar) : acos(b_bar);


	if (i == orient) {
		I_i = 4 * Pi * a_u * a_e * a_e * abs((value - d_bar / b_bar)) / pow((a_l * a_l - a_s * a_s), 1.5);
	}
	else {
		I_i = 2 * Pi * a_u * a_e * a_e * abs((b_bar * d_bar - value)) / pow((a_l * a_l - a_s * a_s), 1.5);
	}

	return I_i;

}


double Potential3DElliptical::I_i_j(int i, int j, double* r_i, int location, double lambda, Ref<VectorXd> a)
{
	double Lambda_j, A_i, I_i_j;

	if (location == 0) {
		Lambda_j = 0;
		A_i = 1;
	}
	else {
		double delta = sqrt((a(0) * a(0) + lambda) * (a(1) * a(1) + lambda) * (a(2) * a(2) + lambda));
		A_i = -2 * Pi * a(0) * a(1) * a(2) / ((a(i) * a(i) + lambda) * delta);
		double F_1 = 2 * r_i[j] / (a(j) * a(j) + lambda);
		double F_2 = r_i[0] * r_i[0] / pow((a(0) * a(0) + lambda), 2) + r_i[1] * r_i[1] / pow((a(1) * a(1) + lambda), 2) + r_i[2] * r_i[2] / pow((a(2) * a(2) + lambda), 2);
		Lambda_j = F_1 / F_2;

	}

	I_i_j = A_i * Lambda_j;

	return I_i_j;
}


double  Potential3DElliptical::I_i_jk(int i, int j, int k, double* r_i, int location, double lambda, Ref<VectorXd> a)
{
	double Lambda_j, Lambda_k, Lambda_jk, A_i, Z_i, I_i_jk;

	if (location == 0) {
		Lambda_j = 0;
		Lambda_k = 0;
		Lambda_jk = 0;
		A_i = 1;
		Z_i = 1;
	}
	else
	{

		double delta = sqrt(((a(0) * a(0) + lambda) * (a(1) * a(1) + lambda) * (a(2) * a(2) + lambda)));

		A_i = -2 * Pi * a(0) * a(1) * a(2) / ((a(i) * a(i) + lambda) * delta);
		Z_i = 1.0 / (a(i) * a(i) + lambda) + 0.5 * (1.0 / (a(0) * a(0) + lambda) + 1.0 / (a(1) * a(1) + lambda) + 1.0 / (a(2) * a(2) + lambda));
		double F_j = 2 * r_i[j] / (a(j) * a(j) + lambda);
		double F_k = 2 * r_i[k] / (a(k) * a(k) + lambda);
		double C_1 = r_i[0] * r_i[0] / pow((a(0) * a(0) + lambda), 2) + r_i[1] * r_i[1] / pow((a(1) * a(1) + lambda), 2) + r_i[2] * r_i[2] / pow((a(2) * a(2) + lambda), 2);
		double C_2 = r_i[0] * r_i[0] / pow((a(0) * a(0) + lambda), 3) + r_i[1] * r_i[1] / pow((a(1) * a(1) + lambda), 3) + r_i[2] * r_i[2] / pow((a(2) * a(2) + lambda), 3);
		Lambda_j = F_j / C_1;
		Lambda_k = F_k / C_1;

		double F_jk = (2 * d[j][k] - F_j * Lambda_k) / (a(j) * a(j) + lambda);

		double C_1k = F_k / (a(k) * a(k) + lambda) - 2 * C_2 * Lambda_k;
		Lambda_jk = (F_jk - Lambda_j * C_1k) / C_1;

	}

	return I_i_jk = A_i * (Lambda_jk - Lambda_j * Lambda_k * Z_i);

}


double Potential3DElliptical::I_i_jkl(int i, int p, int q, int r, double* r_i, int location, double lambda, Ref<VectorXd>a)
{

	double I_i_jkl;
	if (location == 0)
	{
		I_i_jkl = 0;
	}
	else {
		double	delta = sqrt((a(0) * a(0) + lambda) * (a(1) * a(1) + lambda) * (a(2) * a(2) + lambda));
		double 	A_i = -2 * Pi * a(0) * a(1) * a(2) / ((a(i) * a(i) + lambda) * delta);
		double 	F_p = 2 * r_i[p] / (a(p) * a(p) + lambda);
		double 	F_q = 2 * r_i[q] / (a(q) * a(q) + lambda);
		double 	F_r = 2 * r_i[r] / (a(r) * a(r) + lambda);
		double 	C_1 = r_i[0] * r_i[0] / pow((a(0) * a(0) + lambda), 2) + r_i[1] * r_i[1] / pow((a(1) * a(1) + lambda), 2) + r_i[2] * r_i[2] / pow((a(2) * a(2) + lambda), 2);

		double 	Lambda_p = F_p / C_1;
		double 	Lambda_q = F_q / C_1;
		double 	Lambda_r = F_r / C_1;

		double 	C_2 = r_i[0] * r_i[0] / pow((a(0) * a(0) + lambda), 3) + r_i[1] * r_i[1] / pow((a(1) * a(1) + lambda), 3) + r_i[2] * r_i[2] / pow((a(2) * a(2) + lambda), 3);
		double 	F_pq = (2 * d[p][q] - F_p * Lambda_q) / (a(p) * a(p) + lambda);
		double 	F_qr = (2 * d[q][r] - F_q * Lambda_r) / (a(q) * a(q) + lambda);
		double 	F_pr = (2 * d[p][r] - F_p * Lambda_r) / (a(p) * a(p) + lambda);

		double 	C_1q = F_q / (a(q) * a(q) + lambda) - 2 * C_2 * Lambda_q;
		double 	C_1p = F_p / (a(p) * a(p) + lambda) - 2 * C_2 * Lambda_p;
		double	C_1r = F_r / (a(r) * a(r) + lambda) - 2 * C_2 * Lambda_r;

		double 	Lambda_pq = (F_pq - Lambda_p * C_1q) / C_1;
		double	Lambda_pr = (F_pr - Lambda_p * C_1r) / C_1;
		double	Lambda_qr = (F_qr - Lambda_q * C_1r) / C_1;

		double	C_3 = r_i[0] * r_i[0] / pow((a(0) * a(0) + lambda), 4) + r_i[1] * r_i[1] / pow((a(1) * a(1) + lambda), 4) + r_i[2] * r_i[2] / pow((a(2) * a(2) + lambda), 4);
		double	F_p_qr = -(F_pq * Lambda_r + F_pr * Lambda_q + F_p * Lambda_qr) / (a(p) * a(p) + lambda);

		double	C_2p = F_p / pow((a(p) * a(p) + lambda), 2) - 3 * C_3 * Lambda_p;
		double	C_2q = F_q / pow((a(q) * a(q) + lambda), 2) - 3 * C_3 * Lambda_q;
		double	C_2r = F_r / pow((a(r) * a(r) + lambda), 2) - 3 * C_3 * Lambda_r;

		double	C_1pq = (F_pq - F_p * Lambda_q / (a(p) * a(p) + lambda)) / (a(p) * a(p) + lambda) - 2 * (C_2q * Lambda_p + C_2 * Lambda_pq);
		double	C_1pr = (F_pr - F_p * Lambda_r / (a(p) * a(p) + lambda)) / (a(p) * a(p) + lambda) - 2 * (C_2r * Lambda_p + C_2 * Lambda_pr);
		double	C_1qr = (F_qr - F_q * Lambda_r / (a(q) * a(q) + lambda)) / (a(q) * a(q) + lambda) - 2 * (C_2r * Lambda_q + C_2 * Lambda_qr);

		double	Lambda_pqr = -(Lambda_pq * C_1r - F_p_qr + Lambda_pr * C_1q + Lambda_p * C_1qr) / C_1;

		double	Z_1 = 1 / (a(i) * a(i) + lambda) + 0.5 * (1.0 / (a(0) * a(0) + lambda) + 1.0 / (a(1) * a(1) + lambda) + 1.0 / (a(2) * a(2) + lambda));

		double	Z_2 = 1 / pow((a(i) * a(i) + lambda), 2) + 0.5 * (1.0 / pow((a(0) * a(0) + lambda), 2) + 1.0 / pow((a(1) * a(1) + lambda), 2) + 1.0 / pow((a(2) * a(2) + lambda), 2));

		I_i_jkl = A_i * (Lambda_pqr - (Lambda_pq * Lambda_r + Lambda_pr * Lambda_q + Lambda_qr * Lambda_p) * Z_1 + Lambda_p * Lambda_q * Lambda_r * (Z_2 + Z_1 * Z_1));

	}

	return I_i_jkl;
}



double Potential3DElliptical::I_ij(int i, int j, int location, double lambda, Ref<VectorXd> a, int orient, int type)
{
	double a_l = a.maxCoeff();
	double a_s = a.minCoeff();
	double a_u = a(orient);
	double a_e = a(dd[orient + 1]);

	int s;
	double I_si, I_ij;
	double delta = sqrt((a(0) * a(0) + lambda) * (a(1) * a(1) + lambda) * (a(2) * a(2) + lambda));


	if (i == j)  // For I_ij = I_11, I_22 or I_33
	{
		if (i == orient)
		{
			for (int z = 0; z < 3; z++)
			{// Pick random number s from[1, 2, 3], but s ~= "orient #"
				if (z != orient)
				{
					s = z;
				}

			}
			I_si = abs(I_i(s, location, lambda, a, orient, type) - I_i(i, location, lambda, a, orient, type)) / ((a_l * a_l - a_s * a_s));
			I_ij = ((4 * Pi * a_u * a_e * a_e) / ((a_u * a_u + lambda) * delta) - 2 * I_si) / 3;
		}

		else
		{
			s = orient;
			I_si = abs(I_i(s, location, lambda, a, orient, type) - I_i(i, location, lambda, a, orient, type)) / ((a_l * a_l - a_s * a_s));
			I_ij = ((4.0 * Pi * a_u * a_e * a_e) / ((a_e * a_e + lambda) * delta) - I_si) / 4.0;
		}

	}
	else  // For I_ij = I_12, I_13 or I_23
	{
		if (i == orient || j == orient) {
			I_ij = abs(I_i(i, location, lambda, a, orient, type) - I_i(j, location, lambda, a, orient, type)) / ((a_l * a_l - a_s * a_s));
		}
		else {
			s = orient;
			I_si = abs(I_i(s, location, lambda, a, orient, type) - I_i(i, location, lambda, a, orient, type)) / ((a_l * a_l - a_s * a_s));
			I_ij = ((4.0 * Pi * a_u * a_e * a_e) / ((a_e * a_e + lambda) * delta) - I_si) / 4;
		}

	}


	return I_ij;
}


double Potential3DElliptical::I_ij_k(int i, int j, int k, double* r_i, int location, double lambda, Ref<VectorXd> a)
{
	double A_ij, Lambda_k;
	if (location == 0) {
		Lambda_k = 0;
		A_ij = 1;
	}
	else {
		double delta = sqrt((a(0) * a(0) + lambda) * (a(1) * a(1) + lambda) * (a(2) * a(2) + lambda));
		A_ij = -2 * Pi * a(0) * a(1) * a(2) / ((a(i) * a(i) + lambda) * (a(j) * a(j) + lambda) * delta);
		double F_1 = 2.0 * r_i[k] / (a(k) * a(k) + lambda);
		double F_2 = r_i[0] * r_i[0] / pow((a(0) * a(0) + lambda), 2) + r_i[1] * r_i[1] / pow((a(1) * a(1) + lambda), 2) + r_i[2] * r_i[2] / pow((a(2) * a(2) + lambda), 2);
		Lambda_k = F_1 / F_2;
	}

	double I_ij_k = A_ij * Lambda_k;

	return I_ij_k;

}

double  Potential3DElliptical::I_ij_kl(int i, int j, int k, int l, double* r_i, int location, double lambda, Ref<VectorXd> a)
{
	double Lambda_l, I_ij_kl;
	if (location == 0)
	{
		Lambda_l = 0;
	}
	else {
		double F_1 = 2.0 * r_i[l] / (a(l) * a(l) + lambda);
		double F_2 = r_i[0] * r_i[0] / pow((a(0) * a(0) + lambda), 2) + r_i[1] * r_i[1] / pow((a(1) * a(1) + lambda), 2) + r_i[2] * r_i[2] / pow((a(2) * a(2) + lambda), 2);
		Lambda_l = F_1 / F_2;
	}

	I_ij_kl = -Lambda_l * I_i_j(i, k, r_i, location, lambda, a) / pow((a(j) * a(j) + lambda), 2) + I_i_jk(i, k, l, r_i, location, lambda, a) / (a(j) * a(j) + lambda);

	return I_ij_kl;

}


double Potential3DElliptical::I_ij_kls(int i, int j, int p, int q, int r, double* r_i, int location, double lambda, Ref<VectorXd> a)
{
	double I_ij_kls;
	if (location == 0)
	{
		I_ij_kls = 0;
	}
	else {
		double delta = sqrt((a(0) * a(0) + lambda) * (a(1) * a(1) + lambda) * (a(2) * a(2) + lambda));
		double A_ij = -2.0 * Pi * a(0) * a(1) * a(2) / ((a(i) * a(i) + lambda) * (a(j) * a(j) + lambda) * delta);
		double 	F_p = 2 * r_i[p] / (a(p) * a(p) + lambda);
		double 	F_q = 2 * r_i[q] / (a(q) * a(q) + lambda);
		double 	F_r = 2 * r_i[r] / (a(r) * a(r) + lambda);

		double C_1 = r_i[0] * r_i[0] / pow((a(0) * a(0) + lambda), 2) + r_i[1] * r_i[1] / pow((a(1) * a(1) + lambda), 2) + r_i[2] * r_i[2] / pow((a(2) * a(2) + lambda), 2);
		double Lambda_p = F_p / C_1;
		double Lambda_q = F_q / C_1;
		double Lambda_r = F_r / C_1;

		double C_2 = r_i[0] * r_i[0] / pow((a(0) * a(0) + lambda), 3) + r_i[1] * r_i[1] / pow((a(1) * a(1) + lambda), 3) + r_i[2] * r_i[2] / pow((a(2) * a(2) + lambda), 3);
		double F_pq = (2 * d[p][q] - F_p * Lambda_q) / (a(p) * a(p) + lambda);
		double F_qr = (2 * d[q][r] - F_q * Lambda_r) / (a(q) * a(q) + lambda);
		double F_pr = (2 * d[p][r] - F_p * Lambda_r) / (a(p) * a(p) + lambda);

		double C_1q = F_q / (a(q) * a(q) + lambda) - 2 * C_2 * Lambda_q;
		double C_1p = F_p / (a(p) * a(p) + lambda) - 2 * C_2 * Lambda_p;
		double C_1r = F_r / (a(r) * a(r) + lambda) - 2 * C_2 * Lambda_r;

		double Lambda_pq = (F_pq - Lambda_p * C_1q) / C_1;
		double Lambda_pr = (F_pr - Lambda_p * C_1r) / C_1;
		double Lambda_qr = (F_qr - Lambda_q * C_1r) / C_1;

		double C_3 = r_i[0] * r_i[0] / pow((a(0) * a(0) + lambda), 4) + r_i[1] * r_i[1] / pow((a(1) * a(1) + lambda), 4) + r_i[2] * r_i[2] / pow((a(2) * a(2) + lambda), 4);
		double F_p_qr = -(F_pq * Lambda_r + F_pr * Lambda_q + F_p * Lambda_qr) / (a(p) * a(p) + lambda);

		double C_2p = F_p / pow((a(p) * a(p) + lambda), 2) - 3 * C_3 * Lambda_p;
		double C_2q = F_q / pow((a(q) * a(q) + lambda), 2) - 3 * C_3 * Lambda_q;
		double C_2r = F_r / pow((a(r) * a(r) + lambda), 2) - 3 * C_3 * Lambda_r;

		double C_1pq = (F_pq - F_p * Lambda_q / (a(p) * a(p) + lambda)) / (a(p) * a(p) + lambda) - 2 * (C_2q * Lambda_p + C_2 * Lambda_pq);
		double C_1pr = (F_pr - F_p * Lambda_r / (a(p) * a(p) + lambda)) / (a(p) * a(p) + lambda) - 2 * (C_2r * Lambda_p + C_2 * Lambda_pr);
		double C_1qr = (F_qr - F_q * Lambda_r / (a(q) * a(q) + lambda)) / (a(q) * a(q) + lambda) - 2 * (C_2r * Lambda_q + C_2 * Lambda_qr);

		double Lambda_pqr = -(Lambda_pq * C_1r - F_p_qr + Lambda_pr * C_1q + Lambda_p * C_1qr) / C_1;

		double Z_1 = 1.0 / (a(i) * a(i) + lambda) + 1.0 / (a(j) * a(j) + lambda) + 0.5 * (1.0 / (a(0) * a(0) + lambda) + 1.0 / (a(1) * a(1) + lambda) + 1.0 / (a(2) * a(2) + lambda));

		double Z_2 = 1.0 / pow((a(i) * a(i) + lambda), 2) + 1.0 / pow((a(j) * a(j) + lambda), 2) + 0.5 * (1.0 / pow((a(0) * a(0) + lambda), 2) + 1.0 / pow((a(1) * a(1) + lambda), 2) + 1.0 / pow((a(2) * a(2) + lambda), 2));

		I_ij_kls = A_ij * (Lambda_pqr - (Lambda_pq * Lambda_r + Lambda_pr * Lambda_q + Lambda_qr * Lambda_p) * Z_1 + Lambda_p * Lambda_q * Lambda_r * (Z_2 + (Z_1) * (Z_1)));

	}

	return I_ij_kls;

}

double Potential3DElliptical::I_ijk(int i, int j, int k, int location, double lambda, Ref<VectorXd>a, int orient, int type)
{

	double a_l = a.maxCoeff();
	double a_s = a.minCoeff();
	double a_u = a(orient);
	double a_e = a(dd[orient + 1]);
	double I_ijk;
	int s;

	double delta = sqrt((a(0) * a(0) + lambda) * (a(1) * a(1) + lambda) * (a(2) * a(2) + lambda));



	if (i == j && j == k) // For I_ijk = I_111, I_222 or I_333
	{
		if (i == orient) {
			for (int z = 0; z < 3; z++) { // Pick random number s from[1, 2, 3], but s ~= "orient #"
				if (z != orient) {
					s = z;
				}
			}
			double I_si = abs(I_ij(s, i, location, lambda, a, orient, type) - I_ij(i, i, location, lambda, a, orient, type)) / ((a_l * a_l - a_s * a_s));
			I_ijk = ((4 * Pi * a_u * a_e * a_e) / (pow((a_u * a_u + lambda), 2) * delta) - 2 * I_si) / 5;
		}
		else {
			s = orient;
			double I_si = abs(I_ij(s, i, location, lambda, a, orient, type) - I_ij(i, i, location, lambda, a, orient, type)) / ((a_l * a_l - a_s * a_s));
			I_ijk = ((4 * Pi * a_u * a_e * a_e) / (pow((a_e * a_e + lambda), 2) * delta) - I_si) / 6;
		}
	}


	else if ((i == j) && (j != k)) { // % For I_ijk = I_112, I_113, I_221, , I_223, I_331 and I_332
		if ((i == orient) || (k == orient)) {
			I_ijk = abs(I_ij(i, k, location, lambda, a, orient, type) - I_ij(i, i, location, lambda, a, orient, type)) / ((a_l * a_l - a_s * a_s));
		}
		else {
			s = orient;
			double I_si = abs(I_ij(s, i, location, lambda, a, orient, type) - I_ij(i, i, location, lambda, a, orient, type)) / ((a_l * a_l - a_s * a_s));
			I_ijk = ((4 * Pi * a_u * a_e * a_e) / (pow((a_e * a_e + lambda), 2) * delta) - I_si) / 6;
		}
	}



	else if ((j == k) && (k != i)) {   // % For I_ijk = I_211, I_311, I_122, , I_322, I_133 and I_233
		if ((j == orient) || (i == orient)) {
			I_ijk = abs(I_ij(j, i, location, lambda, a, orient, type) - I_ij(j, j, location, lambda, a, orient, type)) / ((a_l * a_l - a_s * a_s));
		}
		else {
			s = orient;
			double I_si = abs(I_ij(s, j, location, lambda, a, orient, type) - I_ij(j, j, location, lambda, a, orient, type)) / ((a_l * a_l - a_s * a_s));
			I_ijk = ((4 * Pi * a_u * a_e * a_e) / (pow((a_e * a_e + lambda), 2) * delta) - I_si) / 6;
		}
	}

	else if ((k == i) && (i != j)) { // % For I_ijk = I_121, I_131, I_212, , I_232, I_313 and I_323
		if ((k == orient) || (j == orient)) {
			I_ijk = abs(I_ij(k, j, location, lambda, a, orient, type) - I_ij(k, k, location, lambda, a, orient, type)) / ((a_l * a_l - a_s * a_s));
		}
		else {
			s = orient;
			double I_si = abs(I_ij(s, k, location, lambda, a, orient, type) - I_ij(k, k, location, lambda, a, orient, type)) / ((a_l * a_l - a_s * a_s));
			I_ijk = ((4 * Pi * a_u * a_e * a_e) / (pow((a_e * a_e + lambda), 2) * delta) - I_si) / 6;
		}
	}

	else {                                          //    % For I_ijk = I_123
		for (int z = 0; z < 3; z++) { // % Pick random number s from[1, 2, 3], but s ~= "orient #"
			if (z != orient) {
				s = z;
			}
		}

		int t = orient;
		I_ijk = abs(I_ij(s, s, location, lambda, a, orient, type) - I_ij(s, t, location, lambda, a, orient, type)) / ((a_l * a_l - a_s * a_s));

	}

	return I_ijk;

}

double Potential3DElliptical::I_ijk_l(int i, int j, int k, int l, double* r_i, int location, double lambda, Ref<VectorXd> a)
{
	double A_ijk, Lambda_l, I_ijk_l;
	if (location == 0) {
		A_ijk = 1;
		Lambda_l = 0;
	}
	else {
		double delta = sqrt((a(0) * a(0) + lambda) * (a(1) * a(1) + lambda) * (a(2) * a(2) + lambda));
		A_ijk = -2 * Pi * a(0) * a(1) * a(2) / ((a(i) * a(i) + lambda) * (a(j) * a(j) + lambda) * (a(k) * a(k) + lambda) * delta);

		double F_1 = 2 * r_i[l] / (a(l) * a(l) + lambda);
		double F_2 = r_i[0] * r_i[0] / pow((a(0) * a(0) + lambda), 2) + r_i[1] * r_i[1] / pow((a(1) * a(1) + lambda), 2) + r_i[2] * r_i[2] / pow((a(2) * a(2) + lambda), 2);
		Lambda_l = F_1 / F_2;
	}

	I_ijk_l = A_ijk * Lambda_l;
	return I_ijk_l;
}



double  Potential3DElliptical::I_ijk_mn(int i, int j, int k, int p, int q, double* r_i, int location, double lambda, Ref<VectorXd> a)
{
	double A_ijk, Lambda_pq, Lambda_p, Lambda_q, Z_ijk, I_ijk_mn;
	if (location == 0) {
		A_ijk = 1;
		Lambda_pq = 0;
		Lambda_p = 0;
		Lambda_q = 0;
		Z_ijk = 1;
	}
	else {
		double delta = sqrt((a(0) * a(0) + lambda) * (a(1) * a(1) + lambda) * (a(2) * a(2) + lambda));
		A_ijk = -2 * Pi * a(0) * a(1) * a(2) / ((a(i) * a(i) + lambda) * (a(j) * a(j) + lambda) * (a(k) * a(k) + lambda) * delta);
		Z_ijk = 1.0 / (a(i) * a(i) + lambda) + 1.0 / (a(j) * a(j) + lambda) + 1.0 / (a(k) * a(k) + lambda) + 0.5 * (1.0 / (a(0) * a(0) + lambda) + 1.0 / (a(1) * a(1) + lambda) + 1.0 / (a(2) * a(2) + lambda));

		double 	F_p = 2 * r_i[p] / (a(p) * a(p) + lambda);
		double 	F_q = 2 * r_i[q] / (a(q) * a(q) + lambda);
		double C_1 = r_i[0] * r_i[0] / pow((a(0) * a(0) + lambda), 2) + r_i[1] * r_i[1] / pow((a(1) * a(1) + lambda), 2) + r_i[2] * r_i[2] / pow((a(2) * a(2) + lambda), 2);
		double C_2 = r_i[0] * r_i[0] / pow((a(0) * a(0) + lambda), 3) + r_i[1] * r_i[1] / pow((a(1) * a(1) + lambda), 3) + r_i[2] * r_i[2] / pow((a(2) * a(2) + lambda), 3);

		Lambda_p = F_p / C_1;
		Lambda_q = F_q / C_1;

		double F_pq = (2 * d[p][q] - F_p * Lambda_q) / (a(p) * a(p) + lambda);
		double C_1q = F_q / (a(q) * a(q) + lambda) - 2 * C_2 * Lambda_q;
		Lambda_pq = (F_pq - Lambda_p * C_1q) / C_1;
	}


	I_ijk_mn = A_ijk * (Lambda_pq - Lambda_p * Lambda_q * Z_ijk);

	return I_ijk_mn;
}


double  Potential3DElliptical::I_ijk_mnp(int i, int j, int k, int p, int q, int r, double* r_i, int location, double lambda, Ref<VectorXd> a)
{
	double I_ijk_mnp, A_ijk;
	if (location == 0) {
		I_ijk_mnp = 0;
	}

	else {
		double delta = sqrt((a(0) * a(0) + lambda) * (a(1) * a(1) + lambda) * (a(2) * a(2) + lambda));
		A_ijk = -2 * Pi * a(0) * a(1) * a(2) / ((a(i) * a(i) + lambda) * (a(j) * a(j) + lambda) * (a(k) * a(k) + lambda) * delta);
		double 	F_p = 2 * r_i[p] / (a(p) * a(p) + lambda);
		double 	F_q = 2 * r_i[q] / (a(q) * a(q) + lambda);
		double 	F_r = 2 * r_i[r] / (a(r) * a(r) + lambda);

		double C_1 = r_i[0] * r_i[0] / pow((a(0) * a(0) + lambda), 2) + r_i[1] * r_i[1] / pow((a(1) * a(1) + lambda), 2) + r_i[2] * r_i[2] / pow((a(2) * a(2) + lambda), 2);
		double Lambda_p = F_p / C_1;
		double Lambda_q = F_q / C_1;
		double Lambda_r = F_r / C_1;

		double C_2 = r_i[0] * r_i[0] / pow((a(0) * a(0) + lambda), 3) + r_i[1] * r_i[1] / pow((a(1) * a(1) + lambda), 3) + r_i[2] * r_i[2] / pow((a(2) * a(2) + lambda), 3);
		double F_pq = (2 * d[p][q] - F_p * Lambda_q) / (a(p) * a(p) + lambda);
		double F_qr = (2 * d[q][r] - F_q * Lambda_r) / (a(q) * a(q) + lambda);
		double F_pr = (2 * d[p][r] - F_p * Lambda_r) / (a(p) * a(p) + lambda);

		double C_1q = F_q / (a(q) * a(q) + lambda) - 2 * C_2 * Lambda_q;
		double C_1p = F_p / (a(p) * a(p) + lambda) - 2 * C_2 * Lambda_p;
		double C_1r = F_r / (a(r) * a(r) + lambda) - 2 * C_2 * Lambda_r;

		double Lambda_pq = (F_pq - Lambda_p * C_1q) / C_1;
		double Lambda_pr = (F_pr - Lambda_p * C_1r) / C_1;
		double Lambda_qr = (F_qr - Lambda_q * C_1r) / C_1;

		double C_3 = r_i[0] * r_i[0] / pow((a(0) * a(0) + lambda), 4) + r_i[1] * r_i[1] / pow((a(1) * a(1) + lambda), 4) + r_i[2] * r_i[2] / pow((a(2) * a(2) + lambda), 4);
		double F_p_qr = -(F_pq * Lambda_r + F_pr * Lambda_q + F_p * Lambda_qr) / (a(p) * a(p) + lambda);
		double C_2p = F_p / pow((a(p) * a(p) + lambda), 2) - 3 * C_3 * Lambda_p;
		double C_2q = F_q / pow((a(q) * a(q) + lambda), 2) - 3 * C_3 * Lambda_q;
		double C_2r = F_r / pow((a(r) * a(r) + lambda), 2) - 3 * C_3 * Lambda_r;


		double C_1pq = (F_pq - F_p * Lambda_q / (a(p) * a(p) + lambda)) / (a(p) * a(p) + lambda) - 2 * (C_2q * Lambda_p + C_2 * Lambda_pq);
		double C_1pr = (F_pr - F_p * Lambda_r / (a(p) * a(p) + lambda)) / (a(p) * a(p) + lambda) - 2 * (C_2r * Lambda_p + C_2 * Lambda_pr);
		double C_1qr = (F_qr - F_q * Lambda_r / (a(q) * a(q) + lambda)) / (a(q) * a(q) + lambda) - 2 * (C_2r * Lambda_q + C_2 * Lambda_qr);

		double Lambda_pqr = -(Lambda_pq * C_1r - F_p_qr + Lambda_pr * C_1q + Lambda_p * C_1qr) / C_1;

		double Z_1 = 1.0 / (a(i) * a(i) + lambda) + 1.0 / (a(j) * a(j) + lambda) + 1.0 / (a(k) * a(k) + lambda) + 0.5 * (1.0 / (a(0) * a(0) + lambda) + 1.0 / (a(1) * a(1) + lambda) + 1.0 / (a(2) * a(2) + lambda));

		double Z_2 = 1.0 / pow((a(i) * a(i) + lambda), 2) + 1.0 / pow((a(j) * a(j) + lambda), 2) + 1.0 / pow((a(k) * a(k) + lambda), 2) + 0.5 * (1.0 / pow((a(0) * a(0) + lambda), 2) + 1.0 / pow((a(1) * a(1) + lambda), 2) + 1.0 / pow((a(2) * a(2) + lambda), 2));


		I_ijk_mnp = A_ijk * (Lambda_pqr - (Lambda_pq * Lambda_r + Lambda_pr * Lambda_q + Lambda_qr * Lambda_p) * Z_1 + Lambda_p * Lambda_q * Lambda_r * (Z_2 + (Z_1) * (Z_1)));

	}

	return I_ijk_mnp;
}


double  Potential3DElliptical::V(double* r_i, int location, double lambda, Ref<VectorXd>a, int orient, int type)
{
	double VV = 0.5 * (I(lambda, a, orient, type) - (r_i[0] * r_i[0] * I_i(0, location, lambda, a, orient, type) + r_i[1] * r_i[1] * I_i(1, location, lambda, a, orient, type) + r_i[2] * r_i[2] * I_i(2, location, lambda, a, orient, type)));
	return VV;
}

double  Potential3DElliptical::V_0_i(int i, double* r_i, int location, double lambda, Ref<VectorXd>a, int orient, int type)
{
	double VV = -r_i[i] * I_i(i, location, lambda, a, orient, type);
	return VV;
}

double Potential3DElliptical::V_0_ij(int i, int j, double* r_i, int location, double lambda, Ref<VectorXd>a, int orient, int type)
{
	double VV = -d[i][j] * I_i(i, location, lambda, a, orient, type) - r_i[i] * I_i_j(i, j, r_i, location, lambda, a);
	return  VV;
}

double Potential3DElliptical::V_0_ijk(int i, int j, int k, double* r_i, int location, double  lambda, Ref<VectorXd>a, int orient, int type)
{
	double	VV = -d[i][j] * I_i_j(i, k, r_i, location, lambda, a) - d[i][k] * I_i_j(i, j, r_i, location, lambda, a) - r_i[i] * I_i_jk(i, j, k, r_i, location, lambda, a);
	return VV;
}

double Potential3DElliptical::V_0_ijkl(int i, int j, int k, int l, double* r_i, int location, double lambda, Ref<VectorXd>a, int orient, int type)
{
	double VV = -d[i][j] * I_i_jk(i, k, l, r_i, location, lambda, a) - d[i][k] * I_i_jk(i, j, l, r_i, location, lambda, a) - d[i][l] * I_i_jk(i, j, k, r_i, location, lambda, a) - r_i[i] * I_i_jkl(i, j, k, l, r_i, location, lambda, a);
	return VV;
}


double Potential3DElliptical::V_i(int i, double* r_i, int location, double lambda, Ref<VectorXd> a, int orient, int type)
{
	double	V_i = 0.5 * (I_i(i, location, lambda, a, orient, type) - (r_i[0] * r_i[0] * I_ij(0, i, location, lambda, a, orient, type) + r_i[1] * r_i[1] * I_ij(1, i, location, lambda, a, orient, type) + r_i[2] * r_i[2] * I_ij(2, i, location, lambda, a, orient, type)));
	return V_i;

}

double  Potential3DElliptical::V_i_j(int i, int p, double* r_i, int location, double lambda, Ref<VectorXd> a, int orient, int type)
{
	double V_i_j = -r_i[p] * I_ij(p, i, location, lambda, a, orient, type);
	return V_i_j;
}

double  Potential3DElliptical::V_i_jk(int i, int p, int q, double* r_i, int location, double lambda, Ref<VectorXd> a, int orient, int type)
{
	double	V_i_jk = -(d[p][q] * I_ij(p, i, location, lambda, a, orient, type) + r_i[p] * I_ij_k(p, i, q, r_i, location, lambda, a));

	return V_i_jk;
}

double  Potential3DElliptical::V_i_jkl(int i, int p, int q, int r, double* r_i, int location, double lambda, Ref<VectorXd>a, int orient, int type)
{
	double V_i_jkl = -(d[p][q] * I_ij_k(p, i, r, r_i, location, lambda, a) + d[p][r] * I_ij_k(p, i, q, r_i, location, lambda, a) + r_i[p] * I_ij_kl(p, i, q, r, r_i, location, lambda, a));
	return V_i_jkl;

}

double Potential3DElliptical::V_i_jkls(int i, int p, int q, int r, int t, double* r_i, int location, double lambda, Ref<VectorXd>a, int orient, int type)
{

	double	V_i_jkls = -(d[p][q] * I_ij_kl(p, i, r, t, r_i, location, lambda, a) + d[p][r] * I_ij_kl(p, i, q, t, r_i, location, lambda, a) + d[p][t] * I_ij_kl(p, i, q, r, r_i, location, lambda, a) + r_i[p] * I_ij_kls(p, i, q, r, t, r_i, location, lambda, a));
	return  V_i_jkls;
}

double Potential3DElliptical::V_ij(int i, int j, double* r_i, int location, double lambda, Ref<VectorXd> a, int orient, int type)
{
	double V_ij = 0.5 * (I_ij(i, j, location, lambda, a, orient, type) - (r_i[0] * r_i[0] * I_ijk(0, i, j, location, lambda, a, orient, type) + r_i[1] * r_i[1] * I_ijk(1, i, j, location, lambda, a, orient, type) + r_i[2] * r_i[2] * I_ijk(2, i, j, location, lambda, a, orient, type)));
	return V_ij;

}

double  Potential3DElliptical::V_ij_k(int i, int j, int p, double* r_i, int location, double lambda, Ref<VectorXd>a, int orient, int type)
{

	double V_ij_k = -r_i[p] * I_ijk(p, i, j, location, lambda, a, orient, type);
	return V_ij_k;
}

double  Potential3DElliptical::V_ij_kl(int i, int j, int p, int q, double* r_i, int location, double lambda, Ref<VectorXd>a, int orient, int type)
{

	double	V_ij_kl = -(d[p][q] * I_ijk(p, i, j, location, lambda, a, orient, type) + r_i[p] * I_ijk_l(p, i, j, q, r_i, location, lambda, a));
	return V_ij_kl;
}

double Potential3DElliptical::V_ij_klmn(int i, int j, int p, int q, int r, int t, double* r_i, int location, double lambda, Ref<VectorXd>a, int orient, int type)
{
	double	V_ij_klmn = -(d[p][q] * I_ijk_mn(p, i, j, r, t, r_i, location, lambda, a) + d[p][r] * I_ijk_mn(p, i, j, q, t, r_i, location, lambda, a)
		+ d[p][t] * I_ijk_mn(p, i, j, q, r, r_i, location, lambda, a) + r_i[p] * I_ijk_mnp(p, i, j, q, r, t, r_i, location, lambda, a));

	return V_ij_klmn;
}

double   Potential3DElliptical::V_ij_kls(int i, int j, int p, int q, int r, double* r_i, int location, double lambda, Ref<VectorXd> a, int orient, int type)
{
	double	V_ij_kls = -(d[p][q] * I_ijk_l(p, i, j, r, r_i, location, lambda, a) + d[p][r] * I_ijk_l(p, i, j, q, r_i, location, lambda, a) + r_i[p] * I_ijk_mn(p, i, j, q, r, r_i, location, lambda, a));
	return  V_ij_kls;
}


///////////////////////////
///////////////////////////


double   Potential3DElliptical::Phi_00p(double* r_p, Ref<VectorXd> a, int orient, int type)
{

	double rho_p = r_p[0] * r_p[0] / (a(0) * a(0)) + r_p[1] * r_p[1] / (a(1) * a(1)) + r_p[2] * r_p[2] / (a(2) * a(2));
	int location;
	double lambda;

	if (rho_p <= 1) {
		location = 0;
		lambda = 0;
	}
	else {
		location = 1;
		double absr_p = sqrt(r_p[0] * r_p[0] + r_p[1] * r_p[1] + r_p[2] * r_p[2]);
		double a_l = a.maxCoeff();
		double a_s = a.minCoeff();
		double a_u = a(orient);
		double a_e = a(dd[orient + 1]);

		lambda = 0.5 * (absr_p * absr_p - a_u * a_u - a_e * a_e + sqrt(pow((absr_p * absr_p + a_u * a_u - a_e * a_e), 2) - 4 * (a_u * a_u - a_e * a_e) * r_p[orient] * r_p[orient]));

	}


	double Phi = V(r_p, location, lambda, a, orient, type);
	return Phi;
}

double  Potential3DElliptical::Phi_01p(double* r_p, Ref<VectorXd> a, int orient, int i, int type)
{

	double rho_p = r_p[0] * r_p[0] / (a(0) * a(0)) + r_p[1] * r_p[1] / (a(1) * a(1)) + r_p[2] * r_p[2] / (a(2) * a(2));
	double a0 = a(0);
	double a1 = a(1);
	double a2 = a(2);
	int location;
	double lambda;
	if (rho_p <= 1) {
		location = 0;
		lambda = 0;
	}
	else {
		location = 1;
		double absr_p = sqrt(r_p[0] * r_p[0] + r_p[1] * r_p[1] + r_p[2] * r_p[2]);
		double a_l = a.maxCoeff();
		double a_s = a.minCoeff();
		double a_u = a(orient);
		double a_e = a(dd[orient + 1]);

		lambda = 0.5 * (absr_p * absr_p - a_u * a_u - a_e * a_e + sqrt(pow((absr_p * absr_p + a_u * a_u - a_e * a_e), 2) - 4 * (a_u * a_u - a_e * a_e) * r_p[orient] * r_p[orient]));

	}

	double Phi = V_0_i(i, r_p, location, lambda, a, orient, type);
	return Phi;
}

double  Potential3DElliptical::Phi_02p(double* r_i, Ref<VectorXd>a, int orient, int i, int j, int type)
{
	double rho_p = r_i[0] * r_i[0] / (a(0) * a(0)) + r_i[1] * r_i[1] / (a(1) * a(1)) + r_i[2] * r_i[2] / (a(2) * a(2));
	int location;
	double lambda;
	if (rho_p <= 1) {
		location = 0;
		lambda = 0;
	}
	else {
		location = 1;
		double absr_p = sqrt(r_i[0] * r_i[0] + r_i[1] * r_i[1] + r_i[2] * r_i[2]);
		double a_l = a.maxCoeff();
		double a_s = a.minCoeff();
		double a_u = a(orient);
		double a_e = a(dd[orient + 1]);

		lambda = 0.5 * (absr_p * absr_p - a_u * a_u - a_e * a_e + sqrt(pow((absr_p * absr_p + a_u * a_u - a_e * a_e), 2) - 4 * (a_u * a_u - a_e * a_e) * r_i[orient] * r_i[orient]));

	}

	double Phi = V_0_ij(i, j, r_i, location, lambda, a, orient, type);
	return Phi;
}

double Potential3DElliptical::Phi_03p(double* r_i, Ref<VectorXd>a, int orient, int i, int j, int  r, int type)
{
	double rho_p = r_i[0] * r_i[0] / (a(0) * a(0)) + r_i[1] * r_i[1] / (a(1) * a(1)) + r_i[2] * r_i[2] / (a(2) * a(2));
	int location;
	double lambda;
	if (rho_p <= 1) {
		location = 0;
		lambda = 0;
	}
	else {
		location = 1;
		double absr_p = sqrt(r_i[0] * r_i[0] + r_i[1] * r_i[1] + r_i[2] * r_i[2]);
		double a_l = a.maxCoeff();
		double a_s = a.minCoeff();
		double a_u = a(orient);
		double a_e = a(dd[orient + 1]);

		lambda = 0.5 * (absr_p * absr_p - a_u * a_u - a_e * a_e + sqrt(pow((absr_p * absr_p + a_u * a_u - a_e * a_e), 2) - 4 * (a_u * a_u - a_e * a_e) * r_i[orient] * r_i[orient]));

	}

	double Phi = V_0_ijk(i, j, r, r_i, location, lambda, a, orient, type);
	return Phi;
}

double Potential3DElliptical::Phi_11p(double* r_i, Ref<VectorXd>a, int orient, int n, int i, int type)
{
	double rho_p = r_i[0] * r_i[0] / (a(0) * a(0)) + r_i[1] * r_i[1] / (a(1) * a(1)) + r_i[2] * r_i[2] / (a(2) * a(2));
	int location;
	double lambda;
	if (rho_p <= 1) {
		location = 0;
		lambda = 0;
	}
	else {
		location = 1;
		double absr_p = sqrt(r_i[0] * r_i[0] + r_i[1] * r_i[1] + r_i[2] * r_i[2]);
		double a_l = a.maxCoeff();
		double a_s = a.minCoeff();
		double a_u = a(orient);
		double a_e = a(dd[orient + 1]);

		lambda = 0.5 * (absr_p * absr_p - a_u * a_u - a_e * a_e + sqrt(pow((absr_p * absr_p + a_u * a_u - a_e * a_e), 2) - 4 * (a_u * a_u - a_e * a_e) * r_i[orient] * r_i[orient]));

	}


	double Phi = a(n) * a(n) * (d[i][n] * V_i(n, r_i, location, lambda, a, orient, type) + r_i[n] * V_i_j(n, i, r_i, location, lambda, a, orient, type));
	return Phi;

}

double Potential3DElliptical::Phi_12p(double* r_i, Ref<VectorXd>a, int orient, int n, int i, int j, int type)
{

	double rho_p = r_i[0] * r_i[0] / (a(0) * a(0)) + r_i[1] * r_i[1] / (a(1) * a(1)) + r_i[2] * r_i[2] / (a(2) * a(2));
	int location;
	double lambda;
	if (rho_p <= 1) {
		location = 0;
		lambda = 0;
	}
	else {
		location = 1;
		double absr_p = sqrt(r_i[0] * r_i[0] + r_i[1] * r_i[1] + r_i[2] * r_i[2]);
		double a_l = a.maxCoeff();
		double a_s = a.minCoeff();
		double a_u = a(orient);
		double a_e = a(dd[orient + 1]);

		lambda = 0.5 * (absr_p * absr_p - a_u * a_u - a_e * a_e + sqrt(pow((absr_p * absr_p + a_u * a_u - a_e * a_e), 2) - 4 * (a_u * a_u - a_e * a_e) * r_i[orient] * r_i[orient]));

	}


	double Phi = a(n) * a(n) * (d[i][n] * V_i_j(n, j, r_i, location, lambda, a, orient, type) + d[j][n] * V_i_j(n, i, r_i, location, lambda, a, orient, type) + r_i[n] * V_i_jk(n, i, j, r_i, location, lambda, a, orient, type));

	return Phi;

}

double Potential3DElliptical::Phi_13p(double* r_i, Ref<VectorXd> a, int orient, int n, int i, int j, int k, int type)
{

	double rho_p = r_i[0] * r_i[0] / (a(0) * a(0)) + r_i[1] * r_i[1] / (a(1) * a(1)) + r_i[2] * r_i[2] / (a(2) * a(2));
	int location;
	double lambda;
	if (rho_p <= 1) {
		location = 0;
		lambda = 0;
	}
	else {
		location = 1;
		double absr_p = sqrt(r_i[0] * r_i[0] + r_i[1] * r_i[1] + r_i[2] * r_i[2]);
		double a_l = a.maxCoeff();
		double a_s = a.minCoeff();
		double a_u = a(orient);
		double a_e = a(dd[orient + 1]);

		lambda = 0.5 * (absr_p * absr_p - a_u * a_u - a_e * a_e + sqrt(pow((absr_p * absr_p + a_u * a_u - a_e * a_e), 2) - 4 * (a_u * a_u - a_e * a_e) * r_i[orient] * r_i[orient]));

	}



	double Phi = a(n) * a(n) * (d[i][n] * V_i_jk(n, j, k, r_i, location, lambda, a, orient, type) + d[j][n] * V_i_jk(n, i, k, r_i, location, lambda, a, orient, type) + d[n][k] * V_i_jk(n, i, j, r_i, location, lambda, a, orient, type)
		+ r_i[n] * V_i_jkl(n, i, j, k, r_i, location, lambda, a, orient, type));

	return Phi;
}

double Potential3DElliptical::Psi_02p(double* r_i, Ref<VectorXd>a, int orient, int i, int j, int type)
{
	double rho_p = r_i[0] * r_i[0] / (a(0) * a(0)) + r_i[1] * r_i[1] / (a(1) * a(1)) + r_i[2] * r_i[2] / (a(2) * a(2));
	int location;
	double lambda;
	if (rho_p <= 1) {
		location = 0;
		lambda = 0;
	}
	else {
		location = 1;
		double absr_p = sqrt(r_i[0] * r_i[0] + r_i[1] * r_i[1] + r_i[2] * r_i[2]);
		double a_l = a.maxCoeff();
		double a_s = a.minCoeff();
		double a_u = a(orient);
		double a_e = a(dd[orient + 1]);

		lambda = 0.5 * (absr_p * absr_p - a_u * a_u - a_e * a_e + sqrt(pow((absr_p * absr_p + a_u * a_u - a_e * a_e), 2) - 4 * (a_u * a_u - a_e * a_e) * r_i[orient] * r_i[orient]));

	}



	double Psi = d[i][j] * (V(r_i, location, lambda, a, orient, type) - a(i) * a(i) * V_i(i, r_i, location, lambda, a, orient, type))
		+ r_i[i] * (V_0_i(j, r_i, location, lambda, a, orient, type) - a(i) * a(i) * V_i_j(i, j, r_i, location, lambda, a, orient, type));
	return Psi;
}


double  Potential3DElliptical::Psi_03p(double* r_i, Ref<VectorXd>a, int  orient, int i, int j, int k, int type)
{

	double rho_p = r_i[0] * r_i[0] / (a(0) * a(0)) + r_i[1] * r_i[1] / (a(1) * a(1)) + r_i[2] * r_i[2] / (a(2) * a(2));
	int location;
	double lambda;
	if (rho_p <= 1) {
		location = 0;
		lambda = 0;
	}
	else {
		location = 1;
		double absr_p = sqrt(r_i[0] * r_i[0] + r_i[1] * r_i[1] + r_i[2] * r_i[2]);
		double a_l = a.maxCoeff();
		double a_s = a.minCoeff();
		double a_u = a(orient);
		double a_e = a(dd[orient + 1]);

		lambda = 0.5 * (absr_p * absr_p - a_u * a_u - a_e * a_e + sqrt(pow((absr_p * absr_p + a_u * a_u - a_e * a_e), 2) - 4 * (a_u * a_u - a_e * a_e) * r_i[orient] * r_i[orient]));

	}

	double Psi = d[i][j] * (V_0_i(k, r_i, location, lambda, a, orient, type) - a(i) * a(i) * V_i_j(i, k, r_i, location, lambda, a, orient, type)) + d[i][k] * (V_0_i(j, r_i, location, lambda, a, orient, type) - a(i) * a(i) * V_i_j(i, j, r_i, location, lambda, a, orient, type))
		+ r_i[i] * (V_0_ij(j, k, r_i, location, lambda, a, orient, type) - a(i) * a(i) * V_i_jk(i, j, k, r_i, location, lambda, a, orient, type));

	return Psi;
}


double Potential3DElliptical::Psi_04p(double* r_i, Ref<VectorXd> a, int orient, int i, int j, int k, int l, int type)
{

	double rho_p = r_i[0] * r_i[0] / (a(0) * a(0)) + r_i[1] * r_i[1] / (a(1) * a(1)) + r_i[2] * r_i[2] / (a(2) * a(2));
	int location;
	double lambda;
	if (rho_p <= 1) {
		location = 0;
		lambda = 0;
	}
	else {
		location = 1;
		double absr_p = sqrt(r_i[0] * r_i[0] + r_i[1] * r_i[1] + r_i[2] * r_i[2]);
		double a_l = a.maxCoeff();
		double a_s = a.minCoeff();
		double a_u = a(orient);
		double a_e = a(dd[orient + 1]);

		lambda = 0.5 * (absr_p * absr_p - a_u * a_u - a_e * a_e + sqrt(pow((absr_p * absr_p + a_u * a_u - a_e * a_e), 2) - 4 * (a_u * a_u - a_e * a_e) * r_i[orient] * r_i[orient]));

	}

	double Psi = d[i][j] * (V_0_ij(k, l, r_i, location, lambda, a, orient, type) - a(i) * a(i) * V_i_jk(i, k, l, r_i, location, lambda, a, orient, type)) + d[i][k] * (V_0_ij(j, l, r_i, location, lambda, a, orient, type) - a(i) * a(i) * V_i_jk(i, j, l, r_i, location, lambda, a, orient, type)) + d[i][l] * (V_0_ij(j, k, r_i, location, lambda, a, orient, type) - a(i) * a(i) * V_i_jk(i, j, k, r_i, location, lambda, a, orient, type))
		+ r_i[i] * (V_0_ijk(j, k, l, r_i, location, lambda, a, orient, type) - a(i) * a(i) * V_i_jkl(i, j, k, l, r_i, location, lambda, a, orient, type));
	return Psi;
}

double  Potential3DElliptical::Psi_05p(double* r_i, Ref<VectorXd>a, int orient, int i, int j, int  k, int l, int s, int type)
{

	double rho_p = r_i[0] * r_i[0] / (a(0) * a(0)) + r_i[1] * r_i[1] / (a(1) * a(1)) + r_i[2] * r_i[2] / (a(2) * a(2));
	int location;
	double lambda;
	if (rho_p <= 1) {
		location = 0;
		lambda = 0;
	}
	else {
		location = 1;
		double absr_p = sqrt(r_i[0] * r_i[0] + r_i[1] * r_i[1] + r_i[2] * r_i[2]);
		double a_l = a.maxCoeff();
		double a_s = a.minCoeff();
		double a_u = a(orient);
		double a_e = a(dd[orient + 1]);

		lambda = 0.5 * (absr_p * absr_p - a_u * a_u - a_e * a_e + sqrt(pow((absr_p * absr_p + a_u * a_u - a_e * a_e), 2) - 4 * (a_u * a_u - a_e * a_e) * r_i[orient] * r_i[orient]));

	}

	double Psi = d[i][j] * (V_0_ijk(k, l, s, r_i, location, lambda, a, orient, type) - a(i) * a(i) * V_i_jkl(i, k, l, s, r_i, location, lambda, a, orient, type)) + d[i][k] * (V_0_ijk(j, l, s, r_i, location, lambda, a, orient, type) - a(i) * a(i) * V_i_jkl(i, j, l, s, r_i, location, lambda, a, orient, type)) + d[i][l] * (V_0_ijk(j, k, s, r_i, location, lambda, a, orient, type) - a(i) * a(i) * V_i_jkl(i, j, k, s, r_i, location, lambda, a, orient, type))
		+ d[i][s] * (V_0_ijk(j, k, l, r_i, location, lambda, a, orient, type) - a(i) * a(i) * V_i_jkl(i, j, k, l, r_i, location, lambda, a, orient, type)) + r_i[i] * (V_0_ijkl(j, k, l, s, r_i, location, lambda, a, orient, type) - a(i) * a(i) * V_i_jkls(i, j, k, l, s, r_i, location, lambda, a, orient, type));
	return Psi;
}

/*
double Psi_12p(double*r_i, Ref<VectorXd>a, int orient, int n, int i, int j, int k, int l)
{
	double rho_p = r_i[0] * r_i[0] / a(0)*a(0) + r_i[1] * r_i[1] / a(1)*a(1) + r_i[2] * r_i[2] / a(2)*a(2);
	int location;
	double lambda;
	if (rho_p <= 1){
		location = 0;
		lambda = 0;
	}
	else{
		int location = 1;
		double absr_p = sqrt(r_i[0] * r_i[0] + r_i[1] * r_i[1] + r_i[2] * r_i[2]);
		double a_l = a.maxCoeff();
		double a_s = a.minCoeff();
		double lambda = 0.5 * (absr_p*absr_p - a_l*a_l - a_s*a_s + sqrt(pow((absr_p*absr_p + a_l*a_l - a_s*a_s), 2) - 4 * (a_l*a_l - a_s*a_s) * r_i[orient] * r_i[orient]));
	}

	double	Psi = -a(n)*a(n)* (d[n][i]*(d[k][l]*d[j][m] + d[j][l]*d[m][k])*(I_ij(i, k, location, lambda, a, orient) - a(n)*a(n) * I_ijk(j, k, n, location, lambda, a, orient))
		+ (d[i][j]*(d[k][l]*d[n][m] + d[n][l]*d[k][m]) + d[n][j]*(d[k][l]*d[i][m] + d[i][l])*d[k][m]))*(I_ij(i, k, location, lambda, a, orient)
		- a(n)*a(n) * I_ijk(n, i, k, location, lambda, a, orient)) + (d[n][k]*(d[i][l]*d[j][m] + d[j][l]*d[i][m]) + d[i][k]*(d[n][l]*d[j][m] + d[n][m]*d[j][l])
		+ d[j][k]*(d[n][m]*d[i][l] + d[n][l]*d[i][m]))*(I_ij(i, j, location, lambda, a, orient) - a(n)*a(n)* I_ijk(n, i, j, location, lambda, a, orient))
		+ d[l][m]*(d[i][j]*d[n][k] + d[n][j]*d[i][k])*(I_ij(i, l, location, lambda, a, orient) - a(n)*a(n) * I_ijk(n, i, l, location, lambda, a, orient))
		+ d[n][i]*d[j][k]*d[l][m]*(I_ij(j, l, location, lambda, a, orient) - a(n)*a(n) * I_ijk(n, j, l, location, lambda, a, orient)))*r_i[m]
		- a(n)*a(n) * r_i[k]*r_i[j]*d[n][i]*(I_ij_k(j, k, l, r_i, location, lambda, a) - a(n)*a(n) * I_ijk_l(n, j, k, l, r_i, location, lambda, a))
		- a(n)*a(n) * r_i[k]*(d[n][j]*r_i[i] + d[i][j]*r_i[n])*(I_ij_k(i, k, l, r_i, location, lambda, a) - a(n)*a(n) * I_ijk_l(n, i, k, l, r_i, location, lambda, a))
		- a(n)*a(n) * (d[n][k]*r_i[i]*r_i[j] + d[i][k]*r_i[n]*r_i[j] + d[k][j]*r_i[n]*r_i[i])*(I_ij_k(i, j, l, r_i, location, lambda, a) - a(n)*a(n) * I_ijk_l(n, i, j, l, r_i, location, lambda, a))
		- a(n)*a(n) * (d[n][l]*r_i[i]*r_i[j] + d[i][l]*r_i[j]*r_i[n] + d[j][n]*r_i[n]*r_i[i])*(I_ij_k(i, j, k, r_i, location, lambda, a) - a(n)*a(n) * I_ijk_l(n, i, j, k, r_i, location, lambda, a))
		- a(n)*a(n) * r_i[i]*r_i[j]*r_i[n]*(I_ij_kl(j, i, k, l, r_i, location, lambda, a) - a(n)*a(n) * I_ijk_mn(i, j, n, k, l, r_i, location, lambda, a));
	return Psi;
}
*/


double Potential3DElliptical::Psi_13p(double* r_i, Ref<VectorXd>a, int orient, int n, int i, int j, int k, int type)
{
	double rho_p = r_i[0] * r_i[0] / (a(0) * a(0)) + r_i[1] * r_i[1] / (a(1) * a(1)) + r_i[2] * r_i[2] / (a(2) * a(2));
	int location;
	double lambda;
	if (rho_p <= 1) {
		location = 0;
		lambda = 0;
	}
	else {
		location = 1;
		double absr_p = sqrt(r_i[0] * r_i[0] + r_i[1] * r_i[1] + r_i[2] * r_i[2]);
		double a_l = a.maxCoeff();
		double a_s = a.minCoeff();
		double a_u = a(orient);
		double a_e = a(dd[orient + 1]);

		lambda = 0.5 * (absr_p * absr_p - a_u * a_u - a_e * a_e + sqrt(pow((absr_p * absr_p + a_u * a_u - a_e * a_e), 2) - 4 * (a_u * a_u - a_e * a_e) * r_i[orient] * r_i[orient]));

	}

	double Psi = a(n) * a(n) * (d[i][n] * d[j][k] * (V_i(j, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij(j, n, r_i, location, lambda, a, orient, type)) + d[i][n] * r_i[j] * (V_i_j(j, k, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_k(j, n, k, r_i, location, lambda, a, orient, type))
		+ (d[n][j] * d[i][k] + d[i][j] * d[n][k]) * (V_i(i, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij(i, n, r_i, location, lambda, a, orient, type))
		+ (d[n][j] * r_i[i] + d[i][j] * r_i[n]) * (V_i_j(i, k, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_k(i, n, k, r_i, location, lambda, a, orient, type))
		+ (d[n][k] * r_i[i] + d[i][k] * r_i[n]) * (V_i_j(i, j, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_k(i, n, j, r_i, location, lambda, a, orient, type))
		+ r_i[n] * r_i[i] * (V_i_jk(i, j, k, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_kl(i, n, j, k, r_i, location, lambda, a, orient, type)));
	return Psi;
}


double Potential3DElliptical::Psi_14p(double* r_i, Ref<VectorXd>a, int orient, int n, int i, int j, int k, int l, int type)

{
	double rho_p = r_i[0] * r_i[0] / (a(0) * a(0)) + r_i[1] * r_i[1] / (a(1) * a(1)) + r_i[2] * r_i[2] / (a(2) * a(2));
	int location;
	double lambda;
	if (rho_p <= 1) {
		location = 0;
		lambda = 0;
	}
	else {
		location = 1;
		double absr_p = sqrt(r_i[0] * r_i[0] + r_i[1] * r_i[1] + r_i[2] * r_i[2]);
		double a_l = a.maxCoeff();
		double a_s = a.minCoeff();
		double a_u = a(orient);
		double a_e = a(dd[orient + 1]);

		lambda = 0.5 * (absr_p * absr_p - a_u * a_u - a_e * a_e + sqrt(pow((absr_p * absr_p + a_u * a_u - a_e * a_e), 2) - 4 * (a_u * a_u - a_e * a_e) * r_i[orient] * r_i[orient]));

	}

	double Psi = a(n) * a(n) * (d[i][n] * d[j][k] * (V_i_j(j, l, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_k(j, n, l, r_i, location, lambda, a, orient, type)) + d[i][n] * d[j][l] * (V_i_j(j, k, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_k(j, n, k, r_i, location, lambda, a, orient, type))
		+ d[i][n] * r_i[j] * (V_i_jk(j, k, l, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_kl(j, n, k, l, r_i, location, lambda, a, orient, type)) + (d[n][j] * d[i][k] + d[i][j] * d[n][k]) * (V_i_j(i, l, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_k(i, n, l, r_i, location, lambda, a, orient, type))
		+ (d[n][j] * d[i][l] + d[i][j] * d[n][l]) * (V_i_j(i, k, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_k(i, n, k, r_i, location, lambda, a, orient, type)) + (d[n][j] * r_i[i] + d[i][j] * r_i[n]) * (V_i_jk(i, k, l, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_kl(i, n, k, l, r_i, location, lambda, a, orient, type))
		+ (d[n][k] * d[i][l] + d[i][k] * d[n][l]) * (V_i_j(i, j, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_k(i, n, j, r_i, location, lambda, a, orient, type))
		+ (d[n][k] * r_i[i] + d[i][k] * r_i[n]) * (V_i_jk(i, j, l, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_kl(i, n, j, l, r_i, location, lambda, a, orient, type))
		+ (d[n][l] * r_i[i] + d[i][l] * r_i[n]) * (V_i_jk(i, j, k, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_kl(i, n, j, k, r_i, location, lambda, a, orient, type)) + r_i[n] * r_i[i] * (V_i_jkl(i, j, k, l, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_kls(i, n, j, k, l, r_i, location, lambda, a, orient, type)));
	return Psi;
}


double Potential3DElliptical::Psi_15p(double* r_i, Ref<VectorXd>a, int orient, int n, int i, int j, int k, int l, int s, int type)
{
	double rho_p = r_i[0] * r_i[0] / (a(0) * a(0)) + r_i[1] * r_i[1] / (a(1) * a(1)) + r_i[2] * r_i[2] / (a(2) * a(2));
	int location;
	double lambda;
	if (rho_p <= 1) {
		location = 0;
		lambda = 0;
	}
	else {
		location = 1;
		double absr_p = sqrt(r_i[0] * r_i[0] + r_i[1] * r_i[1] + r_i[2] * r_i[2]);
		double a_l = a.maxCoeff();
		double a_s = a.minCoeff();
		double a_u = a(orient);
		double a_e = a(dd[orient + 1]);

		lambda = 0.5 * (absr_p * absr_p - a_u * a_u - a_e * a_e + sqrt(pow((absr_p * absr_p + a_u * a_u - a_e * a_e), 2) - 4 * (a_u * a_u - a_e * a_e) * r_i[orient] * r_i[orient]));

	}



	double Psi = a(n) * a(n) * (d[i][n] * d[j][k] * (V_i_jk(j, l, s, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_kl(j, n, l, s, r_i, location, lambda, a, orient, type))
		+ d[i][n] * d[j][l] * (V_i_jk(j, k, s, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_kl(j, n, k, s, r_i, location, lambda, a, orient, type))
		+ d[i][n] * d[j][s] * (V_i_jk(j, k, l, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_kl(j, n, k, l, r_i, location, lambda, a, orient, type))
		+ d[i][n] * r_i[j] * (V_i_jkl(j, k, l, s, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_kls(j, n, k, l, s, r_i, location, lambda, a, orient, type))
		+ (d[n][j] * d[i][k] + d[i][j] * d[n][k]) * (V_i_jk(i, l, s, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_kl(i, n, l, s, r_i, location, lambda, a, orient, type))
		+ (d[n][j] * d[i][l] + d[i][j] * d[n][l]) * (V_i_jk(i, k, s, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_kl(i, n, k, s, r_i, location, lambda, a, orient, type))
		+ (d[n][j] * d[i][s] + d[i][j] * d[n][s]) * (V_i_jk(i, k, l, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_kl(i, n, k, l, r_i, location, lambda, a, orient, type))
		+ (d[n][j] * r_i[i] + d[i][j] * r_i[n]) * (V_i_jkl(i, k, l, s, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_kls(i, n, k, l, s, r_i, location, lambda, a, orient, type))
		+ (d[n][k] * d[i][l] + d[i][k] * d[n][l]) * (V_i_jk(i, j, s, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_kl(i, n, j, s, r_i, location, lambda, a, orient, type))
		+ (d[n][k] * d[i][s] + d[i][k] * d[n][s]) * (V_i_jk(i, j, l, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_kl(i, n, j, l, r_i, location, lambda, a, orient, type))
		+ (d[n][k] * r_i[i] + d[i][k] * r_i[n]) * (V_i_jkl(i, j, l, s, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_kls(i, n, j, l, s, r_i, location, lambda, a, orient, type))
		+ (d[n][l] * d[i][s] + d[i][l] * d[n][s]) * (V_i_jk(i, j, k, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_kl(i, n, j, k, r_i, location, lambda, a, orient, type))
		+ (d[n][l] * r_i[i] + d[i][l] * r_i[n]) * (V_i_jkl(i, j, k, s, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_kls(i, n, j, k, s, r_i, location, lambda, a, orient, type))
		+ (d[n][s] * r_i[i] + d[i][s] * r_i[n]) * (V_i_jkl(i, j, k, l, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_kls(i, n, j, k, l, r_i, location, lambda, a, orient, type))
		+ r_i[n] * r_i[i] * (V_i_jkls(i, j, k, l, s, r_i, location, lambda, a, orient, type) - a(n) * a(n) * V_ij_klmn(i, n, j, k, l, s, r_i, location, lambda, a, orient, type)));
	return Psi;
}

double Potential3DElliptical::IGreen_01(double* x, Ref<VectorXd> a, int orient, int i, int type) {

	double IGreen_ij_l = -1.0 / (4.0 * Pi * k_0) * (Phi_01p(x, a, orient, i, type));

	return IGreen_ij_l;
}

double Potential3DElliptical::IGreen_11(double* x, Ref<VectorXd> a, int orient, int p, int l, int type) {

	double IGreen_ijp_l = -1.0 / (4.0 * Pi * k_0) * (Phi_11p(x, a, orient, p, l, type));

	return IGreen_ijp_l;
}

double Potential3DElliptical::Df_40(int i, int m, int s, int h, Ref<MatrixXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p, Ref<VectorXi> orient, Ref<MatrixXd> bb, Ref<VectorXi> type)
{
	double x1[3];
	double x1p[3];

	x1[0] = x_p(s, 0) - x_o(h, 0);  x1p[0] = 0.0;
	x1[1] = x_p(s, 1) - x_o(h, 1);  x1p[1] = 0.0;
	x1[2] = x_p(s, 2) - x_o(h, 2);  x1p[2] = 0.0;

	for (int ii = 0; ii < 3; ii++) {
		for (int jj = 0; jj < 3; jj++) {
			x1p[ii] = x1p[ii] + bb(ii, jj) * x1[jj];
		}
	}

	VectorXd a = radius.row(h);
	int orientp = orient(h);
	int typep = type(h);

	double D_4 = -(-k_0) * (-1.0 / (4 * Pi * k_0) * (Phi_02p(x1p, a, orientp, i, m, typep)));

	//double D_4 = 1 / (8 * pi*(1 - nu0))*(Psi_14(x1, a, p, k, l, i, j) - 2 * nu0*d[k][l] * Phi_12(x1, a, p, i, j) - (1 - nu0)*(Phi_12(x1, a, p, k, j)*d[i][l] + Phi_12(x1, a, p, k, i)*d[j][l] + Phi_12(x1, a, p, l, j)*d[i][k] + Phi_12(x1, a, p, l, i)*d[j][k]));

	return D_4;
}

double Potential3DElliptical::Df_50(int i, int m, int p, int s, int h, Ref<MatrixXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p, Ref<VectorXi> orient, Ref<MatrixXd> bb, Ref<VectorXi> type)
{
	double x1[3];
	double x1p[3];

	x1[0] = x_p(s, 0) - x_o(h, 0); x1p[0] = 0.0;
	x1[1] = x_p(s, 1) - x_o(h, 1); x1p[1] = 0.0;
	x1[2] = x_p(s, 2) - x_o(h, 2); x1p[2] = 0.0;

	for (int ii = 0; ii < 3; ii++) {
		for (int jj = 0; jj < 3; jj++) {
			x1p[ii] = x1p[ii] + bb(ii, jj) * x1[jj];
		}
	}

	VectorXd a = radius.row(h);
	int orientp = orient(h);
	int typep = type(h);

	double D_5 = -(-k_0) * (-1.0 / (4 * Pi * k_0) * (Phi_12p(x1p, a, orientp, p, i, m, typep)));

	return D_5;

}

double Potential3DElliptical::Df_41(int i, int m, int r, int s, int h, Ref<MatrixXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p, Ref<VectorXi> orient, Ref<MatrixXd>bb, Ref<VectorXi> type)
{
	double x1[3];
	double x1p[3];

	x1[0] = x_p(s, 0) - x_o(h, 0); x1p[0] = 0.0;
	x1[1] = x_p(s, 1) - x_o(h, 1); x1p[1] = 0.0;
	x1[2] = x_p(s, 2) - x_o(h, 2); x1p[2] = 0.0;

	for (int ii = 0; ii < 3; ii++) {
		for (int jj = 0; jj < 3; jj++) {
			x1p[ii] = x1p[ii] + bb(ii, jj) * x1[jj];
		}
	}

	VectorXd a = radius.row(h);
	int orientp = orient(h);
	int typep = type(h);

	double D_4_r = -(-k_0) * (-1.0 / (4 * Pi * k_0) * (Phi_03p(x1p, a, orientp, i, m, r, typep)));

	return D_4_r;
}

double Potential3DElliptical::Df_51(int i, int m, int p, int r, int s, int h, Ref<MatrixXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p, Ref<VectorXi> orient, Ref<MatrixXd>bb, Ref<VectorXi> type)
{
	double x1[3];
	double x1p[3];

	x1[0] = x_p(s, 0) - x_o(h, 0); x1p[0] = 0.0;
	x1[1] = x_p(s, 1) - x_o(h, 1); x1p[1] = 0.0;
	x1[2] = x_p(s, 2) - x_o(h, 2); x1p[2] = 0.0;

	for (int ii = 0; ii < 3; ii++) {
		for (int jj = 0; jj < 3; jj++) {
			x1p[ii] = x1p[ii] + bb(ii, jj) * x1[jj];
		}
	}

	VectorXd a = radius.row(h);
	int orientp = orient(h);
	int typep = type(h);

	double D_5_r = -(-k_0) * (-1.0 / (4 * Pi * k_0) * (Phi_13p(x1p, a, orientp, p, i, m, r, typep)));

	return D_5_r;
}

double Potential3DElliptical::I_40(int i, int j, int m, int n)
{
	//int chi = (m+3)/3-1;
	//m = m - 3*chi;
	double I_ijmn = d[i][m] * d[j][n];
	return I_ijmn;
}

double  Potential3DElliptical::Lf_00(int s, int h, Ref<MatrixXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p, Ref<MatrixXd> aa, Ref<VectorXi> orient, Ref<VectorXi> type)
{
	double* x1, L_ik, * x1p;

	x1 = new double[3];
	x1p = new double[3];


	double x = x_p(s, 2);
	x1[0] = x_p(s, 0) - x_o(h, 0);
	x1[1] = x_p(s, 1) - x_o(h, 1);
	x1[2] = x_p(s, 2) - x_o(h, 2);
	x1p[0] = 0.0;
	x1p[1] = 0.0;
	x1p[2] = 0.0;


	for (int ii = 0; ii < 3; ii++) {
		for (int jj = 0; jj < 3; jj++) {
			x1p[ii] = x1p[ii] + aa(ii, jj) * x1[jj];
		}
	}

	VectorXd a = radius.row(h);
	int orientp = orient(h);
	int typep = type(h);




	L_ik = -1.0 / (4.0 * Pi * k_0) * (Phi_00p(x1p, a, orientp, typep));


	delete[] x1;
	delete[] x1p;

	return L_ik;




}

double  Potential3DElliptical::Mf_30(int m, int s, int h, Ref<MatrixXd>radius, Ref<MatrixXd>x_o, Ref<MatrixXd> x_p, Ref<MatrixXd> aa, Ref<VectorXi> orient, Ref<VectorXi> type)
{
	double* x1, * x1p;
	x1 = new double[3];
	x1p = new double[3];

	double x = x_p(s, 2);
	x1[0] = x_p(s, 0) - x_o(h, 0);
	x1[1] = x_p(s, 1) - x_o(h, 1);
	x1[2] = x_p(s, 2) - x_o(h, 2);

	x1p[0] = 0.0;
	x1p[1] = 0.0;
	x1p[2] = 0.0;

	for (int ii = 0; ii < 3; ii++) {
		for (int jj = 0; jj < 3; jj++) {
			x1p[ii] = x1p[ii] + aa(ii, jj) * x1[jj];
		}
	}

	VectorXd a = radius.row(h);
	int orientp = orient(h);
	int typep = type(h);

	double M_imn = -(-k_0) * (-1.0 / (4.0 * Pi * k_0)) * (Phi_01p(x1p, a, orientp, m, typep));


	delete[] x1;
	delete[] x1p;

	return M_imn;


}

double  Potential3DElliptical::Nf_40(int m, int p, int s, int h, Ref<MatrixXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p, Ref<MatrixXd>aa, Ref<VectorXi> orient, Ref<VectorXi> type)
{
	double* x1, * x1p;
	x1 = new double[3];
	x1p = new double[3];

	double x = x_p(s, 2);
	x1[0] = x_p(s, 0) - x_o(h, 0);
	x1[1] = x_p(s, 1) - x_o(h, 1);
	x1[2] = x_p(s, 2) - x_o(h, 2);

	x1p[0] = 0.0;
	x1p[1] = 0.0;
	x1p[2] = 0.0;

	for (int ii = 0; ii < 3; ii++) {
		for (int jj = 0; jj < 3; jj++) {
			x1p[ii] = x1p[ii] + aa(ii, jj) * x1[jj];
		}
	}

	VectorXd a = radius.row(h);
	int orientp = orient(h);
	int typep = type(h);

	double N_imnp = -(-k_0) * (-1.0 / (4.0 * Pi * k_0)) * (Phi_11p(x1p, a, orientp, p, m, typep));


	delete[] x1;
	delete[] x1p;

	return N_imnp;


}