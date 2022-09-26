#pragma once
#include "GlobalConstant.h"
#include "Eigen/Dense"

using namespace Eigen;

namespace Potential3DElliptical {
	double I(double lambda, Ref<VectorXd> a, int orient, int type);


	double   I_i(int i, int location, double lambda, Ref<VectorXd> a, int orient, int type);


	double  I_i_j(int i, int j, double* r_i, int location, double lambda, Ref<VectorXd> a);


	double   I_i_jk(int i, int j, int k, double* r_i, int location, double lambda, Ref<VectorXd> a);


	double  I_i_jkl(int i, int p, int q, int r, double* r_i, int location, double lambda, Ref<VectorXd>a);



	double  I_ij(int i, int j, int location, double lambda, Ref<VectorXd> a, int orient, int type);


	double  I_ij_k(int i, int j, int k, double* r_i, int location, double lambda, Ref<VectorXd> a);

	double   I_ij_kl(int i, int j, int k, int l, double* r_i, int location, double lambda, Ref<VectorXd> a);


	double  I_ij_kls(int i, int j, int p, int q, int r, double* r_i, int location, double lambda, Ref<VectorXd> a);

	double  I_ijk(int i, int j, int k, int location, double lambda, Ref<VectorXd>a, int orient, int type);

	double  I_ijk_l(int i, int j, int k, int l, double* r_i, int location, double lambda, Ref<VectorXd> a);


	double   I_ijk_mn(int i, int j, int k, int p, int q, double* r_i, int location, double lambda, Ref<VectorXd> a);


	double   I_ijk_mnp(int i, int j, int k, int p, int q, int r, double* r_i, int location, double lambda, Ref<VectorXd> a);


	double   V(double* r_i, int location, double lambda, Ref<VectorXd>a, int orient, int type);

	double   V_0_i(int i, double* r_i, int location, double lambda, Ref<VectorXd>a, int orient, int type);

	double  V_0_ij(int i, int j, double* r_i, int location, double lambda, Ref<VectorXd>a, int orient, int type);

	double  V_0_ijk(int i, int j, int k, double* r_i, int location, double  lambda, Ref<VectorXd>a, int orient, int type);

	double  V_0_ijkl(int i, int j, int k, int l, double* r_i, int location, double lambda, Ref<VectorXd>a, int orient, int type);


	double  V_i(int i, double* r_i, int location, double lambda, Ref<VectorXd> a, int orient, int type);

	double   V_i_j(int i, int p, double* r_i, int location, double lambda, Ref<VectorXd> a, int orient, int type);

	double   V_i_jk(int i, int p, int q, double* r_i, int location, double lambda, Ref<VectorXd> a, int orient, int type);

	double   V_i_jkl(int i, int p, int q, int r, double* r_i, int location, double lambda, Ref<VectorXd>a, int orient, int type);

	double  V_i_jkls(int i, int p, int q, int r, int t, double* r_i, int location, double lambda, Ref<VectorXd>a, int orient, int type);

	double  V_ij(int i, int j, double* r_i, int location, double lambda, Ref<VectorXd> a, int orient, int type);

	double   V_ij_k(int i, int j, int p, double* r_i, int location, double lambda, Ref<VectorXd>a, int orient, int type);

	double   V_ij_kl(int i, int j, int p, int q, double* r_i, int location, double lambda, Ref<VectorXd>a, int orient, int type);

	double  V_ij_klmn(int i, int j, int p, int q, int r, int t, double* r_i, int location, double lambda, Ref<VectorXd>a, int orient, int type);

	double    V_ij_kls(int i, int j, int p, int q, int r, double* r_i, int location, double lambda, Ref<VectorXd> a, int orient, int type);


	///////////////////////////
	///////////////////////////


	double    Phi_00p(double* r_p, Ref<VectorXd> a, int orient, int type);

	double   Phi_01p(double* r_p, Ref<VectorXd> a, int orient, int i, int type);

	double   Phi_02p(double* r_i, Ref<VectorXd>a, int orient, int i, int j, int type);

	double  Phi_03p(double* r_i, Ref<VectorXd>a, int orient, int i, int j, int  r, int type);

	double  Phi_11p(double* r_i, Ref<VectorXd>a, int orient, int n, int i, int type);

	double  Phi_12p(double* r_i, Ref<VectorXd>a, int orient, int n, int i, int j, int type);

	double  Phi_13p(double* r_i, Ref<VectorXd> a, int orient, int n, int i, int j, int k, int type);

	double  Psi_02p(double* r_i, Ref<VectorXd>a, int orient, int i, int j, int type);


	double   Psi_03p(double* r_i, Ref<VectorXd>a, int  orient, int i, int j, int k, int type);


	double  Psi_04p(double* r_i, Ref<VectorXd> a, int orient, int i, int j, int k, int l, int type);

	double   Psi_05p(double* r_i, Ref<VectorXd>a, int orient, int i, int j, int  k, int l, int s, int type);

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



	double  Psi_13p(double* r_i, Ref<VectorXd>a, int orient, int n, int i, int j, int k, int type);


	double  Psi_14p(double* r_i, Ref<VectorXd>a, int orient, int n, int i, int j, int k, int l, int type);


	double  Psi_15p(double* r_i, Ref<VectorXd>a, int orient, int n, int i, int j, int k, int l, int s, int type);

	double IGreen_01(double* x, Ref<VectorXd> a, int orient, int i, int type);

	double IGreen_11(double* x, Ref<VectorXd> a, int orient, int p, int l, int type);

	double Df_40(int i, int m, int s, int h, Ref<MatrixXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p, Ref<VectorXi> orient, Ref<MatrixXd> bb, Ref<VectorXi> type);

	double Df_50(int i, int m, int p, int s, int h, Ref<MatrixXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p, Ref<VectorXi> orient, Ref<MatrixXd> bb, Ref<VectorXi> type);

	double Df_41(int i, int m, int r, int s, int h, Ref<MatrixXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p, Ref<VectorXi> orient, Ref<MatrixXd>bb, Ref<VectorXi> type);

	double Df_51(int i, int m, int p, int r, int s, int h, Ref<MatrixXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p, Ref<VectorXi> orient, Ref<MatrixXd>bb, Ref<VectorXi> type);

	double I_40(int i, int j, int m, int n);

	double Lf_00(int s, int h, Ref<MatrixXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p, Ref<MatrixXd> aa, Ref<VectorXi> orient, Ref<VectorXi> type);

	double Mf_30(int m, int s, int h, Ref<MatrixXd>radius, Ref<MatrixXd>x_o, Ref<MatrixXd> x_p, Ref<MatrixXd> aa, Ref<VectorXi> orient, Ref<VectorXi> type);

	double Nf_40(int m, int p, int s, int h, Ref<MatrixXd> radius, Ref<MatrixXd> x_o, Ref<MatrixXd> x_p, Ref<MatrixXd>aa, Ref<VectorXi> orient, Ref<VectorXi> type);
}