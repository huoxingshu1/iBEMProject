# include "BEMPotentialHelper3D.h"
# include "Eigen/Dense"
using namespace BEMPotential3D;

void BEMPotential3D::Triangle_Coord(int ntr, int l_pi, int** cor_tri)
{
	int n, temp;

	// double *xsii , *etai

	int nod_tri[3];
	int xsii[8];
	int etai[8];

	if (ntr == 0) {

		if (l_pi == 0) {

			nod_tri[0] = 1;
			nod_tri[1] = 2;
			nod_tri[2] = 0;
		}

		if (l_pi == 1) {

			nod_tri[0] = 2;
			nod_tri[1] = 3;
			nod_tri[2] = 1;
		}

		if (l_pi == 2) {

			nod_tri[0] = 0;
			nod_tri[1] = 1;
			nod_tri[2] = 2;
		}

		if (l_pi == 3) {

			nod_tri[0] = 0;
			nod_tri[1] = 1;
			nod_tri[2] = 3;
		}
	}



	if (ntr == 1) {

		if (l_pi == 0) {

			nod_tri[0] = 2;
			nod_tri[1] = 3;
			nod_tri[2] = 0;
		}

		if (l_pi == 1) {

			nod_tri[0] = 3;
			nod_tri[1] = 0;
			nod_tri[2] = 1;
		}

		if (l_pi == 2) {

			nod_tri[0] = 3;
			nod_tri[1] = 0;
			nod_tri[2] = 2;
		}

		if (l_pi == 3) {

			nod_tri[0] = 1;
			nod_tri[1] = 2;
			nod_tri[2] = 3;
		}
	}


	xsii[0] = -1;  xsii[1] = 1;  xsii[2] = 1; xsii[3] = -1;  xsii[4] = 0;  xsii[5] = 1; xsii[6] = 0; xsii[7] = -1;
	etai[0] = -1;  etai[1] = -1;  etai[2] = 1; etai[3] = 1;  etai[4] = -1;  etai[5] = 0; etai[6] = 1; etai[7] = 0;


	for (n = 0; n < 3; n++) {
		temp = nod_tri[n];
		cor_tri[0][n] = xsii[temp];
		cor_tri[1][n] = etai[temp];
	}


}

void BEMPotential3D::Trans_Tri(int ntr, int l_pi, double xsib, double etab, double& xsi, double& eta, double& Jacb)
{

	int i, j, k, n;

	double Nb[3];
	double dNbdxb[3];
	double dNbdeb[3];


	int** cor_tri;

	double Jacobm[2][2];
	double Gradiant[2][3];
	cor_tri = new int* [2];




	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
			Jacobm[i][j] = 0.0;


	for (i = 0; i < 2; i++)
		cor_tri[i] = new int[3];



	Triangle_Coord(ntr, l_pi, cor_tri);

	xsi = 0;
	eta = 0;

	Nb[0] = 0.25 * (1 + xsib) * (1 - etab);
	Nb[1] = 0.25 * (1 + xsib) * (1 + etab);
	Nb[2] = 0.5 * (1 - xsib);

	for (n = 0; n < 3; n++) {
		xsi = xsi + Nb[n] * cor_tri[0][n];
		eta = eta + Nb[n] * cor_tri[1][n];
	}


	dNbdxb[0] = 0.25 * (1 - etab);
	dNbdxb[1] = 0.25 * (1 + etab);
	dNbdxb[2] = -0.5;

	dNbdeb[0] = -0.25 * (1 + xsib);
	dNbdeb[1] = 0.25 * (1 + xsib);
	dNbdeb[2] = 0;

	Gradiant[0][0] = dNbdxb[0];
	Gradiant[0][1] = dNbdxb[1];
	Gradiant[0][2] = dNbdxb[2];
	Gradiant[1][0] = dNbdeb[0];
	Gradiant[1][1] = dNbdeb[1];
	Gradiant[1][2] = dNbdeb[2];

	for (i = 0; i < 2; i++) {
		for (j = 0; j < 2; j++) {
			for (k = 0; k < 3; k++) {
				Jacobm[i][j] = Jacobm[i][j] + Gradiant[i][k] * cor_tri[j][k];
			}
		}
	}

	Jacb = Jacobm[0][0] * Jacobm[1][1] - Jacobm[0][1] * Jacobm[1][0];



}

void BEMPotential3D::Serendip_func(double xsi, double eta, double& N1, double& N2, double& N3, double& N4, double* N, double** Ni)

{

	for (int i = 0; i < 1; i++) {
		for (int j = 0; j < 4; j++) {
			Ni[i][j] = 0.0;
		}
	}

	N1 = 0.25 * (1 - xsi) * (1 - eta);
	N2 = 0.25 * (1 + xsi) * (1 - eta);
	N3 = 0.25 * (1 + xsi) * (1 + eta);
	N4 = 0.25 * (1 - xsi) * (1 + eta);

	N[0] = N1;
	N[1] = N2;
	N[2] = N3;
	N[3] = N4;

	Ni[0][0] = N1;
	Ni[0][1] = N2;
	Ni[0][2] = N3;
	Ni[0][3] = N4;

}

void BEMPotential3D::Normal_Jac(double xsi, double eta, double** Elcor, double& Jac, double* Vnorm)

{
	double dN1_dxsi, dN2_dxsi, dN3_dxsi, dN4_dxsi, dN1_deta, dN2_deta, dN3_deta, dN4_deta;

	double V_xsi[3];
	double V_eta[3];


	int i;


	dN1_dxsi = -0.25 * (1 - eta);
	dN2_dxsi = 0.25 * (1 - eta);
	dN3_dxsi = 0.25 * (1 + eta);
	dN4_dxsi = -0.25 * (1 + eta);

	dN1_deta = -0.25 * (1 - xsi);
	dN2_deta = -0.25 * (1 + xsi);
	dN3_deta = 0.25 * (1 + xsi);
	dN4_deta = 0.25 * (1 - xsi);

	for (i = 0; i < 3; i++) {
		V_xsi[i] = dN1_dxsi * Elcor[0][i] + dN2_dxsi * Elcor[1][i] + dN3_dxsi * Elcor[2][i] + dN4_dxsi * Elcor[3][i];
	}


	for (i = 0; i < 3; i++) {
		V_eta[i] = dN1_deta * Elcor[0][i] + dN2_deta * Elcor[1][i] + dN3_deta * Elcor[2][i] + dN4_deta * Elcor[3][i];
	}



	Vnorm[0] = V_xsi[1] * V_eta[2] - V_eta[1] * V_xsi[2];
	Vnorm[1] = V_xsi[2] * V_eta[0] - V_xsi[0] * V_eta[2];
	Vnorm[2] = V_xsi[0] * V_eta[1] - V_xsi[1] * V_eta[0];

	Jac = sqrt(Vnorm[0] * Vnorm[0] + Vnorm[1] * Vnorm[1] + Vnorm[2] * Vnorm[2]);

	Vnorm[0] = Vnorm[0] / Jac;
	Vnorm[1] = Vnorm[1] / Jac;
	Vnorm[2] = Vnorm[2] / Jac;


}

void BEMPotential3D::Cartesian(double** Elcor, double& N1, double& N2, double& N3, double& N4, double* GCcor)
{

	int i;

	for (i = 0; i < 3; i++) {

		GCcor[i] = N1 * Elcor[0][i] + N2 * Elcor[1][i] + N3 * Elcor[2][i] + N4 * Elcor[3][i];
	}

}

double BEMPotential3D::R(double* x)
{
	double R = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];

	R = sqrt(R);
	return R;
}

double BEMPotential3D::R_01(double* x, int i)
{
	double R_i = x[i] / R(x);
	return R_i;
}

double BEMPotential3D::R_02(double* x, int i, int j)
{
	double R_ij = (d[i][j] - (x[i] / R(x)) * (x[j] / R(x))) / R(x);
	return R_ij;
}

double BEMPotential3D::R_03(double* x, int i, int j, int k)
{
	double R_ijk = (3 * (x[i] / R(x)) * (x[j] / R(x)) * (x[k] / R(x)) - d[i][j] * x[k] / R(x) - d[i][k] * x[j] / R(x) - d[j][k] * x[i] / R(x)) / (R(x) * R(x));
	return R_ijk;
}

double BEMPotential3D::fR(double* x, int sign)
{
	double R1 = sign * x[2] + R(x);
	double fR = 1 / R1;
	return fR;
}

double BEMPotential3D::fR_01(double* x, int i, int sign)
{
	double R1 = sign * x[2] + R(x);
	double fR_i = -(R_01(x, i) + sign * d[2][i]) / (R1 * R1);
	return fR_i;
}

double BEMPotential3D::fR_02(double* x, int i, int j, int sign)
{
	double R1 = R(x) + sign * x[2];
	double fR_ij = -1.00 / (R1 * R1 * R1) * (R1 * R_02(x, i, j) - 2 * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]));
	return fR_ij;
}

double BEMPotential3D::fR_03(double* x, int i, int j, int k, int sign)
{
	double R1 = R(x) + sign * x[2];
	double fR_ijk = -1.00 / (R1 * R1 * R1 * R1) * (R1 * R1 * R_03(x, i, j, k) - 2 * R1 * ((R_01(x, i) + sign * d[2][i]) * R_02(x, j, k) + (R_01(x, j) + sign * d[2][j]) * R_02(x, i, k) + (R_01(x, k) + sign * d[2][k]) * R_02(x, i, j))
		+ 6 * (R_01(x, i) + sign * d[2][i]) * (R_01(x, j) + sign * d[2][j]) * (R_01(x, k) + sign * d[2][k]));
	return fR_ijk;
}

double BEMPotential3D::Green_01(double* x, int i) {

	double G_ijk = -1.0 / (4 * Pi * k_0) * fR_01(x, i, 0);

	return G_ijk;
}

double BEMPotential3D::Green_02(double* x, int i, int j) {

	double G_ijkl = -1.0 / (4 * Pi * k_0) * fR_02(x, i, j, 0);

	return G_ijkl;
}

double BEMPotential3D::Green_03(double* x, int i, int j, int k) {

	double G_ijkls = -1.0 / (4 * Pi * k_0) * fR_03(x, i, j, k, 0);

	return G_ijkls;
}

void BEMPotential3D::UK(double* x, double r, double& UP)
{
	UP = -1.0 / (4 * Pi * k_0) * (1.0 / r);
}

void BEMPotential3D::TK(double* x, double r, double& TP, double* Vnorm)
{


	TP = 0.0;

	for (int j = 0; j < 3; j++) {

		TP = TP - k_0 * Vnorm[j] * Green_01(x, j);

	}


}