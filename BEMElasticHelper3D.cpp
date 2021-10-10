# include "BEMElasticHelper3D.h"
# include "Eigen/Dense"


void Trans_Tri(int ntr, int l_pi, double xsib, double etab, double& xsi, double& eta, double& Jacb)
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

void Triangle_Coord(int ntr, int l_pi, int** cor_tri)
{
	int n, temp;

	// double *xsii , *etai

	int nod_tri[3] = { 1,2,3 };
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


void Serendip_func(double xsi, double eta, double& N1, double& N2, double& N3, double& N4, double* N, double** Ni)

{

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 12; j++) {
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

	Ni[0][0] = Ni[1][1] = Ni[2][2] = N1;
	Ni[0][3] = Ni[1][4] = Ni[2][5] = N2;
	Ni[0][6] = Ni[1][7] = Ni[2][8] = N3;
	Ni[0][9] = Ni[1][10] = Ni[2][11] = N4;

}

void Normal_Jac(double xsi, double eta, double** Elcor, double& Jac, double* Vnorm)

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

void Cartesian(double** Elcor, double& N1, double& N2, double& N3, double& N4, double* GCcor)
{

	int i;

	for (i = 0; i < 3; i++) {

		GCcor[i] = N1 * Elcor[0][i] + N2 * Elcor[1][i] + N3 * Elcor[2][i] + N4 * Elcor[3][i];
	}

}

void UK(double* dxr, double r, double E, double nu, double** UP)
{
	double G, c1, c, conr;

	G = E / (2 * (1 + nu));
	c1 = 3 - 4 * nu;
	c = 1 / (16 * pi * G * (1 - nu));
	conr = c / r;

	UP[0][0] = conr * (c1 + dxr[0] * dxr[0]);
	UP[1][1] = conr * (c1 + dxr[1] * dxr[1]);
	UP[2][2] = conr * (c1 + dxr[2] * dxr[2]);
	UP[0][1] = conr * dxr[0] * dxr[1];
	UP[0][2] = conr * dxr[0] * dxr[2];
	UP[1][0] = UP[0][1];
	UP[1][2] = conr * dxr[1] * dxr[2];
	UP[2][0] = conr * dxr[2] * dxr[0];
	UP[2][1] = conr * dxr[2] * dxr[1];
}


void TK(double* dxr, double r, double E, double nu, double** TP, double* Vnorm)
{
	double c3, c2, Costh, Conr;
	int i, j;

	c3 = 1 - 2 * nu;
	c2 = 1 / (8 * pi * (1 - nu));
	Costh = Vnorm[0] * dxr[0] + Vnorm[1] * dxr[1] + Vnorm[2] * dxr[2];
	Conr = c2 / (r * r);



	TP[0][0] = -Conr * (c3 + 3 * dxr[0] * dxr[0]) * Costh;
	TP[1][1] = -Conr * (c3 + 3 * dxr[1] * dxr[1]) * Costh;
	TP[2][2] = -Conr * (c3 + 3 * dxr[2] * dxr[2]) * Costh;

	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			if (i != j) {
				TP[i][j] = -Conr * (3 * dxr[i] * dxr[j] * Costh - c3 * (Vnorm[j] * dxr[i] - Vnorm[i] * dxr[j]));
			}
		}
	}

}
