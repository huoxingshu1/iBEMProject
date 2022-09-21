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

#include "configElastic3D_arbitrary.h"
#include "integratorElastic3Dhelper.h"
#include "BEMBuilderElastic3D_arbitrary.h"
#include "BEMElasticHelper3D.h"

/*
    Private members
*/

using namespace Elastic3DSpherical;

int BEMBuilderElastic3D_arbitrary::ffloord(int i, int j)
{

	int d;

	if (i == 0 || i == 1 || i == 2) {
		d = 0;
	}

	if (i == 3 || i == 4 || i == 5) {
		d = 1;
	}

	if (i == 6 || i == 7 || i == 8) {
		d = 2;
	}

	if (i == 9 || i == 10 || i == 11) {
		d = 3;
	}

	return d;

}

int BEMBuilderElastic3D_arbitrary::Find(int h, int e, int i, double** record, Ref<MatrixXi> NConnect, int NE)
{


	int t, p;
	p = 3;

	t = ffloord(i, p);

	for (int j = 0; j < 12 * NE; j++) {
		if ((NConnect(e, t) + 0.1 * (i % 3 + 1)) == record[j][1]) {

			h = j;
			goto pao;
		}

	}

pao:
	;

	return h;


}

void BEMBuilderElastic3D_arbitrary::Integ3(double E, double nu, int* sctrb, double** Elcor, double** dUe, double** dTe, Ref<MatrixXd> xP, int Ncol) {

	int i, n, m, k, ii, jj, judge, l_pi, iD, nD, ntr, Ntri, ngp;

	double xsi, eta, W, r, Jac, xsib, etab, N1, N2, N3, N4, Jacb;


	double** UP, ** TP, ** Ni;

	int Ndof = 3;

	double* dxr = new double[3];
	double* N = new double[4];
	double* Vnorm = new double[3];
	double* GCcor = new double[3];


	double w[8];
	double gp[8];


	UP = new double* [3];
	TP = new double* [3];
	Ni = new double* [3];

	for (i = 0; i < 3; i++) {
		UP[i] = new double[3];
		TP[i] = new double[3];
		Ni[i] = new double[12];
	}


	ngp = 8;
	gp[0] = 0.960289856497536;
	gp[1] = -0.960289856497536;
	gp[2] = 0.796666477413627;
	gp[3] = -0.796666477413627;
	gp[4] = 0.525532409916329;
	gp[5] = -0.525532409916329;
	gp[6] = 0.183434642495650;
	gp[7] = -0.183434642495650;





	w[0] = 0.101228536290376;
	w[1] = 0.101228536290376;
	w[2] = 0.222381034453374;
	w[3] = 0.222381034453374;
	w[4] = 0.313706645877887;
	w[5] = 0.313706645877887;
	w[6] = 0.362683783378362;
	w[7] = 0.362683783378362;


	Ntri = 2;


	for (i = 0; i < Ncol; i++) {
		judge = 0;
		for (n = 0; n < 4; n++) {
			if (sctrb[n] == i + 1) {
				judge = 1;
				l_pi = n;
			}
		}



		if (judge == 0) {
			for (m = 0; m < ngp; m++) {
				xsi = gp[m];

				for (k = 0; k < ngp; k++) {
					eta = gp[k];
					W = w[m] * w[k];

					Serendip_func(xsi, eta, N1, N2, N3, N4, N, Ni);
					Normal_Jac(xsi, eta, Elcor, Jac, Vnorm);
					Cartesian(Elcor, N1, N2, N3, N4, GCcor);

					r = sqrt((GCcor[0] - xP(i, 0)) * (GCcor[0] - xP(i, 0)) + (GCcor[1] - xP(i, 1)) * (GCcor[1] - xP(i, 1)) + (GCcor[2] - xP(i, 2)) * (GCcor[2] - xP(i, 2)));


					dxr[0] = (GCcor[0] - xP(i, 0)) / r;
					dxr[1] = (GCcor[1] - xP(i, 1)) / r;
					dxr[2] = (GCcor[2] - xP(i, 2)) / r;





					UK(dxr, r, E, nu, UP);
					TK(dxr, r, E, nu, TP, Vnorm);



					for (ii = 0; ii < Ndof; ii++) {
						iD = 3 * i + ii;

						for (jj = 0; jj < Ndof; jj++) {

							for (n = 0; n < 4; n++) {
								nD = Ndof * n + jj;

								dUe[iD][nD] = dUe[iD][nD] + N[n] * UP[ii][jj] * Jac * W;
								dTe[iD][nD] = dTe[iD][nD] + N[n] * TP[ii][jj] * Jac * W;
							}
						}
					}
				}
			}
		}


		else {
			for (ntr = 0; ntr < Ntri; ntr++) {
				for (m = 0; m < ngp; m++) {
					xsib = gp[m];
					for (k = 0; k < ngp; k++) {
						etab = gp[k];
						W = w[m] * w[k];

						Trans_Tri(ntr, l_pi, xsib, etab, xsi, eta, Jacb);

						Serendip_func(xsi, eta, N1, N2, N3, N4, N, Ni);

						Normal_Jac(xsi, eta, Elcor, Jac, Vnorm);

						Jac = Jac * Jacb;

						Cartesian(Elcor, N1, N2, N3, N4, GCcor);

						r = sqrt((GCcor[0] - xP(i, 0)) * (GCcor[0] - xP(i, 0)) + (GCcor[1] - xP(i, 1)) * (GCcor[1] - xP(i, 1)) + (GCcor[2] - xP(i, 2)) * (GCcor[2] - xP(i, 2)));


						dxr[0] = (GCcor[0] - xP(i, 0)) / r;
						dxr[1] = (GCcor[1] - xP(i, 1)) / r;
						dxr[2] = (GCcor[2] - xP(i, 2)) / r;

						UK(dxr, r, E, nu, UP);
						TK(dxr, r, E, nu, TP, Vnorm);

						for (ii = 0; ii < Ndof; ii++) {

							iD = Ndof * i + ii;
							for (jj = 0; jj < Ndof; jj++) {
								for (n = 0; n < 4; n++) {
									nD = Ndof * n + jj;
									dUe[iD][nD] = dUe[iD][nD] + N[n] * UP[ii][jj] * Jac * W;

									if (sctrb[n] != i + 1) {

										dTe[iD][nD] = dTe[iD][nD] + N[n] * TP[ii][jj] * Jac * W;
									}
								}
							}
						}

					}

				}
			}
		}

	}


	delete[] N;
	delete[] Vnorm;
	delete[] GCcor;


	for (int i = 0; i < 3; i++) {
		delete UP[i];
		delete TP[i];
		delete Ni[i];
	}

	delete[] UP;
	delete[] TP;
	delete[] Ni;

}

void BEMBuilderElastic3D_arbitrary::get_for_stress_BEM(double* Vnorm, double* point, double* GCcor, double** Dst, double** Dsu)
{


	int t = 0;

	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 3; j++) {
			Dst[i][j] = 0.0;
			Dsu[i][j] = 0.0;
		}
	}

	double x[3];

	x[0] = GCcor[0] - point[0];
	x[1] = GCcor[1] - point[1];
	x[2] = GCcor[2] - point[2];


	double TG1[3];
	double TG[9];
	double sgm[9];


	for (int i = 0; i < 3; i++) { TG1[i] = 0.0; }
	for (int i = 0; i < 3 * 3; i++) { TG[i] = 0.0; }
	for (int i = 0; i < 3 * 3; i++) { sgm[i] = 0.0; }

	for (int n = 0; n < 3; n++) {
		for (int i = 0; i < 3; i++) {
			TG1[n] = TG1[n] + Green_01(x, n, i, i);
		}
	}


	int tt = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int l = 0; l < 3; l++) {
				TG[tt] = TG[tt] + Green_02(x, l, i, l, j);
			}
			tt = tt + 1;
		}
	}

	tt = 0;
	for (int m = 0; m < 3; m++) {
		for (int n = 0; n < 3; n++) {
			for (int i = 0; i < 3; i++) {
				sgm[tt] = sgm[tt] + 2.0 * mu_0 * nu0 / (1 - 2 * nu0) * d[m][n] * TG[3 * i + i] + 2.0 * mu_0 * Green_02(x, m, i, n, i);
			}
			tt = tt + 1;
		}
	}



	for (int k = 0; k < 3; k++) {
		for (int l = k; l < 3; l++) {
			for (int n = 0; n < 3; n++) {

				for (int m = 0; m < 3; m++) {

					Dsu[t][n] = Dsu[t][n] + ((2.0 * mu_0 * nu0 / (1 - 2 * nu0) - 2.0 * mu_1 * nu1 / (1 - 2 * nu1)) * d[k][l] * sgm[3 * m + n] + (mu_0 - mu_1) * (2.0 * mu_0 * nu0 / (1 - 2 * nu0) * d[m][n] * (TG[3 * k + l] + TG[3 * l + k]) + mu_0 * Green_02(x, m, k, n, l) + mu_0 * Green_02(x, m, l, n, k) + mu_0 * Green_02(x, n, k, m, l) + mu_0 * Green_02(x, n, l, m, k))) * Vnorm[m];
				}


				Dst[t][n] = Dst[t][n] - ((2.0 * mu_0 * nu0 / (1 - 2 * nu0) - 2.0 * mu_1 * nu1 / (1 - 2 * nu1)) * d[k][l] * TG1[n] + (mu_0 - mu_1) * Green_01(x, n, k, l) + (mu_0 - mu_1) * Green_01(x, n, l, k));

			}

			t = t + 1;
		}
	}



}

void BEMBuilderElastic3D_arbitrary::Integ3_stress(double** Elcor, double** DDsu, double** DDst, double* point)
{
	int i, m, k, ngp;

	double xsi, eta, W, Jac, N1, N2, N3, N4;



	double** Dsu, ** Dst, ** Ni, ** TDDsu, ** TDDst;

	int Ndof = 3;


	double* N = new double[4];
	double* Vnorm = new double[3];
	double* GCcor = new double[3];

	double gp[8];
	double w[8];

	Dsu = new double* [6];
	Dst = new double* [6];
	Ni = new double* [3];
	TDDsu = new double* [6];
	TDDst = new double* [6];

	for (i = 0; i < 6; i++) { Dsu[i] = new double[3];	Dst[i] = new double[3]; TDDsu[i] = new double[12]; TDDst[i] = new double[12]; }
	for (i = 0; i < 3; i++) { Ni[i] = new double[12]; }

	for (i = 0; i < 6; i++) {
		for (int j = 0; j < 12; j++) {
			TDDsu[i][j] = 0.0;
			TDDst[i][j] = 0.0;
		}
	}






	ngp = 8;

	gp[0] = 0.960289856497536;
	gp[1] = -0.960289856497536;
	gp[2] = 0.796666477413627;
	gp[3] = -0.796666477413627;
	gp[4] = 0.525532409916329;
	gp[5] = -0.525532409916329;
	gp[6] = 0.183434642495650;
	gp[7] = -0.183434642495650;
	w[0] = 0.101228536290376;
	w[1] = 0.101228536290376;
	w[2] = 0.222381034453374;
	w[3] = 0.222381034453374;
	w[4] = 0.313706645877887;
	w[5] = 0.313706645877887;
	w[6] = 0.362683783378362;
	w[7] = 0.362683783378362;


	for (m = 0; m < ngp; m++) {
		xsi = gp[m];

		for (k = 0; k < ngp; k++) {
			eta = gp[k];
			W = w[m] * w[k];

			Serendip_func(xsi, eta, N1, N2, N3, N4, N, Ni);
			Normal_Jac(xsi, eta, Elcor, Jac, Vnorm);
			Cartesian(Elcor, N1, N2, N3, N4, GCcor);

			get_for_stress_BEM(Vnorm, point, GCcor, Dst, Dsu);

			for (int mm = 0; mm < 6; mm++) {
				for (int ll = 0; ll < 12; ll++) {
					for (int nn = 0; nn < 3; nn++) {
						DDsu[mm][ll] = DDsu[mm][ll] + Dsu[mm][nn] * Ni[nn][ll] * W * Jac;
						DDst[mm][ll] = DDst[mm][ll] + Dst[mm][nn] * Ni[nn][ll] * W * Jac;
					}
				}
			}

	

		}
	}



	delete[] N;
	delete[] Vnorm;
	delete[] GCcor;



	for (int i = 0; i < 6; i++) {
		delete Dsu[i];
		delete Dst[i];
	}

	for (int i = 0; i < 3; i++) {
		delete Ni[i];
	}

	delete[] Dsu;
	delete[] Dst;
	delete[] Ni;



}


/*
	Public members
*/
// two empty functions
void BEMBuilderElastic3D_arbitrary::addFluxEquivalentFirstOrderBEM(Config& config_)
{}

void BEMBuilderElastic3D_arbitrary::addFluxEquivalentSecondOrderBEM(Config& config_)
{}


void BEMBuilderElastic3D_arbitrary::addBEM(Config& config_) {


	configElastic3D_arbitrary& config = dynamic_cast<configElastic3D_arbitrary&> (config_);

	double& EE = config.EE;
	double& XNU = config.XNU;
	int& NE = config.NE;
	int& NN = config.NN;
	Ref<MatrixXi> NConnect = config.NConnect;
	Ref<MatrixXd> NODES = config.NODES;
	Ref<MatrixXd> HMAT = config.HMAT;
	Ref<MatrixXd> GMAT = config.GMAT;
	int Dimension = config.Dimension;

	//EE = 2 * mu_0 * (1 + nu0);
	//XNU = nu0;

	MatrixXd Diag = MatrixX3d::Zero(Dimension * NN, Dimension);

#pragma omp parallel shared(HMAT,GMAT) 
	{

		double** dUe;
		dUe = new double* [Dimension * NN];
		for (int j = 0; j < Dimension * NN; j++) {
			dUe[j] = new double[12];
		}


		double** dTe;
		dTe = new double* [Dimension * NN];
		for (int j = 0; j < Dimension * NN; j++) {
			dTe[j] = new double[12];
		}

		int* sctrb = new int[4];

		double** Elcor;
		Elcor = new double* [4];
		for (int j = 0; j < 4; j++)
		{
			Elcor[j] = new double[Dimension];
		}

#pragma omp for 

		for (int i = 0; i < NE; i++) {

			sctrb[0] = NConnect(i, 0);
			sctrb[1] = NConnect(i, 1);
			sctrb[2] = NConnect(i, 2);
			sctrb[Dimension] = NConnect(i, Dimension);


			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < Dimension; k++) {
					Elcor[j][k] = NODES(sctrb[j] - 1, k);
				}
			}



			for (int ii = 0; ii < Dimension * NN; ii++) {
				for (int jj = 0; jj < 12; jj++) {
					dUe[ii][jj] = 0.0;
					dTe[ii][jj] = 0.0;
				}
			}



			Integ3(EE, XNU, sctrb, Elcor, dUe, dTe, NODES, NN);




			for (int l = 0; l < Dimension * NN; l++) {
				for (int k = 0; k < 12; k++) {
					GMAT(l, 12 * i + k) = dUe[l][k];
				}
			}

# pragma omp critical
			{
				for (int jj = 0; jj < 4; jj++)
				{
					int pp = 0;
					pp = Dimension * (sctrb[jj] - 1);
					for (int kk = 0; kk < Dimension; kk++) {
						for (int ll = 0; ll < Dimension * NN; ll++) {

							HMAT(ll, pp + kk) = HMAT(ll, pp + kk) + dTe[ll][Dimension * jj + kk];
							//	GMAT[l][p + k] = GMAT[l][p + k] + dUe[l][Dimension * j + k];
						}
					}
				}


			}

			//	printf("%d", i);




		}


		for (int i = 0; i < Dimension * NN; i++) {
			delete dTe[i];
			delete dUe[i];
		}

		delete[] dTe;
		delete[] dUe;



	}



	for (int i = 0; i < NN; i++) {


		for (int j = 0; j < Dimension; j++) {

			int n = i * Dimension + j;

			for (int k = 0; k < Dimension * NN; k++) {

				Diag(k, j) = Diag(k, j) + HMAT(k, n);

			}
		}
	}

	for (int m = 0; m < NN; m++) {

		int posi = m * Dimension;
		for (int n = posi; n < posi + Dimension; n++) {
			for (int l = posi; l < posi + Dimension; l++) {
				HMAT(n, l) = -Diag(n, l - posi);
			}
		}
	}

	//cout << HMAT;
	//cout << "\n";
	//cout << GMAT;
	//cout << "\n";
	Diag.resize(0, 0);

}

void BEMBuilderElastic3D_arbitrary::addFluxEquivalentBEM(Config& config_)
{
	configElastic3D_arbitrary& config = dynamic_cast<configElastic3D_arbitrary&> (config_);

	int& num_NN = config.num_NN;

	int*& count_number = config.count_number;
	double***& changed_coord = config.changed_coord;
	double**& number_weight = config.number_weight;


	double& EE = config.EE;
	double& XNU = config.XNU;
	int& NE = config.NE;
	int& NN = config.NN;
	Ref<MatrixXi> NConnect = config.NConnect;
	Ref<MatrixXd> NODES = config.NODES;
	Ref<MatrixXd> HMAT = config.HMAT;
	Ref<MatrixXd> GMAT = config.GMAT;
	int Dimension = config.Dimension;

# pragma omp parallel shared (HMAT,GMAT)

	{
		double** DDsu, ** DDst;

		double* point = new double[3];

		DDsu = new double* [6];
		DDst = new double* [6];


		double** Elcor;
		Elcor = new double* [4];
		for (int j = 0; j < 4; j++) Elcor[j] = new double[3];


		for (int i = 0; i < 6; i++) {
			DDsu[i] = new double[12];  DDst[i] = new double[12];
		}

		int* sctrb = new int[4];



#pragma omp for schedule(dynamic)

		for (int ee = 0; ee < NE; ee++) {

			for (int KK = 0; KK < num_NN; KK++) {

				for (int s = 0; s < count_number[KK]; s++) {
					for (int mm = 0; mm < 6; mm++) {
						for (int ll = 0; ll < 12; ll++) {
							DDsu[mm][ll] = 0.0;
							DDst[mm][ll] = 0.0;
						}
					}

					point[0] = changed_coord[KK][s][0];  point[1] = changed_coord[KK][s][1]; point[2] = changed_coord[KK][s][2];

					sctrb[0] = NConnect(ee, 0);
					sctrb[1] = NConnect(ee, 1);
					sctrb[2] = NConnect(ee, 2);
					sctrb[3] = NConnect(ee, 3);

					for (int j = 0; j < 4; j++) {
						for (int k = 0; k < 3; k++) {
							Elcor[j][k] = NODES(sctrb[j] - 1, k);
						}
					}

					Integ3_stress(Elcor, DDsu, DDst, point);




# pragma omp  critical 
					{

						for (int j = 0; j < 4; j++)
						{
							int p = 0;
							p = 3 * (sctrb[j] - 1);

							for (int k = 0; k < 3; k++) {
								for (int l = 0; l < 6; l++) {

									HMAT(3 * NN + 6 * KK + l, p + k) = HMAT(3 * NN + 6 * KK + l, p + k) - DDsu[l][3 * j + k] * number_weight[KK][s];
								}
							}
						}

						for (int l = 0; l < 6; l++) {
							for (int k = 0; k < 12; k++) {
								GMAT(3 * NN + 6 * KK + l, 12 * ee + k) += DDst[l][k] * number_weight[KK][s];
							}
						}

					}
				}




			}

		}

		for (int i = 0; i < 6; i++) {
			delete DDsu[i];
			delete DDst[i];
		}

		for (int i = 0; i < 4; i++) {
			delete Elcor[i];
		}

		delete[] DDsu;
		delete[] DDst;
		delete[] sctrb;
		delete[] point;
		delete[] Elcor;
	}

}

void BEMBuilderElastic3D_arbitrary::ApplyBC(Config& config_)
{
	configElastic3D_arbitrary& config = dynamic_cast<configElastic3D_arbitrary&> (config_);

	int& NE = config.NE;
	int& NN = config.NN;
	int& num_NN = config.num_NN;
	Ref<MatrixXi> NConnect = config.NConnect;
	Ref<MatrixXd> NODES = config.NODES;
	Ref<MatrixXd> BC = config.BC;
	Ref<MatrixXd> t_bc = config.t_bc;
	MatrixXd& HMAT = config.HMAT;
	MatrixXd& GMAT = config.GMAT;
	MatrixXd& AA = config.A;
	VectorXd& BB = config.b;

	record = new double* [12 * NE];

	for (int i = 0; i < 12 * NE; i++) {
		record[i] = new double[2];
	}


	for (int i = 0; i < 12 * NE; i++) {
		for (int j = 0; j < 2; j++) {
			record[i][j] = 0.0;
		}
	}

	int t = 0;
	int h;

	for (int j = 0; j < NN * 3; j++) {
		if (BC(j, 0) == 1) {

			for (int i = 0; i < NN * 3 + 6 * num_NN; i++) {
				AA(i, t) = HMAT(i, j);
			}
			t = t + 1;
		}
	}

	for (int j = 0; j < NN * 3; j++) {
		if (BC(j, 0) == 0) {
			for (int i = 0; i < NN * 3 + 6 * num_NN; i++) {
				BB(i) = BB(i) - HMAT(i, j) * BC(j, 1);
			}
		}
	}


	int mm = 0; int p = 3;

	for (int e = 0; e < NE; e++) {

		for (int i = 0; i < 12; i++) {

			if (t_bc(e, i) == 1) {
				for (int row = 0; row < NN * 3 + 6 * num_NN; row++) { BB(row) = BB(row) + t_bc(e, 12 + i) * GMAT(row, 12 * e + i); }
			}

			else if (t_bc(e, i) == 0) {

				h = -1;
				h = Find(h, e, i, record, NConnect, NE);

				if (h == -1) {

					record[mm][0] = t;  int tt = ffloord(i, p);
					record[mm][1] = NConnect(e, tt) + 0.1 * (i % 3 + 1);
					mm = mm + 1;

					for (int ii = 0; ii < NN * 3 + 6 * num_NN; ii++)
					{
						AA(ii, t) = -GMAT(ii, 12 * e + i);
					}
					t = t + 1;
				}

				else {
					int cc = record[h][0];
					for (int ii = 0; ii < NN * 3 + 6 * num_NN; ii++) { AA(ii, cc) = AA(ii, cc) - GMAT(ii, 12 * e + i); }
				}
			}
		}
	}



	for (int i = 0; i < 3 * NN + 6 * num_NN; i++) {
		for (int j = 3 * NN; j < 3 * NN + 6 * num_NN; j++) {
			AA(i, j) = HMAT(i, j);
		}
	}
	config.HMAT.resize(0, 0);
	config.GMAT.resize(0, 0);

}

void BEMBuilderElastic3D_arbitrary::getUTVector(Config& config_)
{
	configElastic3D_arbitrary& config = dynamic_cast<configElastic3D_arbitrary&> (config_);

	VectorXd& U = config.U;
	VectorXd& T = config.T;
	MatrixXd& BC = config.BC;
	MatrixXd& t_bc = config.t_bc;
	VectorXd& XX = config.X;
	MatrixXi& NConnect = config.NConnect;

	int& NN = config.NN;
	int& num_NN = config.num_NN;
	int& NE = config.NE;
	int t, h;
	t = 0;
	h = 0;
	int tt = 0;

	for (int i = 0; i < NN * 3; i++) {
		if (BC(i, 0) == 0) {
			U[i] = BC(i, 1);
		}
		else
		{
			U[i] = XX(t);
			t = t + 1;
		}
	}



	for (int i = 3 * NN; i < 3 * NN + 6 * num_NN; i++) {
		U[i] = XX(i);
	}

	int m = 0;

	for (int e = 0; e < NE; e++) {
		for (int i = 0; i < 12; i++) {

			if (t_bc(e, i) == 1) {
				T[m] = t_bc(e, 12 + i);
				m = m + 1;
			}


			else if (t_bc(e, i) == 0) {
				h = Find(h, e, i, record, NConnect, NE);
				tt = record[h][0];
				T[m] = XX(tt);
				m = m + 1;
			}

		}
	}

	    /*
    To simplify the files in folder, this part has been set as comments

    Users can freely modify the function to obtain U and T as you wish

    Below, we present an example, one can use it by uncomment

    Ufile -> displacements, Tfile -> Tractions, Eigen -> Eigen strains on nodes by sequence
    */

    /*
    configElastic3D_arbitrary& config = dynamic_cast<configElastic3D_arbitrary&> (config_);

    int& NE = config.NE;
    int& NN = config.NN;
    int& num_NN = config.num_NN;
    int& Dimension = config.Dimension;

    // here, U and T vectors have been filled with data (putback from boundary)
    Ref<VectorXd> U = config.U; Ref<VectorXd> T = config.T;

    ofstream myfile_U, myfile_T, myfile_Eigen;

    myfile_U.open("Ufile_3D_arbitrary.txt"); myfile_T.open("Tfile_3D_arbitrary.txt"); myfile_Eigen.open("Eigen_3D_arbitrary.txt");

    for (int i = 0; i < Dimension * NN; i++) {
        myfile_U << U(i) << endl;
    }

    for (int i = 0; i < Dimension * 4 * NE; i++) {
        myfile_T << T(i) << endl;
    }

    for (int i = 0; i < 6 * num_NN; i++) {
        myfile_Eigen << U(Dimension * NN + i) << endl;
    }

    myfile_U.close(); myfile_T.close(); myfile_Eigen.close();
    */


}