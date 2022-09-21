#include "configPotential3D_elliptical.h"
#include "integratorhelperPotential3D_elliptical.h"
#include "BEMBuilderPotential3D.h"
#include "BEMPotentialHelper3D.h"

using namespace BEMPotential3D;

void BEMBuilderPotential3D::Integ3(int* sctrb, double** Elcor, double** dUe, double** dTe, Ref<MatrixXd> xP, int Ncol) {
	
	int i, n, m, k, ii, jj, judge, l_pi, iD, nD, ntr, Ntri, ngp;

	double xsi, eta, W, r, Jac, xsib, etab, N1, N2, N3, N4, Jacb;


	double** Ni;

	int Ndof = 1;

	double* x = new double[3];
	double* N = new double[4];
	double* Vnorm = new double[3];
	double* GCcor = new double[3];


	double w[8];
	double gp[8];


	double UP = 0;
	double TP = 0;

	Ni = new double* [1];

	for (i = 0; i < 1; i++) {
		Ni[i] = new double[4];
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


					x[0] = (GCcor[0] - xP(i, 0));
					x[1] = (GCcor[1] - xP(i, 1));
					x[2] = (GCcor[2] - xP(i, 2));


					UK(x, r, UP);
					TK(x, r, TP, Vnorm);



					for (ii = 0; ii < Ndof; ii++) {
						iD = Ndof * i + ii;

						for (jj = 0; jj < Ndof; jj++) {

							for (n = 0; n < 4; n++) {
								nD = Ndof * n + jj;

								dUe[iD][nD] = dUe[iD][nD] + N[n] * UP * Jac * W;
								dTe[iD][nD] = dTe[iD][nD] + N[n] * TP * Jac * W;
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


						x[0] = (GCcor[0] - xP(i, 0));
						x[1] = (GCcor[1] - xP(i, 1));
						x[2] = (GCcor[2] - xP(i, 2));

						UK(x, r, UP);
						TK(x, r, TP, Vnorm);

						for (ii = 0; ii < Ndof; ii++) {
							iD = Ndof * i + ii;

							for (jj = 0; jj < Ndof; jj++) {
								for (n = 0; n < 4; n++) {

									nD = Ndof * n + jj;
									dUe[iD][nD] = dUe[iD][nD] + N[n] * UP * Jac * W;

									if (sctrb[n] != i + 1) {

										dTe[iD][nD] = dTe[iD][nD] + N[n] * TP * Jac * W;
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


	for (int i = 0; i < 1; i++) {
		delete Ni[i];
	}

	delete[] Ni;
}

void BEMBuilderPotential3D::addBEM(Config& config_) {

	configPotential3D_elliptical& config = dynamic_cast<configPotential3D_elliptical&> (config_);

	int& NE = config.NE;
	int& NN = config.NN;
	Ref<MatrixXi> NConnect = config.NConnect;
	Ref<MatrixXd> NODES = config.NODES;
	Ref<MatrixXd> HMAT = config.HMAT;
	Ref<MatrixXd> GMAT = config.GMAT;
	int Dimension = config.Dimension;


	VectorXd Diag = VectorXd::Zero(NN);

#pragma omp parallel shared(HMAT,GMAT) 
	{

		double** dUe;
		dUe = new double* [NN];
		for (int j = 0; j < NN; j++) {
			dUe[j] = new double[4];
		}


		double** dTe;
		dTe = new double* [NN];
		for (int j = 0; j < NN; j++) {
			dTe[j] = new double[4];
		}

		int* sctrb = new int[4];

		double** Elcor;
		Elcor = new double* [4];
		for (int j = 0; j < 4; j++)
		{
			Elcor[j] = new double[3];
		}

#pragma omp for 

		for (int i = 0; i < NE; i++) {

			sctrb[0] = NConnect(i, 0);
			sctrb[1] = NConnect(i, 1);
			sctrb[2] = NConnect(i, 2);
			sctrb[3] = NConnect(i, 3);


			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 3; k++) {
					Elcor[j][k] = NODES(sctrb[j] - 1, k);
				}
			}



			for (int ii = 0; ii < NN; ii++) {
				for (int jj = 0; jj < 4; jj++) {
					dUe[ii][jj] = 0.0;
					dTe[ii][jj] = 0.0;
				}
			}



			Integ3(sctrb, Elcor, dUe, dTe, NODES, NN);

			for (int l = 0; l < NN; l++) {
				for (int k = 0; k < 4; k++) {
					GMAT(l, 4 * i + k) = dUe[l][k];
				}
			}

# pragma omp critical
			{
				for (int jj = 0; jj < 4; jj++)
				{
					int pp = 0;
					pp = (sctrb[jj] - 1);

					for (int ll = 0; ll < NN; ll++) {
						HMAT(ll, pp) = HMAT(ll, pp) + dTe[ll][jj];
						//	GMAT[l][p + k] = GMAT[l][p + k] + dUe[l][3 * j + k];
					}

				}


			}

			//	printf("%d", i);

		}


		for (int i = 0; i < NN; i++) {
			delete dTe[i];
			delete dUe[i];
		}

		delete[] dTe;
		delete[] dUe;



	}



	for (int i = 0; i < NN; i++) {
		for (int k = 0; k < NN; k++) {
			Diag(k) = Diag(k) + HMAT(k, i);
		}
	}

	for (int m = 0; m < NN; m++) {
		HMAT(m, m) = -Diag(m);
	}
	

	Diag.resize(0);

}

void BEMBuilderPotential3D::addFluxEquivalentBEM(Config& config_) {};

void BEMBuilderPotential3D::addFluxEquivalentFirstOrderBEM(Config& config_) {};
void BEMBuilderPotential3D::addFluxEquivalentSecondOrderBEM(Config& config_) {};
void BEMBuilderPotential3D::ApplyBC(Config& config_) {};
void BEMBuilderPotential3D::getUTVector(Config & config_) {};

int BEMBuilderPotential3D::ffloord(int i, int j) {

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

};
int BEMBuilderPotential3D::Find(int h, int e, int i, double** record, Ref<MatrixXi> NConnect, int NE) {

	for (int j = 0; j < 4 * NE; j++) {
		if ((NConnect(e, i)) == record[j][1]) {

			h = j;
			goto pao;
		}

	}

pao:
	;
	return h;

};