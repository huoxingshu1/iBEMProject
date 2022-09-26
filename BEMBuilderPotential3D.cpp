#include "configPotential3D_elliptical.h"
#include "integratorhelperPotential3D_elliptical.h"
#include "BEMBuilderPotential3D.h"
#include "BEMPotentialHelper3D.h"

using namespace std;
using namespace BEMPotential3D;

void get_for_stress_BEM1(double* Vnorm, double* point, double* GCcor, double** Dstw, double** Dsuw, Ref<MatrixXd> aa)
{



	int t = 0;

	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 1; j++) {
			Dstw[i][j] = 0.0;
			Dsuw[i][j] = 0.0;
		}
	}

	double x[3];

	x[0] = GCcor[0] - point[0];
	x[1] = GCcor[1] - point[1];
	x[2] = GCcor[2] - point[2];



	double*** T_Strain_w;
	double*** Tp_Strain_w;

	T_Strain_w = new double** [3];
	Tp_Strain_w = new double** [3];


	for (int i = 0; i < 3; i++) {
		T_Strain_w[i] = new double* [3];
		Tp_Strain_w[i] = new double* [3];
		for (int j = 0; j < 3; j++) {
			T_Strain_w[i][j] = new double[1];
			Tp_Strain_w[i][j] = new double[1];
		}
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int w = 0; w < 1; w++) {
				T_Strain_w[i][j][w] = 0.0;
				Tp_Strain_w[i][j][w] = 0.0;
			}
		}
	}



	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int w = 0; w < 1; w++) {
				for (int mm = 0; mm < 3; mm++) {
					T_Strain_w[i][j][w] = T_Strain_w[i][j][w] - (-k_0 * Green_03(x, mm, i, j) * Vnorm[mm]);
				}
			}
		}
	}


	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int w = 0; w < 1; w++) {


				for (int ii = 0; ii < 3; ii++) {
					for (int jj = 0; jj < 3; jj++) {

						Tp_Strain_w[i][j][w] = Tp_Strain_w[i][j][w] + aa(i, ii) * aa(j, jj) * T_Strain_w[ii][jj][0];
					}
				}

			}
		}
	}



	t = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			Dsuw[t][0] = Tp_Strain_w[i][j][0] * (-k_0 - (-k_1));
			t = t + 1;
		}
	}


	/////////////////////////////////////


	double*** U_Strain_w;
	double*** Up_Strain_w;

	U_Strain_w = new double** [3];
	Up_Strain_w = new double** [3];

	for (int i = 0; i < 3; i++) {
		U_Strain_w[i] = new double* [3];
		Up_Strain_w[i] = new double* [3];
		for (int j = 0; j < 3; j++) {
			U_Strain_w[i][j] = new double[1];
			Up_Strain_w[i][j] = new double[1];
		}
	}


	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int w = 0; w < 1; w++) {
				U_Strain_w[i][j][w] = 0.0;
				Up_Strain_w[i][j][w] = 0.0;
			}
		}
	}




	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int w = 0; w < 1; w++) {
				U_Strain_w[i][j][w] = Green_02(x, i, j);
			}
		}
	}



	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int w = 0; w < 1; w++) {

				for (int ii = 0; ii < 3; ii++) {
					for (int jj = 0; jj < 3; jj++) {
						Up_Strain_w[i][j][w] = Up_Strain_w[i][j][w] + aa(i, ii) * aa(j, jj) * U_Strain_w[ii][jj][0];

					}
				}
			}
		}
	}


	t = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			Dstw[t][0] = Up_Strain_w[i][j][0] * (-k_0 - (-k_1));
			t = t + 1;
		}
	}




	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			delete[] U_Strain_w[i][j];
			delete[] Up_Strain_w[i][j];
			delete[] T_Strain_w[i][j];
			delete[] Tp_Strain_w[i][j];

		}
		delete[] U_Strain_w[i];
		delete[] Up_Strain_w[i];
		delete[] T_Strain_w[i];
		delete[] Tp_Strain_w[i];
	}


	delete[] U_Strain_w;
	delete[] Up_Strain_w;
	delete[] T_Strain_w;
	delete[] Tp_Strain_w;




}

void Integ3_stress1(double** Elcor, double** DDsuw, double** DDstw, double* point, Ref<MatrixXd> aa)
{
	int i, m, k, ngp;

	double xsi, eta, W, Jac, N1, N2, N3, N4;



	double** Dsuw, ** Dstw, ** Ni;

	int Ndof = 1;


	double* N = new double[4];
	double* Vnorm = new double[3];
	double* GCcor = new double[3];

	double gp[8];
	double w[8];

	Dsuw = new double* [9];
	Dstw = new double* [9];
	Ni = new double* [1];


	for (i = 0; i < 9; i++) { Dsuw[i] = new double[1];	Dstw[i] = new double[1]; }
	for (i = 0; i < 1; i++) { Ni[i] = new double[4]; }




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

			get_for_stress_BEM1(Vnorm, point, GCcor, Dstw, Dsuw, aa);


			for (int mm = 0; mm < 9; mm++) {
				for (int ll = 0; ll < 4; ll++) {
					for (int nn = 0; nn < 1; nn++) {
						DDsuw[mm][ll] = DDsuw[mm][ll] + Dsuw[mm][nn] * Ni[nn][ll] * W * Jac;
						DDstw[mm][ll] = DDstw[mm][ll] + Dstw[mm][nn] * Ni[nn][ll] * W * Jac;
					}
				}
			}


		}
	}





	delete[] N;
	delete[] Vnorm;
	delete[] GCcor;


	for (int i = 0; i < 9; i++) {
		delete Dsuw[i];
		delete Dstw[i];
	}

	for (int i = 0; i < 1; i++) {
		delete Ni[i];
	}

	delete[] Dsuw;
	delete[] Dstw;
	delete[] Ni;






}

void get_for_stress_BEM(double* Vnorm, double* point, double* GCcor, double** Dst, double** Dsu, Ref<MatrixXd> aa)
{

	int t = 0;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 1; j++) {
			Dst[i][j] = 0.0;
			Dsu[i][j] = 0.0;
		}
	}

	double x[3];

	x[0] = GCcor[0] - point[0];
	x[1] = GCcor[1] - point[1];
	x[2] = GCcor[2] - point[2];



	double** T_Strain;
	double** Tp_Strain;

	T_Strain = new double* [3];
	Tp_Strain = new double* [3];

	for (int i = 0; i < 3; i++) {
		T_Strain[i] = new double[1]; Tp_Strain[i] = new double[1];
	}


	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 1; j++) {
			Tp_Strain[i][j] = 0.0;
			T_Strain[i][j] = 0.0;
		}
	}




	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			T_Strain[i][0] = T_Strain[i][0] - k_0 * Green_02(x, j, i) * Vnorm[j];

		}
	}


	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 1; j++) {

			for (int ii = 0; ii < 3; ii++) {
				Tp_Strain[i][j] = Tp_Strain[i][j] + aa(i, ii) * T_Strain[ii][0];

			}
		}
	}


	t = 0;
	for (int i = 0; i < 3; i++) {
		Dsu[t][0] = Tp_Strain[i][0] * (-k_0 - (-k_1));
		t = t + 1;
	}




	/////////////////////////////////////////////
	/////////////////////////////////////////////

	double** U_Strain;
	double** Up_Strain;

	U_Strain = new double* [3];
	Up_Strain = new double* [3];

	for (int i = 0; i < 3; i++) {
		U_Strain[i] = new double[1];
		Up_Strain[i] = new double[1];
	}


	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 1; j++) {
			Up_Strain[i][j] = 0.0;
			Up_Strain[i][j] = 0.0;
		}
	}


	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 1; j++) {
			U_Strain[i][j] = -(Green_01(x, i));
		}
	}


	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 1; j++) {
			for (int ii = 0; ii < 3; ii++) {
				Up_Strain[i][j] = Up_Strain[i][j] + aa(i, ii) * U_Strain[ii][0];
			}
		}
	}



	t = 0;
	for (int i = 0; i < 3; i++) {
		Dst[t][0] = Up_Strain[i][0] * (-k_0 - (-k_1));
		t = t + 1;
	}




	for (int i = 0; i < 3; i++) {
		delete[] U_Strain[i];
		delete[] Up_Strain[i];
		delete[] T_Strain[i];
		delete[] Tp_Strain[i];
	}

	delete[] U_Strain;
	delete[] Up_Strain;
	delete[] T_Strain;
	delete[] Tp_Strain;





}

void Integ3_stress(double** Elcor, double** DDsu, double** DDst, double* point, Ref<MatrixXd> aa)
{
	int i, m, k, ngp;

	double xsi, eta, W, Jac, N1, N2, N3, N4;

	double** Dsu, ** Dst, ** Ni;

	int Ndof = 1;


	double* N = new double[4];
	double* Vnorm = new double[3];
	double* GCcor = new double[3];

	double gp[8];
	double w[8];

	Dsu = new double* [3];
	Dst = new double* [3];
	Ni = new double* [1];


	for (i = 0; i < 3; i++) { Dsu[i] = new double[1];	Dst[i] = new double[1]; }
	for (i = 0; i < 1; i++) { Ni[i] = new double[4]; }



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

			get_for_stress_BEM(Vnorm, point, GCcor, Dst, Dsu, aa);


			for (int mm = 0; mm < 3; mm++) {
				for (int ll = 0; ll < 4; ll++) {
					for (int nn = 0; nn < 1; nn++) {
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



	for (int i = 0; i < 3; i++) {
		delete Dsu[i];
		delete Dst[i];
	}

	for (int i = 0; i < 1; i++) {
		delete Ni[i];
	}

	delete[] Dsu;
	delete[] Dst;
	delete[] Ni;


}

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

void BEMBuilderPotential3D::addFluxEquivalentBEM(Config& config_) {

	configPotential3D_elliptical& config = dynamic_cast<configPotential3D_elliptical&> (config_);

	int& NE = config.NE;
	int& NN = config.NN;
	int& num = config.num;
	Ref<MatrixXi> NConnect = config.NConnect;
	Ref<MatrixXd> NODES = config.NODES;
	Ref<MatrixXd> HMAT = config.HMAT;
	Ref<MatrixXd> GMAT = config.GMAT;
	Ref<MatrixXd> eigen_point = config.eigen_point;
	Ref<VectorXi> orient = config.orient;
	Ref<MatrixXd> direc = config.direc;

	int Dimension = config.Dimension;

# pragma omp parallel shared (HMAT,GMAT)

	{

		double** DDsu, ** DDst;

		double* point = new double[3];

		DDsu = new double* [3];
		DDst = new double* [3];

		double** Elcor;
		Elcor = new double* [4];
		for (int j = 0; j < 4; j++) Elcor[j] = new double[3];


		for (int i = 0; i < 3; i++) {
			DDsu[i] = new double[4];  DDst[i] = new double[4];
		}

		int* sctrb = new int[4];

		MatrixXd aa = MatrixXd::Zero(3, 3);


#pragma omp for

		for (int ee = 0; ee < NE; ee++) {

			for (int KK = 0; KK < num; KK++) {

				for (int mm = 0; mm < 3; mm++) {
					for (int ll = 0; ll < 4; ll++) {
						DDsu[mm][ll] = 0.0;
						DDst[mm][ll] = 0.0;
					}
				}

				point[0] = eigen_point(KK, 0);  point[1] = eigen_point(KK, 1); point[2] = eigen_point(KK, 2);

				for (int sss = 0; sss < 3; sss++) {
					for (int ggg = 0; ggg < 3; ggg++) {
						aa(sss, ggg) = direc(3 * KK + sss, ggg);
					}
				}

				sctrb[0] = NConnect(ee, 0);
				sctrb[1] = NConnect(ee, 1);
				sctrb[2] = NConnect(ee, 2);
				sctrb[3] = NConnect(ee, 3);


				for (int j = 0; j < 4; j++) {
					for (int k = 0; k < 3; k++) {
						Elcor[j][k] = NODES(sctrb[j] - 1, k);
					}
				}

				Integ3_stress(Elcor, DDsu, DDst, point, aa);


# pragma omp  critical 
				{

					for (int j = 0; j < 4; j++)
					{
						int p = 0;
						p = 1 * (sctrb[j] - 1);

						for (int k = 0; k < 1; k++) {

							for (int l = 0; l < 3; l++) {

								HMAT(NN + 12 * (KK)+l, p + k) = HMAT(NN + 12 * (KK)+l, p + k) - DDsu[l][1 * j + k];
							}
						}
					}

				}


				for (int l = 0; l < 3; l++) {
					for (int k = 0; k < 4; k++) {
						GMAT(NN + 12 * (KK)+l, 4 * ee + k) = DDst[l][k];
					}
				}


			}


		}


		for (int i = 0; i < 3; i++) {
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


};

void BEMBuilderPotential3D::addFluxEquivalentFirstOrderBEM(Config& config_) {

	configPotential3D_elliptical& config = dynamic_cast<configPotential3D_elliptical&> (config_);

	int& NE = config.NE;
	int& NN = config.NN;
	int& num = config.num;
	Ref<MatrixXi> NConnect = config.NConnect;
	Ref<MatrixXd> NODES = config.NODES;
	Ref<MatrixXd> HMAT = config.HMAT;
	Ref<MatrixXd> GMAT = config.GMAT;
	Ref<MatrixXd> eigen_point = config.eigen_point;
	Ref<VectorXi> orient = config.orient;
	Ref<MatrixXd> direc = config.direc;

	int Dimension = config.Dimension;

#pragma omp parallel shared (HMAT,GMAT)
	{
		int ttt = 0;

		double** DDsuw, ** DDstw;

		double* point = new double[3];

		DDsuw = new double* [9];
		DDstw = new double* [9];


		double** Elcor;
		Elcor = new double* [4];
		for (int j = 0; j < 4; j++) Elcor[j] = new double[3];

		for (int i = 0; i < 9; i++) {
			DDsuw[i] = new double[4];  DDstw[i] = new double[4];
		}

		int* sctrb = new int[4];
		MatrixXd aa = MatrixXd::Zero(3, 3);

# pragma omp for 

		for (int ee = 0; ee < NE; ee++) {

			for (int KK = 0; KK < num; KK++) {

				for (int mm = 0; mm < 9; mm++) {
					for (int ll = 0; ll < 4; ll++) {
						DDsuw[mm][ll] = 0.0;
						DDstw[mm][ll] = 0.0;
					}
				}

				point[0] = eigen_point(KK, 0);  point[1] = eigen_point(KK, 1); point[2] = eigen_point(KK, 2);
				sctrb[0] = NConnect(ee, 0);     sctrb[1] = NConnect(ee, 1);   sctrb[2] = NConnect(ee, 2);   sctrb[3] = NConnect(ee, 3);

				for (int sss = 0; sss < 3; sss++) {
					for (int ggg = 0; ggg < 3; ggg++) {
						aa(sss, ggg) = direc(3 * KK + sss, ggg);
					}
				}

				for (int j = 0; j < 4; j++) {
					for (int k = 0; k < 3; k++) {
						Elcor[j][k] = NODES(sctrb[j] - 1, k);
					}
				}

				Integ3_stress1(Elcor, DDsuw, DDstw, point, aa);

#pragma omp critical 
				{

					for (int j = 0; j < 4; j++)
					{
						int p = 0;
						p = 1 * (sctrb[j] - 1);

						for (int k = 0; k < 1; k++) {

							for (int l = 0; l < 9; l++) {

								HMAT(NN + 12 * (KK)+l + 3, p + k) = HMAT(NN + 12 * (KK)+l + 3, p + k) - DDsuw[l][1 * j + k];
							}
						}
					}

				}

				for (int l = 0; l < 9; l++) {
					for (int k = 0; k < 4; k++) {
						GMAT(NN + 12 * (KK)+l + 3, 4 * ee + k) = DDstw[l][k];
					}
				}

			}

		}


		for (int i = 0; i < 9; i++) {
			delete DDsuw[i];
			delete DDstw[i];
		}

		for (int i = 0; i < 4; i++) {
			delete Elcor[i];
		}
		delete[] DDsuw;
		delete[] DDstw;
		delete[] sctrb;
		delete[] Elcor;
		delete[] point;

	}

};

void BEMBuilderPotential3D::addFluxEquivalentSecondOrderBEM(Config& config_) {
};
void BEMBuilderPotential3D::ApplyBC(Config& config_) {

	configPotential3D_elliptical& config = dynamic_cast<configPotential3D_elliptical&> (config_);

	int& NE = config.NE;
	int& NN = config.NN;
	int& num = config.num;
	Ref<MatrixXi> NConnect = config.NConnect;
	Ref<MatrixXd> NODES = config.NODES;
	Ref<MatrixXd> BC = config.BC;
	Ref<MatrixXd> t_bc = config.t_bc;
	MatrixXd& HMAT = config.HMAT;
	MatrixXd& GMAT = config.GMAT;
	VectorXd& BB = config.b;

	record = new double* [4 * NE];

	for (int i = 0; i < 4 * NE; i++) {
		record[i] = new double[2];
	}


	for (int i = 0; i < 4 * NE; i++) {
		for (int j = 0; j < 2; j++) {
			record[i][j] = 0.0;
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////// apply bc ////////////////////////////////////////////////////////

	int t = 0;
	int h;

	for (int j = 0; j < NN; j++) {
		if (BC(j, 0) == 1) {
			HMAT.col(t) = HMAT.col(j);
			t = t + 1;
		}
		if (BC(j, 0) == 0) {
			BB = BB - HMAT.col(j) * BC(j, 1);
		}

	}


	int mm = 0; int p = 3;

	for (int e = 0; e < NE; e++) {

		for (int i = 0; i < 4; i++) {

			if (t_bc(e, i) == 1) {
				BB = BB + t_bc(e, 4 + i) * GMAT.col(4 * e + i);
			}

			else if (t_bc(e, i) == 0) {

				h = -1;
				h = Find(h, e, i, record, NConnect, NE);

				if (h == -1) {

					record[mm][0] = t;
					record[mm][1] = NConnect(e, i);
					mm = mm + 1;

					HMAT.col(t) = -GMAT.col(4 * e + i);
					t = t + 1;
				}

				else {
					int cc = record[h][0];
					HMAT.col(cc) = HMAT.col(cc) - GMAT.col(4 * e + i);
				}
			}
		}
	}

	config.GMAT.resize(0, 0);


};
void BEMBuilderPotential3D::getUTVector(Config & config_) {

	configPotential3D_elliptical& config = dynamic_cast<configPotential3D_elliptical&> (config_);

	VectorXd& U = config.U;
	VectorXd& T = config.T;
	MatrixXd& BC = config.BC;
	MatrixXd& t_bc = config.t_bc;
	VectorXd& XX = config.X;
	MatrixXi& NConnect = config.NConnect;

	int& NN = config.NN;
	int& num = config.num;
	int& NE = config.NE;

	int t = 0;
	int h = 0;
	int tt = 0;

	for (int i = 0; i < NN; i++) {
		if (BC(i, 0) == 0) {
			U[i] = BC(i, 1);
		}
		else
		{
			U[i] = XX(t);
			t = t + 1;
		}
	}


	for (int i = NN; i < NN + 12 * num; i++) {
		U[i] = XX(i);

	}

	int m = 0;

	for (int e = 0; e < NE; e++) {
		for (int i = 0; i < 4; i++) {

			if (t_bc(e, i) == 1) {
				T[m] = t_bc(e, 4 + i);
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
	cout << "get U" << endl;
	cout << config.U << endl;
	cout << "get T" << endl;
	cout << config.T << endl;
	*/

	//ofstream testH1, testH2, testH3, testH4, testH5, testH6, testH7, testH8, testH9;
	//ofstream testG1, testG2, testG3, testG4, testG5, testG6, testG7, testG8, testG9;

	//testH9.open("testH9.txt");
	//testG9.open("testG9.txt");

	//for (int i = 0; i < NN + 12 * num; i++) {
	//	testH9 << U(i);
	//	testH9 << "\n";
	//}

	//for (int i = 0; i < 4 * NE; i++) {
	//	testG9 << T(i);
	//	testG9 << "\n";
	//}

	//testH9.close();
	//testG9.close();

};

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