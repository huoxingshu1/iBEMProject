

#include "postProcessorElastic3D_arbitrary.h"
#include "integratorElastic3DHelper.h"
#include <cmath>
#include <iostream>
#include <fstream>

using namespace Eigen;

void postProcessorElastic3D_arbitrary::get_for_stress_BEM_post(double* Vnorm, double* point, double* GCcor, double** Dst, double** Dsu)
{

	int t = 0;

	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 3; j++) {
			Dst[i][j] = 0.0;
			Dsu[i][j] = 0.0;
		}
	}

	double* x = new double[3];

	x[0] = GCcor[0] - point[0];
	x[1] = GCcor[1] - point[1];
	x[2] = GCcor[2] - point[2];


	double* TG1 = new double[3];
	double* TG = new double[3 * 3];
	double* sgm = new double[3 * 3];


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

					Dsu[t][n] = Dsu[t][n] + ((2.0 * mu_0 * nu0 / (1 - 2 * nu0)) * d[k][l] * sgm[3 * m + n] + (mu_0) * (2.0 * mu_0 * nu0 / (1 - 2 * nu0) * d[m][n] * (TG[3 * k + l] + TG[3 * l + k]) + mu_0 * Green_02(x, m, k, n, l) + mu_0 * Green_02(x, m, l, n, k) + mu_0 * Green_02(x, n, k, m, l) + mu_0 * Green_02(x, n, l, m, k))) * Vnorm[m];
				}


				Dst[t][n] = Dst[t][n] - ((2.0 * mu_0 * nu0 / (1 - 2 * nu0)) * d[k][l] * TG1[n] + (mu_0)*Green_01(x, n, k, l) + (mu_0)*Green_01(x, n, l, k));

			}

			t = t + 1;
		}
	}



	delete[] TG;
	delete[] TG1;
	delete[] sgm;
	delete[] x;

}

void postProcessorElastic3D_arbitrary::Integ3_stress_post(double** Elcor, double** DDsup, double** DDstp, double* point)
{
	int i, m, k, ngp;

	double xsi, eta, W, Jac, N1, N2, N3, N4;

	//double *dxr, **UP ,**TP, *N1, *N2, *N3, *N4, *N, *Vnorm, *GCcor;

	double** Dsu, ** Dst, ** Ni;

	//  double  **TDDsu, **TDDst;

	int Ndof = 3;


	double* N = new double[4];
	double* Vnorm = new double[3];
	double* GCcor = new double[3];
	double* gp = new double[8];
	double* w = new double[8];

	Dsu = new double* [6];
	Dst = new double* [6];
	Ni = new double* [3];


	for (i = 0; i < 6; i++) { Dsu[i] = new double[3];	Dst[i] = new double[3]; }


	for (i = 0; i < 3; i++) { Ni[i] = new double[12]; }



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

			get_for_stress_BEM_post(Vnorm, point, GCcor, Dst, Dsu);


			for (int mm = 0; mm < 6; mm++) {
				for (int ll = 0; ll < 12; ll++) {
					for (int nn = 0; nn < 3; nn++) {
						DDsup[mm][ll] = DDsup[mm][ll] + Dsu[mm][nn] * Ni[nn][ll] * W * Jac;
						DDstp[mm][ll] = DDstp[mm][ll] + Dst[mm][nn] * Ni[nn][ll] * W * Jac;
					}
				}
			}


		}
	}



	delete[] N;
	delete[] Vnorm;
	delete[] GCcor;
	delete[] gp;
	delete[] w;


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

void postProcessorElastic3D_arbitrary::post_stress_BEM(int nump, int NE, int NN, Ref<MatrixXd> x_p, Ref<MatrixXi> NConnect, Ref<MatrixXd> NODES,\
	double** s_pb, Ref<VectorXd> T, Ref<VectorXd> U)
{



# pragma omp parallel shared(nump,NE,NN,x_p,NConnect,NODES,s_pb,T,U)
	{

		double** DDsup, ** DDstp;
		double** GSUe, ** GSTe;

		double* point = new double[3];

		DDsup = new double* [6];
		DDstp = new double* [6];
		GSUe = new double* [6];
		GSTe = new double* [6];



		double** Elcor;
		Elcor = new double* [4];
		for (int j = 0; j < 4; j++) { Elcor[j] = new double[3]; }


		for (int i = 0; i < 6; i++) {
			DDsup[i] = new double[12];  DDstp[i] = new double[12];
			GSUe[i] = new double[12 * NE]; GSTe[i] = new double[3 * NN];
		}


# pragma omp for

		for (int KK = 0; KK < nump; KK++) {


			for (int i = 0; i < 6; i++) {
				for (int j = 0; j < 3 * NN; j++) {
					GSTe[i][j] = 0.0;
				}

				for (int j = 0; j < 12 * NE; j++) {
					GSUe[i][j] = 0.0;
				}
			}


			for (int ee = 0; ee < NE; ee++) {

				for (int mm = 0; mm < 6; mm++) {
					for (int ll = 0; ll < 12; ll++) {
						DDsup[mm][ll] = 0.0;
						DDstp[mm][ll] = 0.0;
					}
				}

				point[0] = x_p(KK, 0);  point[1] = x_p(KK, 1); point[2] = x_p(KK, 2);


				int* sctrb = new int[4];
				sctrb[0] = NConnect(ee, 0);
				sctrb[1] = NConnect(ee, 1);
				sctrb[2] = NConnect(ee, 2);
				sctrb[3] = NConnect(ee, 3);



				for (int j = 0; j < 4; j++) {
					for (int k = 0; k < 3; k++) {
						Elcor[j][k] = NODES((sctrb[j] - 1), k);
					}
				}




				Integ3_stress_post(Elcor, DDsup, DDstp, point);



# pragma omp critical
				{
					for (int j = 0; j < 4; j++)
					{
						int p = 0;
						p = 3 * (sctrb[j] - 1);

						for (int k = 0; k < 3; k++) {
							for (int l = 0; l < 6; l++) {
								GSTe[l][p + k] = GSTe[l][p + k] + DDsup[l][3 * j + k];
							}
						}
					}

				}



				for (int l = 0; l < 6; l++) {
					for (int k = 0; k < 12; k++) {
						GSUe[l][12 * ee + k] = DDstp[l][k];
					}
				}

			}


			for (int ii = 0; ii < 6; ii++) {
				for (int mm = 0; mm < 3 * NN; mm++) {
					s_pb[KK][ii] = s_pb[KK][ii] + GSTe[ii][mm] * U(mm);
				}

				for (int mm = 0; mm < 12 * NE; mm++) {
					s_pb[KK][ii] = s_pb[KK][ii] + GSUe[ii][mm] * T(mm);
				}

			}

			//printf("KK= %d", KK);



		}


		for (int j = 0; j < 4; j++) { delete[] Elcor[j]; }

		for (int i = 0; i < 6; i++) {
			delete[] DDsup[i];  delete[] DDstp[i];
			delete[] GSUe[i];   delete[] GSTe[i];
		}

		delete[] Elcor;
		delete[] DDsup;
		delete[] DDstp;
		delete[] GSUe;
		delete[] GSTe;
	}



}


void postProcessorElastic3D_arbitrary::postprocessor(int nump, int num, int num_NE, double** s_p, double** s_pb, double** e_p, \
	int** index_E_ij, Ref<MatrixXd> x_o, Ref<MatrixXi> eigen_connect, Ref<MatrixXd> x_p, double* Eigen, \
	Ref<VectorXd> U, Ref<VectorXd> T, int NE, int NN, Ref<MatrixXi> NConnect, Ref<MatrixXd> NODES)
{

	post_stress_BEM(nump, NE, NN, x_p, NConnect, NODES, s_pb, T, U);		// obtain the BEM part


	int* test = new int[nump];
	for (int i = 0; i < nump; i++) {
		test[i] = -1;
	}

	// create a point to judge
	double xxx[3];
	double** vertex = new double* [4];

	for (int mm = 0; mm < 4; mm++) {
		vertex[mm] = new double[3];
	}

	int flag_tet;

	for (int s = 0; s < nump; s++)
	{
		xxx[0] = x_p(s, 0); xxx[1] = x_p(s, 1); xxx[2] = x_p(s, 2);
		flag_tet = -1;
		for (int h = 0; h < num_NE; h++)
		{

			for (int mm = 0; mm < 4; mm++) {
				vertex[mm][0] = x_o(eigen_connect(h, mm) - 1, 0);
				vertex[mm][1] = x_o(eigen_connect(h, mm) - 1, 1);
				vertex[mm][2] = x_o(eigen_connect(h, mm) - 1, 2);
			}

			// assign the flag based on its inside or outside
			flag_tet = INT_HELP.judge_tet_point(xxx, vertex);

			if (flag_tet == 0) {	// inside
				test[s] = h;
			}

		}
	}


# pragma omp parallel shared(e_p)

	{

		double* e_p0 = new double[6];

		double* e_pe0 = new double[6];


		int i, j, p, s, h, x;
		int m, n;

		int judge_flag, flag_ppp;
		double** coord = new double* [10];
		for (int pp = 0; pp < 10; pp++) {
			coord[pp] = new double[3];
		}

		double* xx = new double[3];
		double* SH = new double[4];

		int xsf, xsq, pp;

# pragma omp for schedule(dynamic)

		for (s = 0; s < nump; s++) {
			xx[0] = x_p(s, 0); xx[1] = x_p(s, 1); xx[2] = x_p(s, 2);
			for (x = 0; x < 6; x++) {
				e_pe0[x] = 0.0;
				e_p0[x] = 0.0;
			}


			if (test[s] != -1)
			{
				// assign the coordinate
				for (xsf = 0; xsf < 10; xsf++) {
					coord[xsf][0] = x_o(eigen_connect(test[s], xsf) - 1, 0);
					coord[xsf][1] = x_o(eigen_connect(test[s], xsf) - 1, 1);                          // assign the coordinates for one triangle element
					coord[xsf][2] = x_o(eigen_connect(test[s], xsf) - 1, 2);
				}

				for (i = 0; i < 3; i++)
				{
					for (j = i; j < 3; j++)
					{
						for (pp = 0; pp < 10; pp++) {

							double shape_L = INT_HELP.tet_shape_10(coord, pp, xx);
							e_pe0[index_E_ij[i][j]] += Eigen[index_E_ij[3 * (eigen_connect(test[s], pp) - 1) + i][j]] * shape_L;

						}

					}
				}
			}



			for (h = 0; h < num_NE; h++) {
				for (xsf = 0; xsf < 10; xsf++) {
					coord[xsf][0] = x_o(eigen_connect(h, xsf) - 1, 0);
					coord[xsf][1] = x_o(eigen_connect(h, xsf) - 1, 1);                          // assign the coordinates for one triangle element
					coord[xsf][2] = x_o(eigen_connect(h, xsf) - 1, 2);
				}
				judge_flag = 0;
				for (int ppp = 0; ppp < 10; ppp++) {
					if (x_p(s, 0) == x_o(eigen_connect(h, ppp) - 1, 0) && x_p(s, 1) == x_o(eigen_connect(h, ppp) - 1, 1) && x_p(s, 2) == x_o(eigen_connect(h, ppp) - 1, 2)) {
						judge_flag = 1;
						flag_ppp = ppp;
					}
				}

				if (judge_flag == 0 || judge_flag != 0) {


					for (i = 0; i < 3; i++) {
						for (j = i; j < 3; j++) {
							double result[10] = {};
							for (m = 0; m < 3; m++) {

								for (n = m; n < 3; n++) {

									INT_A.Polyhedral_Eshelby_strain_all_10(i, j, m, n, coord, xx, result);
									for (int nnn = 0; nnn < 10; nnn++)
										e_p0[index_E_ij[i][j]] = e_p0[index_E_ij[i][j]] + map[m][n] * Eigen[index_E_ij[3 * (eigen_connect(h, nnn) - 1) + m][n]] * result[nnn];

								}
							}
						}
					}




				} // end of if judge_flag

			}


			for (i = 0; i < 3; i++) {
				for (j = i; j < 3; j++) {
					e_p[s][index_E_ij[i][j]] = e_p0[index_E_ij[i][j]] - e_pe0[index_E_ij[i][j]];  //strain of interested point//
				}
			}



			printf("s= %d", s);

		}


		delete[] e_p0;
		delete[] SH;
		delete[] xx;
		delete[] e_pe0;


	}


	for (int s = 0; s < nump; s++) {

		for (int i = 0; i < 3; i++) {
			for (int j = i; j < 3; j++) {
				for (int m = 0; m < 3; m++) {
					for (int n = m; n < 3; n++) {
						s_p[s][index_E_ij[i][j]] = s_p[s][index_E_ij[i][j]] + Constitutive(i, j, m, n) * map[m][n] * e_p[s][index_E_ij[m][n]];
					}
				}
			}
		}


	}




}

void postProcessorElastic3D_arbitrary::getField(Config& config_) {}

void postProcessorElastic3D_arbitrary::getFlux(Config& config_)
{
	configElastic3D_arbitrary& config = dynamic_cast<configElastic3D_arbitrary&> (config_);
	int& nump = config.nump;
	int& num_NN = config.num_NN; int& num_NE = config.num_NE;
	int& NE = config.NE;
	int& NN = config.NN;
	Ref<VectorXd> U = config.U;
	Ref<VectorXd> T = config.T;
	int**& index_E_ij = config.index_E_ij;
	Ref<MatrixXd>NODES = config.NODES;
	Ref<MatrixXd>eigen_point = config.eigen_point;
	Ref<MatrixXi>eigen_connect = config.eigen_connect;
	Ref<MatrixXd> x_p = config.x_p;
	Ref<MatrixXi> NConnect = config.NConnect;


	double** s_p, ** s_pb, ** s_pt, ** e_p;
	s_p = new double* [nump];
	s_pb = new double* [nump];
	s_pt = new double* [nump];
	e_p = new double* [nump];

	double* Eigen = new double[6 * num_NN];
	for (int i = 0; i < nump; i++) {
		s_p[i] = new double[6];  s_pb[i] = new double[6];  
		s_pt[i] = new double[6];  e_p[i] = new double[6];
	}

	for (int i = 0; i < nump; i++) {
		for (int j = 0; j < 6; j++) {
			s_p[i][j] = 0.0; e_p[i][j] = 0.0; 
			s_pb[i][j] = 0.0; s_pt[i][j] = 0.0; e_p[i][j] = 0.0; 
		}
	}

	for (int i = 0; i < 6 * num_NN; i++) {
		Eigen[i] = U[i + 3 * NN];
	}

	postprocessor(nump, num_NN, num_NE, s_p, s_pb, e_p, index_E_ij, eigen_point, eigen_connect, x_p, Eigen, U, T, NE, NN, NConnect, NODES);

	for (int KK = 0; KK < nump; KK++) {
		for (int ii = 0; ii < 6; ii++) {
			s_pt[KK][ii] = s_p[KK][ii] + s_pb[KK][ii];
		}
	}

	config.s_pt = s_pt;

	ofstream myfile7;
	myfile7.open("stress_3D.txt");

	for (int KK = 0; KK < nump; KK++) {
		myfile7 << std::scientific << config.s_pt[KK][0] << " " << config.s_pt[KK][1] << " " << config.s_pt[KK][2] << config.s_pt[KK][3] << config.s_pt[KK][4] << config.s_pt[KK][5];
		myfile7 << "\n";
	}
	myfile7.close();

}
