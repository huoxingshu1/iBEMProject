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

# include "postProcessorElastic3D.h"
# include "GlobalConstant.h"
# include "integratorElastic3DHelper.h"
# include "BEMElasticHelper3D.h"

using namespace Elastic3DSpherical;

inline void get_for_stress_BEM_post(double* Vnorm, double* point, double* GCcor, double** Dst, double** Dsu)
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

void Integ3_stress_post(double** Elcor, double** DDsup, double** DDstp, double* point)
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



inline void post_stress_BEM(int nump, int NE, int NN, Ref<MatrixXd> x_p, Ref<MatrixXi> NConnect, Ref<MatrixXd> NODES, double** s_pb, Ref<VectorXd> T, Ref<VectorXd> U)
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

			printf("KK= %d", KK);



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

void postprocessor(int nump, int num, double** s_p, double** s_pb, double** e_p, int** index_E_ij, int*** index_E_ijk, int**** index_E_ijkl, Ref<MatrixXd> x_o, Ref<VectorXd> radius, Ref<MatrixXd> x_p, double* Eigen, Ref<VectorXd> U, Ref<VectorXd> T, int NE, int NN, Ref<MatrixXi> NConnect, Ref<MatrixXd> NODES)
{


	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////


	post_stress_BEM(nump, NE, NN, x_p, NConnect, NODES, s_pb, T, U);


	int* test = new int[nump];
	for (int i = 0; i < nump; i++) {
		test[i] = -1;
	}

	for (int s = 0; s < nump; s++)
	{
		for (int h = 0; h < num; h++)
		{
			if (((x_p(s, 0) - x_o(h, 0)) * (x_p(s, 0) - x_o(h, 0)) + (x_p(s, 1) - x_o(h, 1)) * (x_p(s, 1) - x_o(h, 1)) + (x_p(s, 2) - x_o(h, 2)) * (x_p(s, 2) - x_o(h, 2)) - radius[h] * radius[h]) <= 0)
				test[s] = h;
		}
	}


# pragma omp parallel shared(x_p, x_o, test, e_p, s_p, index_E_ij, index_E_ijk , index_E_ijkl , Eigen, nump , radius)

	{

		double* e_p0 = new double[6];
		double* e_p1 = new double[6];
		double* e_p2 = new double[6];
		double* e_pe0 = new double[6];
		double* e_pe1 = new double[6];
		double* e_pe2 = new double[6];

		int i, j, p, q, s, h, x;
		int m, n;

# pragma omp for 

		for (s = 0; s < nump; s++) {

			for (x = 0; x < 6; x++) {
				e_pe0[x] = 0.0;
				e_pe1[x] = 0.0;
				e_pe2[x] = 0.0;
				e_p0[x] = 0.0;
				e_p1[x] = 0.0;
				e_p2[x] = 0.0;
			}

			if (test[s] != -1)
			{
				for (i = 0; i < 3; i++)
				{
					for (j = i; j < 3; j++)
					{
						e_pe0[index_E_ij[i][j]] = Eigen[index_E_ij[3 * test[s] + i][j]];
						for (int k = 0; k < 3; k++)
						{
							e_pe1[index_E_ij[i][j]] = e_pe1[index_E_ij[i][j]] + Eigen[index_E_ijk[3 * test[s] + i][j][k]] * (x_p(s, k) - x_o(test[s], k));
						}
						for (int k = 0; k < 3; k++)
						{
							for (int l = k; l < 3; l++)
							{
								e_pe2[index_E_ij[i][j]] = e_pe2[index_E_ij[i][j]] + map[k][l] * Eigen[index_E_ijkl[3 * test[s] + i][j][k][l]] * (x_p(s, k) - x_o(test[s], k)) * (x_p(s, l) - x_o(test[s], l));
							}
						}
					}
				}
			}



			for (h = 0; h < num; h++) {
				for (i = 0; i < 3; i++) {
					for (j = i; j < 3; j++) {
						for (m = 0; m < 3; m++) {

							for (n = m; n < 3; n++) {


								e_p0[index_E_ij[i][j]] = e_p0[index_E_ij[i][j]] + D_40(i, j, m, n, s, h, radius, x_o, x_p) * map[m][n] * Eigen[index_E_ij[3 * h + m][n]];


								for (p = 0; p < 3; p++) {

									e_p1[index_E_ij[i][j]] = e_p1[index_E_ij[i][j]] + D_50(i, j, m, n, p, s, h, radius, x_o, x_p) * map[m][n] * Eigen[index_E_ijk[3 * h + m][n][p]];

								}

								for (p = 0; p < 3; p++) {
									for (q = 0; q < 3; q++) {

										e_p2[index_E_ij[i][j]] = e_p2[index_E_ij[i][j]] + D_60(i, j, m, n, p, q, s, h, radius, x_o, x_p) * map[m][n] * Eigen[index_E_ijkl[3 * h + m][n][p][q]];

									}
								}

							}
						}
					}
				}
			}


			for (i = 0; i < 3; i++) {
				for (j = i; j < 3; j++) {
					e_p[s][index_E_ij[i][j]] = e_p0[index_E_ij[i][j]] + e_p1[index_E_ij[i][j]] + e_p2[index_E_ij[i][j]]
						- e_pe0[index_E_ij[i][j]] - e_pe1[index_E_ij[i][j]] - e_pe2[index_E_ij[i][j]];  //strain of interested point//
				}
			}



			printf("s= %d", s);

		}


		delete[] e_p0;
		delete[] e_p1;
		delete[] e_p2;
		delete[] e_pe0;
		delete[] e_pe1;
		delete[] e_pe2;

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




	//	for (j = 0; j < 6; j++){
	//	     printf("%e\n", s_p[0][j]);
	//	}

}


void postProcessorElastic3D::getFlux (Config& config_) {
	
	configElastic3D& config = dynamic_cast<configElastic3D&> (config_);
	int& nump = config.nump;
	int& num = config.num;
	int& NE = config.NE;
	int& NN = config.NN;
	Ref<VectorXd> U = config.U;
	Ref<VectorXd> T = config.T;
	int**& index_E_ij = config.index_E_ij;
	int*** &index_E_ijk = config.index_E_ijk;
	int****& index_E_ijkl = config.index_E_ijkl;
	Ref<MatrixXd>NODES = config.NODES;
	Ref<MatrixXd>eigen_point = config.eigen_point;
	Ref<VectorXd> radius = config.radius;
	Ref<MatrixXd> x_p = config.x_p;
	Ref<MatrixXi> NConnect = config.NConnect;


	double** s_p, ** s_pb, ** s_pt, ** e_p;
	s_p = new double* [nump];
	s_pb = new double* [nump];
	s_pt = new double* [nump];
	e_p = new double* [nump];


	double* Eigen = new double[78 * num];


	for (int i = 0; i < nump; i++) {
		s_p[i] = new double[6];  s_pb[i] = new double[6];  s_pt[i] = new double[6];  e_p[i] = new double[6];
	}

	for (int i = 0; i < nump; i++) {
		for (int j = 0; j < 6; j++) {
			s_p[i][j] = 0.0; s_pb[i][j] = 0.0; s_pt[i][j] = 0.0; e_p[i][j] = 0.0;
		}
	}

	for (int i = 0; i < 78 * num; i++) {
		Eigen[i] = U[i + 3 * NN];
	}

	postprocessor(nump, num, s_p, s_pb, e_p, index_E_ij, index_E_ijk, index_E_ijkl, eigen_point, radius, x_p, Eigen, U, T, NE, NN, NConnect, NODES);

	for (int KK = 0; KK < nump; KK++) {
		for (int ii = 0; ii < 6; ii++) {
			s_pt[KK][ii] = s_p[KK][ii] + s_pb[KK][ii];
		}
	}

	config.s_pt = s_pt;

	ofstream myfile7;
	myfile7.open("stress.txt");

	for (int KK = 0; KK < nump; KK++) {
		myfile7 << std::scientific << config.s_pt[KK][0] << " " << config.s_pt[KK][1] << " " << config.s_pt[KK][2] << " " << config.s_pt[KK][3] << " " << config.s_pt[KK][4] << " " << config.s_pt[KK][5];
		myfile7 << "\n";
	}
	myfile7.close();

};


void postProcessorElastic3D::getField(Config& config_) {

};