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

# include "postProcessorElastic2D.h"
# include "GlobalConstant.h"
# include "integratorElastic2DHelper.h"
# include "BEMElasticHelper2D.h"

void postProcessorElastic2D::getFlux(Config& config_) {

	configElastic2D& config = dynamic_cast<configElastic2D&> (config_);
	int& nump = config.nump;
	int& num = config.num;
	int& NE = config.NE;
	int& NN = config.NN;
	Ref<VectorXd> U = config.U;
	Ref<VectorXd> T = config.T;
	int**& index_E_ij = config.index_E_ij;
	int***& index_E_ijk = config.index_E_ijk;
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


	double* Eigen = new double[21 * num];


	for (int i = 0; i < nump; i++) {
		s_p[i] = new double[3];  s_pb[i] = new double[3];  s_pt[i] = new double[3];  e_p[i] = new double[3];
	}

	for (int i = 0; i < nump; i++) {
		for (int j = 0; j < 3; j++) {
			s_p[i][j] = 0.0; s_pb[i][j] = 0.0; s_pt[i][j] = 0.0; e_p[i][j] = 0.0;
		}
	}

	for (int i = 0; i < 21 * num; i++) {
		Eigen[i] = U[i + 2 * NN];
	}

	postprocessor(nump, num, s_p, s_pb, e_p, index_E_ij, index_E_ijk, index_E_ijkl, eigen_point, radius, x_p, Eigen, U, T, NE, NN, NConnect, NODES);

	for (int KK = 0; KK < nump; KK++) {
		for (int ii = 0; ii < 3; ii++) {
			s_pt[KK][ii] = s_p[KK][ii] + s_pb[KK][ii];
		}
	}

	config.s_pt = s_pt;

	ofstream myfile7;
	myfile7.open("stress_2D.txt");

	for (int KK = 0; KK < nump; KK++) {
        myfile7 << std::scientific << config.s_pt[KK][0] << " " << config.s_pt[KK][1] << " " << config.s_pt[KK][2];
		myfile7 << "\n";
	}
	myfile7.close();

};


void postProcessorElastic2D::getField(Config& config_) {

};

void postProcessorElastic2D::postprocessor(int nump, int num, double** s_p, double** s_pb, double** e_p, int** index_E_ij, int*** index_E_ijk, int**** index_E_ijkl, Ref<MatrixXd> x_o, Ref<VectorXd> radius, Ref<MatrixXd> x_p, double* Eigen, Ref<VectorXd> U \
    , Ref<VectorXd> T, int NE, int NN, Ref<MatrixXi> NConnect, Ref<MatrixXd> NODES)
{
    post_bem_stress(nump, NE, NN, x_p, NConnect, NODES, s_pb, T, U);				// obtain the stress/strain field through the BEM solution part
                                                                        // s_pb is the stress by BEM part
    int* test = new int[nump];
    for (int i = 0; i < nump; i++) {
        test[i] = -1;
    }

    for (int s = 0; s < nump; s++) {
        for (int h = 0; h < num; h++) {

            if (((x_p(s, 0) - x_o(h, 0)) * (x_p(s, 0) - x_o(h, 0)) + (x_p(s, 1) - x_o(h, 1)) * (x_p(s, 1) - x_o(h, 1)) - radius[h] * radius[h]) <= 0)
                test[s] = h;                        // make sure whether it is inside a particle
        }
    }

# pragma omp parallel shared(x_p, x_o, test, e_p, s_p, index_E_ij, Eigen, nump , radius)
    {

        double* e_p0 = new double[3];
        double* e_p1 = new double[3];                                           // results from the lineat eigen field
        double* e_p2 = new double[3];
        double* e_pe0 = new double[3];
        double* e_pe1 = new double[3];                                          // assign the linear eigen field
        double* e_pe2 = new double[3];
        int i, j, p, q, s, h, x, m, n;

# pragma omp for

        for (s = 0; s < nump; s++) {

            for (x = 0; x < 3; x++) {
                e_pe0[x] = 0.0;
                e_pe1[x] = 0.0;
                e_p0[x] = 0.0;
                e_p1[x] = 0.0;
                e_pe2[x] = 0.0;
                e_p2[x] = 0.0;
            }

            if (test[s] != -1) {                   // interior points
                for (i = 0; i < 2; i++) {
                    for (j = i; j < 2; j++) {
                        e_pe0[index_E_ij[i][j]] = Eigen[index_E_ij[2 * test[s] + i][j]];        // assign the eigen strain field to e_pe0, minus the non-mechanical eigen strain
                        for (int k = 0; k < 2; k++)
                        {
                            e_pe1[index_E_ij[i][j]] = e_pe1[index_E_ij[i][j]] + Eigen[index_E_ijk[2 * test[s] + i][j][k]] * (x_p(s, k) - x_o(test[s], k));
                        }

                        for (int k = 0; k < 2; k++)
                        {
                            for (int l = k; l < 2; l++)
                            {
                                e_pe2[index_E_ij[i][j]] = e_pe2[index_E_ij[i][j]] + map[k][l] * Eigen[index_E_ijkl[2 * test[s] + i][j][k][l]] * (x_p(s, k) - x_o(test[s], k)) * (x_p(s, l) - x_o(test[s], l));
                            }
                        }


                    }

                }

            }                                       // end if test[s]!=-1

            for (h = 0; h < num; h++) {
                for (i = 0; i < 2; i++) {
                    for (j = i; j < 2; j++) {
                        for (m = 0; m < 2; m++) {
                            for (n = m; n < 2; n++) {
                                e_p0[index_E_ij[i][j]] = e_p0[index_E_ij[i][j]] + map[m][n] * Eigen[index_E_ij[2 * h + m][n]] * Eshelby_strain_integration(i, j, m, n, s, h, radius, x_p, x_o);
                                // the Eshelby's tensor times eigen strain obtain the disturbed strain

                                for (p = 0; p < 2; p++) {

                                    e_p1[index_E_ij[i][j]] = e_p1[index_E_ij[i][j]] + Eshelby_strain_integration_linear(i, j, m, n, p, s, h, radius, x_p, x_o) * map[m][n] * Eigen[index_E_ijk[2 * h + m][n][p]];

                                }

                                for (p = 0; p < 2; p++) {
                                    for (q = 0; q < 2; q++) {

                                        e_p2[index_E_ij[i][j]] = e_p2[index_E_ij[i][j]] + Eshelby_strain_integration_quadratic(i, j, m, n, p, q, s, h, radius, x_p, x_o) * map[m][n] * Eigen[index_E_ijkl[2 * h + m][n][p][q]];

                                    }
                                }

                            }
                        }
                    }
                }
            }

            for (i = 0; i < 2; i++) {
                for (j = i; j < 2; j++) {
                    e_p[s][index_E_ij[i][j]] = e_p0[index_E_ij[i][j]] + e_p1[index_E_ij[i][j]] + e_p2[index_E_ij[i][j]] - e_pe0[index_E_ij[i][j]] - e_pe1[index_E_ij[i][j]] - e_pe2[index_E_ij[i][j]];
                }
            }                                       // math: epsilon_{kl} - epsilon_{kl}^*, mechanical strain caused by the eigen strain field (haven't plus the BEM strain part)



            printf("s=%d", s);                      // print the results
        }


        delete[] e_p0;
        delete[] e_pe0;
        delete[] e_p1;
        delete[] e_p2;
        delete[] e_pe1;
        delete[] e_pe2;

    } // end of the omp parallel

    for (int s = 0; s < nump; s++) {
        for (int i = 0; i < 2; i++) {
            for (int j = i; j < 2; j++) {
                for (int m = 0; m < 2; m++) {
                    for (int n = m; n < 2; n++) {
                        s_p[s][index_E_ij[i][j]] = s_p[s][index_E_ij[i][j]] + cons_law(nu0, mu_0, i, j, m, n) * map[m][n] * e_p[s][index_E_ij[m][n]];
                    }
                }

            }
        }
    }
    // s_p is the final stress, e_p is the final strain (superposition iwht BEM and Inclusions)


}

void postProcessorElastic2D::get_for_stress_BEM_post(double* Vnorm, double* point, double* GCcor, double** Dst, double** Dsu)
{
    int t = 0;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            Dst[i][j] = 0.0;
            Dsu[i][j] = 0.0;
        }
    }

    double* x = new double[2];
    x[0] = GCcor[0] - point[0];
    x[1] = GCcor[1] - point[1];

    double* TG1 = new double[2];
    double* TG = new double[2 * 2];
    double* sgm = new double[2 * 2];

    for (int i = 0; i < 2; i++) { TG1[i] = 0.0; }
    for (int i = 0; i < 2 * 2; i++) { TG[i] = 0.0; }
    for (int i = 0; i < 2 * 2; i++) { sgm[i] = 0.0; }

    for (int n = 0; n < 2; n++) {
        for (int i = 0; i < 2; i++) {
            TG1[n] = TG1[n] + B2D.Green_01(x, n, i, i);
        }
    }

    int tt = 0;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int l = 0; l < 2; l++) {
                TG[tt] = TG[tt] + B2D.Green_02(x, l, i, l, j);
            }
            tt = tt + 1;
        }
    }

    tt = 0;
    for (int m = 0; m < 2; m++) {
        for (int n = 0; n < 2; n++) {
            for (int i = 0; i < 2; i++) {
                sgm[tt] = sgm[tt] + 2.0 * mu_0 * nu0 / (1 - 2 * nu0) * d[m][n] * TG[2 * i + i] + 2.0 * mu_0 * B2D.Green_02(x, m, i, n, i);
            }
            tt = tt + 1;
        }
    }

    for (int k = 0; k < 2; k++) {
        for (int l = k; l < 2; l++) {
            for (int n = 0; n < 2; n++) {

                for (int m = 0; m < 2; m++) {

                    Dst[t][n] = Dst[t][n] + ((2.0 * mu_0 * nu0 / (1 - 2 * nu0)) * d[k][l] * sgm[2 * m + n] + (mu_0) * (2.0 * mu_0 * nu0 / (1 - 2 * nu0) * d[m][n] * (TG[2 * k + l] + TG[2 * l + k]) + mu_0 * B2D.Green_02(x, m, k, n, l) + mu_0 * B2D.Green_02(x, m, l, n, k) + mu_0 * B2D.Green_02(x, n, k, m, l) + mu_0 * B2D.Green_02(x, n, l, m, k))) * Vnorm[m];
                }


                Dsu[t][n] = Dsu[t][n] - ((2.0 * mu_0 * nu0 / (1 - 2 * nu0)) * d[k][l] * TG1[n] + (mu_0)*B2D.Green_01(x, n, k, l) + (mu_0)*B2D.Green_01(x, n, l, k));

            }

            t = t + 1;
        }
    }
    delete[] TG;
    delete[] TG1;
    delete[] sgm;
    delete[] x;

}

void postProcessorElastic2D::Integ2_stress_post(double** Elcor, double** DDsup, double** DDstp, double* point)
{
      int i, m, k, ngp;
    double xsi, W, Jac, N1, N2, N3;
    double** Dsu, ** Dst, ** Ni;

    int Ndof = 2;

    double* N = new double[3];
    double* Vnorm = new double[2];
    double* GCcor = new double[2];

    ngp = 8;
    double w[8] = { 0.1012285363, 0.1012285363, 0.2223810345, 0.2223810345, 0.3137066459, 0.3137066459, 0.3626837834, 0.3626837834 };
    double gp[8] = { 0.9602898565, -0.9602898565, 0.7966664774, -0.7966664774, 0.5255324099, -0.5255324099, 0.1834346425, -0.1834346425 };

    Dsu = new double* [3];                          // 3 variables in 2D
    Dst = new double* [3];
    Ni = new double* [2];

    for (i = 0; i < 3; i++) { Dsu[i] = new double[2]; Dst[i] = new double[2]; }

    for (i = 0; i < 2; i++) { Ni[i] = new double[6]; }


	for (m = 0; m < ngp; m++) {
		xsi = gp[m]; W = w[m];

		B2D.Serendip_func(xsi, N1, N2, N3, Ni, N);
        B2D.Normal_Jac(xsi, Elcor, Jac, Vnorm);
        B2D.Cartesian(Elcor, N1, N2, N3, GCcor);

		//get_for_stress_BEM_post(Vnorm, point, GCcor, Dst, Dsu);
        get_for_stress_BEM_post(Vnorm, point, GCcor, Dst, Dsu);

		int ll = 0;
		for (int mm = 0; mm < 3; mm++) {
			DDsup[mm][0] = DDsup[mm][0] + Dsu[mm][0] * N[0] * W * Jac;
			DDsup[mm][1] = DDsup[mm][1] + Dsu[mm][1] * N[0] * W * Jac;
			DDsup[mm][2] = DDsup[mm][2] + Dsu[mm][0] * N[1] * W * Jac;
			DDsup[mm][3] = DDsup[mm][3] + Dsu[mm][1] * N[1] * W * Jac;
			DDsup[mm][4] = DDsup[mm][4] + Dsu[mm][0] * N[2] * W * Jac;
			DDsup[mm][5] = DDsup[mm][5] + Dsu[mm][1] * N[2] * W * Jac;


			DDstp[mm][0] = DDstp[mm][0] + Dst[mm][0] * N[0] * W * Jac;
			DDstp[mm][1] = DDstp[mm][1] + Dst[mm][1] * N[0] * W * Jac;
			DDstp[mm][2] = DDstp[mm][2] + Dst[mm][0] * N[1] * W * Jac;
			DDstp[mm][3] = DDstp[mm][3] + Dst[mm][1] * N[1] * W * Jac;
			DDstp[mm][4] = DDstp[mm][4] + Dst[mm][0] * N[2] * W * Jac;
			DDstp[mm][5] = DDstp[mm][5] + Dst[mm][1] * N[2] * W * Jac;

		}
	}
    delete[] N;
    delete[] Vnorm;
    delete[] GCcor;

    for (int i = 0; i < 3; i++) {
        delete Dsu[i];
        delete Dst[i];
    }

    for (int i = 0; i < 2; i++) {
        delete Ni[i];
    }

    delete[] Dsu;
    delete[] Dst;
    delete[] Ni;


}


void postProcessorElastic2D::post_bem_stress(int nump, int NE, int NN, Ref<MatrixXd> x_p, Ref<MatrixXi> NConnect, Ref<MatrixXd> NODES, double** s_pb, Ref<VectorXd> T, Ref<VectorXd> U)
{
# pragma omp parallel shared(nump,NE,NN,x_p,NConnect,NODES,s_pb,T,U)
    {
        double** DDsup, ** DDstp;
        double** GSUe, ** GSTe;

        double* point = new double[2];
        DDsup = new double* [3];
        DDstp = new double* [3];
        GSUe = new double* [3];
        GSTe = new double* [3];

        double** Elcor;
        Elcor = new double* [3];
        for (int j = 0; j < 3; j++) {
            Elcor[j] = new double[2];                           // two-d and 3 nodes elements
        }

        for (int i = 0; i < 3; i++) {
            DDsup[i] = new double[6];
            DDstp[i] = new double[6];
            GSUe[i] = new double[6 * NE];
            GSTe[i] = new double[2 * NN];
        }
        double lambda = (2.0 * mu_0 * nu0) / (1.0 - 2.0 * nu0);

        // end initialization
# pragma omp for

        for (int KK = 0; KK < nump; KK++) {

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 2 * NN; j++) {
                    GSTe[i][j] = 0.0;
                }

                for (int j = 0; j < 6 * NE; j++) {
                    GSUe[i][j] = 0.0;
                }
            }

            for (int ee = 0; ee < NE; ee++) {

                for (int mm = 0; mm < 3; mm++) {
                    for (int ll = 0; ll < 6; ll++) {
                        DDsup[mm][ll] = 0.0;
                        DDstp[mm][ll] = 0.0;
                    }
                }

                point[0] = x_p(KK, 0);  point[1] = x_p(KK, 1);

                int* sctrb = new int[3];
                sctrb[0] = NConnect(ee, 0);
                sctrb[1] = NConnect(ee, 1);
                sctrb[2] = NConnect(ee, 2);

                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 2; k++) {
                        Elcor[j][k] = NODES((sctrb[j] - 1), k);
                    }
                }                                                       // record the element coordinate


                Integ2_stress_post(Elcor, DDsup, DDstp, point);


# pragma omp critical
                {
                    for (int j = 0; j < 3; j++) {

                        int p = 0;
                        p = 2 * (sctrb[j] - 1);

                        for (int k = 0; k < 2; k++) {
                            for (int l = 0; l < 3; l++) {
                                GSTe[l][p + k] = GSTe[l][p + k] + DDstp[l][2 * j + k];
                            }
                        }
                    }

                }                                                           // end critical


                for (int l = 0; l < 3; l++) {
                    for (int k = 0; k < 6; k++) {
                        GSUe[l][6 * ee + k] = DDsup[l][k];
                    }
                }


            }                                                               // end of the ee 

            for (int ii = 0; ii < 3; ii++) {
                for (int mm = 0; mm < 2 * NN; mm++) {
                    s_pb[KK][ii] = s_pb[KK][ii] + GSTe[ii][mm] * U(mm);
                }

                for (int mm = 0; mm < 6 * NE; mm++) {
                    s_pb[KK][ii] = s_pb[KK][ii] + GSUe[ii][mm] * T(mm);
                }
            }                                                              // apply BC to get stress
    

        }

        for (int j = 0; j < 3; j++) { delete[] Elcor[j]; }

        for (int i = 0; i < 3; i++) {
            delete[] DDsup[i];
            delete[] DDstp[i];
            delete[] GSUe[i];
            delete[] GSTe[i];
        }
        delete[] Elcor;
        delete[] DDsup;
        delete[] DDstp;
        delete[] GSUe;
        delete[] GSTe;


    }

}