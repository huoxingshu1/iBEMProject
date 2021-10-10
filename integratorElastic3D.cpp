#include "configElastic3D.h"
#include "integratorElastic3D.h"
#include "integratorElastic3DHelper.h"

void get_addu(double* x, double a, double** DD_u3, int** index)

{

	for (int ii = 0; ii < 3; ii++) {
		for (int jj = 0; jj < 6; jj++) {
			DD_u3[ii][jj] = 0.0;
		}
	}

	for (int i = 0; i < 3; i++) {

		for (int m = 0; m < 3; m++) {
			for (int n = m; n < 3; n++) {

				for (int jj = 0; jj < 3; jj++) {
					for (int ll = 0; ll < 3; ll++) {

						//		DD_u3[i][index[m][n]] = DD_u3[i][index[m][n]] - map[m][n] * Constitutive(jj, ll, m, n)*( 1/(4*pi*mu_0)* d[i][jj]*Phi_01( x, a,ll) - 1/(16*pi*mu_0*(1-nu0))*Psi_03(x,a, i , jj, ll) );
						DD_u3[i][index[m][n]] = DD_u3[i][index[m][n]] - map[m][n] * Constitutive(jj, ll, m, n) * IGreen_01(x, a, i, jj, ll);


					}
				}

				//printf("%.12lf\n", DD_u3[i][index[m][n]]);
			}
		}

	}


}

void get_addu1(double* x, double a, double** DD_u4, int*** indexp)
{

	for (int ii = 0; ii < 3; ii++) {
		for (int jj = 0; jj < 18; jj++) {
			DD_u4[ii][jj] = 0.0;
		}
	}

	for (int i = 0; i < 3; i++) {

		for (int m = 0; m < 3; m++) {
			for (int n = m; n < 3; n++) {
				for (int p = 0; p < 3; p++) {


					for (int jj = 0; jj < 3; jj++) {
						for (int ll = 0; ll < 3; ll++) {

							DD_u4[i][indexp[m][n][p]] = DD_u4[i][indexp[m][n][p]] - map[m][n] * Constitutive(jj, ll, m, n) * IGreen_11(x, a, i, jj, p, ll);

						}
					}


				}
			}
		}

	}


}

void get_addu2(double* x, double a, double** DD_u5, int**** indexpp)

{

	for (int ii = 0; ii < 3; ii++) {
		for (int jj = 0; jj < 54; jj++) {
			DD_u5[ii][jj] = 0.0;
		}
	}

	for (int i = 0; i < 3; i++) {

		for (int m = 0; m < 3; m++) {
			for (int n = m; n < 3; n++) {
				for (int p = 0; p < 3; p++) {
					for (int q = 0; q < 3; q++) {


						for (int jj = 0; jj < 3; jj++) {
							for (int ll = 0; ll < 3; ll++) {

								DD_u5[i][indexpp[m][n][p][q]] = DD_u5[i][indexpp[m][n][p][q]] - map[m][n] * Constitutive(jj, ll, m, n) * IGreen_21(x, a, i, jj, p, q, ll);


							}
						}

					}
				}
			}
		}

	}


}

void integratorElastic3D::addFieldToBEM(Config& config_) {
	configElastic3D& config = dynamic_cast<configElastic3D&> (config_);
	Ref<MatrixXd> HMAT = config.HMAT;
	Ref<MatrixXd> eigen_point = config.eigen_point;
	Ref<VectorXd> radius = config.radius;
	Ref<MatrixXd> NODES = config.NODES;
	int num = config.num;
	int NN = config.NN;

# pragma omp parallel shared (HMAT)
	{
		double* x = new double[3];
		double** DD_u3, ** DD_u4, ** DD_u5;

		int** index;
		int*** indexp;
		int**** indexpp;
		int id = 0;
		int pp;

		index = new int* [3];
		for (int ii = 0; ii < 3; ii++)
			index[ii] = new int[3];


		indexp = new int** [3];
		for (int ii = 0; ii < 3; ii++) {
			indexp[ii] = new int* [3];
			for (int jj = 0; jj < 3; jj++) {
				indexp[ii][jj] = new int[3];
			}
		}

		indexpp = new int*** [3];

		for (int ii = 0; ii < 3; ii++) {
			indexpp[ii] = new int** [3];
			for (int jj = 0; jj < 3; jj++) {
				indexpp[ii][jj] = new int* [3];
				for (int kk = 0; kk < 3; kk++) {
					indexpp[ii][jj][kk] = new int[3];
				}
			}
		}


		for (int ii = 0; ii < 3; ii++) {
			for (int jj = 0; jj < 3; jj++) {
				index[ii][jj] = 0;
			}
		}


		for (int ii = 0; ii < 3; ii++) {
			for (int jj = 0; jj < 3; jj++) {
				for (int kk = 0; kk < 3; kk++) {
					indexp[ii][jj][kk] = 0;
				}
			}
		}


		for (int ii = 0; ii < 3; ii++) {
			for (int jj = 0; jj < 3; jj++) {
				for (int kk = 0; kk < 3; kk++) {
					for (int ll = 0; ll < 3; ll++) {
						indexpp[ii][jj][kk][ll] = 0;
					}
				}
			}
		}




		for (int ii = 0; ii < 3; ii++) {
			for (int jj = ii; jj < 3; jj++) {
				index[ii][jj] = id;
				id = id + 1;
			}
		}

		id = 0;
		for (int ii = 0; ii < 3; ii++) {
			for (int jj = ii; jj < 3; jj++) {
				for (int kk = 0; kk < 3; kk++) {
					indexp[ii][jj][kk] = id;
					id = id + 1;
				}
			}
		}

		id = 0;
		for (int ii = 0; ii < 3; ii++) {
			for (int jj = ii; jj < 3; jj++) {
				for (int kk = 0; kk < 3; kk++) {
					for (int ll = 0; ll < 3; ll++) {
						indexpp[ii][jj][kk][ll] = id;
						id = id + 1;
					}
				}
			}
		}


		DD_u3 = new double* [3];
		DD_u4 = new double* [3];
		DD_u5 = new double* [3];

		for (int ii = 0; ii < 3; ii++) {
			DD_u3[ii] = new double[6];
			DD_u4[ii] = new double[18];
			DD_u5[ii] = new double[54];

		}


		for (int ii = 0; ii < 3; ii++) {
			for (int jj = 0; jj < 6; jj++) {
				DD_u3[ii][jj] = 0.0;
			}

			for (int kk = 0; kk < 18; kk++) {
				DD_u4[ii][kk] = 0.0;
			}

			for (int ll = 0; ll < 54; ll++) {
				DD_u5[ii][ll] = 0.0;
			}
		}


#pragma omp for

		for (int KK = 0; KK < num; KK++) {

			for (int jj = 0; jj < NN; jj++) {

				x[0] = NODES(jj, 0) - eigen_point(KK, 0);
				x[1] = NODES(jj, 1) - eigen_point(KK, 1);
				x[2] = NODES(jj, 2) - eigen_point(KK, 2);

				get_addu(x, radius[KK], DD_u3, index);
				get_addu1(x, radius[KK], DD_u4, indexp);
				get_addu2(x, radius[KK], DD_u5, indexpp);


				pp = 3 * jj;
				int  tt = 78 * KK;


				for (int qq = pp; qq < pp + 3; qq++) {

					for (int mm = 0; mm < 6; mm++) {
						HMAT(qq, 3 * NN + mm + tt) = -DD_u3[qq - pp][mm];
					}


					for (int mm = 0; mm < 18; mm++) {
						HMAT(qq, 3 * NN + mm + tt + 6) = -DD_u4[qq - pp][mm];
					}

					for (int mm = 0; mm < 54; mm++) {
						HMAT(qq, 3 * NN + mm + tt + 24) = -DD_u5[qq - pp][mm];
					}
				}

			}

		}


		for (int i = 0; i < 3; i++) {
			delete DD_u3[i];
			delete DD_u4[i];
			delete DD_u5[i];
		}

		delete[] DD_u3;
		delete[] DD_u4;
		delete[] DD_u5;


		for (int ii = 0; ii < 3; ii++)
		{
			for (int jj = 0; jj < 3; jj++)
			{
				for (int kk = 0; kk < 3; kk++)
				{
					delete[] indexpp[ii][jj][kk];
				}

				delete[] indexp[ii][jj];
			}
			delete[] index[ii];

		}


		delete[] index;
		delete[] indexp;
		delete[] indexpp;
	}


};
void integratorElastic3D::addFluxEquivalentInclusion(Config& config_) {

	configElastic3D& config = dynamic_cast<configElastic3D&> (config_);
	Ref<MatrixXd> HMAT = config.HMAT;
	Ref<MatrixXd> x_o = config.eigen_point;
	Ref<VectorXd> radius = config.radius;
	Ref<MatrixXd> NODES = config.NODES;
	int**& index_E_ij = config.index_E_ij;
	int***& index_E_ijk = config.index_E_ijk;
	int****& index_E_ijkl = config.index_E_ijkl;
	double** A;

	int num = config.num;
	int NN = config.NN;
	A = new double* [78 * num];

	for (int i = 0; i < 78 * num; i++)
		A[i] = new double[78 * num];

	for (int i = 0; i < 78 * num; i++)
		for (int j = 0; j < 78 * num; j++)
			A[i][j] = 0.0;


# pragma omp parallel shared(A,num,radius,x_o,index_E_ij,index_E_ijk,index_E_ijkl)

	{
		int t = 0;
		int m, n, p, q, h;
		int l;
		double sym = 0.0;

		////Assembly Matrix Tensor///////

#	pragma omp for

		for (int s = 0; s < num; s++) {
			for (int i = 0; i < 3; i++) {
				for (int j = i; j < 3; j++) {
					for (h = 0; h < num; h++) {
						for (m = 0; m < 3; m++) {
							for (n = m; n < 3; n++) {
								double D_ijmn = 0.0;
								for (int k = 0; k < 3; k++) {
									for (l = 0; l < 3; l++) {

										D_ijmn = D_ijmn - (Constitutive(i, j, k, l) - Constitutive1(i, j, k, l)) * D_40(k, l, m, n, s, h, radius, x_o, x_o);

									}
								}

								if (h == s)
									sym = 1.0;
								else
									sym = 0.0;

								A[index_E_ij[3 * s + i][j]][index_E_ij[3 * h + m][n]] = map[m][n] * (D_ijmn + sym * Constitutive(i, j, m, n));

								for (p = 0; p < 3; p++) {
									double D_ijmnp = 0.0;
									for (int k = 0; k < 3; k++) {
										for (l = 0; l < 3; l++) {

											D_ijmnp = D_ijmnp - (Constitutive(i, j, k, l) - Constitutive1(i, j, k, l)) * D_50(k, l, m, n, p, s, h, radius, x_o, x_o);

										}
									}
									A[index_E_ij[3 * s + i][j]][index_E_ijk[3 * h + m][n][p]] = map[m][n] * D_ijmnp;
								}

								for (p = 0; p < 3; p++) {
									for (q = 0; q < 3; q++) {
										double D_ijmnpq = 0.0;
										for (int k = 0; k < 3; k++) {
											for (l = 0; l < 3; l++) {

												D_ijmnpq = D_ijmnpq - (Constitutive(i, j, k, l) - Constitutive1(i, j, k, l)) * D_60(k, l, m, n, p, q, s, h, radius, x_o, x_o);

											}
										}
										A[index_E_ij[3 * s + i][j]][index_E_ijkl[3 * h + m][n][p][q]] = map[m][n] * D_ijmnpq;

									}
								}
							}
						}
					}

					//	t = t + 1;

				}
			}


		}

	}

	for (int i = 0; i < 78 * num; i++)
	{
		for (int j = 0; j < 78 * num; j++)
		{
			HMAT(3 * NN + i, 3 * NN + j) = A[i][j];
		}
	}

	for (int i = 0; i < 78 * num; i++) {
		delete A[i];
	}
	delete[] A;



};
void integratorElastic3D::addFluxEquivalentFirstOrderInclusion(Config& config_) {

	configElastic3D& config = dynamic_cast<configElastic3D&> (config_);
	Ref<MatrixXd> HMAT = config.HMAT;
	Ref<MatrixXd> x_o = config.eigen_point;
	Ref<VectorXd> radius = config.radius;
	Ref<MatrixXd> NODES = config.NODES;
	int**& index_E_ij = config.index_E_ij;
	int***& index_E_ijk = config.index_E_ijk;
	int****& index_E_ijkl = config.index_E_ijkl;
	int num = config.num;
	int NN = config.NN;
	double** A;

	A = new double* [78 * num];

	for (int i = 0; i < 78 * num; i++)
		A[i] = new double[78 * num];

	for (int i = 0; i < 78 * num; i++)
		for (int j = 0; j < 78 * num; j++)
			A[i][j] = 0.0;

# pragma omp parallel shared(A,num,radius,x_o,index_E_ij,index_E_ijk,index_E_ijkl)

	{
		int t = 0;
		int m, n, p, q, r, h;
		int l;
		double sym = 0.0;

		////Assembly Matrix Tensor///////

#	pragma omp for

		for (int s = 0; s < num; s++) {

			for (int i = 0; i < 3; i++) {
				for (int j = i; j < 3; j++) {
					for (r = 0; r < 3; r++) {
						for (h = 0; h < num; h++) {
							for (m = 0; m < 3; m++) {
								for (n = m; n < 3; n++) {
									double D_ijmn_r = 0.0;
									for (int k = 0; k < 3; k++) {
										for (l = 0; l < 3; l++) {

											D_ijmn_r = D_ijmn_r - (Constitutive(i, j, k, l) - Constitutive1(i, j, k, l)) * D_41(k, l, m, n, r, s, h, radius, x_o, x_o);

										}
									}

									A[index_E_ijk[3 * s + i][j][r]][index_E_ij[3 * h + m][n]] = map[m][n] * D_ijmn_r;

									if (h == s)
										sym = 1.0;
									else
										sym = 0.0;

									for (p = 0; p < 3; p++) {
										double D_ijmnp_r = 0.0;
										for (int k = 0; k < 3; k++) {
											for (l = 0; l < 3; l++) {

												D_ijmnp_r = D_ijmnp_r - (Constitutive(i, j, k, l) - Constitutive1(i, j, k, l)) * D_51(k, l, m, n, p, r, s, h, radius, x_o, x_o);

											}
										}
										A[index_E_ijk[3 * s + i][j][r]][index_E_ijk[3 * h + m][n][p]] = map[m][n] * (D_ijmnp_r + sym * Constitutive(i, j, m, n) * d[p][r]);
									}

									for (p = 0; p < 3; p++) {
										for (q = 0; q < 3; q++) {
											double D_ijmnpq_r = 0.0;
											for (int k = 0; k < 3; k++) {
												for (l = 0; l < 3; l++) {

													D_ijmnpq_r = D_ijmnpq_r - (Constitutive(i, j, k, l) - Constitutive1(i, j, k, l)) * D_61(k, l, m, n, p, q, r, s, h, radius, x_o, x_o);

												}
											}
											A[index_E_ijk[3 * s + i][j][r]][index_E_ijkl[3 * h + m][n][p][q]] = map[m][n] * D_ijmnpq_r;
										}
									}


								}
							}
						}

						//	t = t + 1;
					}
				}
			}

		}

	}

	for (int i = 0; i < 78 * num; i++)
	{
		for (int j = 0; j < 78 * num; j++)
		{
			HMAT(3 * NN + i, 3 * NN + j) = A[i][j];
		}
	}

	for (int i = 0; i < 78 * num; i++) {
		delete A[i];
	}
	delete[] A;

};
void integratorElastic3D::addFluxEquivalentSecondOrderInclusion(Config& config_) {

	configElastic3D& config = dynamic_cast<configElastic3D&> (config_);
	Ref<MatrixXd> HMAT = config.HMAT;
	Ref<MatrixXd> x_o = config.eigen_point;
	Ref<VectorXd> radius = config.radius;
	Ref<MatrixXd> NODES = config.NODES;
	int**& index_E_ij = config.index_E_ij;
	int***& index_E_ijk = config.index_E_ijk;
	int****& index_E_ijkl = config.index_E_ijkl;
	int num = config.num;
	int NN = config.NN;
	double** A;

	A = new double* [78 * num];

	for (int i = 0; i < 78 * num; i++)
		A[i] = new double[78 * num];

	for (int i = 0; i < 78 * num; i++)
		for (int j = 0; j < 78 * num; j++)
			A[i][j] = 0.0;

#	pragma omp parallel shared(A, num, radius, x_o, index_E_ij, index_E_ijk, index_E_ijkl)

	{
		int t = 0;
		int m, n, p, q, r, z, h;
		int l;
		double sym = 0.0;

		////Assembly Matrix Tensor///////

#	pragma omp for

		for (int s = 0; s < num; s++) {

			for (int i = 0; i < 3; i++) {
				for (int j = i; j < 3; j++) {
					for (r = 0; r < 3; r++) {
						for (z = 0; z < 3; z++) {
							for (h = 0; h < num; h++) {
								for (m = 0; m < 3; m++) {
									for (n = m; n < 3; n++) {
										double D_ijmn_rz = 0.0;
										for (int k = 0; k < 3; k++) {
											for (l = 0; l < 3; l++) {

												D_ijmn_rz = D_ijmn_rz - 0.5 * (Constitutive(i, j, k, l) - Constitutive1(i, j, k, l)) * D_42(k, l, m, n, r, z, s, h, radius, x_o, x_o);

											}
										}
										A[index_E_ijkl[3 * s + i][j][r][z]][index_E_ij[3 * h + m][n]] = map[m][n] * D_ijmn_rz;

										for (p = 0; p < 3; p++) {
											double D_ijmnp_rz = 0.0;
											for (int k = 0; k < 3; k++) {
												for (l = 0; l < 3; l++) {

													D_ijmnp_rz = D_ijmnp_rz - 0.5 * (Constitutive(i, j, k, l) - Constitutive1(i, j, k, l)) * D_52(k, l, m, n, p, r, z, s, h, radius, x_o, x_o);

												}
											}
											A[index_E_ijkl[3 * s + i][j][r][z]][index_E_ijk[3 * h + m][n][p]] = map[m][n] * D_ijmnp_rz;
										}

										if (h == s)
											sym = 1.0;
										else
											sym = 0.0;

										for (p = 0; p < 3; p++) {
											for (q = 0; q < 3; q++) {
												double D_ijmnpq_rz = 0.0;
												for (int k = 0; k < 3; k++) {
													for (l = 0; l < 3; l++) {

														D_ijmnpq_rz = D_ijmnpq_rz - 0.5 * (Constitutive(i, j, k, l) - Constitutive1(i, j, k, l)) * D_62(k, l, m, n, p, q, r, z, s, h, radius, x_o, x_o);

													}
												}

												A[index_E_ijkl[3 * s + i][j][r][z]][index_E_ijkl[3 * h + m][n][p][q]] = map[m][n] * (D_ijmnpq_rz + sym * Constitutive(i, j, m, n) * d[p][r] * d[q][z]);

											}
										}

									}
								}
							}

							//	t = t + 1;
						}
					}
				}
			}

		}

	}

	for (int i = 0; i < 78 * num; i++)
	{
		for (int j = 0; j < 78 * num; j++)
		{
			HMAT(3 * NN + i, 3 * NN + j) = A[i][j];
		}
	}

	for (int i = 0; i < 78 * num; i++) {
		delete A[i];
	}
	delete[] A;

};