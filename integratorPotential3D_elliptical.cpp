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

#include "configPotential3D_elliptical.h"
#include "integratorPotential3D_elliptical.h"
#include "integratorhelperPotential3D_elliptical.h"

using namespace Potential3DElliptical;

void get_addu(double* x, Ref<VectorXd> a, double** DD_u3, int* index, int orient, int type)

{

	for (int ii = 0; ii < 1; ii++) {
		for (int jj = 0; jj < 3; jj++) {
			DD_u3[ii][jj] = 0.0;
		}
	}

	for (int i = 0; i < 1; i++) {

		for (int m = 0; m < 3; m++) {
			DD_u3[i][index[m]] = -(-k_0) * IGreen_01(x, a, orient, m, type);
		}

	}


}

void get_addu1(double* x, Ref<VectorXd> a, double** DD_u4, int** indexp, int orient, int type)
{

	for (int ii = 0; ii < 1; ii++) {
		for (int jj = 0; jj < 9; jj++) {
			DD_u4[ii][jj] = 0.0;
		}
	}

	for (int i = 0; i < 1; i++) {

		for (int m = 0; m < 3; m++) {
			for (int p = 0; p < 3; p++) {

				DD_u4[i][indexp[m][p]] = -(-k_0) * IGreen_11(x, a, orient, m, p, type);


			}
		}

	}



}

void integratorPotential3D_elliptical::addFieldToBEM(Config& config_) {

	configPotential3D_elliptical& config = dynamic_cast<configPotential3D_elliptical&> (config_);
	Ref<MatrixXd> HMAT = config.HMAT;
	Ref<MatrixXd> eigen_point = config.eigen_point;
	Ref<MatrixXd> radius = config.radius;
	Ref<MatrixXd> NODES = config.NODES;
	Ref<VectorXi> orient = config.orient;
	Ref<MatrixXd> direc = config.direc;
	Ref<VectorXi> type = config.type;

	//cout << NODES.size() << endl;
	//cout << NODES.rows() << endl;
	//cout << NODES.cols() << endl;

	//cout << eigen_point.size() << endl;
	//cout << eigen_point.rows() << endl;
	//cout << eigen_point.cols() << endl;

	int num = config.num;
	int NN = config.NN;


# pragma omp parallel shared (HMAT)
	{

		double* x = new double[3];
		double* xp = new double[3];

		double** DD_u3, ** DD_u4;

		MatrixXd glo_direc = MatrixXd::Zero(3, 3);
		MatrixXd aa = MatrixXd::Zero(3, 3);
		MatrixXd aar = MatrixXd::Zero(3, 3);

		glo_direc(0, 0) = 1.0;
		glo_direc(1, 1) = 1.0;
		glo_direc(2, 2) = 1.0;


		int* index;
		int** indexp;

		int id = 0;
		int pp;

		index = new int[3];
		indexp = new int* [3];

		for (int ii = 0; ii < 3; ii++) {
			indexp[ii] = new int[3];
		}



		for (int ii = 0; ii < 3; ii++) {
			index[ii] = 0;
		}


		for (int ii = 0; ii < 3; ii++) {
			for (int jj = 0; jj < 3; jj++) {
				indexp[ii][jj] = 0;
			}
		}



		id = 0;
		for (int ii = 0; ii < 3; ii++) {
			index[ii] = id;
			id = id + 1;
		}


		id = 0;
		for (int ii = 0; ii < 3; ii++) {
			for (int jj = 0; jj < 3; jj++) {
				indexp[ii][jj] = id;
				id = id + 1;
			}
		}



		DD_u3 = new double* [1];
		DD_u4 = new double* [1];



		for (int ii = 0; ii < 1; ii++) {
			DD_u3[ii] = new double[3];
			DD_u4[ii] = new double[9];

		}


		for (int ii = 0; ii < 1; ii++) {
			for (int jj = 0; jj < 3; jj++) {
				DD_u3[ii][jj] = 0.0;
			}

			for (int kk = 0; kk < 9; kk++) {
				DD_u4[ii][kk] = 0.0;
			}
		}


#pragma omp for

		for (int KK = 0; KK < num; KK++) {

			for (int jj = 0; jj < NN; jj++) {

				//cout << NODES.size() << endl;
				//cout << NODES.rows() << endl;
				//cout << NODES.cols() << endl;

				//cout << eigen_point.size() << endl;
				//cout << eigen_point.rows() << endl;
				//cout << eigen_point.cols() << endl;

				x[0] = NODES(jj, 0) - eigen_point(KK, 0);
				x[1] = NODES(jj, 1) - eigen_point(KK, 1);
				x[2] = NODES(jj, 2) - eigen_point(KK, 2);


				xp[0] = 0.0; xp[1] = 0.0; xp[2] = 0.0;

				for (int sss = 0; sss < 3; sss++) {
					for (int ggg = 0; ggg < 3; ggg++) {
						aa(sss, ggg) = direc(3 * KK + sss, ggg);
					}
				}

				for (int sss = 0; sss < 3; sss++) {
					for (int ggg = 0; ggg < 3; ggg++) {
						aar(sss, ggg) = glo_direc.row(sss).dot(direc.row(3 * KK + ggg));
					}
				}


				for (int sss = 0; sss < 3; sss++) {
					for (int ggg = 0; ggg < 3; ggg++) {
						xp[sss] = xp[sss] + aa(sss, ggg) * x[ggg];
					}
				}



				VectorXd a = radius.row(KK);
				get_addu(xp, a, DD_u3, index, orient(KK), type(KK));
				get_addu1(xp, a, DD_u4, indexp, orient(KK), type(KK));



				pp = jj;
				int  tt = 12 * (KK);


				for (int qq = pp; qq < pp + 1; qq++) {

					for (int mm = 0; mm < 3; mm++) {
						HMAT(qq, NN + mm + tt) = -DD_u3[0][mm];
					}

					for (int mm = 0; mm < 9; mm++) {
						HMAT(qq, NN + mm + tt + 3) = -DD_u4[0][mm];
					}

				}


			}

		}


		for (int i = 0; i < 1; i++) {
			delete DD_u3[i];
			delete DD_u4[i];
		}

		delete[] DD_u3;
		delete[] DD_u4;


		for (int ii = 0; ii < 3; ii++)
		{
			delete[] indexp[ii];
		}

		delete[] index;
		delete[] indexp;

	}

};

void integratorPotential3D_elliptical::addFluxEquivalentInclusion(Config& config_) {

	configPotential3D_elliptical& config = dynamic_cast<configPotential3D_elliptical&> (config_);
	Ref<MatrixXd> HMAT = config.HMAT;
	Ref<MatrixXd> eigen_point = config.eigen_point;
	Ref<MatrixXd> x_o = config.eigen_point;
	Ref<MatrixXd> radius = config.radius;
	Ref<MatrixXd> NODES = config.NODES;
	Ref<VectorXi> orient = config.orient;
	Ref<MatrixXd> direc = config.direc;
	Ref<VectorXi> type = config.type;
	
	int*& index_E_i = config.index_E_i;
	int**& index_E_ij = config.index_E_ij;


	int num = config.num;
	int NN = config.NN;

	MatrixXd A = MatrixXd::Zero(12 * num, 12 * num);
    //////////////////////////////////////////////////////

# pragma omp parallel shared(A,num,radius,x_o,index_E_i,index_E_ij,orient, direc)

	{
		//ofstream testH3;
		//testH3.open("testH3.txt");


		int t = 0;
		int m, p, r, h;

		double sym = 0.0;
		MatrixXd aa = MatrixXd::Zero(3, 3);
		MatrixXd bb = MatrixXd::Zero(3, 3);
		////Assembly Matrix Tensor///////

#	pragma omp for

		for (int s = 0; s < num; s++) {


			for (int i = 0; i < 3; i++) {

				for (h = 0; h < num; h++) {

					if (h == s) {
						sym = 1.0;
					}
					else {
						sym = 0.0;
					}
					for (int ii = 0; ii < 3; ii++) {
						for (int jj = 0; jj < 3; jj++) {
							aa(ii, jj) = direc.row(3 * s + ii).dot(direc.row(3 * h + jj));
							bb(ii, jj) = direc(3 * h + ii, jj);
						}
					}




					for (m = 0; m < 3; m++) {

						double D_ijmn = 0.0;

						for (int cc = 0; cc < 3; cc++) {
							D_ijmn = D_ijmn - aa(i, cc) * (-k_0 - (-k_1)) * Df_40(cc, m, s, h, radius, x_o, x_o, orient, bb, type);
						}



						A(index_E_i[3 * s + i], index_E_i[3 * h + m]) = D_ijmn + sym * (-k_0) * d[i][m];


						for (p = 0; p < 3; p++) {
							double D_ijmnp = 0.0;

							for (int cc = 0; cc < 3; cc++) {
								D_ijmnp = D_ijmnp - aa(i, cc) * (-k_0 - (-k_1)) * Df_50(cc, m, p, s, h, radius, x_o, x_o, orient, bb, type);
							}


							A(index_E_i[3 * s + i], index_E_ij[3 * h + m][p]) = D_ijmnp;

						}


					}
				}
			}
			
		}
	}


	//////////////////////////////////////////////////////
	for (int i = 0; i < 12 * num; i++) {
		for (int j = 0; j < 12 * num; j++) {
			HMAT(NN + i, NN + j) += A(i, j);
		}
	}

	A.resize(0, 0);
};

void integratorPotential3D_elliptical::addFluxEquivalentFirstOrderInclusion(Config& config_) {
	configPotential3D_elliptical& config = dynamic_cast<configPotential3D_elliptical&> (config_);
	Ref<MatrixXd> HMAT = config.HMAT;
	Ref<MatrixXd> eigen_point = config.eigen_point;
	Ref<MatrixXd> x_o = config.eigen_point;
	Ref<MatrixXd> radius = config.radius;
	Ref<MatrixXd> NODES = config.NODES;
	Ref<VectorXi> orient = config.orient;
	Ref<MatrixXd> direc = config.direc;
	Ref<VectorXi> type = config.type;

	int*& index_E_i = config.index_E_i;
	int**& index_E_ij = config.index_E_ij;


	int num = config.num;
	int NN = config.NN;

	MatrixXd A = MatrixXd::Zero(12 * num, 12 * num);
	//////////////////////////////////////////////////////

# pragma omp parallel shared(A,num,radius,x_o,index_E_i,index_E_ij,orient, direc)

	{
		//ofstream testH3;
		//testH3.open("testH3.txt");


		int t = 0;
		int m, p, r, h;

		double sym = 0.0;
		MatrixXd aa = MatrixXd::Zero(3, 3);
		MatrixXd bb = MatrixXd::Zero(3, 3);
		////Assembly Matrix Tensor///////

#	pragma omp for

		for (int s = 0; s < num; s++) {
			

			for (int i = 0; i < 3; i++) {
				for (r = 0; r < 3; r++) {
					for (h = 0; h < num; h++) {

						for (int ii = 0; ii < 3; ii++) {
							for (int jj = 0; jj < 3; jj++) {
								aa(ii, jj) = direc.row(3 * s + ii).dot(direc.row(3 * h + jj));
								bb(ii, jj) = direc(3 * h + ii, jj);
							}
						}


						for (m = 0; m < 3; m++) {
							double D_ijmn_r = 0.0;

							for (int cc = 0; cc < 3; cc++) {
								for (int hh = 0; hh < 3; hh++) {
									D_ijmn_r = D_ijmn_r - aa(i, cc) * aa(r, hh) * (-k_0 - (-k_1)) * Df_41(cc, m, hh, s, h, radius, x_o, x_o, orient, bb, type);
								}
							}

							//	if (s == 1 && h == 0){
							//		testH3 << i << j << r << m << n << " " << D_ijmn_r << endl;
							//	}

							A(index_E_ij[3 * s + i][r], index_E_i[3 * h + m]) = D_ijmn_r;

							if (h == s) {
								sym = 1.0;
							}
							else {
								sym = 0.0;
							}

							for (p = 0; p < 3; p++) {
								double D_ijmnp_r = 0.0;

								for (int cc = 0; cc < 3; cc++) {
									for (int hh = 0; hh < 3; hh++) {
										D_ijmnp_r = D_ijmnp_r - aa(i, cc) * aa(r, hh) * (-k_0 - k_1) * Df_51(cc, m, p, hh, s, h, radius, x_o, x_o, orient, bb, type);
									}
								}

								//	if (s == 0 && h == 1){
								//		testH3 << i << j << r << m << n << p << " " << D_ijmnp_r << endl;
								//	}

								A(index_E_ij[3 * s + i][r], index_E_ij[3 * h + m][p]) = D_ijmnp_r + sym * (-k_0) * I_40(i, r, m, p);
							}

						}
					}

				}


			}
		}
	}


	//////////////////////////////////////////////////////
	for (int i = 0; i < 12 * num; i++) {
		for (int j = 0; j < 12 * num; j++) {
			HMAT(NN + i, NN + j) += A(i, j);
		}
	}

	A.resize(0, 0);

};

void integratorPotential3D_elliptical::addFluxEquivalentSecondOrderInclusion(Config& config_) {};