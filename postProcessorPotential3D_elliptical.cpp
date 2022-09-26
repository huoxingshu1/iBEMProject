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

# include "postProcessorPotential3D_elliptical.h"
# include "GlobalConstant.h"
# include "integratorhelperPotential3D_elliptical.h"
# include "BEMPotentialHelper3D.h"

using namespace Potential3DElliptical;
using namespace BEMPotential3D;

void Integ3_velocity_post(double** Elcor, double** dUe, double** dTe, double* point)
{
	int i, m, k, ngp;

	double xsi, eta, W, Jac, N1, N2, N3, N4;

	double** Ni;



	int Ndof = 1;


	double* N = new double[4];
	double* Vnorm = new double[3];
	double* GCcor = new double[3];
	double* gp = new double[8];
	double* w = new double[8];
	double* x = new double[3];
	double r = 0.0;
	double UP = 0.0;
	double TP = 0.0;
	Ni = new double* [1];


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


			r = sqrt((GCcor[0] - point[0]) * (GCcor[0] - point[0]) + (GCcor[1] - point[1]) * (GCcor[1] - point[1]) + (GCcor[2] - point[2]) * (GCcor[2] - point[2]));
			x[0] = (GCcor[0] - point[0]);
			x[1] = (GCcor[1] - point[1]);
			x[2] = (GCcor[2] - point[2]);

			UK(x, r, UP);
			TK(x, r, TP, Vnorm);

			for (int mm = 0; mm < 1; mm++) {
				for (int ll = 0; ll < 4; ll++) {
					for (int nn = 0; nn < 1; nn++) {
						dUe[mm][ll] = dUe[mm][ll] + UP * Ni[nn][ll] * W * Jac;
						dTe[mm][ll] = dTe[mm][ll] + TP * Ni[nn][ll] * W * Jac;
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




	for (int i = 0; i < 1; i++) {
		delete Ni[i];
	}

	delete[] Ni;



}

inline void post_velocity_BEM(int nump, int NE, int NN, Ref<MatrixXd> x_p, Ref<MatrixXi> NConnect, Ref<MatrixXd> NODES, double** s_pb, Ref<VectorXd> T, Ref<VectorXd> U)
{



# pragma omp parallel shared(nump,NE,NN,x_p,NConnect,NODES,s_pb,T,U)
	{

		double** dUe, ** dTe;
		double** HMAT, ** GMAT;

		double* point = new double[3];

		dUe = new double* [1];
		dTe = new double* [1];
		HMAT = new double* [1];
		GMAT = new double* [1];



		double** Elcor;
		Elcor = new double* [4];
		for (int j = 0; j < 4; j++) { Elcor[j] = new double[3]; }


		for (int i = 0; i < 1; i++) {
			dUe[i] = new double[4];  dTe[i] = new double[4];
			HMAT[i] = new double[NN]; GMAT[i] = new double[4 * NE];
		}


# pragma omp for

		for (int KK = 0; KK < nump; KK++) {


			for (int i = 0; i < 1; i++) {
				for (int j = 0; j < NN; j++) {
					HMAT[i][j] = 0.0;
				}

				for (int j = 0; j < 4 * NE; j++) {
					GMAT[i][j] = 0.0;
				}
			}


			for (int ee = 0; ee < NE; ee++) {

				for (int mm = 0; mm < 1; mm++) {
					for (int ll = 0; ll < 4; ll++) {
						dUe[mm][ll] = 0.0;
						dTe[mm][ll] = 0.0;
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




				Integ3_velocity_post(Elcor, dUe, dTe, point);



# pragma omp critical
				{
					for (int j = 0; j < 4; j++)
					{
						int p = 0;
						p = 1 * (sctrb[j] - 1);

						for (int k = 0; k < 1; k++) {
							for (int l = 0; l < 1; l++) {
								HMAT[l][p + k] = HMAT[l][p + k] + dTe[l][1 * j + k];
							}
						}
					}

				}



				for (int l = 0; l < 1; l++) {
					for (int k = 0; k < 4; k++) {
						GMAT[l][4 * ee + k] = dUe[l][k];
					}
				}

			}


			for (int ii = 0; ii < 1; ii++) {
				for (int mm = 0; mm < NN; mm++) {
					s_pb[KK][ii] = s_pb[KK][ii] - HMAT[ii][mm] * U(mm);
				}

				for (int mm = 0; mm < 4 * NE; mm++) {
					s_pb[KK][ii] = s_pb[KK][ii] + GMAT[ii][mm] * T(mm);
				}

			}

			//printf("KK= %d", KK);

		}


		for (int j = 0; j < 4; j++) { delete[] Elcor[j]; }

		for (int i = 0; i < 1; i++) {
			delete[] dUe[i];  delete[] dTe[i];
			delete[] HMAT[i];   delete[] GMAT[i];
		}

		delete[] Elcor;
		delete[] dUe;
		delete[] dTe;
		delete[] HMAT;
		delete[] GMAT;
	}



}

void postprocessor(int nump, int num, double** s_p, double** s_pb, int* index_E_i, int** index_E_ij, Ref<MatrixXd> x_o, Ref<MatrixXd> radius, Ref<MatrixXd> x_p, double* Eigen, Ref<VectorXd> U, Ref<VectorXd> T, int NE, int NN, Ref<MatrixXi> NConnect, Ref<MatrixXd> NODES, Ref<VectorXd> f, Ref<VectorXi>orient, Ref<MatrixXd> direc, Ref<VectorXi> type)
{


	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////


	post_velocity_BEM(nump, NE, NN, x_p, NConnect, NODES, s_pb, T, U);


# pragma omp parallel shared ( num, index_E_i, index_E_ij, radius, f, Eigen, s_p, x_o, x_p , direc, orient)

	{
		double* v_f, * v_e0, * v_e1;
		int i, h;
		v_f = new double[1];
		v_e0 = new double[1];
		v_e1 = new double[1];


		int s, p, m, x;

		MatrixXd aa = MatrixXd::Zero(3, 3);
		MatrixXd bb = MatrixXd::Zero(3, 3);
		MatrixXd glo_dire = MatrixXd::Zero(3, 3);
		glo_dire(0, 0) = 1.0;
		glo_dire(1, 1) = 1.0;
		glo_dire(2, 2) = 1.0;


# pragma omp for

		for (s = 0; s < nump; s++) {

			v_f[0] = 0.0;
			v_e0[0] = 0.0;
			v_e1[0] = 0.0;
			//Initializing///////////////////////

			////////////////////////////////////////////
			for (h = 0; h < num; h++) {

				for (int sss = 0; sss < 3; sss++) {
					for (int ggg = 0; ggg < 3; ggg++) {
						aa(sss, ggg) = glo_dire.row(sss).dot(direc.row(3 * h + ggg));
						bb(sss, ggg) = direc(3 * h + sss, ggg);
					}
				}


				for (i = 0; i < 1; i++) {


					v_f[i] = v_f[i] + Lf_00(s, h, radius, x_o, x_p, bb, orient, type) * f(h); /// f is set as local coordinate


					for (m = 0; m < 3; m++) {


						v_e0[i] = v_e0[i] + Mf_30(m, s, h, radius, x_o, x_p, bb, orient, type) * Eigen[index_E_i[3 * h + m]];


						for (p = 0; p < 3; p++) {

							v_e1[i] = v_e1[i] + Nf_40(m, p, s, h, radius, x_o, x_p, bb, orient, type) * Eigen[index_E_ij[3 * h + m][p]];

						}
					}
				}

			}


			for (x = 0; x < 1; x++) {
				s_p[s][x] = (v_f[x] + v_e0[x] + v_e1[x]);

			}


		}


		delete[] v_f;
		delete[] v_e0;
		delete[] v_e1;

	}


}


void postProcessorPotential3D_elliptical::getFlux(Config& config_) {
};

void postProcessorPotential3D_elliptical::getField(Config& config_) {

	configPotential3D_elliptical& config = dynamic_cast<configPotential3D_elliptical&> (config_);
	
	int& nump = config.nump;
	int& num = config.num;
	int& NE = config.NE;
	int& NN = config.NN;

	Ref<VectorXd> U = config.U;
	Ref<VectorXd> T = config.T;
	int*& index_E_i = config.index_E_i;
	int**& index_E_ij = config.index_E_ij;

	Ref<MatrixXd>NODES = config.NODES;
	Ref<MatrixXd>eigen_point = config.eigen_point;
	Ref<MatrixXd> radius = config.radius;
	Ref<MatrixXd> x_p = config.x_p;
	Ref<MatrixXi> NConnect = config.NConnect;
	Ref<VectorXd> f = config.f;
	Ref<VectorXi> orient = config.orient;
	Ref<MatrixXd> direc = config.direc;
	Ref<VectorXi> type = config.type;

	double* Eigen = new double[12 * num];
	for (int i = 0; i < 12 * num; i++) {
		Eigen[i] = U[i + NN];
		//cout << "Eigen"<<Eigen[i] << endl;
	}
	double** s_pf, ** s_pbf, * s_ptf;
	double** str_pf, ** str_pbf, * str_ptf;

	////////////////////////////////////////////////////////


	s_pf = new double* [nump];
	s_pbf = new double* [nump];
	s_ptf = new double[nump];

	for (int i = 0; i < nump; i++) {
		s_pf[i] = new double[1];  s_pbf[i] = new double[1];
	}

	for (int i = 0; i < nump; i++) {
		for (int j = 0; j < 1; j++) {
			s_pf[i][j] = 0.0;  
			s_pbf[i][j] = 0.0;
		}
	}

	postprocessor(nump, num, s_pf, s_pbf, index_E_i, index_E_ij, eigen_point, radius, x_p, Eigen, U, T, NE, NN, NConnect, NODES, f, orient, direc, type);
	//  postprocessor_strain(numpf, num, str_pf, str_pbf, index_E_ij, index_E_ijk, eigen_point, radius, x_pf, Eigen, U, T, NE, NN, NConnect, NODES, f, taskid);



	for (int KK = 0; KK < nump; KK++) {
		for (int ii = 0; ii < 1; ii++) {
			s_ptf[KK] = s_pf[KK][ii] + s_pbf[KK][ii];
		}
	}

	config.s_pt = s_ptf;

	/*
	ttt = 0;
	for (int KK = 0; KK < nump; KK++) {
		for (int ii = 0; ii < 6; ii++) {
			str_ptf[ttt] = str_pf[KK][ii] + str_pbf[KK][ii];
			ttt = ttt + 1;
		}
	}
	*/

	ofstream myfile7;
	myfile7.open("Temperature.txt");

	for (int KK = 0; KK < nump; KK++) {
			myfile7 << s_ptf[KK];
			myfile7 << "\n";
	}
	myfile7.close();

/*
	for (int KK = 0; KK < nump; KK++) {
			cout << config.s_pt[KK]<<endl;
	}
*/

};