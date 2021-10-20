#pragma once
#include <string>
#include <fstream>
#include "Eigen/Dense"
#include <iostream>
#include "Config.h"

using namespace std;
using namespace Eigen;

class configElastic3D : public Config {

public:
	double nu0 = 0.0;
	double mu_0 = 0.0;
	double nu1 = 0.0;
	double mu_1 = 0.0;
	int NN = 0;
	int NE = 0;
	int num = 0;
	int nump = 0;
	double XNU = 0.0;
	double EE = 0.0;
	int Dimension = 3;
	double** s_pt;

	string particleInput = "position.txt";
	string bemInput = "BEM.txt";
	string postProcessInput = "postprocess.txt";
	MatrixXd eigen_point;
	VectorXd radius;
	int** index_E_ij;
	int*** index_E_ijk;
	int**** index_E_ijkl;

	MatrixXd NODES;
	MatrixXi NConnect;
	MatrixXd BC;
	MatrixXd t_bc;

	MatrixXd HMAT;
	MatrixXd GMAT;
	VectorXd U;
	VectorXd T;

	MatrixXd A;
	VectorXd b;

	void Readparticle_num() {

		ifstream myfile;
		myfile.open(particleInput);
		myfile >> num;
		myfile.close();
	}

	void Readpartical_info() {

		ifstream myfile8;
		myfile8.open(particleInput);

		myfile8 >> num;

		eigen_point = MatrixXd::Zero(num, 3);     //particle position
		radius = VectorXd::Zero(num);             //particle radius
		index_E_ij = new int* [3 * num];
		index_E_ijk = new int** [3 * num];
		index_E_ijkl = new int*** [3 * num];

		for (int KK = 0; KK < num; KK++) {
			myfile8 >> eigen_point(KK, 0) >> eigen_point(KK, 1) >> eigen_point(KK, 2) >> radius(KK);
		}
		myfile8.close();


		for (int i = 0; i < 3 * num; i++) {
			index_E_ij[i] = new int[3];
		}

		for (int i = 0; i < 3 * num; i++)
		{
			index_E_ijk[i] = new int* [3];
			for (int j = 0; j < 3; j++)
			{
				index_E_ijk[i][j] = new int[3];
			}
		}


		for (int i = 0; i < 3 * num; i++)
		{
			index_E_ijkl[i] = new int** [3];
			for (int j = 0; j < 3; j++)
			{
				index_E_ijkl[i][j] = new int* [3];
				for (int k = 0; k < 3; k++) {
					index_E_ijkl[i][j][k] = new int[3];
				}
			}
		}


		int id = 0;

		for (int h = 0; h < num; h++)
		{
			for (int i = 0; i < 3; i++) {
				for (int j = i; j < 3; j++)
				{
					index_E_ij[3 * h + i][j] = id;
					id = id + 1;
				}
			}

			for (int i = 0; i < 3; i++)
			{
				for (int j = i; j < 3; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						index_E_ijk[3 * h + i][j][k] = id;
						id = id + 1;
					}
				}
			}

			for (int i = 0; i < 3; i++) {
				for (int j = i; j < 3; j++) {
					for (int k = 0; k < 3; k++) {
						for (int l = 0; l < 3; l++) {
							index_E_ijkl[3 * h + i][j][k][l] = id;
							id = id + 1;
						}
					}
				}
			}

		}
	}

	void ReadBEM_num() {

		ifstream myfile;

		myfile.open(bemInput, ios::in);

		string T;
		string Title;
		getline(myfile, T);
		Title = T;
		cout << "Title is: " << T << "\n";
		getline(myfile, T);
		EE = stoi(T);
		getline(myfile, T);
		XNU = stoi(T);
		getline(myfile, T);
		NN = stoi(T);
		getline(myfile, T);
		NE = stoi(T);
		myfile.close();


		NODES = MatrixXd(NN, 3);
		NConnect = MatrixXi(NE, 4);
		BC = MatrixXd(3 * NN, 2);
		t_bc = MatrixXd(NE, 24);
	}

	void ReadBEM_info() {

		ifstream myfile;

		myfile.open(bemInput, ios::in);

		string TT;
		string Title;

		getline(myfile, TT);
		Title = TT;
		cout << "Title is: " << TT << "\n";

		getline(myfile, TT);
		EE = stod(TT);

		printf("EE = %f \n", EE);		

		getline(myfile, TT);
		XNU = stod(TT);
		printf("mu_0 = %f \n", XNU);

		getline(myfile, TT);
		NN = stoi(TT);

		getline(myfile, TT);
		NE = stoi(TT);


		for (int i = 0; i < NN; i++) {
			myfile >> NODES(i, 0) >> NODES(i, 1) >> NODES(i, 2);
		}
		cout << "nodes\n";



		for (int i = 0; i < NE; i++) {
			for (int j = 0; j < 4; j++) {
				myfile >> NConnect(i, j);
			}
		}

		cout << "nconnect\n";


		cout << "BC\n";
		for (int i = 0; i < NN * 3; i++) {
			for (int j = 0; j < 2; j++) {
				myfile >> BC(i, j);
			}
		}

		for (int i = 0; i < NE; i++) {
			for (int j = 0; j < 24; j++) {
				myfile >> t_bc(i, j);
				//	printf("%f", t_bc[i][j]);
			}
			//	printf("\n");
		}

		myfile.close();

		HMAT = MatrixXd::Zero(3 * NN + 78 * num, 3 * NN + 78 * num);
		GMAT = MatrixXd::Zero(3 * NN + 78 * num, 12 * NE);
		U = VectorXd::Zero(3 * NN + 78 * num);
		T = VectorXd::Zero(12 * NE);

		A = MatrixXd::Zero(3 * NN + 78 * num, 3 * NN + 78 * num);
		b = VectorXd::Zero(3 * NN + 78 * num);

	}

	void ReadPostProcess_info() {
		ifstream myfile;
		myfile.open(postProcessInput);
		myfile >> nump;

		x_p = MatrixXd::Zero(nump, 3);

		for (int KK = 0; KK < nump; KK++) {
			myfile >> x_p(KK, 0) >> x_p(KK, 1) >> x_p(KK, 2);
		}

		myfile.close();

	}


	string getbembuilder() {
		return "BEMBuilderElastic3D";
	};

	string getintegrator() {
		return "integratorElastic3D";
	};

	string getsystemsolver() {
		return "systemSolverElastic3D";
	};

	string getpostprocessor() {
		return "postProcessorElastic3D";
	}

};
