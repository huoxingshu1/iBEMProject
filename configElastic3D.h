#pragma once
#include <string>
#include <fstream>
#include "Eigen/Dense"
#include <iostream>
#include "Config.h"

using namespace std;
using namespace Eigen;

extern double mu_0;
extern double nu0;
extern double mu_1;
extern double nu1;


class configElastic3D : public Config {

public:

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

	void Readpartical_info();

	void ReadBEM_num() {

		ifstream myfile;

		myfile.open(bemInput, ios::in);

		string T;
		string Title;
		getline(myfile, T);
		Title = T;
		cout << "Title is: " << T << "\n";
		getline(myfile, T);
		EE = stod(T);
		getline(myfile, T);
		XNU = stod(T);
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

	void ReadBEM_info();

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
