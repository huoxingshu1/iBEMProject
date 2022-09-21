#pragma once
#include <string>
#include <fstream>
#include "Eigen/Dense"
#include <iostream>
#include "Config.h"

using namespace std;
using namespace Eigen;

extern double k_0;
extern double k_1;


class configPotential3D_elliptical : public Config {

public:

	int NN = 0;
	int NE = 0;
	int num = 0;
	int nump = 0;
	double k_0_ = 0.0;
	double heat_src = 0.0;
	int Dimension = 3;
	double** s_pt;

	string particleInput = "position.txt";
	string bemInput = "BEM.txt";
	string postProcessInput = "postprocess.txt";
	MatrixXd eigen_point;
	MatrixXd radius;
	MatrixXd direc;
	VectorXi orient;
	VectorXi type;


	int* index_E_i;
	int** index_E_ij;


	MatrixXd NODES;
	MatrixXi NConnect;
	MatrixXd BC;
	MatrixXd t_bc;

	MatrixXd HMAT;
	MatrixXd GMAT;
	VectorXd U;
	VectorXd T;

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
		k_0_ = stod(T);
		getline(myfile, T);
		NN = stoi(T);
		getline(myfile, T);
		NE = stoi(T);
		myfile.close();


		NODES = MatrixXd(NN, 3);
		NConnect = MatrixXi(NE, 4);
		BC = MatrixXd(NN, 2);
		t_bc = MatrixXd(NE, 8);
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
		return "BEMBuilderPotential3D";
	};

	string getintegrator() {
		return "integratorPotential3D_elliptical";
	};

	string getsystemsolver() {
		return "systemSolverPotential3D_elliptical";
	};

	string getpostprocessor() {
		return "postProcessorPotential3D_elliptical";
	}

};