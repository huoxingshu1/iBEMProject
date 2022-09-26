#include "configPotential3D_elliptical.h"

double k_0 = 0.0;
double k_1 = 0.0;

void configPotential3D_elliptical::ReadBEM_info() {

	ifstream myfile;

	myfile.open(bemInput, ios::in);

	string TT;
	string Title;

	getline(myfile, TT);
	Title = TT;
	cout << "Title is: " << TT << "\n";

	getline(myfile, TT);
	k_0_ = stod(TT);
	k_0 = k_0_;

	printf("k_0 = %f \n", k_0);

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
	for (int i = 0; i < NN; i++) {
		for (int j = 0; j < 2; j++) {
			myfile >> BC(i, j);
		}
	}

	for (int i = 0; i < NE; i++) {
		for (int j = 0; j < 8; j++) {
			myfile >> t_bc(i, j);
			//	printf("%f", t_bc[i][j]);
		}
		//	printf("\n");
	}

	HMAT = MatrixXd::Zero(NN + 12 * num, NN + 12 * num);
	GMAT = MatrixXd::Zero(NN + 12 * num, 4 * NE);
	U = VectorXd::Zero(3 * NN + 12 * num);
	T = VectorXd::Zero(4 * NE);
	Ff = VectorXd::Zero(9 * num);
	Fb = VectorXd::Zero(NN);

	b = VectorXd::Zero(NN + 12 * num);
	myfile.close();

}
void configPotential3D_elliptical::Readpartical_info() {

	ifstream myfile8;
	myfile8.open(particleInput);

	heat_src = 0.0;
	myfile8 >> num;
	myfile8 >> k_1;
	cout << "k_1 = " << k_1 << endl;

	eigen_point = MatrixXd::Zero(num, 3);
	radius = MatrixXd::Zero(num, 3);
	direc = MatrixXd::Zero(3 * num, 3);
	orient = VectorXi::Zero(num);
	type = VectorXi::Zero(num);
	
	index_E_i = new int [3 * num];
	index_E_ij = new int* [3 * num];

	for (int KK = 0; KK < num; KK++) {
		myfile8 >> eigen_point(KK, 0) >> eigen_point(KK, 1) >> eigen_point(KK, 2) >> radius(KK, 0) >> radius(KK, 1) >> radius(KK, 2) >> orient(KK) >> type(KK);
	}

	for (int KK = 0; KK < num; KK++) {
		myfile8 >> direc(3 * KK, 0) >> direc(3 * KK, 1) >> direc(3 * KK, 2);
		myfile8 >> direc(3 * KK + 1, 0) >> direc(3 * KK + 1, 1) >> direc(3 * KK + 1, 2);
		myfile8 >> direc(3 * KK + 2, 0) >> direc(3 * KK + 2, 1) >> direc(3 * KK + 2, 2);
	}

	myfile8.close();


	for (int i = 0; i < 3 * num; i++) {
		index_E_ij[i] = new int[3];
	}


	int id = 0;

	for (int h = 0; h < num; h++)
	{
		for (int i = 0; i < 3; i++) {
			index_E_i[3 * h + i] = id;
			id = id + 1;
		}

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++)
			{
				index_E_ij[3 * h + i][j] = id;
				id = id + 1;
			}
		}

	}

	f = VectorXd::Zero(num);
	
	for (int h = 0; h < num; h++) {
		f(h) = heat_src;
	}
	cout << "partical position" << endl;
	cout << eigen_point << endl;
	cout << "radius" << endl;
	cout << radius << endl;
	cout << "direc" << endl;
	cout << direc << endl;
	cout << "orient" << endl;
	cout << orient << endl;
	cout << "type" << endl;
	cout << type << endl;

}