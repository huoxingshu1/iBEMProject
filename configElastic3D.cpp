#include "configElastic3D.h"

double mu_0 = 0.0;
double nu0 = 0.0;
double mu_1 = 0.0;
double nu1 = 0.0;

void configElastic3D::ReadBEM_info() {

	ifstream myfile;

	myfile.open(bemInput, ios::in);

	string TT;
	string Title;

	getline(myfile, TT);
	Title = TT;
	cout << "Title is: " << TT << "\n";

	getline(myfile, TT);
	EE = stod(TT);

	printf("E0 = %f \n", EE);

	getline(myfile, TT);
	nu0 = stod(TT);
	XNU = nu0;
	printf("nu0 = %f \n", nu0);
	mu_0 = EE / (2.0 * (1.0 + nu0));


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

}
void configElastic3D::Readpartical_info() {

	ifstream myfile8;
	myfile8.open(particleInput);

	double EE1 = 0;
	myfile8 >> num;
	myfile8 >> EE1;
	myfile8 >> nu1;

	printf("E1 = %f \n", EE1);
	printf("nu1 = %f \n", nu1);

	mu_1 = EE1 / (2.0 * (1.0 + nu1));

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

	HMAT = MatrixXd::Zero(3 * NN + 78 * num, 3 * NN + 78 * num);
	GMAT = MatrixXd::Zero(3 * NN + 78 * num, 12 * NE);
	U = VectorXd::Zero(3 * NN + 78 * num);
	T = VectorXd::Zero(12 * NE);

	A = MatrixXd::Zero(3 * NN + 78 * num, 3 * NN + 78 * num);
	b = VectorXd::Zero(3 * NN + 78 * num);
}