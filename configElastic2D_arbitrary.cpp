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

#include "configElastic2D_arbitrary.h"
#include <string>
#include <fstream>
#include "Eigen/Dense"
#include "Config.h"

using namespace Eigen;

void configElastic2D_arbitrary::Initialization_record()
{
    number_in_array = 500;
    changed_coord = new double** [num_NN];
    number_weight = new double* [num_NN];
    count_number = new int[num_NN];

    for (int i = 0; i < num_NN; i++) {
        changed_coord[i] = new double* [number_in_array];
        number_weight[i] = new double[number_in_array];
    }

    for (int i = 0; i < num_NN; i++) {
        for (int j = 0; j < number_in_array; j++) {
            changed_coord[i][j] = new double[2];
        }

    }
}


void configElastic2D_arbitrary::Readparticle_num()
{
    ifstream myfile; myfile.open(particleInput);

    int line_num = 0; string line;

    while (!myfile.eof() && line_num < 2) {
        // only EE1, XNU1, num_NE and num_NN are needed here
        getline(myfile, line);
        // judging the # symbol and blank line
        if (line[0] != '#' && line.length() != 0) {

            if (line_num == 0) {
                EE1 = stod(line);
                line_num++;
            }

            else if (line_num == 1) {
                XNU1 = stod(line);
                mu_1 = EE1 * 0.5 / (1.0 + XNU1); nu1 = XNU1;
                line_num++;
            }

            else if (line_num == 2) {
                num_NN = stoi(line); line_num++;  
            }

            else if (line_num == 3) {
                num_NE = stoi(line); line_num++;
            }

        }

    }

    myfile.close();
}

void configElastic2D_arbitrary::Readpartical_info()
{
    ifstream myfile; myfile.open(particleInput);

    int line_num = 0; string line; stringstream ss;

    while (!myfile.eof()) {
        // only num_NE and num_NN are needed here
        getline(myfile, line);
        // judging the # symbol and blank line
        if (line[0] != '#' && line.length() != 0) {

            if (line_num < 2) {
                line_num++;
            }

            else if (line_num == 2) {
                num_NN = stoi(line); line_num++;
                eigen_point = MatrixXd::Zero(num_NN, Dimension);
                Initialization_record();
            }

            else if (line_num == 3) {
                num_NE = stoi(line); line_num++;
                // assume 6-node quadratic triangular elements
                eigen_connect = MatrixXi::Zero(num_NE, 6);

                // create arrays to store eigenstrain
                index_E_ij = new int* [Dimension * num_NN];

                // cout details about the particle
                cout << "number of 2D quadratic triangular elements is" << " " << num_NE << endl;
            }

            // input the positions of nodes and elements
            else if(line_num < 4 + num_NN){
                ss << line;
                ss >> eigen_point(line_num - 4, 0) >> eigen_point(line_num - 4, 1);
                line_num++; ss.clear();
            }

            else if (line_num < 4 + num_NN + num_NE) {
                ss << line; 
                for (int ii = 0; ii < 6; ii++) {
                    ss >> eigen_connect(line_num - 4 - num_NN, ii);
                }
                line_num++; ss.clear();
            }

        }

        else if (line[0] == '#') {
            cout << line << endl;
        }

    }

    myfile.close();

    // assign values to array index_E_ij
    for (int i = 0; i < Dimension * num_NN; i++) {
        index_E_ij[i] = new int[2];
    }

    int id = 0;

    for (int h = 0; h < num_NN; h++)
    {
        for (int i = 0; i < 2; i++) {
            for (int j = i; j < 2; j++)
            {
                index_E_ij[2 * h + i][j] = id;
                id = id + 1;
            }
        }
    }
}

void configElastic2D_arbitrary::ReadBEM_num() {

    ifstream myfile;

    myfile.open(bemInput, ios::in);

    string line; int line_num = 0;

    while (!myfile.eof() && line_num < 5) {
        getline(myfile, line);

        if (line[0] != '#' && line.length() != 0) {

            switch (line_num) {

            case 0:
                cout << "Title file is: " << line << endl;
                line_num++;
                break;

            case 1:
                EE = stod(line);
                cout << "Young's moduli of matrix is: " << EE << endl;
                line_num++;
                break;

            case 2:
                XNU = stod(line);
                cout << "Poisson's ratio of matrix is: " << XNU << endl;

                nu0 = XNU; mu_0 = EE * 0.5 / (1.0 + nu0);

                line_num++;
                break;

            case 3:
                NN = stoi(line);
                cout << "Number of Nodes is: " << NN << endl;
                line_num++;
                break;

            case 4:
                NE = stoi(line);
                cout << "Number of Elements is: " << NE << endl;
                line_num++;
                break;

            default:
                cout << "No speficic value or missing value has been defined for the BEM part!" << endl;


            }



        }

    }

    myfile.close();

    // define the matrix for NODES positions, Element connections, displacements and traction boundary conditions

    NODES = MatrixXd(NN, Dimension);

    // Connection depends on the types of the elements, here, 3-node quadratic is used
    NConnect = MatrixXi(NE, 3);

    // Bc is boundary conditions about displacements
    BC = MatrixXd(Dimension * NN, 2);

    // t_BC is boundary conditions about tractions, (Dimension * number of nodes in one element * define and number)
    t_bc = MatrixXd(NE, 2 * 3 * Dimension);
}


void configElastic2D_arbitrary::ReadBEM_info() {

    ifstream myfile; string line; stringstream ss; int line_num = 0;

    myfile.open(bemInput, ios::in);

    while (!myfile.eof()) {
        getline(myfile, line);

        if (line[0] != '#' && line.length() != 0) {
            // cases with line_num < 5 has been dealed with in ReadBEM() as matrix constructions
            if (line_num < 5) {
                line_num++;
            }

            else if (line_num < 5 + NN) {
                // read positions of BEM boundary nodes
                ss << line; ss >> NODES(line_num - 5, 0) >> NODES(line_num - 5, 1);
                line_num++; ss.clear();
            }

            else if (line_num < 5 + NN + NE) {
                // input the connection of the Elements, currently equal 3
                ss << line; ss >> NConnect(line_num - 5 - NN, 0) >> NConnect(line_num - 5 - NN, 1) >> NConnect(line_num - 5 - NN, 2);
                line_num++; ss.clear();
            }

            else if (line_num < 5 + NN + NE + Dimension * NN) {
                // input the BCs of displacements
                ss << line; ss >> BC(line_num - 5 - NN - NE, 0) >> BC(line_num - 5 - NN - NE, 1);
                line_num++; ss.clear();
            }

            else if (line_num < 5 + NN + NE + Dimension * NN + NE) {
                // input the BCs of tractions
                ss << line;
                for (int ii = 0; ii < Dimension * 3 * 2; ii++) {
                    ss >> t_bc(line_num - 5 - NN - NE - Dimension * NN, ii);
                }
                line_num++;
                //cout << line << endl;
                ss.clear();
            }

        }

        else if (line[0] == '#') {

            cout << line << endl;
        }
    }

    myfile.close();
    // Assign the memory of iBEM matrix
    // Each particle node has 3 different eigenstrain

    HMAT = MatrixXd::Zero(Dimension * NN + 3 * num_NN, Dimension * NN + 3 * num_NN);
    GMAT = MatrixXd::Zero(Dimension * NN + 3 * num_NN, 3 * Dimension * NE);

    // Boundary solution vectors and eigenstrain in each particles
    U = VectorXd::Zero(Dimension * NN + 3 * num_NN);

    // tractions = number of nodes in one element * Dimension * Number of elements = 6 * NE
    T = VectorXd(Dimension * 3 * NE);

    // Define the Re-organized matrix
    A = MatrixXd::Zero(Dimension * NN + 3 * num_NN, Dimension * NN + 3 * num_NN);
    b = VectorXd::Zero(Dimension * NN + 3 * num_NN);

}

// Input the observing points
void configElastic2D_arbitrary::ReadPostProcess_info() {
    ifstream myfile;
    myfile.open(postProcessInput); string line; stringstream ss;
    int line_num = 0;

    while (!myfile.eof() && line_num <= nump) {
        getline(myfile, line);

        if (line[0] != '#' && line.length() != 0) {
            if (line_num == 0) {
                nump = stoi(line);

                // x_p is the position of postprocess observing points
                x_p = MatrixXd::Zero(nump, Dimension);
                line_num++;
            }
            else {
                ss << line; ss >> x_p(line_num - 1, 0) >> x_p(line_num - 1, 1);
                line_num++; ss.clear();
            }
        }

        else if (line[0] == '#') {
            cout << line << endl;
        }

    }

    myfile.close();

}

/*
String functions to initialize services in ClassFinder.h
*/

string configElastic2D_arbitrary::getbembuilder()
{
    // although BEM part will not change, but function ApplyBC is different.
    return "BEMBuilderElastic2D_arbitrary";
}

string configElastic2D_arbitrary::getintegrator()
{
    return "integratorElastic2D_arbitrary";
}

string configElastic2D_arbitrary::getsystemsolver()
{
    return "systemSolverElastic2D_arbitrary";
}

string configElastic2D_arbitrary::getpostprocessor()
{
    return "postProcessorElastic2D_arbitrary";
}