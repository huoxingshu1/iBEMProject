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


/*
 This header file defines the input file kinds in iBEM
 
 such as reading the boundary elements mesh, boundary conditions and particle assignments
 
 added the comments line based on skipping symbol "#", you can added more if you like by changing line[0] != "$"
 */
#pragma once
#include <string>
#include <fstream>
#include "Eigen/Dense"
#include <iostream>
#include <sstream>
#include "Config.h"

using namespace std;
using namespace Eigen;

extern double mu_0;
extern double nu0;
extern double mu_1;
extern double nu1;

class configElastic2D : public Config {

public:
    
    int NN = 0;
    int NE = 0;
    int num = 0;
    int nump = 0;
    double XNU = 0.0; double XNU1;
    double EE = 0.0;  double EE1;
    int Dimension = 2;
    double** s_pt;

    string particleInput = "position_2D.txt";
    string bemInput = "BEM_2D.txt";
    string postProcessInput = "postprocess_2D.txt";
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

        // the comment symbol is "#"

        ifstream myfile; int line_num = 0;
        myfile.open(particleInput); string line;
        
        while(!myfile.eof() && line_num < 3){
            // obtain the information in 1 line
            getline(myfile, line);
            // juding the # symbol of blank line
            if(line[0] != '#' && line.length() != 0){
                if (line_num == 0) {
                    num = stoi(line); line_num++;
                }
                else if (line_num == 1) {
                    EE1 = stod(line); line_num++;
                }
                
                else if (line_num == 2) {
                    XNU1 = stod(line);

                    mu_1 = EE1 * 0.5 / (1.0 + XNU1); nu1 = XNU1;

                    line_num++;
                }
            }
        }
        
        myfile.close();
    }

    
    void Readpartical_info() {

        ifstream myfile8; string line;
        myfile8.open(particleInput);
        stringstream ss; int line_num = 0;
        // 1 by num + num * positions
        while(!myfile8.eof() && line_num < 3 + num){
            getline(myfile8, line);
            if(line[0] != '#' && line.length() != 0){
                
                if(line_num == 0){
                    num = stoi(line); line_num++;
                    
                    // create the position matrix and radius
                    eigen_point = MatrixXd::Zero(num, Dimension);
                    radius = VectorXd::Zero(num);
                    
                    // create the index matrix
                    index_E_ij = new int* [Dimension * num];
                    index_E_ijk = new int** [Dimension * num];
                    index_E_ijkl = new int*** [Dimension * num];
                    
                    // output details
                    cout << "number of 2D elliptical particles =" << " " << num << endl;
                }
                
                else if (line_num < 3) {
                    line_num++;
                }

                // input the positions and radius
                else{
                    // assign text string to stringstream ss
                    ss << line;
                    
                    // input the positions, be careful "line_num - 1"
                    ss >> eigen_point(line_num - 3, 0) >>eigen_point(line_num - 3, 1) >> radius(line_num - 3);
                    
                    line_num++; ss.clear();
                }
            }

            // print the comment as author required
            else if(line[0] == '#'){

                cout << line << endl;

            }
            
            
        }
        
        myfile8.close();

        // create the index matrix, 2D problem with 2
        for (int i = 0; i < 2 * num; i++) {
            index_E_ij[i] = new int[2];
        }

        for (int i = 0; i < 2 * num; i++)
        {
            index_E_ijk[i] = new int* [2];
            for (int j = 0; j < 2; j++)
            {
                index_E_ijk[i][j] = new int[2];
            }
        }

        for (int i = 0; i < 2 * num; i++)
        {
            index_E_ijkl[i] = new int** [2];
            for (int j = 0; j < 2; j++)
            {
                index_E_ijkl[i][j] = new int* [2];
                for (int k = 0; k < 2; k++) {
                    index_E_ijkl[i][j][k] = new int[2];
                }
            }
        }

        int id = 0;

        for (int h = 0; h < num; h++)
        {
            for (int i = 0; i < 2; i++) {
                for (int j = i; j < 2; j++)
                {
                    index_E_ij[2 * h + i][j] = id;
                    id = id + 1;
                }
            }

            for (int i = 0; i < 2; i++)
            {
                for (int j = i; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        index_E_ijk[2 * h + i][j][k] = id;
                        id = id + 1;
                    }
                }
            }


            for (int i = 0; i < 2; i++) {
                for (int j = i; j < 2; j++) {
                    for (int k = 0; k < 2; k++) {
                        for (int l = 0; l < 2; l++) {
                            index_E_ijkl[2 * h + i][j][k][l] = id;
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
        
        string line; int line_num = 0;
        
        while (!myfile.eof() && line_num < 5){
            getline(myfile, line);
            
            if(line[0] != '#' && line.length() != 0){
                
                switch (line_num){
                        
                    case 0:
                        cout << "Title file is: " << line << endl;
                        line_num++;
                        break;
                        
                    case 1:
                        EE = stoi(line);
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
        t_bc = MatrixXd(NE, Dimension * 3 * Dimension);
    }

    
    void ReadBEM_info() {

        ifstream myfile; string line; stringstream ss; int line_num = 0;
        
        myfile.open(bemInput, ios::in);

        while(!myfile.eof()){
            getline(myfile, line);
            
            if(line[0] != '#' && line.length() != 0){
                // cases with line_num < 5 has been dealed with in ReadBEM() as matrix constructions
                if(line_num < 5){
                    line_num++;
                }
                
                else if (line_num < 5 + NN){
                    // read positions of BEM boundary nodes
                    ss << line; ss >> NODES(line_num - 5, 0) >> NODES(line_num - 5, 1);
                    line_num++; ss.clear();
                }
                
                else if (line_num < 5 + NN + NE){
                    // input the connection of the Elements, currently equal 3
                    ss << line; ss >> NConnect(line_num - 5 - NN, 0) >> NConnect(line_num - 5 - NN, 1) >> NConnect(line_num - 5 - NN, 2);
                    line_num++; ss.clear();
                }
                
                else if (line_num < 5 + NN + NE + Dimension * NN){
                    // input the BCs of displacements
                    ss << line; ss >> BC(line_num - 5 - NN - NE, 0) >> BC(line_num - 5 - NN - NE, 1);
                    line_num++; ss.clear();
                }
                
                else if (line_num < 5 + NN + NE + Dimension * NN + NE){
                    // input the BCs of tractions
                    ss << line;
                    for(int ii = 0; ii < Dimension * 3 * 2;ii++){
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
        // 21 = 3 + 3 * 2 + 3 * 2^2, based on the Taylor series of expansion of elliptical particles
        
        HMAT = MatrixXd::Zero(Dimension * NN + 21 * num, Dimension * NN + 21 * num);
        GMAT = MatrixXd::Zero(Dimension * NN + 21 * num, 3 * Dimension * NE);
        
        // Boundary solution vectors and eigenstrain in each particles
        U = VectorXd::Zero(Dimension * NN + 21 * num);
        
        // tractions = number of nodes in one element * Dimension * Number of elements = 6 * NE
        T = VectorXd(Dimension * 3 * NE);

        // Define the Re-organized matrix
        A = MatrixXd::Zero(Dimension * NN + 21 * num, Dimension * NN + 21 * num);
        b = VectorXd::Zero(Dimension * NN + 21 * num);

    }

    // Input the observing points
    void ReadPostProcess_info() {
        ifstream myfile;
        myfile.open(postProcessInput); string line; stringstream ss;
        int line_num = 0;
        
        while (!myfile.eof()){
            getline(myfile, line);
            
            if(line[0] != '#' && line.length() != 0){
                if(line_num == 0){
                    nump = stoi(line);
                    
                    // x_p is the position of postprocess observing points
                    x_p = MatrixXd::Zero(nump, Dimension);
                    line_num++;
                }
                else{
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


    string getbembuilder() {
        return "BEMBuilderElastic2D";
    };

    string getintegrator() {
        return "integratorElastic2D";
    };

    string getsystemsolver() {
        return "systemSolverElastic2D";
    };

    string getpostprocessor() {
        return "postProcessorElastic2D";
    }

};
