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

#pragma once
#include <string>
#include "Eigen/Dense"
using namespace std;
using namespace Eigen;

class Config {
public:
	MatrixXd A;  // AX = b;
	VectorXd b;
	VectorXd X;
	MatrixXd x_p;
	virtual void setUp() {		
		this->ReadBEM_num();       // read BEM node and element numbers
		this->ReadBEM_info();      // read BEM geometry
		this->Readparticle_num();  // read partical number
		this->Readpartical_info(); // read partical positons and radius
		this->ReadPostProcess_info();
	};

	virtual void Readparticle_num() {};
	virtual void Readpartical_info() {};
	virtual void ReadBEM_num() {};
	virtual void ReadBEM_info() {};
	virtual void ReadPostProcess_info() {};
	virtual string getbembuilder() { return " "; };  // return corresponding bembuilder class name
	virtual string getintegrator() { return " "; };  // return corresponding integrator class name
	virtual string getpostprocessor() { return " "; };  // return corresponding postprocessor class name
	virtual string getsystemsolver() { return " "; };// return corresponding systemsolver class name	
	virtual void solve() {};
};
