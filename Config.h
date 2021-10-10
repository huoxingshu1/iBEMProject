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
		this->Readparticle_num();  // read partical number
		this->Readpartical_info(); // read partical positons and radius
		this->ReadBEM_num();       // read BEM node and element numbers
		this->ReadBEM_info();      // read BEM geometry
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
