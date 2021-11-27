
#pragma once
# include "systemSolver.h"
# include "Eigen/Dense"
# include "configElastic2D.h"

using namespace Eigen;

class systemSolverElastic2D : public systemSolver {

    void solve(Config& config_);  // solve the syetem of equations
    template <typename T> void _solve(Config& config_);
};

void systemSolverElastic2D::solve( Config& config_ ) {

    _solve<configElastic2D>(config_);

}

template <typename T>
void systemSolverElastic2D::_solve(Config& config_) {

    T& config = dynamic_cast<T&> (config_);
    MatrixXd& AA = config.A;
    VectorXd& BB = config.b;

//    cout << config.A;
//    cout << "\n";
//    cout << config.b;

    config.X = AA.fullPivLu().solve(BB);

//    cout << config.X;
    double r = (AA * config.X - BB).norm();

}
