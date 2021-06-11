//
// Created by jingche on 2020/8/21.
//
#include "mpc_solver.h"

void MPC_Solver::operator()(ADvector &fg, const ADvector &x) {
    assert(fg.size() == 5);
    assert(x.size() == 4);
    // variables
    AD<double> x1 = x[0];
    AD<double> x2 = x[1];
    AD<double> x3 = x[2];
    AD<double> x4 = x[3];
    // f(x) objective function
    fg[0] = x1 * x4 * (x1 + x2 + x3) + x3;
    // constraints
    fg[1] = x1 * x2 * x3 * x4;
    fg[2] = x1 * x1 + x2 * x2 + x3 * x3 + x4 * x4;
    return;
}