//
// Created by jingche on 2020/8/21.
//

#ifndef MPC_MPC_SOLVER_H
#define MPC_MPC_SOLVER_H
#include <iostream>
#include <cppad/ipopt/solve.hpp>

using namespace std;

using namespace CppAD;
class MPC_Solver {
public:
    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
    void operator()(ADvector& fg, const ADvector& x);
    bool start(void);
};
#endif //MPC_MPC_SOLVER_H
