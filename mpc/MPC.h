//
// Created by jingche on 2020/8/25.
//

#ifndef MPC_MPC_H
#define MPC_MPC_H

#include <math.h>
#include <iostream>
#include <algorithm>
#include <eigen3/Eigen/Core>
#include <vector>
#include <cppad/ipopt/solve.hpp>
#include "eigen3/Eigen/QR"

using namespace std;
class MPC {
    typedef CPPAD_TESTVECTOR(double) Dvector;
private:
    int N = 10; // how many states we "lookahead" in the future
    double dt = 0.1; // how much time we expect environment changes

    double Lf = 2.67; // this is the length from front of vehicle to Center-of-Gravity
    double VELOCITY_MAX = 100.0; // this is what we ideally want our speed to always be

    int NUMBER_OF_STATES = 6; // px, py, psi, v, cte, epsi
    int NUMBER_OF_ACTUATIONS = 2; // steering angle, acceleration
    int NX =  N * NUMBER_OF_STATES + (N - 1) * NUMBER_OF_ACTUATIONS; // number of state + actuation variables
    int NG = N * NUMBER_OF_STATES; // number of constraints

// where the first element of each state variable is stored in the vector to be feeded the optimization algorithm
    int ID_FIRST_px = 0;
    int ID_FIRST_py = ID_FIRST_px + N;
    int ID_FIRST_psi = ID_FIRST_py + N;
    int ID_FIRST_v = ID_FIRST_psi + N;
    int ID_FIRST_cte = ID_FIRST_v + N;
    int ID_FIRST_epsi = ID_FIRST_cte + N;
    int ID_FIRST_delta = ID_FIRST_epsi + N;
    int ID_FIRST_a = ID_FIRST_delta + N - 1;

// weights for cost computations
    double W_cte = 1100.0;
    double W_epsi = 1200.0;
    double W_v = 1.0;
    double W_delta = 1000.0;
    double W_a = 10;
    double W_ddelta = 2000; // weight cost for high difference between consecutive steering actuations
    double W_da = 15.0; // weight cost for high difference between consecutive acceleration actuations
    double MAX_STEER = (45.0)/180*M_PI;

public:
    std::map<string, double> param_map;

    bool DELAY_MODE = true;
    bool DEBUG_INFO = false;

    double steer;
    double throttle;
    double ref_speed;

    Eigen::VectorXd init_state;

    Dvector x; // where all the state and actuation variables will be stored
    Dvector x_lowerbound; //lower limit for each corresponding variable in x
    Dvector x_upperbound; //upper limit for each corresponding variable in x
    Dvector g_lowerbound; // value constraint for each corresponding constraint expression
    Dvector g_upperbound; // value constraint for each corresponding constraint expression

    std::vector<double> future_xs;
    std::vector<double> future_ys;

    // errors
    double e_cte;
    double e_psi;

    MPC();
    virtual ~MPC();

    // this function solves the model given the current state and road curve coefficients.
    vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd K, vector<double>& speed);
    void updateDelayState(Eigen::VectorXd& predState, double v, double cte, double epsi ,double steering, double throttle);
    double polyeval(Eigen::VectorXd coeffs, double x);
    Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals, int order);
    vector<double> controlCB(vector<double>& cx, vector<double>& cy, vector<double> &v_ref, double steering, double throttle);
    void LoadParams (const std::map<string, double> &params);
};

#endif /* MPC_H */