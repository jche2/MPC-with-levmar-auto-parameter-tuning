#include "MPC.h"

using CppAD::AD;
using namespace std;

class FG_eval {
private:
    int N = 10; // how many states we "lookahead" in the future
    double dt = 0.1; // how much time we expect environment changes

    double Lf = 0.1; // this is the length from front of vehicle to Center-of-Gravity
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

    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
    Eigen::VectorXd K; // Fitted road curve polynomial coefficients
    vector<double> TARGET_SPEED;
    void setSpeed(vector<double>& speed) {
        this->TARGET_SPEED = speed;
    }

    FG_eval(Eigen::VectorXd Kin) : K(Kin) {}

    void LoadParams (const std::map<string, double> &params) {
        N = params.find("N") != params.end() ? params.at("N") : N;
        dt = params.find("dt") != params.end() ? params.at("dt") : dt;
        Lf = params.find("Lf") != params.end() ? params.at("Lf") : Lf;
        VELOCITY_MAX = params.find("VELOCITY_MAX") != params.end() ? params.at("VELOCITY_MAX") : VELOCITY_MAX;
        W_cte = params.find("W_cte") != params.end() ? params.at("W_cte") : W_cte;
        W_epsi = params.find("W_epsi") != params.end() ? params.at("W_epsi") : W_epsi;
        W_v = params.find("W_v") != params.end() ? params.at("W_v") : W_v;
        W_delta = params.find("W_delta") != params.end() ? params.at("W_delta") : W_delta;
        W_a = params.find("W_a") != params.end() ? params.at("W_a") : W_a;
        W_ddelta = params.find("W_ddelta") != params.end() ? params.at("W_ddelta") : W_ddelta;
        W_da = params.find("W_da") != params.end() ? params.at("W_da") : W_da;
        MAX_STEER = params.find("MAX_STEER") != params.end() ? params.at("MAX_STEER") : MAX_STEER;

        NX =  N * NUMBER_OF_STATES + (N - 1) * NUMBER_OF_ACTUATIONS; // number of state + actuation variables
        NG = N * NUMBER_OF_STATES; // number of constraints

        // where the first element of each state variable is stored in the vector to be feeded the optimization algorithm
        ID_FIRST_px = 0;
        ID_FIRST_py = ID_FIRST_px + N;
        ID_FIRST_psi = ID_FIRST_py + N;
        ID_FIRST_v = ID_FIRST_psi + N;
        ID_FIRST_cte = ID_FIRST_v + N;
        ID_FIRST_epsi = ID_FIRST_cte + N;
        ID_FIRST_delta = ID_FIRST_epsi + N;
        ID_FIRST_a = ID_FIRST_delta + N - 1;
    }

    void operator()(ADvector& fg, const ADvector& x) {
        // fg a vector containing the cost and all constraints
        // x is a vector containing all states and actuations for N "lookahead" states and actuations.

        //*********************************************************
        //* COST DEFINED HERE
        //*********************************************************

        fg[0] = 0.0;
//        cout << "W_V IS : " << W_v << endl;
        for (int i = 0; i < N; ++i) {

            const auto cte = x[ID_FIRST_cte + i];
            const auto epsi = x[ID_FIRST_epsi + i];
            const auto v = x[ID_FIRST_v + i] - TARGET_SPEED[i];

            fg[0] += (W_cte * cte * cte + W_epsi * epsi * epsi + W_v * v * v);
            if (i >= 3*N/4) {
                fg[0] += (W_cte * cte * cte + W_epsi * epsi * epsi + W_v * v * v);
            }
//            cout<< "error epsi is : " << epsi << endl;
        }

        for (int i = 0; i < N - 1; ++i) {

            const auto delta = x[ID_FIRST_delta + i];
            const auto a = x[ID_FIRST_a + i];

            fg[0] += (W_delta * delta * delta + W_a * a * a);
        }

        for (int i = 0; i < N - 2; ++i) {

            const auto ddelta = x[ID_FIRST_delta + i + 1] - x[ID_FIRST_delta + i];
            const auto da = x[ID_FIRST_a + i + 1] - x[ID_FIRST_a + i];

            fg[0] += (W_ddelta * ddelta * ddelta + W_da * da * da);
        }
        //*********************************************************
        //* CONSTRAINTS DEFINED HERE
        //*********************************************************

        // given state does not vary
        fg[ID_FIRST_px + 1] = x[ID_FIRST_px];
        fg[ID_FIRST_py + 1] = x[ID_FIRST_py];
        fg[ID_FIRST_psi + 1] = x[ID_FIRST_psi];
        fg[ID_FIRST_v + 1] = x[ID_FIRST_v];
        fg[ID_FIRST_cte + 1] = x[ID_FIRST_cte];
        fg[ID_FIRST_epsi + 1] = x[ID_FIRST_epsi];

        // constraints based on our kinematic model
        auto pre_epsi = x[ID_FIRST_epsi];

        for (int i = 0; i < N - 1; ++i) {

            // where the current state variables of interest are stored
            // stored for readability
            const int ID_CURRENT_px = ID_FIRST_px + i;
            const int ID_CURRENT_py = ID_FIRST_py + i;
            const int ID_CURRENT_psi = ID_FIRST_psi + i;
            const int ID_CURRENT_v = ID_FIRST_v + i;
            const int ID_CURRENT_cte = ID_FIRST_cte + i;
            const int ID_CURRENT_epsi = ID_FIRST_epsi + i;
            const int ID_CURRENT_delta = ID_FIRST_delta + i;
            const int ID_CURRENT_a = ID_FIRST_a + i;

            //current state and actuations
            const auto px0 = x[ID_CURRENT_px];
            const auto py0 = x[ID_CURRENT_py];
            const auto psi0 = x[ID_CURRENT_psi];
            const auto v0 = x[ID_CURRENT_v];
            const auto cte0 = x[ID_CURRENT_cte];
            const auto epsi0 = x[ID_CURRENT_epsi];
            const auto delta0 = x[ID_CURRENT_delta];
            const auto a0 = x[ID_CURRENT_a];

            // next state
            const auto px1 = x[ID_CURRENT_px + 1];
            const auto py1 = x[ID_CURRENT_py + 1];
            const auto psi1 = x[ID_CURRENT_psi + 1];
            const auto v1 = x[ID_CURRENT_v + 1];
            const auto cte1 = x[ID_CURRENT_cte + 1];
            const auto epsi1 = x[ID_CURRENT_epsi + 1];

            // desired py and psi
            const auto py_desired = K[3] * px0 * px0 * px0 + K[2] * px0 * px0 + K[1] * px0 + K[0];
            const auto psi_desired = CppAD::atan(3.0 * K[3] * px0 * px0 + 2.0 * K[2] * px0 + K[1]);

//            cout<< " psi_desired is : " << psi_desired << endl;
            // cout<< " d_psi is : " << d_psi << endl;


            // relationship of current state + actuations and next state
            // based on our kinematic model
            const auto px1_f = px0 + v0 * CppAD::cos(psi0) * dt;
            const auto py1_f = py0 + v0 * CppAD::sin(psi0) * dt;
            const auto psi1_f = psi0 + v0 * (-delta0) / Lf * dt;
            const auto v1_f = v0 + a0 * dt;
            const auto cte1_f = py_desired - py0 + v0 * CppAD::sin(epsi0) * dt;
            const auto epsi1_f = psi0 - psi_desired + v0 * (-delta0) / Lf * dt;

            // store the constraint expression of two consecutive states
            fg[ID_CURRENT_px + 2] = px1 - px1_f;
            fg[ID_CURRENT_py + 2] = py1 - py1_f;
            fg[ID_CURRENT_psi + 2] = psi1 - psi1_f;
            fg[ID_CURRENT_v + 2] = v1 - v1_f;
            fg[ID_CURRENT_cte + 2] = cte1 - cte1_f;
            fg[ID_CURRENT_epsi + 2] = epsi1 - epsi1_f;
        }
    }
};

MPC::MPC() {

    //**************************************************************
    //* SET INITIAL VALUES OF VARIABLES
    //**************************************************************
    this->x.resize(NX);

    // all states except the ID_FIRST are set to zero
    // the aformentioned states will be initialized when solve() is called

    for (int i = 0; i < NX; ++i) {
        this->x[i] = 0.0;
    }

    //**************************************************************
    //* SET UPPER AND LOWER LIMITS OF VARIABLES
    //**************************************************************

    this->x_lowerbound.resize(NX);
    this->x_upperbound.resize(NX);

    // all other values large values the computer can handle
    for (int i = 0; i < ID_FIRST_delta; ++i) {
        this->x_lowerbound[i] = -1.0e10;
        this->x_upperbound[i] = 1.0e10;
    }

    // all actuation inputs (steering, acceleration) should have values between [-1, 1]
    for (int i = ID_FIRST_delta; i < ID_FIRST_a; ++i) {
        this->x_lowerbound[i] = - MAX_STEER;;
        this->x_upperbound[i] = MAX_STEER;;
    }

    for (int i = ID_FIRST_a; i < NX; ++i) {
        this->x_lowerbound[i] = -0.5;
        this->x_upperbound[i] = 1.0;
    }

    //**************************************************************
    //* SET UPPER AND LOWER LIMITS OF CONSTRAINTS
    //**************************************************************
    this->g_lowerbound.resize(NG);
    this->g_upperbound.resize(NG);

    // the first constraint for each state veriable
    // refer to the initial state conditions
    // this will be initialized when solve() is called
    // the succeeding constraints refer to the relationship
    // between succeeding states based on our kinematic model of the system

    for (int i = 0; i < NG; ++i) {
        this->g_lowerbound[i] = 0.0;
        this->g_upperbound[i] = 0.0;
    }

    this->init_state = Eigen::VectorXd (4);
}

MPC::~MPC() {}

void MPC::LoadParams (const std::map<string, double> &params) {
    this->param_map = params;
    N = params.find("N") != params.end() ? params.at("N") : N;
    dt = params.find("dt") != params.end() ? params.at("dt") : dt;
    Lf = params.find("Lf") != params.end() ? params.at("Lf") : Lf;
    VELOCITY_MAX = params.find("VELOCITY_MAX") != params.end() ? params.at("VELOCITY_MAX") : VELOCITY_MAX;
    W_cte = params.find("W_cte") != params.end() ? params.at("W_cte") : W_cte;
    W_epsi = params.find("W_epsi") != params.end() ? params.at("W_epsi") : W_epsi;
    W_v = params.find("W_v") != params.end() ? params.at("W_v") : W_v;
    W_delta = params.find("W_delta") != params.end() ? params.at("W_delta") : W_delta;
    W_a = params.find("W_a") != params.end() ? params.at("W_a") : W_a;
    W_ddelta = params.find("W_ddelta") != params.end() ? params.at("W_ddelta") : W_ddelta;
    W_da = params.find("W_da") != params.end() ? params.at("W_da") : W_da;
    MAX_STEER = params.find("MAX_STEER") != params.end() ? params.at("MAX_STEER") : MAX_STEER;

    NX =  N * NUMBER_OF_STATES + (N - 1) * NUMBER_OF_ACTUATIONS; // number of state + actuation variables
    NG = N * NUMBER_OF_STATES; // number of constraints

// where the first element of each state variable is stored in the vector to be feeded the optimization algorithm
    ID_FIRST_px = 0;
    ID_FIRST_py = ID_FIRST_px + N;
    ID_FIRST_psi = ID_FIRST_py + N;
    ID_FIRST_v = ID_FIRST_psi + N;
    ID_FIRST_cte = ID_FIRST_v + N;
    ID_FIRST_epsi = ID_FIRST_cte + N;
    ID_FIRST_delta = ID_FIRST_epsi + N;
    ID_FIRST_a = ID_FIRST_delta + N - 1;

    //**************************************************************
    //* SET INITIAL VALUES OF VARIABLES
    //**************************************************************
    this->x.resize(NX);

    // all states except the ID_FIRST are set to zero
    // the aformentioned states will be initialized when solve() is called

    for (int i = 0; i < NX; ++i) {
        this->x[i] = 0.0;
    }

    //**************************************************************
    //* SET UPPER AND LOWER LIMITS OF VARIABLES
    //**************************************************************

    this->x_lowerbound.resize(NX);
    this->x_upperbound.resize(NX);

    // all other values large values the computer can handle
    for (int i = 0; i < ID_FIRST_delta; ++i) {
        this->x_lowerbound[i] = -1.0e10;
        this->x_upperbound[i] = 1.0e10;
    }

    // all actuation inputs (steering, acceleration) should have values between [-1, 1]
    for (int i = ID_FIRST_delta; i < ID_FIRST_a; ++i) {
        this->x_lowerbound[i] = - MAX_STEER;;
        this->x_upperbound[i] = MAX_STEER;;
    }

    for (int i = ID_FIRST_a; i < NX; ++i) {
        this->x_lowerbound[i] = -0.5;
        this->x_upperbound[i] = 1.0;
    }

    //**************************************************************
    //* SET UPPER AND LOWER LIMITS OF CONSTRAINTS
    //**************************************************************
    this->g_lowerbound.resize(NG);
    this->g_upperbound.resize(NG);

    // the first constraint for each state veriable
    // refer to the initial state conditions
    // this will be initialized when solve() is called
    // the succeeding constraints refer to the relationship
    // between succeeding states based on our kinematic model of the system

    for (int i = 0; i < NG; ++i) {
        this->g_lowerbound[i] = 0.0;
        this->g_upperbound[i] = 0.0;
    }

    this->init_state = Eigen::VectorXd (4);
}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd K, vector<double>& speed) {

    const double px = state[0];
    const double py = state[1];
    const double psi = state[2];
    const double v = state[3];
    const double cte = state[4];
    const double epsi = state[5];

    this->x[ID_FIRST_px] = px;
    this->x[ID_FIRST_py] = py;
    this->x[ID_FIRST_psi] = psi;
    this->x[ID_FIRST_v] = v;
    this->x[ID_FIRST_cte] = cte;
    this->x[ID_FIRST_epsi] = epsi;

    this->g_lowerbound[ID_FIRST_px] = px;
    this->g_lowerbound[ID_FIRST_py] = py;
    this->g_lowerbound[ID_FIRST_psi] = psi;
    this->g_lowerbound[ID_FIRST_v] = v;
    this->g_lowerbound[ID_FIRST_cte] = cte;
    this->g_lowerbound[ID_FIRST_epsi] = epsi;

    this->g_upperbound[ID_FIRST_px] = px;
    this->g_upperbound[ID_FIRST_py] = py;
    this->g_upperbound[ID_FIRST_psi] = psi;
    this->g_upperbound[ID_FIRST_v] = v;
    this->g_upperbound[ID_FIRST_cte] = cte;
    this->g_upperbound[ID_FIRST_epsi] = epsi;

    //**************************************************************
    //* SOLVE
    //**************************************************************

    // object that computes objective and constraints
    FG_eval fg_eval(K);
    fg_eval.LoadParams(param_map);
    fg_eval.setSpeed(speed);
    // options for IPOPT solver
    std::string options;
    // Uncomment this if you'd like more print information
    options += "Integer print_level  0\n";
    // NOTE: Setting sparse to true allows the solver to take advantage
    // of sparse routines, this makes the computation MUCH FASTER.
    options += "Sparse  true        forward\n";
    options += "Sparse  true        reverse\n";
    // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
    options += "Numeric max_cpu_time          0.5\n";

    // place to return solution
    CppAD::ipopt::solve_result<Dvector> solution;

    // solve the problem
    CppAD::ipopt::solve<Dvector, FG_eval>(
            options,
            x,
            x_lowerbound,
            x_upperbound,
            g_lowerbound,
            g_upperbound,
            fg_eval,
            solution);

    // comment out the lines below to debug!
    /*
    bool ok = true;
    auto cost = solution.obj_value;
    ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;
    if (ok) {
      std::cout << "OK! Cost:" << cost << std::endl;
    } else {
      std::cout << "SOMETHING IS WRONG!" << cost << std::endl;
    }
    */

    //**************************************************************
    //* STORE RELEVANT INFORMATION FROM SOLUTION
    //**************************************************************

    this->steer = solution.x[ID_FIRST_delta];
    this->throttle = solution.x[ID_FIRST_a];

    this->future_xs = {};
    this->future_ys = {};

    vector<double> res;

    for (int i = 0; i < N; ++i) {

        const double px = solution.x[ID_FIRST_px + i];
        const double py = solution.x[ID_FIRST_py + i];

        this->future_xs.emplace_back(px);
        this->future_ys.emplace_back(py);
    }
    res.push_back(steer);
    res.push_back(throttle);
    return res;
}

void MPC::updateDelayState(Eigen::VectorXd &predState, double v, double cte, double epsi, double steering,
                           double throttle) {
    const double px_act = v * dt;
    const double py_act = 0;
    const double psi_act = - v * tan(steering * MAX_STEER) * dt / Lf;
    const double v_act = v + throttle * dt;
    const double cte_act = cte + v * sin(epsi) * dt;
    const double epsi_act = epsi + psi_act;

    predState[0] = px_act;
    predState[1] = py_act;
    predState[2] = psi_act;
    predState[3] = v_act;
    predState[4] = cte_act;
    predState[5] = epsi_act;
}

vector<double> MPC::controlCB(vector<double> &cx, vector<double> &cy, vector<double> &v_ref, double steering,
                              double throttle) {
    // Waypoints related parameters
    const double x = init_state[0];
    const double y = init_state[1];
    const double psi = init_state[2];
    const double v = init_state[3];

//    TARGET_SPEED = 0.5;

//    while (psi < -M_PI) {
//        psi = psi + M_PI;
//    }
//    while (psi > M_PI) {
//        psi = psi - M_PI;
//    }

    const int N = cx.size(); // Number of waypoints
    const double cospsi = cos(-psi);
    const double sinpsi = sin(-psi);

    // Convert to the vehicle coordinate system
    Eigen::VectorXd x_veh(N);
    Eigen::VectorXd y_veh(N);

    // N is the number of ref way points
    double prex = 0;
    double dir = 1.0;

    for(int i = 0; i < N; i++)
    {
        const double dx = cx[i] - x;
        const double dy = cy[i] - y;
        x_veh[i] = dx * cospsi - dy * sinpsi;
        y_veh[i] = dy * cospsi + dx * sinpsi;
        if (i >= 1) {
            if (x_veh[i] - prex< 0) {
                dir = - 1.0;
            }
        }
        prex = x_veh[i];
//        cout << " x_veh[i] " << x_veh[i] << endl;
//        cout << " y_veh[i] " << y_veh[i] << endl;
    }


    // Fit waypoints
    auto coeffs = polyfit(x_veh, y_veh, 3);
//    cout << " coeffs : " << endl;
//    cout << coeffs << endl;
    const double cte  = polyeval(coeffs, 0.0);
    const double epsi = - atan(coeffs[1]);

//    cout << " cte input : " << cte << endl;
//    cout << " epsi input : " << epsi << endl;
//    cout << coeffs << endl;
    Eigen::VectorXd State(6);
    if (DELAY_MODE)
    {
        // Kinematic model is used to predict vehicle state at the actual
        // moment of control (current time + delay dt)
        updateDelayState(State, v, cte, epsi, steering, throttle);
    } else {
        State << 0, 0, psi, v, cte, epsi;
    }

    this->e_cte = State[4];
    this->e_psi = State[5];

    // Solve MPC Problem
    vector<double> mpc_results = Solve(State, coeffs, v_ref);

    // MPC result (all described in car frame)
    double o_steering = mpc_results[0]; // radian
    double o_throttle = mpc_results[1]; // acceleration
    double o_speed = v + throttle*dt;  // speed

    if(o_speed <= 0.0)
        o_speed = 0.0;

    vector<double> res;
//    res.push_back(o_steering);

    res.push_back(o_steering);
    res.push_back(o_throttle);
//    res.push_back(0.5);
//    res.push_back(0);

    // debug
    if(DEBUG_INFO) {
        cout << "\n\nDEBUG" << endl;
        cout << "psi: " << psi << endl;
        cout << "V: " << v << endl;
//        cout << "coeffs: \n" << coeffs << endl;
        cout << "_steering: \n" << o_steering << endl;
        cout << "_throttle: \n" << o_throttle << endl;
        cout << "_speed: \n" << o_speed << endl;
    }
    return res;
}

// Evaluate a polynomial.
double MPC::polyeval(Eigen::VectorXd coeffs, double x)
{
    double result = 0.0;
    for (int i = 0; i < coeffs.size(); i++)
    {
        result += coeffs[i] * pow(x, i);
    }
    return result;
}

Eigen::VectorXd MPC::polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
    assert(xvals.size() == yvals.size());
    assert(order >= 1 && order <= xvals.size() - 1);
    Eigen::MatrixXd A(xvals.size(), order + 1);

    for (int i = 0; i < xvals.size(); i++) {
        A(i, 0) = 1.0;
    }

    for (int j = 0; j < xvals.size(); j++) {
        for (int i = 0; i < order; i++) {
            A(j, i + 1) = A(j, i) * xvals(j);
        }
    }

    auto Q = A.householderQr();
    auto result = Q.solve(yvals);

    for(int i = 0; i < result.size(); i ++) {
        if (isnan(result[i]) || isinf(result[i])) {
            cout << " path cannot be turn to trajectory : " << endl;
            Eigen::VectorXd res(4);
            res << 0, 0, 0, 0;
            return res;
        }
    }

    return result;
}

