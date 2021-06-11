//
// Created by jingche on 2020/8/26.
//
#include <iostream>
#include "MPC.h"
#include <ceres/ceres.h>

using namespace Eigen;
using namespace std;

double GOAL_DIS = 1.5;
//double DT = 0.1; // DELAY TIME
double DL = 1.0; // course tick
double WB = 2.67;
double STOP_SPEED = 0.5/3.6;
double MAX_TIME = 100;
double DT = 0.1;
int PATH_STEP = 30;
//int N_IND_SEARCH = 10;
double MAX_STEER = (45.0)/180*M_PI;

double MAX_SPEED = (55.0 / 3.6);
double MIN_SPEED = -20.0 / 3.6;

using namespace ceres;

void updateState(VectorXd &state, vector<double> vector);
double deg2rad(double angle) {
    return angle * M_PI / 180.0;
}

static double pi_2_pi(double angle) {
    while (angle > M_PI) {
        angle -= 2.0 * M_PI;
    }
    while (angle < -M_PI) {
        angle += 2.0 * M_PI;
    }
    return angle;
}

/**
 * @param state
 * @param cx
 * @param cy
 * @param cyaw
 * @param pind
 * @return ind & mind (mind not used)
 */


void readTargets(string filename, vector<double>& input) {
    vector<double>::iterator it;
    ifstream data("curve/" + filename);
    double d;
    while (data >> d)
        input.push_back(d);//将数据压入堆栈。//
    data.close();
}

void writeResults(string filename, vector<double>& input) {

    ofstream outfile("trajectory/" + filename, std::ios::app);
    for (int i = 0; i < input.size(); ++i) {
        outfile<< input[i]<<endl;
    }
    outfile.close();
}


bool check_goal(VectorXd state, double x, double y, int tind, int nind) {

    double dx = state[0] - x;
    double dy = state[1] - y;
    double d = hypot(dx, dy);

    bool isGoal = (d <= GOAL_DIS);
    if (abs(tind - nind) >= 7) {
        isGoal = false;
    }

//    bool isStop = (abs(state.v) <= STOP_SPEED);

    if (isGoal) {
        return true;
    }
    return false;
}


/**
 * @param state
 * @param cx
 * @param cy
 * @param cyaw
 * @param pind
 * @return ind & mind (mind not used)
 */
int calc_nearest_index (VectorXd state, vector<double> &cx, vector<double> &cy, int pind) {

    int mind = 0;
    double smallest = (double) INT_MAX;
    double x = state[0];
    double y = state[1];
//    cout << "State X is : "<< x << endl;
//    cout << "State Y is : "<< y << endl;
    for (int i = 0; i < PATH_STEP; i++) {
        double d = pow(x - cx[pind + i], 2) + pow(y - cy[pind + i], 2);
//        cout << "State X is : "<< state.x << endl;
        if (d < smallest) {
            smallest = d;
            mind = i;
        }
    }
//    double angle = pi_2_pi(cyaw.at(ind) - atan2(dyl, dxl));
//    if (angle < 0) {
//        mind *= -1;
//    }
    return mind + pind;
}

int calc_ref_trajectory (VectorXd state,  vector<double> &cx, vector<double> &cy, vector<double> &sp, int pind, vector<double>& x_ref, vector<double>& y_ref, vector<double>& v_ref) {
//    const double x = state[0];
//    const double y = state[1];
//    const double psi = state[2];
//    const double v = state[3];

    int ncourse = cx.size();
    int ind = calc_nearest_index(state, cx, cy, pind);
//    cout << "pre index : " << pind << endl;
//    cout << "nearest index : " << ind << endl;
    if (pind > ind) {
        ind = pind;
    }

    for (int i = 0; i < PATH_STEP; ++i) {
        int dind = int(ind + i);
//        cout<<"ref x val : " << x_ref[i] << endl;

        if ((dind) < ncourse) {
            x_ref[i] = cx[dind];
//            cout<<"ref x val : " << x_ref[i] << endl;
            y_ref[i] = cy[dind];
//            cout<<"ref y val : " << y_ref[i] << endl;
            v_ref[i] = sp[dind];
        } else {
            x_ref[i] = cx[ncourse - 1];
            y_ref[i] = cy[ncourse - 1];
            v_ref[i] = sp[ncourse - 1];
        }
    }
//    cout << "process : " << double(ind)/double(ncourse)*100<<" %"<<endl;
    return ind;
}

vector<double> doSimulation(double *x) {

    MPC mpc;
    vector<double> result;

    // Convert to the vehicle coordinate system
    vector<double> cx;
    vector<double> cy;
    vector<double> cyaw;
    vector<double> ck;
    vector<double> sp;

    readTargets("cxs", cx);
    readTargets("cys", cy);

    vector<double> recorded_t;
    vector<double> recorded_x;
    vector<double> recorded_y;
    vector<double> recorded_v;
    vector<double> recorded_yaw;
    vector<double> recorded_d;
    vector<double> recorded_a;

    double steering = 0.0;
    double throttle = 0.0;
    double t = 0.0;

    map<string, double> _mpc_params;
//    _mpc_params["W_cte"] = 1473.72;
//    _mpc_params["W_epsi"] = 1519.63;
//    _mpc_params["W_v"] = 5.07765;
//    _mpc_params["W_delta"] = 7.06496;
//    _mpc_params["W_a"] = 9.88325;
//    _mpc_params["W_ddelta"] = 12.791;
//    _mpc_params["W_da"] = 16.175;
//    _mpc_params["N"] = 15;
//    _mpc_params["Lf"] = 2.67;

    _mpc_params["W_cte"] = x[0];
    _mpc_params["W_epsi"] = x[1];
    _mpc_params["W_v"] = x[2];
    _mpc_params["W_delta"] = x[3];
    _mpc_params["W_a"] = x[4];
    _mpc_params["W_ddelta"] = x[5];
    _mpc_params["W_da"] = x[6];
    _mpc_params["N"] = x[7];
    _mpc_params["Lf"] = 2.67;

    mpc.LoadParams(_mpc_params);
    mpc.init_state << cx[0], cy[0], deg2rad(0), 15;
    double e_cte_error_sum = 0;
    double e_psi_error_sum = 0;

    vector<double> x_ref(PATH_STEP);
    vector<double> y_ref(PATH_STEP);
    vector<double> v_ref(PATH_STEP);

    // design ref speed
    for (int i = 0; i < cx.size(); ++i) {
        double speed =  10;
        sp.push_back(speed);
    }
    int target_ind = calc_nearest_index(mpc.init_state, cx, cy, 0);
    while (t <= MAX_TIME) {
//        cout << " target_ind " << target_ind << endl;
        target_ind = calc_ref_trajectory(mpc.init_state, cx, cy, sp, target_ind, x_ref, y_ref, v_ref);
        vector<double> res = mpc.controlCB(x_ref, y_ref, v_ref, steering, throttle);
//        cout << "steering : " << res[0] << endl;
        e_cte_error_sum = e_cte_error_sum + std::abs(mpc.e_cte);
        e_psi_error_sum = e_psi_error_sum + std::abs(mpc.e_psi);

        updateState(mpc.init_state, res);
        recorded_x.push_back(mpc.init_state[0]);
        recorded_y.push_back(mpc.init_state[1]);
        recorded_v.push_back(mpc.init_state[3]);
        t += DT;
        if (target_ind + 5 >= cx.size() - 1) {
            cout << " Reached Goal " << endl;
            break;
        }
    }
    cout << "One loop over : " << endl;
    // (VectorXd& init_state, vector<double>& cx, vector<double>& cy, double steering, double throttle, MPC& mpc)

//    system("rm trajectory/recorded_x");
//    system("rm trajectory/recorded_y");
//    system("rm trajectory/recorded_v");
//    system("rm trajectory/sp");
//    writeResults("sp", sp);
//    writeResults("recorded_v", recorded_v);
//    writeResults("recorded_x", recorded_x);
//    writeResults("recorded_y", recorded_y);

    vector<double> res;
    res.push_back(e_cte_error_sum);
    res.push_back(e_psi_error_sum);
    return res;
}

void updateState(VectorXd& state, vector<double> actuator) {
    // vector 0: steering 1: throttle  2: speed

    const double x = state[0];
    const double y = state[1];
    const double psi = state[2];
    const double v = state[3];

//    while (psi < - M_PI) {
//        psi += M_PI;
//    }
//
//    while (psi > M_PI) {
//        psi -= M_PI;
//    }

//    cout << "psi is : " << psi << endl;
//    cout << "v is : " << v << endl;

    double  steering = actuator[0];
    double  throttle = actuator[1];
//    cout << "throttle is : " << throttle << endl;
//    cout << "steering is : " << steering << endl;
//    cout << "MAX_STEER : " << MAX_STEER << endl;

    if (steering >= MAX_STEER) {
        steering = MAX_STEER;
    } else if (steering <= -MAX_STEER) {
        steering = -MAX_STEER;
    }

    double predX = x + v * std::cos(psi) * DT;
    double predY = y + v * std::sin(psi) * DT;
    double predYaw = psi - v * (steering) * DT / WB;
    double predV = v + throttle * DT;

//    cout << "predYaw : " << predYaw << endl;

    if (predV > MAX_SPEED) {
        predV = MAX_SPEED;
    } else if (predV < MIN_SPEED) {
        predV = MIN_SPEED;
    }
    state[0] = predX;
    state[1] = predY;
    state[2] = predYaw;
    state[3] = predV;
}

struct NumericDiffCostFunctor {
    bool operator()(const double* const x, double* residual) const {
        double x_i[8];
        for (int i = 0; i < 8; ++i) {
            x_i[i] =  x[i];
        }
        vector<double> res = doSimulation(x_i);
//      ceres::Jet<double, 8> ass;
//      ass.v[0];
        residual[0] = (res[0]);
        residual[1] = (res[1]);
        std::cout<< " residual[0] : " << residual[0] << std::endl;
        std::cout<< " residual[1] : " << residual[1] << std::endl;
        for(int i = 0; i < 8; i++) {
            std::cout<< " x : " << x[i] << std::endl;
        }
        std::cout<< " one loop end : " << std::endl;
        std::cout<< "********************* " << std::endl;
        return true;
    }
};

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    double x_init[8] = {1500, 1500, 5.0, 10, 10, 10, 15, 15};
    Problem problem;
    CostFunction* cost_function =
            new NumericDiffCostFunction<NumericDiffCostFunctor, CENTRAL, 2, 8>(new NumericDiffCostFunctor); //使用自动求导，将之前的代价函数结构体传入，第一个1是输出维度，即残差的维度，第二个1是输入维度，即待寻优参数x的维度。
    problem.AddResidualBlock(cost_function, NULL, x_init);
//    _mpc_params["W_cte"] = 1500;
//    _mpc_params["W_epsi"] = 1500;
//    _mpc_params["W_v"] = 5.0;
//    _mpc_params["W_delta"] = 10;
//    _mpc_params["W_a"] = 10;
//    _mpc_params["W_ddelta"] = 10;
//    _mpc_params["W_da"] = 15;
//    _mpc_params["N"] = 15;
//    _mpc_params["Lf"] = 2.67;
    problem.SetParameterLowerBound(x_init, 0, 1000);
    problem.SetParameterLowerBound(x_init, 1, 1000);
    problem.SetParameterLowerBound(x_init, 2, 3.0);
    problem.SetParameterLowerBound(x_init, 3, 5);
    problem.SetParameterLowerBound(x_init, 4, 5);
    problem.SetParameterLowerBound(x_init, 5, 5);
    problem.SetParameterLowerBound(x_init, 6, 10);
    problem.SetParameterLowerBound(x_init, 7, 10);

    problem.SetParameterUpperBound(x_init, 0, 2000);
    problem.SetParameterUpperBound(x_init, 1, 2000);
    problem.SetParameterUpperBound(x_init, 2, 10.0);
    problem.SetParameterUpperBound(x_init, 3, 15);
    problem.SetParameterUpperBound(x_init, 4, 15);
    problem.SetParameterUpperBound(x_init, 5, 15);
    problem.SetParameterUpperBound(x_init, 6, 20);
    problem.SetParameterUpperBound(x_init, 7, 20);

    //第三部分： 配置并运行求解器
    Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR; //配置增量方程的解法
    options.minimizer_progress_to_stdout = true;//输出到cout
    Solver::Summary summary;//优化信息
    ceres::Solve(options, &problem, &summary);//求解!!!
    std::cout << summary.BriefReport() << "\n";//输出优化的简要信息
    for(int i = 0; i < 8; i++) {
        std::cout<< " x : " << x_init[i] << std::endl;
    }
    return 0;
}