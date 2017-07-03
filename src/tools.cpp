#include <iostream>
#include "tools.h"
#include <functional>
#include <numeric>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

//static
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
  const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
    std::cout << "Invalid estimation or ground_truth data" << std::endl;
    return rmse;
  }

  auto op2 = [](const VectorXd& a, const VectorXd& b)->VectorXd {
    const VectorXd residual = a - b;
    return residual.array()*residual.array();
  };

  rmse = std::inner_product(estimations.begin(), estimations.end(), ground_truth.begin(), rmse, std::plus<VectorXd>(), op2);
  rmse /= estimations.size();
  rmse = rmse.array().sqrt();

  return rmse;
}

//static
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3, 4);
  //recover state parameters
  const double px = x_state(0);
  const double py = x_state(1);
  const double vx = x_state(2);
  const double vy = x_state(3);

  //pre-compute a set of terms to avoid repeated calculation
  const double c1 = px*px + py*py;
  const double c2 = sqrt(c1);
  const double c3 = (c1*c2);

  //check division by zero
  if (fabs(c1) < 1e-8) {
    std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
    return Hj;
  }

  //compute the Jacobian matrix
  Hj << (px / c2), (py / c2), 0, 0,
    -(py / c1), (px / c1), 0, 0,
    py*(vx*py - vy*px) / c3, px*(px*vy - py*vx) / c3, px / c2, py / c2;

  return Hj;
}
