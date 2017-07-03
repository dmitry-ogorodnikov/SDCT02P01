#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_*x_;
  P_ = F_*P_*F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  const VectorXd y = z - H_*x_;
  /*Si - S inverse*/
  const MatrixXd Si = (H_*P_*H_.transpose() + R_).inverse();
  const MatrixXd K = P_*H_.transpose()*Si;
  
  x_ = x_ + K*y; 
  P_ = (MatrixXd::Identity(x_.size(), x_.size()) - K*H_)*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  VectorXd y = z - cart2Polar(x_);
  y(1) = normalizeAngle(y(1));

  /*Si - S inverse*/
  const MatrixXd Si = (H_*P_*H_.transpose() + R_).inverse();
  const MatrixXd K = P_*H_.transpose()*Si;

  x_ = x_ + K*y;
  P_ = (MatrixXd::Identity(x_.size(), x_.size()) - K*H_)*P_;
}

VectorXd KalmanFilter::cart2Polar(const VectorXd& x_state) {
  VectorXd result(3);
  result << 0, 0, 0;

	const double px = x_state(0);
	const double py = x_state(1);
	const double vx = x_state(2);
	const double vy = x_state(3);
  
  const double rho = sqrt(px*px + py*py);
  const double phi = atan2(py, px);

  //errno == EDOM if px and py are zero.
  if(errno == EDOM) {
    std::cout << "Error in the range of valid values." << std::endl;
    return result;
  }
  const double rhoDot = (px*vx + py*vy) / sqrt(px*px + py*py);
  result << rho, phi, rhoDot;
  return result;
}

double KalmanFilter::normalizeAngle(const double angle) {
  auto result = angle;
  while (result < -M_PI) {
    result += 2 * M_PI;
  }
  while (result > M_PI) {
    result -= 2 * M_PI;
  }
  return result;
}
