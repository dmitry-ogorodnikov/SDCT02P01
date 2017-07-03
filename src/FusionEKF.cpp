#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd::Zero(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  ekf_.x_ = VectorXd::Zero(4);
  ekf_.P_ = MatrixXd::Identity(4, 4) * 1000;
  ekf_.F_ = MatrixXd::Identity(4, 4);
  ekf_.Q_ = MatrixXd::Zero(4, 4);
  H_laser_(0, 0) = 1;
  H_laser_(1, 1) = 1;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    std::cout << "EKF: " << std::endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      const auto rho = measurement_pack.raw_measurements_[0];
      const auto phi = measurement_pack.raw_measurements_[1];
      ekf_.x_(0) = rho*cos(phi);
      ekf_.x_(1) = rho*sin(phi);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_(0) = measurement_pack.raw_measurements_[0];
      ekf_.x_(1) = measurement_pack.raw_measurements_[1];
    }

    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  const double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
  const double dt_2 = dt*dt;
  const double dt_3 = dt_2*dt;
  const double dt_4 = dt_3*dt;

  const double noise_ax = 9;
  const double noise_ay = 9;

  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  ekf_.Q_ << dt_4*noise_ax/4, 0, dt_3*noise_ax/2, 0,
              0, dt_4*noise_ay/4, 0, dt_3*noise_ay/2,
              dt_3*noise_ax/2, 0, dt_2*noise_ax, 0,
              0, dt_3*noise_ay/2, 0, dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.H_ = Tools::CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  std::cout << "x_ = " << ekf_.x_ << std::endl;
  std::cout << "P_ = " << ekf_.P_ << std::endl;
}
