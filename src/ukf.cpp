#include "ukf.h"
#include "measurement_package.h"

#include "Eigen/Dense"
#include <iostream>

//State dimension
const int UKF::N_X   = 5;

//Augmented state dimension
const int UKF::N_AUG = 7;

//Sigma point spreading parameter
const double UKF::LAMBDA = (3.0 - UKF::N_X);

//Process noise values

//Process noise standard deviation longitudinal acceleration in m/s^2
const double UKF::STD_A      = 30.0; //@TODO: need to tune and use a different smaller value!

//Process noise standard deviation yaw acceleration in rad/s^2
const double UKF::STD_YAW_DD = 30.0; //@TODO: need to tune and use a different smaller value!

//Measurement noise values

//Laser measurement noise standard deviation position1 in m
const double UKF::STD_LASER_PX = 0.15;

//Laser measurement noise standard deviation position2 in m
const double UKF::STD_LASER_PY = 0.15;

//Radar measurement noise standard deviation radius in m
const double UKF::STD_RAD_R    = 0.3;

//Radar measurement noise standard deviation angle in rad
const double UKF::STD_RAD_PHI  = 0.03;

//Radar measurement noise standard deviation radius change in m/s
const double UKF::STD_RAD_R_D  = 0.3;


UKF::UKF()
:is_initialized_(false),
 time_us_(0),
 use_laser_(true),
 use_radar_(true)
{
    x_ = VectorXd(5);

    x_ << 0.0, 0.0, 0.0, 0.0, 0.0;

    P_ = MatrixXd(5, 5);

    //@TODO: Initialise weights_, Xsig_pred_ with correct dimensions
}

UKF::~UKF()
{
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage& meas_package)
{
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  std::cout << "Testing Calling ProcessMeasurement!" << std::endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t)
{
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage& meas_package)
{
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const MeasurementPackage& meas_package)
{
  /**
  TODO:


  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}

const VectorXd& UKF::GetStates() const
{
    return x_;
}

