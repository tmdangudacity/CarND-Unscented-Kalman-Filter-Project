#include "ukf.h"
#include "measurement_package.h"
#include "tools.h"

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

UKF::UKF(DataOption in_data)
:data_option_(in_data),
 is_initialized_(false),
 time_us_(0)
{
    x_ = VectorXd(5);
    x_ << 0.0, 0.0, 0.0, 0.0, 0.0;

    P_ = MatrixXd(5, 5);
    P_ << 1.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 1.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 1.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 1.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 1.0;

    //@TODO: Initialise weights_, Xsig_pred_ with correct dimensions
    //Adding matrices R_laser and R_radar, Q?
}

UKF::~UKF()
{
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
bool UKF::ProcessMeasurement(const MeasurementPackage& meas_package)
{
    std::cout << "Test running ProcessMeasurement!" << std::endl;

    bool ret = false;

    //If not initialised yet
    if (!is_initialized_)
    {
        ret = Initialise(meas_package);
    }
    else
    {
        //If already initialised:
        ret = (meas_package.timestamp_ > time_us_);

        if(ret)
        {
            double dt = 1.0e-6 * (meas_package.timestamp_ - time_us_);

            if((data_option_ != LASER_DATA_ONLY) && (meas_package.sensor_type_ == MeasurementPackage::RADAR))
            {
                //Data option is using Radar data or both

                ret = Prediction(dt);
                if(ret)
                {
                    time_us_ = meas_package.timestamp_;
                    UpdateRadar(meas_package);
                }

                std::cout << "x_ =" << std::endl;
                std::cout << x_     << std::endl;

                std::cout << "P_ =" << std::endl;
                std::cout << P_     << std::endl;
            }
            else if((data_option_ != RADAR_DATA_ONLY) && (meas_package.sensor_type_ == MeasurementPackage::LASER))
            {
                //Data option is using Laser data or both

                ret = Prediction(dt);
                if(ret)
                {
                    time_us_ = meas_package.timestamp_;
                    UpdateLidar(meas_package);
                }

                std::cout << "x_ =" << std::endl;
                std::cout << x_     << std::endl;

                std::cout << "P_ =" << std::endl;
                std::cout << P_     << std::endl;
            }
            else
            {
                std::cout << "ERROR, sensor type: " << meas_package.sensor_type_
                          << ", Data option: "      << ToString(data_option_)
                          << std::endl;
            }
        }
        else
        {
            std::cout << "ERROR, Data timestamp: " << meas_package.timestamp_
                      << " < previous timestamp: " << time_us_
                      << std::endl;
        }
    }

    return ret;
}

/**
 * Initialisation of the state with measurement data
 */
bool UKF::Initialise(const MeasurementPackage& meas_package)
{
    //Initialize the state ekf_.x_ with the first measurement.
    //Create the covariance matrix.

    //@TODO: initialise state vector x_(5) instead of (4)
    // first measurement

    is_initialized_ = false;
    VectorXd x_init = VectorXd(4);

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
        //Convert radar from polar to cartesian coordinates and initialize state.
        float ro     = meas_package.raw_measurements_(0);
        float theta  = Tools::UnwrapAngle(meas_package.raw_measurements_(1));
        float ro_dot = meas_package.raw_measurements_(2);

        x_init(0) = ro * cos(theta);
        x_init(1) = ro * sin(theta);
        x_init(2) = ro_dot * cos(theta);
        x_init(3) = ro_dot * sin(theta);
        std::cout << "Fusion data option: " << ToString(data_option_) << std::endl;

        is_initialized_ = true;
        time_us_ = meas_package.timestamp_;

        std::cout << "UKF initialised with Radar data:" << std::endl;
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
        //Initialize state.
        x_init(0) = meas_package.raw_measurements_(0);
        x_init(1) = meas_package.raw_measurements_(1);
        x_init(2) = 1.0;
        x_init(3) = 1.0;

        is_initialized_ = true;
        time_us_ = meas_package.timestamp_;

        std::cout << "UKF initialised with Laser data:" << std::endl;
    }
    else
    {
        std::cout << "Initialisation ERROR. Unknown sensor type!" << std::endl;
    }

    if(is_initialized_)
    {
        std::cout << x_init << std::endl;
        std::cout << "Fusion data option: " << ToString(data_option_) << std::endl;
    }

    return is_initialized_;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
bool UKF::Prediction(double delta_t)
{
    bool ret = true;

    /**
    TODO:
    Complete this function! Estimate the object's location. Modify the state
    vector, x_. Predict sigma points, the state, and the state covariance matrix.
    */

    std::cout << "Test Prediction, dT: " << delta_t << std::endl;

    return ret;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
bool UKF::UpdateLidar(const MeasurementPackage& meas_package)
{
    bool ret = true;

    /**
    TODO:
    Complete this function! Use lidar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.
    You'll also need to calculate the lidar NIS.
    */

    std::cout << "Test Update with Lidar data" << std::endl;

    return ret;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
bool UKF::UpdateRadar(const MeasurementPackage& meas_package)
{
    bool ret = true;

    /**
    TODO:
    Complete this function! Use radar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.
    You'll also need to calculate the radar NIS.
    */

    std::cout << "Test Update with Radar data" << std::endl;

    return ret;
}

const VectorXd& UKF::GetStates() const
{
    return x_;
}

std::string UKF::ToString(DataOption data_option)
{
    std::string ret = "UNKNOWN_DATA_OPTION";

    switch(data_option)
    {
        case LASER_DATA_ONLY:
            ret = "LASER DATA ONLY";
            break;

        case RADAR_DATA_ONLY:
            ret = "RADAR DATA ONLY";
            break;

        case LASER_AND_RADAR_DATA:
            ret = "LASER AND RADAR DATA";
            break;

        default:
            break;
    }

    return ret;
}

