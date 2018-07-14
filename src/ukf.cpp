#include "ukf.h"
#include "measurement_package.h"
#include "tools.h"

#include "Eigen/Dense"
#include <iostream>

//State dimension
const int UKF::N_X   = 5;

//Augmented state dimension
const int UKF::N_AUG = 7;

// 2 * N_AUG + 1
const int UKF::N_AUG_2_PLUS_1 = (2 * UKF::N_AUG + 1);

//Sigma point spreading parameter
const double UKF::LAMBDA = (3.0 - UKF::N_AUG);

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

    weights_ = VectorXd(N_AUG_2_PLUS_1);

    double weight_multiplier = LAMBDA;
    double weight_denom = 1.0 / (LAMBDA + N_AUG);

    for (int i = 0; i < N_AUG_2_PLUS_1 ; ++i)
    {
        if(i > 0)
        {
            weight_multiplier = 0.5;
        }

        weights_(i) = weight_multiplier * weight_denom;
    }

    Xsig_pred_ = MatrixXd(N_X, N_AUG_2_PLUS_1);
    Xsig_pred_.fill(0.0);

    //@TODO: Adding matrices R_laser and R_radar, Q?
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
            double delta_t = 1.0e-6 * (meas_package.timestamp_ - time_us_);

            if((data_option_ != LASER_DATA_ONLY) && (meas_package.sensor_type_ == MeasurementPackage::RADAR))
            {
                //Data option is using Radar data or both

                ret = Prediction(delta_t);
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

                ret = Prediction(delta_t);
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
    is_initialized_ = false;

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
        //Convert radar from polar to cartesian coordinates and initialize state.
        double ro    = meas_package.raw_measurements_(0);
        double theta = Tools::UnwrapAngle(meas_package.raw_measurements_(1));

        x_(0) = ro * cos(theta);
        x_(1) = ro * sin(theta);
        x_(2) = 0.0;
        x_(3) = 0.0;
        x_(4) = 0.0;

        std::cout << "Fusion data option: " << ToString(data_option_) << std::endl;

        is_initialized_ = true;
        time_us_ = meas_package.timestamp_;

        std::cout << "UKF initialised with Radar data:" << std::endl;
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
        //Initialize state.
        x_(0) = meas_package.raw_measurements_(0);
        x_(1) = meas_package.raw_measurements_(1);
        x_(2) = 0.0;
        x_(3) = 0.0;
        x_(4) = 0.0;

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
        std::cout << x_ << std::endl;
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

    std::cout << "Test Prediction, dT: " << delta_t << std::endl;

    MatrixXd x_augmented_sigma;

    GenerateAugmentedSigmaPoints(x_augmented_sigma);
    PredictSigmaPoints(x_augmented_sigma, delta_t);
    PredictMeanAndCovariance();

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

void UKF::GenerateAugmentedSigmaPoints(MatrixXd& Xsig_aug) const
{
    //sigma point matrix
    Xsig_aug = MatrixXd(N_AUG, N_AUG_2_PLUS_1);

    //augmented mean vector
    VectorXd x_aug = VectorXd(N_AUG);

    x_aug.head(N_X) = x_;
    x_aug(N_AUG - 2) = 0.0;
    x_aug(N_AUG - 1) = 0.0;

    //augmented state covariance
    MatrixXd P_aug = MatrixXd(N_AUG, N_AUG);

    P_aug.fill(0.0);
    P_aug.topLeftCorner(N_X, N_X) = P_;
    P_aug((N_AUG - 2), (N_AUG - 2)) = STD_A * STD_A;
    P_aug((N_AUG - 1), (N_AUG - 1)) = STD_YAW_DD * STD_YAW_DD;

    //create square root matrix
    MatrixXd L = P_aug.llt().matrixL();
    L *= sqrt(LAMBDA + N_AUG);

    //create augmented sigma points
    Xsig_aug.col(0)  = x_aug;

    for (int i = 0; i < N_AUG; ++i)
    {
        Xsig_aug.col(i + 1)         = x_aug + L.col(i);
        Xsig_aug.col(i + 1 + N_AUG) = x_aug - L.col(i);
    }

    //print result
    std::cout << "Xsig_aug: " << std::endl << Xsig_aug << std::endl;
}

void UKF::PredictSigmaPoints(const MatrixXd& Xsig_aug, double delta_t)
{
    //predict sigma points
    for (int i = 0; i < N_AUG_2_PLUS_1; ++i)
    {
        //extract values for better readability
        double p_x      = Xsig_aug(0,i);
        double p_y      = Xsig_aug(1,i);
        double v        = Xsig_aug(2,i);
        double yaw      = Xsig_aug(3,i);
        double yawd     = Xsig_aug(4,i);
        double nu_a     = Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);

        //predicted state values
        double px_p, py_p;

        //avoid division by zero
        if( fabs(yawd) > 0.001)
        {
            px_p = p_x + v / yawd * ( sin (yaw + yawd * delta_t) - sin(yaw));
            py_p = p_y + v / yawd * ( cos(yaw) - cos(yaw + yawd * delta_t) );
        }
        else
        {
            //If dyaw is very small
            px_p = p_x + v * delta_t * cos(yaw);
            py_p = p_y + v * delta_t * sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawd * delta_t;
        double yawd_p = yawd;

        //add noise from analytical model of Q wity noise_a and noise_yawdd
        px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
        py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
        v_p  = v_p + nu_a * delta_t;

        yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
        yawd_p = yawd_p + nu_yawdd * delta_t;

        //write predicted sigma point into right column
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        Xsig_pred_(4,i) = yawd_p;
    }

    //print result
    std::cout << "Xsig_pred: " << std::endl << Xsig_pred_ << std::endl;
}

void UKF::PredictMeanAndCovariance()
{
    //predicted state mean
    for(int i = 0; i < N_AUG_2_PLUS_1; ++i)
    {  //iterate over predicted sigma points
        x_ = x_ + weights_(i) * Xsig_pred_.col(i);
    }

    //predicted state covariance matrix
    for(int i = 0; i < N_AUG_2_PLUS_1; ++i)
    { //iterate over predicted sigma points

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;

        double x_diff_3 = x_diff(3);
        x_diff(3) = Tools::UnwrapAngle(x_diff_3);

        P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
    }

    //print result
    std::cout << "Predicted X:" << std::endl << x_ << std::endl;
    std::cout << "Predicted P:" << std::endl << P_ << std::endl;
}


