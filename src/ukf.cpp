#include "ukf.h"
#include "measurement_package.h"
#include "tools.h"

#include "Eigen/Dense"
#include <iostream>

//State dimension
const int UKF::N_X   = 5;

//Laser measurement dimension
const int UKF::N_Z_LASER = 2;

//Radar measurement dimension
const int UKF::N_Z_RADAR = 3;

//Augmented state dimension
const int UKF::N_AUG = 7;

// 2 * N_AUG + 1
const int UKF::N_AUG_2_PLUS_1 = (2 * UKF::N_AUG + 1);

//Sigma point spreading parameter
const double UKF::LAMBDA = (3.0 - UKF::N_AUG);

//Process noise values

//Process noise standard deviation longitudinal acceleration in m/s^2
const double UKF::STD_A      = 0.1; //30.0; //@TODO: need to tune and use a different smaller value!

//Process noise standard deviation yaw acceleration in rad/s^2
const double UKF::STD_YAW_DD = 0.1; //30.0; //@TODO: need to tune and use a different smaller value!

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

    //Laser measurement noise covariance matrix
    R_LASER = MatrixXd(N_Z_LASER, N_Z_LASER);
    R_LASER << (STD_LASER_PX * STD_LASER_PX),  0.0,
                0.0, (STD_LASER_PY * STD_LASER_PY);

    //Radar measurement noise covariance matrix
    R_RADAR = MatrixXd(N_Z_RADAR, N_Z_RADAR);
    R_RADAR << (STD_RAD_R * STD_RAD_R),   0.0,   0.0,
                0.0, (STD_RAD_PHI * STD_RAD_PHI), 0.0,
                0.0, 0.0, (STD_RAD_R_D * STD_RAD_R_D);
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
    std::cout << "ProcessMeasurement" << std::endl;

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
                time_us_ = meas_package.timestamp_;

                Prediction(delta_t);
                UpdateRadar(meas_package);
            }
            else if((data_option_ != RADAR_DATA_ONLY) && (meas_package.sensor_type_ == MeasurementPackage::LASER))
            {
                //Data option is using Laser data or both
                time_us_ = meas_package.timestamp_;

                Prediction(delta_t);
                UpdateLidar(meas_package);
            }
            else
            {
                ret = false;
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

    std::cout << "Update with Lidar data" << std::endl;
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
    std::cout << "- Xsig_aug: " << std::endl << Xsig_aug << std::endl;
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
            px_p = p_x + (v / yawd) * ( sin (yaw + yawd * delta_t) - sin(yaw));
            py_p = p_y + (v / yawd) * ( cos(yaw) - cos(yaw + yawd * delta_t) );
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
    std::cout << "- Xsig_pred: " << std::endl << Xsig_pred_ << std::endl;
}

void UKF::PredictMeanAndCovariance()
{
    //predicted state mean
    x_.fill(0.0);
    for(int i = 0; i < N_AUG_2_PLUS_1; ++i)
    {  //iterate over predicted sigma points
        x_ = x_ + weights_(i) * Xsig_pred_.col(i);
    }

    //predicted state covariance matrix
    P_.fill(0.0);
    for(int i = 0; i < N_AUG_2_PLUS_1; ++i)
    { //iterate over predicted sigma points

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;

        double x_diff_3 = x_diff(3);
        x_diff(3) = Tools::UnwrapAngle(x_diff_3);

        P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
    }

    //print result
    std::cout << "- Predicted x_:" << std::endl << x_ << std::endl;
    std::cout << "- Predicted P_:" << std::endl << P_ << std::endl;
}

void UKF::Prediction(double delta_t)
{
    MatrixXd x_augmented_sigma = MatrixXd(N_AUG, N_AUG_2_PLUS_1);
    x_augmented_sigma.fill(0.0);

    GenerateAugmentedSigmaPoints(x_augmented_sigma);
    PredictSigmaPoints(x_augmented_sigma, delta_t);
    PredictMeanAndCovariance();
}

void UKF::PredictRadarMeasurement(MatrixXd& Zsig, VectorXd& z_pred, MatrixXd& S_inn)
{
    std::cout << "- PredictRadarMeasurement" << std::endl;

    //transform sigma points into measurement space
    Zsig.fill(0.0);
    for (int i = 0; i < N_AUG_2_PLUS_1; ++i)
    {
        // extract values for better readibility
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v   = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);

        double v1  = cos(yaw) * v;
        double v2  = sin(yaw) * v;

        // measurement model
        Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                                //r
        Zsig(1,i) = atan2(p_y,p_x);                                         //phi
        Zsig(2,i) = (p_x * v1 + p_y * v2 ) / sqrt(p_x * p_x + p_y * p_y);   //r_dot
    }

    //mean predicted measurement
    z_pred.fill(0.0);
    for (int i=0; i < N_AUG_2_PLUS_1; ++i)
    {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    //innovation covariance matrix S
    S_inn.fill(0.0);
    for (int i = 0; i < N_AUG_2_PLUS_1; ++i)
    {
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        //angle normalization
        double z_diff_1 = z_diff(1);
        z_diff(1) = Tools::UnwrapAngle(z_diff_1);

        S_inn = S_inn + weights_(i) * z_diff * z_diff.transpose();
    }

    S_inn = S_inn + R_RADAR;

    //print result
    std::cout << "- Predicted Radar measurement z_pred: " << std::endl << z_pred << std::endl;
    std::cout << "- Radar innovation covariance S_inn: "  << std::endl << S_inn  << std::endl;
}

void UKF::UpdateRadarState(const MatrixXd& Zsig,
                           const VectorXd& z_pred,
                           const MatrixXd& S_inn,
                           const MeasurementPackage& meas_package)
{
    std::cout << "- UpdateRadarState" << std::endl;

    //calculate cross correlation matrix
    MatrixXd Tc = MatrixXd(N_X, N_Z_RADAR);
    Tc.fill(0.0);

    for (int i = 0; i < N_AUG_2_PLUS_1; ++i)
    {
        //residual

        VectorXd z_diff = Zsig.col(i) - z_pred;
        z_diff(1) = Tools::UnwrapAngle(z_diff(1));

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        x_diff(3) = Tools::UnwrapAngle(x_diff(3));

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain K;
    MatrixXd K = Tc * S_inn.inverse();

    //Convert radar from polar to cartesian coordinates and initialize state.
    VectorXd z_radar = VectorXd(N_Z_RADAR);
    z_radar(0) = meas_package.raw_measurements_(0);
    z_radar(1) = Tools::UnwrapAngle(meas_package.raw_measurements_(1));
    z_radar(2) = meas_package.raw_measurements_(2);

    std::cout << "- Radar prediction: "  << std::endl << z_pred  << std::endl;
    std::cout << "- Radar measurement: " << std::endl << z_radar << std::endl;

    //residual
    VectorXd z_diff = z_radar - z_pred;
    z_diff(1) = Tools::UnwrapAngle(z_diff(1));

    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S_inn * K.transpose();

    //print result
    std::cout << "- Updated (Radar) state x: "      << std::endl << x_ << std::endl;
    std::cout << "- Updated (Radar) covariance P: " << std::endl << P_ << std::endl;
}

void UKF::UpdateRadar(const MeasurementPackage& meas_package)
{
    /**
    TODO:
    Complete this function! Use radar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.
    You'll also need to calculate the radar NIS.
    */

    std::cout << "Update with Radar data" << std::endl;

    MatrixXd Zsig   = MatrixXd(N_Z_RADAR, N_AUG_2_PLUS_1);
    VectorXd z_pred = VectorXd(N_Z_RADAR);
    MatrixXd S_inn  = MatrixXd(N_Z_RADAR, N_Z_RADAR);

    PredictRadarMeasurement(Zsig, z_pred, S_inn);
    UpdateRadarState(Zsig, z_pred, S_inn, meas_package);
}



