#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class MeasurementPackage;

class UKF
{
    public:

        /**
         * Constructor
         */
        UKF();

        /**
         * Destructor
         */
        virtual ~UKF();

        /**
         * ProcessMeasurement
         * @param meas_package The latest measurement data of either radar or laser
         */
         void ProcessMeasurement(const MeasurementPackage& meas_package);

         /**
          * Get States
          */
         const VectorXd& GetStates() const;

    private:

        /**
         * Prediction Predicts sigma points, the state, and the state covariance
         * matrix
         * @param delta_t Time between k and k+1 in s
         */
        void Prediction(double delta_t);

        /**
         * Updates the state and the state covariance matrix using a laser measurement
         * @param meas_package The measurement at k+1
         */
         void UpdateLidar(const MeasurementPackage& meas_package);

        /**
         * Updates the state and the state covariance matrix using a radar measurement
         * @param meas_package The measurement at k+1
         */
        void UpdateRadar(const MeasurementPackage& meas_package);

        ///* initially set to false, set to true in first call of ProcessMeasurement
        bool is_initialized_;

        ///* time when the state is true, in us
        long long time_us_;

        ///* if this is false, laser measurements will be ignored (except for init)
        bool use_laser_;

        ///* if this is false, radar measurements will be ignored (except for init)
        bool use_radar_;

        ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
        VectorXd x_;

        ///* state covariance matrix
        MatrixXd P_;

        ///* Weights of sigma points
        VectorXd weights_;

        ///* predicted sigma points matrix
        MatrixXd Xsig_pred_;

        ///* State dimension
        static const int N_X;

        ///* Augmented state dimension
        static const int N_AUG;

        ///* Sigma point spreading parameter
        static const double LAMBDA;

        ///* Process noise standard deviation longitudinal acceleration in m/s^2
        static const double STD_A;

        ///* Process noise standard deviation yaw acceleration in rad/s^2
        static const double STD_YAW_DD;

        ///* Laser measurement noise standard deviation position1 in m
        static const double STD_LASER_PX;

        ///* Laser measurement noise standard deviation position2 in m
        static const double STD_LASER_PY;

        ///* Radar measurement noise standard deviation radius in m
        static const double STD_RAD_R;

        ///* Radar measurement noise standard deviation angle in rad
        static const double STD_RAD_PHI;

        ///* Radar measurement noise standard deviation radius change in m/s
        static const double STD_RAD_R_D;
};

#endif /* UKF_H */
