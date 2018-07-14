#ifndef UKF_H
#define UKF_H

#include <string>

#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class MeasurementPackage;

class UKF
{
    public:

        typedef enum
        {
            LASER_DATA_ONLY = 0,
            RADAR_DATA_ONLY,
            LASER_AND_RADAR_DATA
        }
        DataOption;

        /**
        * Constructor.
        */
        UKF(DataOption in_data = LASER_AND_RADAR_DATA);

        /**
         * Destructor
         */
        virtual ~UKF();

        /**
         * ProcessMeasurement
         * @param meas_package The latest measurement data of either radar or laser
         */
         bool ProcessMeasurement(const MeasurementPackage& meas_package);

         /**
          * Get States
          */
         const VectorXd& GetStates() const;

         static std::string ToString(DataOption data_option);

    private:

        // what stream of data is used
        DataOption data_option_;

        ///* initially set to false, set to true in first call of ProcessMeasurement
        bool is_initialized_;

        ///* time when the state is true, in us
        long long time_us_;

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

        /// 2 * N_AUG + 1
        static const int N_AUG_2_PLUS_1;

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

        /**
         * Initialisation of the state with measurement data
         */
        bool Initialise(const MeasurementPackage& meas_package);

        /**
         * Prediction Predicts sigma points, the state, and the state covariance
         * matrix
         * @param delta_t Time between k and k+1 in s
         */
        bool Prediction(double delta_t);

        /**
         * Updates the state and the state covariance matrix using a laser measurement
         * @param meas_package The measurement at k+1
         */
        bool UpdateLidar(const MeasurementPackage& meas_package);

        /**
         * Updates the state and the state covariance matrix using a radar measurement
         * @param meas_package The measurement at k+1
         */
        bool UpdateRadar(const MeasurementPackage& meas_package);

        /**
         * Generate augmented sigma points
         * */
        void GenerateAugmentedSigmaPoints(MatrixXd& Xsig_aug) const;

        /**
         * Sigma point prediction
         * */
        void PredictSigmaPoints(const MatrixXd& Xsig_aug, double delta_t);

        /**
         * Predict mean and covariance
         * */
        void PredictMeanAndCovariance();
};

#endif /* UKF_H */
