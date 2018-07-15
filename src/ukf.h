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

        ///* Laser measurement noise covariance matrix
        MatrixXd R_LASER;

        //* Radar measurement noise covariance matrix
        MatrixXd R_RADAR;

        ///* State dimension
        static const int N_X;

        ///* Laser measurement dimension
        static const int N_Z_LASER;

        ///* Radar measurement dimension
        static const int N_Z_RADAR;

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

        /** Prediction step **/

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

        /**
         * Prediction Predicts sigma points, the state, and the state covariance
         * matrix
         * @param delta_t Time between k and k+1 in s
         */
        void Prediction(double delta_t);

        /** Update with Laser measurment **/

        /**
         * Updates the state and the state covariance matrix using a laser measurement
         * @param meas_package The measurement at k+1
         */
        void UpdateLidar(const MeasurementPackage& meas_package);


        /** Update with Radar measurement **/

        /**
         *  Predict Radar measurement
         **/
        void PredictRadarMeasurement(MatrixXd& Zsig, VectorXd& z_pred, MatrixXd& S_inn);

        /**
         * Update Radar state
         **/
        void UpdateRadarState(const MatrixXd& Zsig,
                              const VectorXd& z_pred,
                              const MatrixXd& S_inn,
                              const MeasurementPackage& meas_package);

        /**
         * Updates the state and the state covariance matrix using a radar measurement
         * @param meas_package The measurement at k+1
         */
        void UpdateRadar(const MeasurementPackage& meas_package);
};

#endif /* UKF_H */
