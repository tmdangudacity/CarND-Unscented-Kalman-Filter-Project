#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools
{
    public:
        /**
        * Constructor.
        */
        Tools();

        /**
        * Destructor.
        */
        virtual ~Tools();

        /**
        * A helper method to calculate RMSE.
        */
        VectorXd CalculateRMSE(const vector<VectorXd>& estimations,
                               const vector<VectorXd>& ground_truth,
                               bool* ok = NULL);

        /**
        * A helper method to calculate Jacobians.
        */
        MatrixXd CalculateJacobian(const VectorXd& x_state,
                                   bool* ok = NULL);

        //Return the value of Pi
        static double Pi();

        //Unwrapping angle to between -Pi and Pi
        static double UnwrapAngle(double in_angle);
};

#endif /* TOOLS_H_ */
