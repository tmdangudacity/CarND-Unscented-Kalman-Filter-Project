#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth, bool* ok)
{
    bool status = estimations.size() && (ground_truth.size() == estimations.size());

    VectorXd rmse(4);
    rmse << 0.0, 0.0, 0.0, 0.0;

    if(status)
    {
        VectorXd temp_error;

        for(unsigned int i=0; i < estimations.size(); ++i)
        {
            temp_error = estimations[i] - ground_truth[i];
            rmse += (temp_error.cwiseProduct(temp_error));
        }

        rmse = rmse / (estimations.size());

        rmse = rmse.array().sqrt();
    }
    else
    {
        //Error
        std::cout << "ERROR! Invalid inputs"
                  << ", Estimation size: "   << estimations.size()
                  << ", Ground truth size: " << ground_truth.size()
                  << std::endl;
    }

    if(ok)
    {
        //Return status if required
        *ok = status;
    }

    return rmse;
}

double Tools::Pi()
{
    static double PI_VALUE = acos(-1.0);

    return PI_VALUE;
}

double Tools::UnwrapAngle(double in_angle)
{
    double ret = in_angle;
    double pi_value = Pi();
    double two_pi_value = 2.0 * pi_value;

    while(ret > pi_value)
    {
        ret -= two_pi_value;
    }

    while(ret < -pi_value)
    {
        ret += two_pi_value;
    }

    return ret;
}


