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

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state, bool* ok)
{
    MatrixXd Hj(3,4);
    Hj << 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0;

    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    float denom = px * px + py * py;
    bool status = (denom > 0.0);

    if(status)
    {
        float denom2 = denom;
        denom = sqrt(denom);

        float denom3 = denom2 * denom;
        float nom = vx * py - vy * px;

        Hj(0, 0) =  px / denom;
        Hj(0, 1) =  py / denom;
        Hj(1, 0) = -py / denom2;
        Hj(1, 1) =  px / denom2;
        Hj(2, 0) =  py * nom / denom3;
        Hj(2, 1) =  px * (0.0 - nom) / denom3;
        Hj(2, 2) =  px / denom;
        Hj(2, 3) =  py / denom;
    }
    else
    {
        std::cout << "ERROR! Denominator value: " << denom << std::endl;
    }

    if(ok)
    {
        //Return status if required
        *ok = status;
    }

    return Hj;
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

