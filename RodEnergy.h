#ifndef RODENERGY_H
#define RODENERGY_H

#include <Eigen/Core>

struct RodParams
{
    double width;
    double thickness;
    double kstretching;
    double kbending;
    double ktwist;
    double rho;
};

double angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2);
void dAngle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2, Eigen::Vector3d &dv1, Eigen::Vector3d &dv2);

double rodEnergy(const Eigen::MatrixXd &centerline, const Eigen::VectorXd &restlens, const RodParams &params, Eigen::MatrixXd &dE);

void masses(const Eigen::VectorXd &restlens, Eigen::VectorXd &M, const RodParams &params);

#endif