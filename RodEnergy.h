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

struct RodState
{
    Eigen::MatrixXd centerline;
    Eigen::MatrixXd directors;
    Eigen::VectorXd thetas;
    Eigen::MatrixXd ceterlineVel;
    Eigen::VectorXd directorAngVel;
};

class Rod
{
public:
    Rod(const RodState &startState, const RodParams &params);

    RodState curState;
    RodState startState;

    Eigen::VectorXd restlens;
    Eigen::VectorXd masses;
    Eigen::VectorXd momInertia;
    RodParams params;

private:
    void initializeRestQuantities();
};

double rodEnergy(const Rod &rod, const RodState &state, Eigen::MatrixXd *dE, Eigen::VectorXd *dtheta);

Eigen::Vector3d parallelTransport(const Eigen::Vector3d &v, const Eigen::Vector3d &e1, const Eigen::Vector3d &e2);

#endif