#include "RodEnergy.h"
#include <Eigen/Geometry>
#include <cassert>
#include <iostream>

double angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2)
{
    return 2.0 * atan2((v1.cross(v2)).norm(), v1.norm()*v2.norm() + v1.dot(v2));
}

void dAngle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2, Eigen::Vector3d &dv1, Eigen::Vector3d &dv2)
{
    Eigen::Vector3d n = v1.cross(v2);
    if (n.norm() == 0.0)
    {
        dv1.setZero();
        dv2.setZero();
        return;
    }
    n /= n.norm();
    dv1 = v1.cross(n) / v1.dot(v1);
    dv2 = -v2.cross(n) / v2.dot(v2);
}

static double rodBendingEnergy(const Eigen::MatrixXd &centerline, const Eigen::VectorXd &restlens, const RodParams &params, Eigen::MatrixXd &dE)
{
    double energy = 0;
    int nverts = centerline.rows();
    assert(restlens.size() == nverts - 1);
    for (int i = 0; i < nverts - 2; i++)
    {
        Eigen::Vector3d v1 = centerline.row(i + 1) - centerline.row(i);
        Eigen::Vector3d v2 = centerline.row(i + 2) - centerline.row(i + 1);
        double theta = angle(v1, v2);
        double k = theta*2.0 / (restlens[i] + restlens[i + 1]);
        energy += params.width * params.thickness * params.thickness * params.thickness * params.kbending * k*k*0.5*(restlens[i] + restlens[i + 1]);
        Eigen::Vector3d d1, d2;
        dAngle(v1, v2, d1, d2);
        dE.row(i) -= params.width * params.thickness * params.thickness * params.thickness * params.kbending * 2.0 * k * d1.transpose();
        dE.row(i + 1) += params.width * params.thickness * params.thickness * params.thickness * params.kbending * 2.0 * k * (d1.transpose() - d2.transpose());
        dE.row(i + 2) += params.width * params.thickness * params.thickness * params.thickness * params.kbending * 2.0 * k * d2.transpose();
    }
    return energy;
}

static double rodStretchingEnergy(const Eigen::MatrixXd &centerline, const Eigen::VectorXd &restlens, const RodParams &params, Eigen::MatrixXd &dE)
{
    double energy = 0;
    int nverts = centerline.rows();
    assert(restlens.size() == nverts - 1);
    for (int i = 0; i < nverts - 1; i++)
    {
        Eigen::Vector3d edgev = centerline.row(i + 1) - centerline.row(i);
        energy += params.width * params.thickness * params.kstretching * (edgev.norm() - restlens[i])*(edgev.norm() - restlens[i]);
        dE.row(i) -= params.width * params.thickness * params.kstretching * 2.0 * (edgev.norm() - restlens[i]) * edgev.transpose() / edgev.norm();
        dE.row(i+1) += params.width * params.thickness * params.kstretching * 2.0 * (edgev.norm() - restlens[i]) * edgev.transpose() / edgev.norm();
    }
    return energy;
}


static double rodTwistingEnergy(const Eigen::MatrixXd &centerline, const Eigen::VectorXd &restlens, const RodParams &params, Eigen::MatrixXd &dE)
{
    double energy = 0;
    int nverts = centerline.rows();
    assert(restlens.size() == nverts - 1);
    for (int i = 0; i < nverts - 3; i++)
    {
        Eigen::Vector3d v1 = centerline.row(i + 1) - centerline.row(i);
        Eigen::Vector3d v2 = centerline.row(i + 2) - centerline.row(i + 1);
        Eigen::Vector3d v3 = centerline.row(i + 3) - centerline.row(i + 2);
        Eigen::Vector3d B1 = v1.cross(v2);
        Eigen::Vector3d B2 = v2.cross(v3);
        double theta = angle(B1, B2);
        double tau = theta/(restlens[i+1]);
        energy += params.width * params.thickness * params.thickness * params.thickness * params.ktwist * tau*tau*restlens[i + 1];
        Eigen::Vector3d d1, d2;
        dAngle(B1, B2, d1, d2);

        Eigen::Vector3d dv1 = v2.cross(d1);
        Eigen::Vector3d dv2 = -v1.cross(d1) + v3.cross(d2);
        Eigen::Vector3d dv3 = -v2.cross(d2);        
        dE.row(i) -= params.width * params.thickness * params.thickness * params.thickness * params.ktwist * 2.0 * tau * dv1.transpose();
        dE.row(i+1) += params.width * params.thickness * params.thickness * params.thickness * params.ktwist * 2.0 * tau * dv1.transpose() - params.width * params.thickness * params.thickness * params.thickness * params.ktwist * 2.0 * tau * dv2.transpose();
        dE.row(i+2) += params.width * params.thickness * params.thickness * params.thickness * params.ktwist * 2.0 * tau * dv2.transpose() - params.width * params.thickness * params.thickness * params.thickness * params.ktwist * 2.0 * tau * dv3.transpose();
        dE.row(i+3) += params.width * params.thickness * params.thickness * params.thickness * params.ktwist * 2.0 * tau * dv3.transpose();
    }
    return energy;
}

double rodEnergy(const Eigen::MatrixXd &centerline, const Eigen::VectorXd &restlens, const RodParams &params, Eigen::MatrixXd &dE)
{
    dE.resize(centerline.rows(), 3);
    dE.setZero();
    double totenergy=0;
    totenergy += rodStretchingEnergy(centerline, restlens, params, dE);
    totenergy += rodBendingEnergy(centerline, restlens, params, dE);
    totenergy += rodTwistingEnergy(centerline, restlens, params, dE);
    return totenergy;
}

void masses(const Eigen::VectorXd &restlens, Eigen::VectorXd &M, const RodParams &params)
{
    int nedges = restlens.size();
    M.resize(nedges + 1);
    M.setZero();
    for (int i = 0; i < nedges; i++)
    {
        M[i] += restlens[i] / 2.0 * params.width * params.thickness * params.rho;
        M[i+1] += restlens[i] / 2.0 * params.width * params.thickness * params.rho;
    }
    //ends are fixed
    M[0] = std::numeric_limits<double>::infinity();
    M[nedges] = std::numeric_limits<double>::infinity();
}

