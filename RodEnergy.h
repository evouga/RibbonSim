#ifndef RODENERGY_H
#define RODENERGY_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include "RodConfig.h"

struct SimParams
{
    double constraintWeight;
    bool allowSliding;
    Eigen::MatrixXd *anchorPoints;
    Eigen::MatrixXd *anchorNormals;
    bool gravityEnabled;
    Eigen::Vector3d gravityDir;
    double floorHeight;
    double floorWeight;
};

void rAndJ(RodConfig &config,
    Eigen::VectorXd &r,
    Eigen::SparseMatrix<double> *Jr,
    double &linearEnergy,
    Eigen::VectorXd &Jlinear,
    const SimParams &params);

Eigen::Vector3d parallelTransport(const Eigen::Vector3d &v, const Eigen::Vector3d &e1, const Eigen::Vector3d &e2);
Eigen::Vector3d perpToVector(const Eigen::Vector3d &v);
#endif
