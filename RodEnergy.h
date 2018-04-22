#ifndef RODENERGY_H
#define RODENERGY_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include "RodConfig.h"

void rAndJ(RodConfig &config, 
    Eigen::VectorXd &r, 
    Eigen::SparseMatrix<double> *Jr, 
    double angleWeight, 
    bool allowSliding, 
    Eigen::MatrixXd *anchorPoints,
    Eigen::MatrixXd *anchorNormals);

Eigen::Vector3d parallelTransport(const Eigen::Vector3d &v, const Eigen::Vector3d &e1, const Eigen::Vector3d &e2);

#endif
