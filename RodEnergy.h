#ifndef RODENERGY_H
#define RODENERGY_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include "RodConfig.h"

void rAndJ(RodConfig &config, Eigen::VectorXd &r, Eigen::SparseMatrix<double> *Jr);

double constraintEnergy(RodConfig &config, std::vector<Eigen::VectorXd> *dEs, std::vector<Eigen::VectorXd> *dthetas);

Eigen::Vector3d parallelTransport(const Eigen::Vector3d &v, const Eigen::Vector3d &e1, const Eigen::Vector3d &e2);

#endif