#ifndef RODCONFIG_H
#define RODCONFIG_H

#include <vector>
#include <Eigen/Core>
#include <igl/AABB.h>
#include "Rod.h"


struct Constraint
{
    int rod1, rod2;
    int seg1, seg2;
    double bary1, bary2;
    double stiffness;
};

class RodConfig
{
public:
    ~RodConfig();

    void addRod(Rod *rod);
    void addConstraint(Constraint c);
    int numRods() const { return (int)rods.size(); }
    void reset();
    bool loadTargetMesh(const std::string &objname);
    void createVisualizationMesh(Eigen::MatrixXd &Q, Eigen::MatrixXi &F);
    void saveRodGeometry(const std::string &prefix);

    std::vector<Rod *> rods;
    std::vector<Constraint> constraints;

    // Target mesh - for optimization 
    Eigen::MatrixXd V_mesh; // |V| x 3
    Eigen::MatrixXi F_mesh; // |F| x 3

    igl::AABB<Eigen::MatrixXd,3> mesh_tree;
    bool hasMesh;
};

#endif
