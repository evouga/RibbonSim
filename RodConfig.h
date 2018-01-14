#ifndef RODCONFIG_H
#define RODCONFIG_H

#include <vector>
#include <Eigen/Core>
#include <igl/AABB.h>

struct RodParams
{
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
    Eigen::MatrixXd centerlineVel;
    Eigen::VectorXd directorAngVel;
    
    // Projection Variables
    Eigen::MatrixXd closestFaceNormals;
    Eigen::MatrixXd closestFaceCentroids; 
};

class Rod
{
public:
    Rod(const RodState &startState, const Eigen::VectorXd &widths, RodParams &params, bool isClosed);

    bool isClosed() const { return isClosed_; }
    int numVertices() const { return (int)curState.centerline.rows(); }
    int numSegments() const { return (int)curState.thetas.size(); }
    RodState curState;
    RodState startState;

    Eigen::VectorXd widths;
    Eigen::VectorXd restlens;
    Eigen::VectorXd masses;
    Eigen::VectorXd momInertia;
    RodParams params;

    igl::AABB<Eigen::MatrixXd,3> *tree;
    void updateProjectionVars(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
     
private:
    bool isClosed_;

    void initializeRestQuantities();
};

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
