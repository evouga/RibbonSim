#ifndef ROD_H
#define ROD_H

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
   
    Eigen::VectorXd widths;

    // Projection Variables
    Eigen::MatrixXd closestFaceNormals;
    Eigen::MatrixXd closestFaceCentroids; 
    double kproject; 
};

Eigen::Vector3d parallelTransport(const Eigen::Vector3d &v, const Eigen::Vector3d &e1, const Eigen::Vector3d &e2);

class Rod
{
public:
    Rod(const RodState &startState, RodParams &params, bool isClosed);

    bool isClosed() const { return isClosed_; }
    int numVertices() const { return (int)curState.centerline.rows(); }
    int numSegments() const { return (int)curState.thetas.size(); }
    RodState curState;
    RodState startState;

    RodParams params;
    igl::AABB<Eigen::MatrixXd,3> *tree;
    void updateProjectionVars(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);

// Material parameters - recompute after reset!
    Eigen::VectorXd restlens;
    Eigen::VectorXd masses;
    Eigen::VectorXd momInertia;

    Eigen::MatrixXd c_points;

    void projectToMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);

private:
    bool isClosed_;

    void initializeRestQuantities();
};

#endif
