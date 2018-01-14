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
    
    // Projection Variables
    Eigen::MatrixXd closestFaceNormals;
    Eigen::MatrixXd closestFaceCentroids; 
    double kproject; 
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

#endif
