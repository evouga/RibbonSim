#ifndef RODCONFIG_H
#define RODCONFIG_H

#include <vector>
#include <Eigen/Core>

#include "utility/simple_svg_1.0.0.hpp"

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

    Eigen::MatrixXd colors;
    bool visible; 

    Eigen::VectorXd widths;
    Eigen::VectorXd restlens;
    Eigen::VectorXd masses;
    Eigen::VectorXd momInertia;
    RodParams params;

    void initializeRestQuantities();
private:
    bool isClosed_;

};

struct Constraint
{
    int rod1, rod2;
    int seg1, seg2;
    double bary1, bary2;
    double stiffness;
    int assignment;
    bool visited;
    int color;
};

class RodConfig
{
public:
    ~RodConfig();

    void addRod(Rod *rod);
    void addConstraint(Constraint c);
    void initWeave(); // Call after all constraints initialized
    int numRods() const { return (int)rods.size(); }
    int numConstraints() const { return (int)constraints.size(); }
    void reset();
    void createVisualizationMesh(Eigen::MatrixXd &Q, Eigen::MatrixXi &F);
    void setVisualizationMeshColors();
    void saveRodGeometry(const std::string &prefix);

    std::vector<Rod *> rods;
    std::vector<Constraint> constraints;

    bool showConstraints;
};

#endif
