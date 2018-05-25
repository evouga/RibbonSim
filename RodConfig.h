#ifndef RODCONFIG_H
#define RODCONFIG_H

#include <vector>
#include <Eigen/Core>

#include "utility/simple_svg_1.0.0.hpp"

const static double rod_colors[][3] = { 
    { .001, .001, .95 },
    { 0.95, 0.001, 0.001 },
    { 0.95, 0.95, 0.001 },
    { 0.01, 0.95, 0.01 },
    { 0.95, 0.7, 0.01 },
    { 0.5, 0.01, 0.5 },
    { 0.001, 0.95, 0.95 }      
};

const static double top_color[3] = { .5, 0, 0 };
const static double bottom_color[3] = { .7, .7, .2 };

const static int num_rod_colors = sizeof(rod_colors)/(3*sizeof(double));


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
    Rod(const RodState &startState, const Eigen::VectorXd &widths, RodParams &params, bool isClosed, int colorID);

    enum RodVisibilityState
    {
        RS_TRANSLUCENT=0,
        RS_VISIBLE,
        RS_HIDDEN
    };

    bool isClosed() const { return isClosed_; }

    void setVisibilityState(RodVisibilityState state) { visState_ = state; }
    RodVisibilityState visibilityState() const { return visState_; }

    Eigen::Vector3d rodColor() const;
    int rodColorID() const { return colorID_; }
    void cycleColor();
    int numVertices() const { return (int)curState.centerline.rows(); }
    int numSegments() const { return (int)curState.thetas.size(); }
    double arclength() const;
    RodState curState;
    RodState startState;
    
    Eigen::VectorXd widths;
    Eigen::VectorXd restlens;
    Eigen::VectorXd masses;
    Eigen::VectorXd momInertia;
    RodParams params;


    void initializeRestQuantities();
private:
    int colorID_;
    bool isClosed_;
    RodVisibilityState visState_;
};

struct Constraint
{
    int rod1, rod2;
    int seg1, seg2;
    double bary1, bary2;
    double stiffness;
    int assignment;
    bool visited;
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
    Eigen::Vector3d shadeRodSegment(int rod, int segment) const;
    void saveRodGeometry(const std::string &prefix);

    std::vector<Rod *> rods;
    std::vector<Constraint> constraints;

    bool showConstraints;
};

#endif
