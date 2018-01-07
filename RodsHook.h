#include "PhysicsHook.h"
#include "RodEnergy.h"

class RodsHook : public PhysicsHook
{
public:
    RodsHook() : PhysicsHook(), iter(0), forceResidual(0.0), dirty(true), config(NULL) {
    }

    virtual void initGUI(igl::viewer::Viewer &viewer);

    virtual void initSimulation();
    
    virtual void updateRenderGeometry()
    {
        if (Q.rows() != renderQ.rows() || F.rows() != renderF.rows())
            dirty = true;
        renderQ = Q;
        renderF = F;        
    }

    virtual bool simulateOneStep();
    
    virtual void renderRenderGeometry(igl::viewer::Viewer &viewer)
    {
        if (dirty)
        {
            viewer.data.clear();
            dirty = false;
        }
        viewer.data.set_mesh(renderQ, renderF);
        if(forceEdges.rows() > 0)
            viewer.data.set_edges(forcePoints, forceEdges, forceColors);
    }

private:    
    void createVisualizationMesh();
    void showForces(int rod, const Eigen::MatrixXd &dE);

    double dt;
    double damp;

    double iter;
    double forceResidual;
    
    RodConfig *config;

    // for visualization
    Eigen::MatrixXd Q;
    Eigen::MatrixXi F;
    Eigen::MatrixXd forcePoints;
    Eigen::MatrixXi forceEdges;
    Eigen::MatrixXd forceColors;

    int segsPerEdge;
    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
    bool dirty;
};