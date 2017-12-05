#include "PhysicsHook.h"
#include "RodEnergy.h"

class RodsHook : public PhysicsHook
{
public:
    RodsHook() : PhysicsHook(), dirty(true) {
    }

    virtual void initGUI(igl::viewer::Viewer &viewer);

    virtual void initSimulation()
    {
        //make helix
        centerline.resize(100, 3);
        vel.resize(100, 3);
        vel.setZero();
        for (int i = 0; i < 100; i++)
        {
            double r = 0.2 + 0.1 * sin(double(i) / 10);
            centerline(i, 0) = r*cos(double(i) / 4.0);
            centerline(i, 1) = r*sin(double(i) / 4.0);
            centerline(i, 2) = double(i) / 100;            
        }
        int nverts = centerline.rows();
        restlens.resize(nverts - 1);
        for (int i = 0; i < nverts - 1; i++)
        {
            restlens[i] = (centerline.row(i + 1) - centerline.row(i)).norm();
        }

        masses(restlens, M, params);
        createVisualizationMesh();
    }

    virtual void updateRenderGeometry()
    {
        if (Q.rows() != renderQ.rows() || F.rows() != renderF.rows())
            dirty = true;
        renderQ = Q;
        renderF = F;
    }

    virtual bool simulateOneStep()
    {
        int nverts = centerline.rows();
        centerline += dt*vel;
        Eigen::MatrixXd dE;
        double energy = rodEnergy(centerline, restlens, params, dE);
        for (int i = 0; i < nverts; i++)
            dE.row(i) /= M[i];
        vel -= dt*dE;
        std::cout << "Energy: " << energy << std::endl;
        createVisualizationMesh();
        return false;
    }

    virtual void renderRenderGeometry(igl::viewer::Viewer &viewer)
    {
        if (dirty)
        {
            viewer.data.clear();
            dirty = false;
        }
        viewer.data.set_mesh(renderQ, renderF);
    }

private:    
    void createVisualizationMesh();

    double dt;
    RodParams params;
    Eigen::MatrixXd centerline;
    Eigen::MatrixXd vel;
    Eigen::VectorXd restlens;
    Eigen::VectorXd M;

    // for visualization
    Eigen::MatrixXd Q;
    Eigen::MatrixXi F;

    int segsPerEdge;
    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
    bool dirty;
};