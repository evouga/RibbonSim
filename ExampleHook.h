#include "PhysicsHook.h"

class ExampleHook : public PhysicsHook
{
public:
    ExampleHook() : PhysicsHook() {}

    virtual void initGUI(igl::viewer::Viewer &viewer)
    {

    }

    virtual void initSimulation()
    {
        origQ.resize(4, 3);
        origQ << -1, -1, 0,
            1, -1, 0,
            -1, 1, 0,
            1, 1, 0;
        Q = origQ*1.1;
        V.resize(4, 3);
        V.setZero();
        F.resize(2, 3);
        F << 0, 1, 2,
            2, 1, 3;

        dt = 1e-3;
        k = 1e-2;
    }

    virtual void updateRenderGeometry()
    {
        renderQ = Q;
        renderF = F;        
    }

    virtual bool simulateOneStep()
    {
        Q += dt * V;
        Eigen::MatrixXd Force = k*(origQ - Q);
        V += dt * Force;
        return false;
    }

    virtual void renderRenderGeometry(igl::viewer::Viewer &viewer)
    {
        viewer.data.set_mesh(renderQ, renderF);
    }

private:
    double k;
    double dt;
    Eigen::MatrixXd origQ;
    Eigen::MatrixXd Q;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
};