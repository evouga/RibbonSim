#include "PhysicsHook.h"
#include "RodEnergy.h"

class RodsHook : public PhysicsHook
{
public:
    RodsHook() : PhysicsHook(), dirty(true), rod(NULL) {
    }

    virtual void initGUI(igl::viewer::Viewer &viewer);

    virtual void initSimulation()
    {
        RodState rs;
        int nverts = 100;
        double PI = 3.1415926525898;
        rs.centerline.resize(nverts, 3);
        rs.centerlineVel.resize(nverts, 3);
        rs.centerlineVel.setZero();

        for (int i = 0; i < nverts; i++)
        {
            double r = 0.2 +0.1 * sin(2.0*PI*double(i) / 50);
            rs.centerline(i, 0) = r*cos(2.0*PI*double(i) / 100.0);
            rs.centerline(i, 1) = r*sin(2.0*PI*double(i) / 100.0);
            rs.centerline(i, 2) = 0; //double(i) / nverts;            
        }
        
        
        rs.directors.resize(nverts, 3);
        rs.directors.setZero();
        rs.directorAngVel.resize(nverts);
        rs.directorAngVel.setZero();

        rs.thetas.resize(nverts);
        rs.thetas.setZero();
        int nsegs = nverts;
        for (int i = 0; i < nverts; i++)
        {
            Eigen::Vector3d v1 = rs.centerline.row(i);
            Eigen::Vector3d v2 = rs.centerline.row((i + 1) % nverts);
            Eigen::Vector3d d = Eigen::Vector3d(0, 0, 1).cross(v2 - v1);
            rs.directors.row(i) = d;
        }
        for (int i = 0; i < nsegs; i++)
        {
            rs.directors.row(i) /= rs.directors.row(i).norm();
        }

        
        if (rod)
            delete rod;
        rod = new Rod(rs, params, true);
        createVisualizationMesh();
        Eigen::MatrixXd dE;
        rodEnergy(*rod, rod->curState, &dE, NULL);
        showForces(dE);

        /*double energy = rodEnergy(*rod, rod->curState, NULL, NULL);
        for (int i = 0; i < nverts; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                RodState cp = rod->curState;
                cp.centerline(i, j) += 1e-6;
                Eigen::MatrixXd dE;
                double newenergy = rodEnergy(*rod, cp, &dE, NULL);
                double findiff = (newenergy - energy) / 1e-6;
                std::cout << findiff << " " << dE(i, j) << std::endl;
            }
        }
        for (int i = 0; i < nverts - 1; i++)
        {
            RodState cp = rod->curState;
            cp.thetas[i] += 1e-6;
            Eigen::VectorXd dtheta;
            double newenergy = rodEnergy(*rod, cp, NULL, &dtheta);
            double findiff = (newenergy - energy) / 1e-6;
            std::cout << findiff << " " << dtheta[i] << std::endl;
        }
        
        while (true);*/
    }

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
        viewer.data.set_edges(forcePoints, forceEdges, forceColors);
    }

private:    
    void createVisualizationMesh();
    void showForces(const Eigen::MatrixXd &dE);

    double dt;
    double damp;
    RodParams params;
    Rod *rod;

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