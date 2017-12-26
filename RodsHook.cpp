#include "RodsHook.h"

void RodsHook::initGUI(igl::viewer::Viewer &viewer)
{
    double Y = 1e9;
    params.kbending = Y;
    params.kstretching = Y;
    params.ktwist = Y / 3.0;
    params.rho = 1.0;
    params.thickness = 1e-4;
    params.width = .01;
    dt = 1e-7;
    damp = 100;

    viewer.ngui->addGroup("Rod Parameters");
    viewer.ngui->addVariable("Thickness", params.thickness);
    viewer.ngui->addVariable("Width", params.width);
    viewer.ngui->addVariable("Stretching k", params.kstretching);
    viewer.ngui->addVariable("Bending k", params.kbending);
    viewer.ngui->addVariable("Twisting k", params.ktwist);

    viewer.ngui->addGroup("Sim Options");
    viewer.ngui->addVariable("Time Step", dt);
    viewer.ngui->addVariable("Damping Factor", damp);
    
}

void RodsHook::showForces(const Eigen::MatrixXd &dE)
{
    int nverts = rod->curState.centerline.rows();
    forcePoints.resize(2*nverts, 3);
    forceEdges.resize(nverts, 2);
    forceColors.resize(nverts, 3);
    for (int i = 0; i < nverts; i++)
    {
        forcePoints.row(2 * i) = rod->curState.centerline.row(i);
        forcePoints.row(2 * i+1) = rod->curState.centerline.row(i) - 1e-1*dE.row(i);
        forceEdges(i, 0) = 2 * i;
        forceEdges(i, 1) = 2 * i + 1;
        forceColors.row(i) = Eigen::Vector3d(1, 0, 0);
    }
}

void RodsHook::createVisualizationMesh()
{
    int nverts = rod->curState.centerline.rows();
    int nsegs = rod->isClosed() ? nverts : nverts - 1;
    Eigen::MatrixXd N(nsegs, 3);
    Eigen::MatrixXd B(nsegs, 3);
    for (int i = 0; i < nsegs; i++)
    {
        Eigen::Vector3d v0 = rod->curState.centerline.row(i);
        Eigen::Vector3d v1 = rod->curState.centerline.row((i + 1) % nverts);
        Eigen::Vector3d e = v1 - v0;
        e /= e.norm();
        Eigen::Vector3d d1 = rod->curState.directors.row(i);
        Eigen::Vector3d d2 = e.cross(d1);
        double theta = rod->curState.thetas[i];
        N.row(i) = d1*cos(theta) + d2*sin(theta);
        B.row(i) = -d1*sin(theta) + d2*cos(theta);
    }
    
    Q.resize(8 * nsegs, 3);
    F.resize(8 * nsegs, 3);
    for (int i = 0; i < nsegs; i++)
    {
        Eigen::Vector3d v0 = rod->curState.centerline.row(i);
        Eigen::Vector3d v1 = rod->curState.centerline.row((i + 1) % nverts);
        Eigen::Vector3d T = v1 - v0;
        T /= T.norm();
        Q.row(8 * i + 0) = (v0.transpose() + params.thickness / 2.0 * N.row(i) - params.width / 2.0 * B.row(i));
        Q.row(8 * i + 1) = (v0.transpose() + params.thickness / 2.0 * N.row(i) + params.width / 2.0 * B.row(i));
        Q.row(8 * i + 2) = (v0.transpose() - params.thickness / 2.0 * N.row(i) + params.width / 2.0 * B.row(i));
        Q.row(8 * i + 3) = (v0.transpose() - params.thickness / 2.0 * N.row(i) - params.width / 2.0 * B.row(i));
        Q.row(8 * i + 4) = (v1.transpose() + params.thickness / 2.0 * N.row(i) - params.width / 2.0 * B.row(i));
        Q.row(8 * i + 5) = (v1.transpose() + params.thickness / 2.0 * N.row(i) + params.width / 2.0 * B.row(i));
        Q.row(8 * i + 6) = (v1.transpose() - params.thickness / 2.0 * N.row(i) + params.width / 2.0 * B.row(i));
        Q.row(8 * i + 7) = (v1.transpose() - params.thickness / 2.0 * N.row(i) - params.width / 2.0 * B.row(i));
        for (int j = 0; j < 4; j++)
        {
            F(8 * i + 2 * j, 0) = 8 * i + j;
            F(8 * i + 2 * j, 2) = 8 * i + 4 + j;
            F(8 * i + 2 * j, 1) = 8 * i + 4 + ((j + 1) % 4);
            F(8 * i + 2 * j + 1, 0) = 8 * i + (j + 1) % 4;
            F(8 * i + 2 * j + 1, 2) = 8 * i + j;
            F(8 * i + 2 * j + 1, 1) = 8 * i + 4 + ((j + 1) % 4);
        }
    }
}

int iter = 0;

bool RodsHook::simulateOneStep()
{
    int nverts = rod->curState.centerline.rows();
    int nsegs = rod->isClosed() ? nverts : nverts - 1;
    for (int i = 0; i < nsegs; i++)
    {
        Eigen::Vector3d oldv1 = rod->curState.centerline.row(i);
        Eigen::Vector3d oldv2 = rod->curState.centerline.row((i + 1) % nverts);
        Eigen::Vector3d v1 = oldv1 + dt*rod->curState.centerlineVel.row(i).transpose();
        Eigen::Vector3d v2 = oldv2 + dt*rod->curState.centerlineVel.row((i + 1) % nverts).transpose();

        rod->curState.directors.row(i) = parallelTransport(rod->curState.directors.row(i), oldv2 - oldv1, v2 - v1);
    }
    rod->curState.centerline += dt*rod->curState.centerlineVel;
    rod->curState.thetas += dt*rod->curState.directorAngVel;

    Eigen::MatrixXd dE;
    Eigen::VectorXd dtheta;
    double energy = rodEnergy(*rod, rod->curState, &dE, &dtheta);
    
    createVisualizationMesh();
    showForces(dE);
    
    for (int i = 0; i < nverts; i++)
    {
        dE.row(i) /= rod->masses[i];
    }
    for (int i = 0; i < nsegs; i++)
    {
        dtheta[i] /= rod->momInertia[i];
    }
        
    rod->curState.centerlineVel -= dt*dE;
    rod->curState.directorAngVel -= dt*dtheta;

    double dampfactor = exp(-dt*damp);
    rod->curState.centerlineVel *= dampfactor;
    rod->curState.directorAngVel *= dampfactor;
    
    std::cout << "Energy: " << energy << std::endl;
    iter++;

    /*if (iter == 100)
    {
        rod->params.kstretching = 0;
        double energy = rodEnergy(*rod, rod->curState, NULL, NULL);
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
        while (true);
    }*/

    return false;
}

