#include "RodsHook.h"
#include "RodParser.h"

void RodsHook::initGUI(igl::viewer::Viewer &viewer)
{
    dt = 1e-7;
    damp = 100;
    savePrefix = "rod_";
    
    loadName = "torus";
    refreshLoadBuffers();

    viewer.ngui->addVariable("Config File", loadName);
    viewer.ngui->addVariable("Show Target", showTarget);

    viewer.ngui->addGroup("Sim Options");
    viewer.ngui->addVariable("Time Step", dt);
    viewer.ngui->addVariable("Damping Factor", damp);
    viewer.ngui->addButton("Save Geometry", std::bind(&RodsHook::saveRods, this));
    viewer.ngui->addVariable("Save Prefix", savePrefix);
    
    viewer.ngui->addGroup("Sim Status");
    viewer.ngui->addVariable("Iteration", iter, false);
    viewer.ngui->addVariable("Force Residual", forceResidual, false);
}

void RodsHook::refreshLoadBuffers()
{
    sprintf( loadBuffer, "../configs/%s.rod", loadName.c_str() );
    sprintf( meshBuffer, "../../relax-field/meshes/%s.obj", loadName.c_str() );
}

void RodsHook::initSimulation()
{
    iter = 0;
    forceResidual = 0;

    if (config)
        delete config;

    refreshLoadBuffers();
    config = readRod(loadBuffer);
    if (!config)
        exit(-1);
    hasMesh = config->loadTargetMesh(meshBuffer);
    config->hasMesh = hasMesh;

    createVisualizationMesh();
    dirty = true;
}


void RodsHook::showForces(int rod, const Eigen::VectorXd &dE)
{
    int nverts = config->rods[rod]->numVertices();
    forcePoints.resize(2*nverts, 3);
    forceEdges.resize(nverts, 2);
    forceColors.resize(nverts, 3);
    for (int i = 0; i < nverts; i++)
    {
        forcePoints.row(2 * i) = config->rods[rod]->curState.centerline.row(i);
        forcePoints.row(2 * i+1) = config->rods[rod]->curState.centerline.row(i) - dE.segment<3>(3*i).transpose();
        forceEdges(i, 0) = 2 * i;
        forceEdges(i, 1) = 2 * i + 1;
        forceColors.row(i) = Eigen::Vector3d(1, 0, 0);
    }
}

void RodsHook::createVisualizationMesh()
{
    config->createVisualizationMesh(Q, F);
}

bool RodsHook::simulateOneStep()
{
    int nrods = config->numRods();
    // update positions
    for (int rod = 0; rod < nrods; rod++)
    {
        int nverts = config->rods[rod]->numVertices();
        int nsegs = config->rods[rod]->numSegments();
        for (int i = 0; i < nsegs; i++)
        {
            Eigen::Vector3d oldv1 = config->rods[rod]->curState.centerline.row(i);
            Eigen::Vector3d oldv2 = config->rods[rod]->curState.centerline.row((i + 1) % nverts);
            Eigen::Vector3d v1 = oldv1 + dt*config->rods[rod]->curState.centerlineVel.row(i).transpose();
            Eigen::Vector3d v2 = oldv2 + dt*config->rods[rod]->curState.centerlineVel.row((i + 1) % nverts).transpose();

            config->rods[rod]->curState.directors.row(i) = parallelTransport(config->rods[rod]->curState.directors.row(i), oldv2 - oldv1, v2 - v1);
        }
        config->rods[rod]->curState.centerline += dt*config->rods[rod]->curState.centerlineVel;
        config->rods[rod]->curState.thetas += dt*config->rods[rod]->curState.directorAngVel;
    }

    createVisualizationMesh();    

    double newresid = 0;

    std::vector<Eigen::VectorXd> constraintdE;
    std::vector<Eigen::VectorXd> constraintdtheta;
    constraintEnergy(*config, &constraintdE, &constraintdtheta);

    // update velocities
    for (int rod = 0; rod < nrods; rod++)
    {
        Eigen::VectorXd dE;
        Eigen::VectorXd dtheta;
        double energy = rodEnergy(*config->rods[rod], config->rods[rod]->curState, &dE, &dtheta, NULL, NULL);
        if (rod == 0)
            showForces(0, dE);
        dE += constraintdE[rod];
        dtheta += constraintdtheta[rod];

        int nverts = config->rods[rod]->numVertices();
        int nsegs = config->rods[rod]->numSegments();
        for (int i = 0; i < nverts; i++)
        {
            newresid += dE.segment<3>(3*i).squaredNorm() / config->rods[rod]->masses[i];
            config->rods[rod]->curState.centerlineVel.row(i) -= dt * dE.segment<3>(3*i).transpose() / config->rods[rod]->masses[i];            
        }
        for (int i = 0; i < nsegs; i++)
        {
            newresid += dtheta[i]*dtheta[i] / config->rods[rod]->momInertia[i];
            config->rods[rod]->curState.directorAngVel[i] -= dt * dtheta[i] / config->rods[rod]->momInertia[i];
        }

        double dampfactor = exp(-dt*damp);
        config->rods[rod]->curState.centerlineVel *= dampfactor;
        config->rods[rod]->curState.directorAngVel *= dampfactor;        
    }
    forceResidual = newresid;
    iter++;
    /*
    if (iter == 10000)
    {
        config->rods[0]->params.kstretching = 0;
        double energy = rodEnergy(*config->rods[0], config->rods[0]->curState, NULL, NULL);
        for (int i = 0; i < config->rods[0]->numVertices(); i++)
        {
            for (int j = 0; j < 3; j++)
            {
                RodState cp = config->rods[0]->curState;
                cp.centerline(i, j) += 1e-6;
                Eigen::MatrixXd dE;
                double newenergy = rodEnergy(*config->rods[0], cp, &dE, NULL);
                double findiff = (newenergy - energy) / 1e-6;
                std::cout << findiff << " " << dE(i, j) << std::endl;
            }
        }
        for (int i = 0; i < config->rods[0]->numSegments(); i++)
        {
            RodState cp = config->rods[0]->curState;
            cp.thetas[i] += 1e-6;
            Eigen::VectorXd dtheta;
            double newenergy = rodEnergy(*config->rods[0], cp, NULL, &dtheta);
            double findiff = (newenergy - energy) / 1e-6;
            std::cout << findiff << " " << dtheta[i] << std::endl;
        }
        while (true);
    }*/

    return false;
}

void RodsHook::saveRods()
{
    config->saveRodGeometry(savePrefix);
}
