#include "RodsHook.h"
#include "RodParser.h"

void RodsHook::initGUI(igl::viewer::Viewer &viewer)
{
    dt = 1e-7;
    damp = 100;

    viewer.ngui->addGroup("Sim Options");
    viewer.ngui->addVariable("Time Step", dt);
    viewer.ngui->addVariable("Damping Factor", damp);
    
    viewer.ngui->addGroup("Sim Status");
    viewer.ngui->addVariable("Iteration", iter, false);
    viewer.ngui->addVariable("Force Residual", forceResidual, false);
}

void RodsHook::initSimulation()
{
    iter = 0;
    forceResidual = 0;

    if (config)
        delete config;

    config = readRod("../configs/example.rod");

    /*RodParams params;
    double Y = 1e8;
    params.kbending = Y;
    params.kstretching = Y;
    params.ktwist = Y / 3.0;
    params.rho = 1.0;
    params.thickness = 1e-4;    

   
    config = new RodConfig;
    double PI = 3.1415926525898;

    int nrods = 4;
    for (int rod = 0; rod < nrods; rod++)
    {
        Eigen::Matrix3d rot;
        double angle = PI*double(rod) / double(nrods);
        rot << 1.0, 0.0, 0.0,
            0.0, cos(angle), sin(angle),
            0.0, -sin(angle), cos(angle);

        RodState rs;
        int nverts = 100;
        rs.centerline.resize(nverts, 3);
        rs.centerlineVel.resize(nverts, 3);
        rs.centerlineVel.setZero();

        for (int i = 0; i < nverts; i++)
        {
            double r = 0.2 + 0.1 * sin(2.0*PI*double(i) / 50);
            rs.centerline(i, 0) = r*cos(2.0*PI*(double(i)-0.5) / 100.0);
            rs.centerline(i, 1) = r*sin(2.0*PI*(double(i)-0.5) / 100.0);
            rs.centerline(i, 2) = 0; //double(i) / nverts;            
            rs.centerline.row(i) = rs.centerline.row(i)*rot.transpose();
        }


        rs.directors.resize(nverts, 3);
        rs.directors.setZero();
        rs.directorAngVel.resize(nverts);
        rs.directorAngVel.setZero();

        rs.thetas.resize(nverts);
        rs.thetas.setZero();
        int nsegs = nverts;
        Eigen::VectorXd widths(nsegs);
        for (int i = 0; i < nverts; i++)
        {
            Eigen::Vector3d v1 = rs.centerline.row(i);
            Eigen::Vector3d v2 = rs.centerline.row((i + 1) % nverts);
            Eigen::Vector3d perp(0, 0, 1);
            perp = rot*perp;
            Eigen::Vector3d d = perp.cross(v2 - v1);
            rs.directors.row(i) = d;

            widths[i] = 0.01 + double(rod) / double(nrods) * 0.02 * (v1[2] + v2[2]);
        }
        for (int i = 0; i < nsegs; i++)
        {
            rs.directors.row(i) /= rs.directors.row(i).norm();            
        }

        Rod *therod = new Rod(rs, widths, params, true);
        config->addRod(therod);
    }

    for (int i = 0; i < nrods; i++)
    {
        for (int j = i + 1; j < nrods; j++)
        {
            Constraint c;
            c.rod1 = i;
            c.rod2 = j;
            c.seg1 = 0;
            c.seg2 = 0;
            c.bary1 = 0.5;
            c.bary2 = 0.5;
            c.stiffness = 100;
            config->addConstraint(c);
        }
    }

    std::ofstream ofs("example.rod");
    ofs << config->numRods() << std::endl;
    ofs << config->constraints.size() << std::endl;
    ofs << params.thickness << std::endl;
    ofs << Y << std::endl;
    ofs << params.rho << std::endl;
    ofs << std::endl << std::endl;
    for (int i = 0; i < config->numRods(); i++)
    {
        ofs << config->rods[i]->numVertices() << std::endl;
        ofs << (config->rods[i]->isClosed() ? '1' : '0') << std::endl;
        for (int j = 0; j < config->rods[i]->numVertices(); j++)
        {
            for (int k = 0; k < 3; k++)
            {
                ofs << config->rods[i]->startState.centerline(j, k) << " ";
            }            
        }
        ofs << std::endl;
        for (int j = 0; j < config->rods[i]->numSegments(); j++)
        {
            for (int k = 0; k < 3; k++)
            {
                ofs << config->rods[i]->startState.directors(j, k) << " ";
            }            
        }
        ofs << std::endl;
        for (int j = 0; j < config->rods[i]->numSegments(); j++)
        {
            ofs << config->rods[i]->widths[j] << " ";
        }
        ofs << std::endl;
        ofs << std::endl;
    }
    for (int i = 0; i < config->constraints.size(); i++)
    {
        ofs << config->constraints[i].rod1 << std::endl;
        ofs << config->constraints[i].rod2 << std::endl;
        ofs << config->constraints[i].seg1 << std::endl;
        ofs << config->constraints[i].seg2 << std::endl;
        ofs << config->constraints[i].bary1 << std::endl;
        ofs << config->constraints[i].bary2 << std::endl;
        ofs << config->constraints[i].stiffness << std::endl;
        ofs << std::endl;
    }*/
    
    createVisualizationMesh();
    dirty = true;
}


void RodsHook::showForces(int rod, const Eigen::MatrixXd &dE)
{
    int nverts = config->rods[rod]->numVertices();
    forcePoints.resize(2*nverts, 3);
    forceEdges.resize(nverts, 2);
    forceColors.resize(nverts, 3);
    for (int i = 0; i < nverts; i++)
    {
        forcePoints.row(2 * i) = config->rods[rod]->curState.centerline.row(i);
        forcePoints.row(2 * i+1) = config->rods[rod]->curState.centerline.row(i) - dE.row(i);
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

    std::vector<Eigen::MatrixXd> constraintdE;
    constraintEnergy(*config, &constraintdE);

    // update velocities
    for (int rod = 0; rod < nrods; rod++)
    {
        Eigen::MatrixXd dE;
        Eigen::VectorXd dtheta;
        double energy = rodEnergy(*config->rods[rod], config->rods[rod]->curState, &dE, &dtheta);

        dE += constraintdE[rod];

        int nverts = config->rods[rod]->numVertices();
        int nsegs = config->rods[rod]->numSegments();
        for (int i = 0; i < nverts; i++)
        {
            newresid += dE.row(i).squaredNorm() / config->rods[rod]->masses[i];
            dE.row(i) /= config->rods[rod]->masses[i];            
        }
        for (int i = 0; i < nsegs; i++)
        {
            newresid += dtheta[i]*dtheta[i] / config->rods[rod]->momInertia[i];
            dtheta[i] /= config->rods[rod]->momInertia[i];
        }

        config->rods[rod]->curState.centerlineVel -= dt*dE;
        config->rods[rod]->curState.directorAngVel -= dt*dtheta;

        double dampfactor = exp(-dt*damp);
        config->rods[rod]->curState.centerlineVel *= dampfactor;
        config->rods[rod]->curState.directorAngVel *= dampfactor;        
    }
    forceResidual = newresid;
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

