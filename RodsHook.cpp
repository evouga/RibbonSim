#include "RodsHook.h"
#include "RodParser.h"

void RodsHook::initGUI(igl::viewer::Viewer &viewer)
{
    savePrefix = "rod_";
    loadName = "../configs/torus.rod";

    viewer.ngui->addVariable("Config File", loadName);

    viewer.ngui->addGroup("Sim Options");
    viewer.ngui->addButton("Save Geometry", std::bind(&RodsHook::saveRods, this));
    viewer.ngui->addVariable("Save Prefix", savePrefix);
    viewer.ngui->addButton("Subdivide", std::bind(&RodsHook::linearSubdivision, this));

    viewer.ngui->addVariable("Orientation Weight", angleWeight);
    
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

    config = readRod(loadName.c_str());
    if (!config)
        exit(-1);
   
    createVisualizationMesh();
    dirty = true;
}

void RodsHook::createVisualizationMesh()
{
    config->createVisualizationMesh(Q, F);
}

double lineSearch(RodConfig &config, const Eigen::VectorXd &update, double angleWeight, bool optimizeWidths)
{
    double t = 1.0;
    double c1 = 0.1;
    double c2 = 0.9;
    double alpha = 0;
    double infinity = 1e6;
    double beta = infinity;

    Eigen::VectorXd r;
    Eigen::SparseMatrix<double> J;
    rAndJ(config, r, &J, angleWeight, optimizeWidths);
    
    Eigen::VectorXd dE;
    Eigen::VectorXd newdE;
    std::vector<RodState> start;
    std::vector<Eigen::VectorXd> startWidths;
    for (int i = 0; i < config.numRods(); i++)
    {
        start.push_back(config.rods[i]->curState);
        startWidths.push_back(config.rods[i]->widths);
    }
    double orig = 0.5 * r.transpose() * r;
    dE = J.transpose() * r;
    double deriv = -dE.dot(update);
    assert(deriv < 0);

    std::cout << "Starting line search, original energy " << orig << ", descent magnitude " << deriv << std::endl;

    while (true)
    {
        int nrods = config.numRods();
        int dofoffset = 0;
        for (int rod = 0; rod < nrods; rod++)
        {
            int nverts = config.rods[rod]->numVertices();
            int nsegs = config.rods[rod]->numSegments();
            if(!optimizeWidths)
            {
                for (int i = 0; i < nsegs; i++)
                {
                    Eigen::Vector3d oldv1 =start[rod].centerline.row(i);
                    Eigen::Vector3d oldv2 = start[rod].centerline.row((i + 1) % nverts);
                    Eigen::Vector3d v1 = oldv1 - t * update.segment<3>(dofoffset + 3 * i);
                    Eigen::Vector3d v2 = oldv2 - t * update.segment<3>(dofoffset + 3 * ((i + 1) % nverts));

                    config.rods[rod]->curState.directors.row(i) = parallelTransport(start[rod].directors.row(i), oldv2 - oldv1, v2 - v1);
                }
                for (int i = 0; i < nverts; i++)
                    config.rods[rod]->curState.centerline.row(i) = start[rod].centerline.row(i) - t * update.segment<3>(dofoffset + 3 * i).transpose();
                for (int i = 0; i < nsegs; i++)
                {
                    config.rods[rod]->curState.thetas[i] = start[rod].thetas[i] - t * update[dofoffset + 3 * nverts + i];
                }
            }
            else
            {
                for (int i = 0; i < nsegs; i++)
                {
                    double newwidth = startWidths[rod][i] - t * update[dofoffset + 3 * nverts + nsegs + i];
                    config.rods[rod]->widths[i] = newwidth;
                }
            }
            dofoffset += 3 * nverts + 2 * nsegs;
        }
        rAndJ(config, r, &J, angleWeight, optimizeWidths);

        double newenergy = 0.5 * r.transpose() * r;
        newdE = J.transpose() * r;

        std::cout << "Trying t = " << t << ", energy now " << newenergy << std::endl;

        if (std::isnan(newenergy) || newenergy > orig + t*deriv*c1)
        {
            beta = t;
            t = 0.5*(alpha + beta);
        }
        else if (-newdE.dot(update) < c2*deriv)
        {
            alpha = t;
            if (beta == infinity)
            {
                t = 2 * alpha;
            }
            else
            {
                t = 0.5*(alpha + beta);
            }

            if (beta - alpha < 1e-8)
            {
                return newdE.squaredNorm();
            }
        }
        else
        {
            return newdE.squaredNorm();
        }
    }
}

bool RodsHook::simulateOneStep()
{
    Eigen::VectorXd r;
    Eigen::SparseMatrix<double> Jr;
    rAndJ(*config, r, &Jr, angleWeight, false);

    std::cout << "Orig energy: " << 0.5 * r.transpose() * r << std::endl;
    Eigen::SparseMatrix<double> mat = Jr.transpose() * Jr;
    Eigen::SparseMatrix<double> I(mat.rows(), mat.cols());
    I.setIdentity();    
    double Treg = 1e-6;
    mat += Treg*I;

    Eigen::VectorXd rhs = Jr.transpose() * r;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver(mat);
    Eigen::VectorXd delta = solver.solve(rhs);
    if (solver.info() != Eigen::Success)
        exit(-1);
    std::cout << "Solver residual: " << (mat*delta - rhs).norm() << std::endl;
    forceResidual = lineSearch(*config, delta, angleWeight, false);
    iter++;

    createVisualizationMesh();    

    return false;
}

void RodsHook::saveRods()
{
    config->saveRodGeometry(savePrefix);
}

void RodsHook::showConstraints()
{
    int nconstraints = config->constraints.size();
    forcePoints.resize(2 *nconstraints,3);
    forceEdges.resize(nconstraints,2);
    forceColors.resize(nconstraints, 3);
    forceColors.col(0).setConstant(1);
    forceColors.col(1).setConstant(0);
    forceColors.col(2).setConstant(0);
    for (int i = 0; i < nconstraints; i++)
    {
        Constraint &c = config->constraints[i];
        Eigen::Vector3d v1 = config->rods[c.rod1]->curState.centerline.row(c.seg1);
        Eigen::Vector3d v2 = config->rods[c.rod1]->curState.centerline.row((c.seg1 + 1) % config->rods[c.rod1]->numVertices());
        Eigen::Vector3d w1 = config->rods[c.rod2]->curState.centerline.row(c.seg2);
        Eigen::Vector3d w2 = config->rods[c.rod2]->curState.centerline.row((c.seg2 + 1) % config->rods[c.rod2]->numVertices());
        forcePoints.row(2 * i) = (1.0 - c.bary1)*v1 + c.bary1 * v2;
        forcePoints.row(2 * i + 1) = (1.0 - c.bary2)*w1 + c.bary2 * w2;
        forceEdges(i, 0) = 2 * i;
        forceEdges(i, 1) = 2 * i + 1;
    }
}

void RodsHook::linearSubdivision()
{
    int nrods = config->numRods();
    for(int i=0; i<nrods; i++)
    {
        int nverts = config->rods[i]->numVertices();        
        RodState newstate;
        newstate.centerline.resize(2*nverts - 1,3);
        for(int j=0; j<nverts; j++)
        {
            newstate.centerline.row(2*j) = config->rods[i]->startState.centerline.row(j);
            if(j != nverts-1)
            {
                newstate.centerline.row(2*j+1) = 0.5 * (config->rods[i]->startState.centerline.row(j) + config->rods[i]->startState.centerline.row(j+1) );
            }
        }
        newstate.centerlineVel.resize(2*nverts - 1,3);
        newstate.centerlineVel.setZero();
        int nsegs = config->rods[i]->numSegments();
        newstate.directors.resize(2*nsegs, 3);
        for(int j=0; j<nsegs; j++)
        {
            newstate.directors.row(2*j) = config->rods[i]->startState.directors.row(j);
            newstate.directors.row(2*j+1) = config->rods[i]->startState.directors.row(j);
        }
        newstate.thetas.resize(2*nsegs);
        newstate.thetas.setZero();
        newstate.directorAngVel.resize(2*nsegs);
        newstate.directorAngVel.setZero();
        config->rods[i]->startState = newstate;
        config->rods[i]->curState = newstate;
        Eigen::VectorXd newwidths(2*nsegs);
        for(int j=0; j<nsegs; j++)
        {
            newwidths[2*j] = config->rods[i]->widths[j];
            newwidths[2*j+1] = config->rods[i]->widths[j];
        }
        config->rods[i]->widths = newwidths;
        config->rods[i]->initializeRestQuantities();
    }

    std::vector<Constraint> newconstraints;
    int nconstraints = config->constraints.size();
    for(int i=0; i<nconstraints; i++)
    {
        Constraint c = config->constraints[i];
        Constraint newc;
        newc.rod1 = c.rod1;
        newc.rod2 = c.rod2;
        newc.stiffness = c.stiffness;
        if(c.bary1 < 0.5)
        {
            newc.seg1 = 2*c.seg1;
            newc.bary1 = c.bary1 * 2.0;
        }
        else
        {
            newc.seg1 = 2*c.seg1 + 1;
            newc.bary1 = c.bary1 * 2.0 - 1.0;
        }
        if(c.bary2 < 0.5)
        {
            newc.seg2 = 2*c.seg2;
            newc.bary2 = c.bary2 * 2.0;
        }
        else
        {
            newc.seg2 = 2*c.seg2 + 1;
            newc.bary2 = c.bary2 * 2.0 - 1.0;
        }
        newconstraints.push_back(newc);
    }
    config->constraints = newconstraints;
    createVisualizationMesh();
    updateRenderGeometry();
}
