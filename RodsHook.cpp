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

double lineSearch(RodConfig &config, const Eigen::VectorXd &update, double angleWeight)
{
    double t = 1.0;
    double c1 = 0.1;
    double c2 = 0.9;
    double alpha = 0;
    double infinity = 1e6;
    double beta = infinity;

    Eigen::VectorXd r;
    Eigen::SparseMatrix<double> J;
    rAndJ(config, r, &J, angleWeight);
    
    Eigen::VectorXd dE;
    Eigen::VectorXd newdE;
    std::vector<RodState> start;
    for (int i = 0; i < config.numRods(); i++)
        start.push_back(config.rods[i]->curState);
    
    double orig = 0.5 * r.squaredNorm();
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
                config.rods[rod]->curState.thetas[i] = start[rod].thetas[i] - t * update[dofoffset + 3 * nverts + i];

            dofoffset += 3 * nverts + 2 * nsegs;
        }
        
        rAndJ(config, r, &J, angleWeight);

        double newenergy = 0.5 * r.squaredNorm();
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
    rAndJ(*config, r, &Jr, angleWeight);

    std::cout << "Orig energy: " << r.squaredNorm() << std::endl;
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
    forceResidual = lineSearch(*config, delta, angleWeight);
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
