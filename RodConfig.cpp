#include "RodConfig.h"

Rod::Rod(const RodState &startState, const Eigen::VectorXd &segwidths, RodParams &params, bool isClosed) : startState(startState), params(params), isClosed_(isClosed)
{
    int nverts = startState.centerline.rows();
    int nsegs = isClosed ? nverts : nverts - 1;
    assert(startState.centerlineVel.rows() == nverts);
    assert(startState.directors.rows() == nsegs);
    assert(startState.directorAngVel.size() == nsegs);
    assert(startState.thetas.size() == nsegs);
    assert(segwidths.size() == nsegs);
    for (int i = 0; i < nsegs; i++)
    {
        Eigen::Vector3d v1 = startState.centerline.row(i);
        Eigen::Vector3d v2 = startState.centerline.row((i+1)%nverts);
        Eigen::Vector3d e = (v2 - v1);
        e /= e.norm();
        assert(fabs(e.dot(startState.directors.row(i))) < 1e-6);
    }
    curState = startState;
    widths = segwidths;
    initializeRestQuantities();
}

void Rod::initializeRestQuantities()
{
    int nverts = startState.centerline.rows();
    int nsegs = isClosed() ? nverts : nverts - 1;
    restlens.resize(nsegs);
    for (int i = 0; i < nsegs; i++)
    {
        Eigen::Vector3d v1 = startState.centerline.row(i).transpose();
        Eigen::Vector3d v2 = startState.centerline.row((i+1)%nverts).transpose();
        double len = (v1 - v2).norm();
        restlens[i] = len;
    }

    masses.resize(nverts);
    masses.setZero();
    for (int i = 0; i < nsegs; i++)
    {
        double len = restlens[i];
        double totmass = widths[i]*params.thickness*len*params.rho;
        masses[i] += totmass / 2.0;
        masses[(i + 1)%nverts] += totmass / 2.0;
    }

    momInertia.resize(nsegs);
    for (int i = 0; i < nsegs; i++)
    {
        double len = restlens[i];
        double mass = widths[i]*params.thickness*len*params.rho;
        momInertia[i] = mass / 12.0 * (widths[i]*widths[i] + params.thickness*params.thickness);
    }
}

RodConfig::~RodConfig()
{
    int nrods = numRods();
    for (int i = 0; i < nrods; i++)
        delete rods[i];
}

void RodConfig::addRod(Rod *rod)
{
    rods.push_back(rod);
}

void RodConfig::reset()
{
    int nrods = (int)rods.size();
    for (int i = 0; i < nrods; i++)
    {
        rods[i]->curState = rods[i]->startState;
    }
}

void RodConfig::createVisualizationMesh(Eigen::MatrixXd &Q, Eigen::MatrixXi &F)
{
    int nrods = (int)rods.size();
    int totalsegs = 0;
    for (int i = 0; i < nrods; i++)
    {
        totalsegs += rods[i]->numSegments();
    }
    Q.resize(8 * totalsegs, 3);
    F.resize(8 * totalsegs, 3);

    int offset = 0;

    for (int rod = 0; rod < nrods; rod++)
    {
        int nverts = rods[rod]->numVertices();
        int nsegs = rods[rod]->numSegments();
        Eigen::MatrixXd N(nsegs, 3);
        Eigen::MatrixXd B(nsegs, 3);
        for (int i = 0; i < nsegs; i++)
        {
            Eigen::Vector3d v0 = rods[rod]->curState.centerline.row(i);
            Eigen::Vector3d v1 = rods[rod]->curState.centerline.row((i + 1) % nverts);
            Eigen::Vector3d e = v1 - v0;
            e /= e.norm();
            Eigen::Vector3d d1 = rods[rod]->curState.directors.row(i);
            Eigen::Vector3d d2 = e.cross(d1);
            double theta = rods[rod]->curState.thetas[i];
            N.row(i) = d1*cos(theta) + d2*sin(theta);
            B.row(i) = -d1*sin(theta) + d2*cos(theta);
        }

        for (int i = 0; i < nsegs; i++)
        {
            Eigen::Vector3d v0 = rods[rod]->curState.centerline.row(i);
            Eigen::Vector3d v1 = rods[rod]->curState.centerline.row((i + 1) % nverts);
            Eigen::Vector3d T = v1 - v0;
            T /= T.norm();
            Q.row(offset + 8 * i + 0) = (v0.transpose() + rods[rod]->params.thickness / 2.0 * N.row(i) - rods[rod]->widths[i] / 2.0 * B.row(i));
            Q.row(offset + 8 * i + 1) = (v0.transpose() + rods[rod]->params.thickness / 2.0 * N.row(i) + rods[rod]->widths[i] / 2.0 * B.row(i));
            Q.row(offset + 8 * i + 2) = (v0.transpose() - rods[rod]->params.thickness / 2.0 * N.row(i) + rods[rod]->widths[i] / 2.0 * B.row(i));
            Q.row(offset + 8 * i + 3) = (v0.transpose() - rods[rod]->params.thickness / 2.0 * N.row(i) - rods[rod]->widths[i] / 2.0 * B.row(i));
            Q.row(offset + 8 * i + 4) = (v1.transpose() + rods[rod]->params.thickness / 2.0 * N.row(i) - rods[rod]->widths[i] / 2.0 * B.row(i));
            Q.row(offset + 8 * i + 5) = (v1.transpose() + rods[rod]->params.thickness / 2.0 * N.row(i) + rods[rod]->widths[i] / 2.0 * B.row(i));
            Q.row(offset + 8 * i + 6) = (v1.transpose() - rods[rod]->params.thickness / 2.0 * N.row(i) + rods[rod]->widths[i] / 2.0 * B.row(i));
            Q.row(offset + 8 * i + 7) = (v1.transpose() - rods[rod]->params.thickness / 2.0 * N.row(i) - rods[rod]->widths[i] / 2.0 * B.row(i));
            for (int j = 0; j < 4; j++)
            {
                F(offset + 8 * i + 2 * j, 0) = offset + 8 * i + j;
                F(offset + 8 * i + 2 * j, 2) = offset + 8 * i + 4 + j;
                F(offset + 8 * i + 2 * j, 1) = offset + 8 * i + 4 + ((j + 1) % 4);
                F(offset + 8 * i + 2 * j + 1, 0) = offset + 8 * i + (j + 1) % 4;
                F(offset + 8 * i + 2 * j + 1, 2) = offset + 8 * i + j;
                F(offset + 8 * i + 2 * j + 1, 1) = offset + 8 * i + 4 + ((j + 1) % 4);
            }
        }
        offset += 8 * rods[rod]->numSegments();
    }
}