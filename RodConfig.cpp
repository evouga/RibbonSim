#include "RodConfig.h"
#include <iostream>
#include <sstream>
#include <igl/writeOBJ.h>
#include <Eigen/Geometry>

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
        double dotprod = e.dot(startState.directors.row(i));
        assert(fabs(dotprod) < 1e-6);
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

void RodConfig::addConstraint(Constraint c)
{
    assert(c.rod1 >= 0 && c.rod1 < numRods());
    assert(c.rod2 >= 0 && c.rod2 < numRods());
    assert(c.seg1 >= 0 && c.seg1 < rods[c.rod1]->numSegments());
    assert(c.seg2 >= 0 && c.seg2 < rods[c.rod2]->numSegments());
    assert(c.bary1 >= 0.0 && c.bary1 <= 1.0);
    assert(c.bary2 >= 0.0 && c.bary2 <= 1.0);
    constraints.push_back(c);
}

struct weaveConstraint
{
    int cId;
    int segId; 
};

bool compWeaveConstraint(const weaveConstraint &a, const weaveConstraint &b)
{
    return a.segId < b.segId;
}

void RodConfig::initWeave()
{
// make a list of constraints per rod 
    std::vector< std::vector < weaveConstraint > > cmap; 
     
    for (int i = 0; i < constraints.size(); i++)
    {
        std::vector< weaveConstraint > m;
        cmap.push_back(m);
    } 
    
    for (int i = 0; i < constraints.size(); i++) 
    { 
        Constraint c = constraints[i];
        weaveConstraint wc1, wc2;
        wc1.cId = i;
        wc2.cId = i;
        wc1.segId = c.seg1;
        wc2.segId = c.seg2;
        
        cmap[c.rod1].push_back(wc1);
        cmap[c.rod2].push_back(wc2);
    }
 
    int startassignment = -1;
    for (int i = 0; i < constraints.size(); i++)
    {
        startassignment *= -1;
        int curassignment = startassignment;
        std::sort(cmap[i].begin(), cmap[i].end(), compWeaveConstraint);
        for (int j = 0; j < cmap[i].size(); j++)
        { 
            if ( constraints[cmap[i][j].cId].visited )
            { 
                curassignment = constraints[cmap[i][j].cId].assignment * -1;   
            } 
            else 
            { 
                constraints[cmap[i][j].cId].assignment = curassignment; 
                curassignment *= -1;
                constraints[cmap[i][j].cId].visited = true;
            } 

        } 
    } 

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
    setVisualizationMeshColors();
}

static double face_colors[9][3] = { { .001, .001, .95 },
                                          { 0.95, 0.001, 0.001 },
                                          { 0.95, 0.95, 0.001 },
                                          { 0.01, 0.95, 0.01 },
                                          { 0.95, 0.7, 0.01 },
                                          { 0.5, 0.01, 0.5 },
                                          { 0.001, 0.95, 0.95 },
                                          { .5, 0, 0},
                                          {.7, .7, .2 } };

void RodConfig::setVisualizationMeshColors()
{
    double minZ = 100000.0;
    double maxZ = -100000.0;
    for (int i = 0; i < numRods(); i++)
    {
        for (int j = 0; j < rods[i]->numSegments(); j++)
        {
            double curPos = rods[i]->startState.centerline(j, 1);
            if ( minZ > curPos )
            {
                minZ = curPos;
            }
            if ( maxZ <  curPos )
            {
                maxZ = curPos;
            }
        }
    }


    int num_colors = 7;
    for (int i = 0; i < numRods(); i++)
    {
        Eigen::MatrixXd col = Eigen::MatrixXd::Zero(rods[i]->numSegments(), 3);
        for (int j = 0; j < rods[i]->numSegments(); j++)
        {
          double zLevel = rods[i]->startState.centerline(j , 1);
          double scale = (maxZ - zLevel) / (maxZ - minZ);

          double tscale = (zLevel - minZ) / (maxZ - minZ);
          double highlight = tscale * tscale * tscale * tscale * .5;

          double* c = face_colors[(i % num_colors)];
          if ( scale > .995 )
          {
            c = face_colors[num_colors];
          }
          else if ( scale < .008 )
          {
            c = face_colors[num_colors + 1];
            scale = 1.;
            tscale = 1.;
          }
          else 
          {
              scale = 1 - scale * scale * scale * scale * scale;
          }

          col(j, 0) = std::min(c[0] * scale + highlight, 1.);
          col(j, 1) = std::min(c[1] * scale + highlight, 1.);
          col(j, 2) = std::min(c[2] * scale + highlight, 1.);
        }
        rods[i]->colors = col;
    }

}

void RodConfig::saveRodGeometry(const std::string &prefix)
{
    int nrods = (int)rods.size();
    
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

        Eigen::MatrixXd Q(4 * nverts, 3);
        Eigen::MatrixXi F(8 * nsegs, 3);

        double negThick = 100.;

        for (int i = 0; i < nsegs; i++)
        {
            Eigen::Vector3d v0 = rods[rod]->curState.centerline.row(i);
            Eigen::Vector3d v1 = rods[rod]->curState.centerline.row((i + 1) % nverts);
            Eigen::Vector3d T = v1 - v0;
            T /= T.norm();
            
            Q.row(4 * i + 0) = (v0.transpose() - rods[rod]->widths[i] / 2.0 * B.row(i) + rods[rod]->params.thickness / 2.0 * N.row(i));
            Q.row(4 * i + 1) = (v0.transpose() + rods[rod]->widths[i] / 2.0 * B.row(i) + rods[rod]->params.thickness / 2.0 * N.row(i));
            Q.row(4 * (i+1) + 0) = (v1.transpose() - rods[rod]->widths[i] / 2.0 * B.row(i) + rods[rod]->params.thickness / 2.0 * N.row(i));
            Q.row(4 * (i+1) + 1) = (v1.transpose() + rods[rod]->widths[i] / 2.0 * B.row(i) + rods[rod]->params.thickness / 2.0 * N.row(i));

            Q.row(4 * i + 2) = (v0.transpose() - rods[rod]->widths[i] / 2.0 * B.row(i) - negThick * rods[rod]->params.thickness / 2.0 * N.row(i));
            Q.row(4 * i + 3) = (v0.transpose() + rods[rod]->widths[i] / 2.0 * B.row(i) - negThick * rods[rod]->params.thickness / 2.0 * N.row(i));
            Q.row(4 * (i+1) + 2) = (v1.transpose() - rods[rod]->widths[i] / 2.0 * B.row(i) -  negThick * rods[rod]->params.thickness / 2.0 * N.row(i));
            Q.row(4 * (i+1) + 3) = (v1.transpose() + rods[rod]->widths[i] / 2.0 * B.row(i) -  negThick * rods[rod]->params.thickness / 2.0 * N.row(i));
            
            F(8 * i + 0, 0) = 4 * i + 0;
            F(8 * i + 0, 1) = 4 * i + 1;
            F(8 * i + 0, 2) = 4 * (i+1) + 0;
            F(8 * i + 1, 0) = 4 * (i+1) + 0;
            F(8 * i + 1, 1) = 4 * i + 1;
            F(8 * i + 1, 2) = 4 * (i+1) + 1;

            F(8 * i + 2, 0) = 4 * i + 3;
            F(8 * i + 2, 1) = 4 * i + 2;
            F(8 * i + 2, 2) = 4 * (i+1) + 2;
            F(8 * i + 3, 0) = 4 * (i+1) + 3;
            F(8 * i + 3, 1) = 4 * i + 3;
            F(8 * i + 3, 2) = 4 * (i+1) + 2;

            F(8 * i + 4, 0) = 4 * i + 0;
            F(8 * i + 4, 1) = 4 * (i+1) + 0;
            F(8 * i + 4, 2) = 4 * i + 2;
            F(8 * i + 5, 0) = 4 * i + 2;
            F(8 * i + 5, 1) = 4 * (i+1) + 0;
            F(8 * i + 5, 2) = 4 * (i+1) + 2;

            F(8 * i + 6, 0) = 4 * i + 3;
            F(8 * i + 6, 1) = 4 * (i+1) + 1;
            F(8 * i + 6, 2) = 4 * i + 1;
            F(8 * i + 7, 0) = 4 * i + 3;
            F(8 * i + 7, 1) = 4 * (i+1) + 3;
            F(8 * i + 7, 2) = 4 * (i+1) + 1;
        }

        std::stringstream ss;
        ss << prefix << rod << ".obj";
        igl::writeOBJ(ss.str(), Q, F);
    }
}
