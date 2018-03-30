#include "RodParser.h"
#include "RodEnergy.h"
#include <fstream>
#include <iostream>

RodConfig *readRod(const char *filename)
{
    RodParams params;
    std::ifstream ifs(filename);
    if (!ifs)
    {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return NULL;
    }

    int nrods, nconstraints;
    ifs >> nrods >> nconstraints;
    ifs >> params.thickness;
    double Y;
    ifs >> Y;
    params.kstretching = Y;
    params.kbending = Y;
    params.ktwist = Y / 3.0;
    ifs >> params.rho;

    if (!ifs)
        return NULL;

    RodConfig *ret = new RodConfig();
    for (int i = 0; i < nrods; i++)
    {
        int nverts;
        ifs >> nverts;
        bool isclosed;
        ifs >> isclosed;
        int nsegs = isclosed ? nverts : nverts - 1;
        RodState rs;
        rs.centerline.resize(nverts, 3);
        for (int j = 0; j < nverts; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                ifs >> rs.centerline(j, k);
            }
        }        
        rs.centerlineVel.resize(nverts, 3);
        rs.centerlineVel.setZero();
        rs.directors.resize(nsegs, 3);
        for (int j = 0; j < nsegs; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                ifs >> rs.directors(j, k);
            }
            Eigen::Vector3d v0 = rs.centerline.row(j);
            Eigen::Vector3d v1 = rs.centerline.row((j + 1) % nverts);
            Eigen::Vector3d e = v1 - v0;
            e /= e.norm();
            double dote = rs.directors.row(j).dot(e);
            if (fabs(dote) > 1e-4)
                std::cout << dote << "Warning: directors not orthogonal to the centerline" << std::endl;
     //       rs.directors.row(j) -= dote * e.transpose();
            rs.directors.row(j) /= rs.directors.row(j).norm();  
            
        }        
        bool quiet = false;
        for (int j = 0; j < nsegs - 1; j++)
        {
            Eigen::Vector3d td = parallelTransport(rs.directors.row(j), rs.centerline.row(j+1) - rs.centerline.row(j), rs.centerline.row(j+2) - rs.centerline.row(j+1));
            
            if (td.dot(rs.directors.row(j + 1)) < 0)
            {
                if(!quiet)
                {
                    quiet=true;
                    std::cerr << "Warning: director angle > pi/2!" << std::endl;
                }
                rs.directors.row(j+1) *= -1;
            }
        }
        rs.thetas.resize(nsegs);
        rs.thetas.setZero();
        rs.directorAngVel.resize(nsegs);
        rs.directorAngVel.setZero();

        Eigen::VectorXd widths(nsegs);
        for (int i = 0; i < nsegs; i++)
        {
            ifs >> widths[i];
        }        
        if(nverts >= 2)
        {
            Rod *r = new Rod(rs, widths, params, isclosed);
            ret->addRod(r);
        }
        else
        {
            std::cerr << "Rod with only " << nverts << " vertices detected" << std::endl;
            exit(-1);
        }
    }
    int assignment = 1;
    for (int i = 0; i < nconstraints; i++)
    {
        Constraint c;
        ifs >> c.rod1 >> c.rod2;
        ifs >> c.seg1 >> c.seg2;
        if(c.rod1 < 0 || c.rod1 >= ret->numRods() || c.rod2 < 0 || c.rod2 >= ret->numRods())
        {
            std::cerr << "Bad rod in constraint " << i << std::endl;
            exit(-1);
        }
        if(c.seg1 < 0 || c.seg1 >= ret->rods[c.rod1]->numSegments() || c.seg2 < 0 || c.seg2 >= ret->rods[c.rod2]->numSegments())
        { 
            std::cerr << "Bad segment in constraint " << i << std::endl;
            std::cerr << "Segment 1: " << c.seg1 << "/" << ret->rods[c.rod1]->numSegments() << std::endl;
            std::cerr << "Segment 2: " << c.seg2 << "/" << ret->rods[c.rod2]->numSegments() << std::endl;
            exit(-1);
        }
        ifs >> c.bary1 >> c.bary2;
        ifs >> c.stiffness;
        c.assignment = assignment;
        c.visited = false;
        assignment *= -1;
        ret->addConstraint(c);
    }
    ret->initWeave(); // initializes constraints more intelligently

    // fix orientations

    Eigen::VectorXd orients(nrods);
    orients.setZero();
    while(true)
    {
        bool changed = true;
        while(changed)
        {
            changed = false;
            for(int i=0; i<nconstraints; i++)
            {
                int rod1 = ret->constraints[i].rod1;
                int rod2 = ret->constraints[i].rod2;
                if(orients[rod1] != 0 && orients[rod2] == 0)
                {
                    Eigen::Vector3d d1 = ret->rods[rod1]->startState.directors.row(ret->constraints[i].seg1);
                    Eigen::Vector3d d2 = ret->rods[rod2]->startState.directors.row(ret->constraints[i].seg2);
                    double factor = (d1.dot(d2) > 0 ? 1.0 : -1.0);
                    orients[rod2] = factor*orients[rod1];
                    changed = true;
                }
                else if(orients[rod2] != 0 && orients[rod1] == 0)
                {
                    Eigen::Vector3d d1 = ret->rods[rod1]->startState.directors.row(ret->constraints[i].seg1);
                    Eigen::Vector3d d2 = ret->rods[rod2]->startState.directors.row(ret->constraints[i].seg2);
                    double factor = (d1.dot(d2) > 0 ? 1.0 : -1.0);
                    orients[rod1] = factor*orients[rod2];
                    changed = true;
                }
            }
        }
        int next = -1;
        for(int i=0; i<nrods; i++)
        {
            if(orients[i] == 0)
                next = i;
        }       
        if(next == -1) break;
        orients[next] = 1.0;
    }

    for(int i=0; i<nrods; i++)
    {
        for(int j=0; j<ret->rods[i]->numSegments(); j++)
        {
            ret->rods[i]->startState.directors.row(j) *= orients[i];
            ret->rods[i]->curState.directors.row(j) *= orients[i];
        }
    }

    if (!ifs)
        exit(-1);
    return ret;
}
