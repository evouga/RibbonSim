#include "RodParser.h"
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
                std::cout << "Warning: directors not orthogonal to the centerline" << std::endl;
            rs.directors.row(j) -= dote * e.transpose();
            rs.directors.row(j) /= rs.directors.row(j).norm();  
            
        }        
        for (int j = 0; j < nsegs - 1; j++)
        {
            if (rs.directors.row(j).dot(rs.directors.row(j + 1)) < 0)
            {
                std::cerr << "Warning: director angle > pi/2!" << std::endl;
                //exit(-1);
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
        ret->addConstraint(c);
    }
    if (!ifs)
        exit(-1);
    return ret;
}
