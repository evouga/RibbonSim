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
            rs.directors.row(j) /= rs.directors.row(j).norm();
        }
        rs.thetas.resize(nsegs);
        rs.thetas.setZero();
        rs.directorAngVel.resize(nsegs);
        rs.directorAngVel.setZero();

        Eigen::VectorXd widths(nsegs);
        for (int i = 0; i < nsegs; i++)
            ifs >> widths[i];
        Rod *r = new Rod(rs, widths, params, isclosed);
        ret->addRod(r);
    }
    for (int i = 0; i < nconstraints; i++)
    {
        Constraint c;
        ifs >> c.rod1 >> c.rod2;
        ifs >> c.seg1 >> c.seg2;
        ifs >> c.bary1 >> c.bary2;
        ifs >> c.stiffness;
        ret->addConstraint(c);
    }
    return ret;
}