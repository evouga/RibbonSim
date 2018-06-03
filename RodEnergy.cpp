#include "RodEnergy.h"
#include <Eigen/Geometry>
#include <iostream>
#include <Eigen/Sparse>
#include <fstream>

double angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2, const Eigen::Vector3d &axis)
{
    double result =  2.0 * atan2((v1.cross(v2)).dot(axis), v1.norm()*v2.norm() + v1.dot(v2));
    return result;
}

Eigen::Vector3d perpToVector(const Eigen::Vector3d &v)
{
    Eigen::Vector3d result;
    int mincoord = 0;
    double minval = std::numeric_limits<double>::infinity();
    for (int i = 0; i < 3; i++)
    {
        if (fabs(v[i]) < minval)
        {
            mincoord = i;
            minval = fabs(v[i]);
        }
    }
    Eigen::Vector3d axis(0, 0, 0);
    axis[mincoord] = 1.0;
    result = axis.cross(v);
    result /= result.norm();
    return result;
}

void rAndJ(RodConfig &config, 
    Eigen::VectorXd &r, 
    Eigen::SparseMatrix<double> *Jr, 
    double angleWeight, 
    bool allowSliding, 
    Eigen::MatrixXd *anchorPoints, 
    Eigen::MatrixXd *anchorNormals)
{
    int nterms = 0;
    int ndofs = 0;
    for (int i = 0; i < config.numRods(); i++)
    {
        ndofs += 3 * config.rods[i]->numVertices() + config.rods[i]->numSegments();
        int njoints = config.rods[i]->isClosed() ? config.rods[i]->numSegments() : (config.rods[i]->numSegments() - 1);
        nterms += config.rods[i]->numSegments(); // stretching
        nterms += 2 * njoints; // bending
        nterms += njoints; // twisting
        if (anchorPoints)
            nterms += 3*config.rods[i]->numVertices();
    } 

    int baryoffset = ndofs;
    ndofs += 2 * config.numConstraints();

    nterms += 6 * config.constraints.size(); // constraint positions
    nterms += config.constraints.size(); // constraint directions
    nterms += 4 * config.constraints.size(); // barycentric coord inequality constraints

    r.resize(nterms);
    r.setConstant(std::numeric_limits<double>::infinity());

    if (Jr)
    {
        Jr->resize(nterms, ndofs);
    }

    std::vector<Eigen::Triplet<double> > J;

    int roffset = 0;
    int dofoffset = 0;
    for (int rodidx = 0; rodidx < config.numRods(); rodidx++)
    {
        RodState &state = config.rods[rodidx]->curState;
        Rod &rod = *config.rods[rodidx];
        int nverts = (int)state.centerline.rows();
        int nsegs = rod.numSegments();
        int thetaoffset = 3 * rod.numVertices();
        
        // stretching terms

        for (int i = 0; i < nsegs; i++)
        {
            Eigen::Vector3d v1 = state.centerline.row(i).transpose();
            Eigen::Vector3d v2 = state.centerline.row((i + 1) % nverts).transpose();
            double len = (v1 - v2).norm();
            double restlen = rod.restlens[i];
            double factor = 0.5 * rod.params.kstretching * rod.params.thickness / restlen;
            double segr = sqrt(factor * rod.widths[i]) * (len - restlen);
            r[roffset + i] = segr;

            if (Jr)
            {
                Eigen::Vector3d dlen = (v1 - v2) / len;
                for (int j = 0; j < 3; j++)
                {
                    J.push_back(Eigen::Triplet<double>(roffset + i, dofoffset + 3 * i + j, sqrt(factor*rod.widths[i])*dlen[j]));
                    J.push_back(Eigen::Triplet<double>(roffset + i, dofoffset + 3 * ((i + 1) % nverts) + j, -sqrt(factor*rod.widths[i])*dlen[j]));
                }
            }
        }

        roffset += rod.numSegments();

        // bending terms
        int startseg = rod.isClosed() ? 0 : 1;
        int idx = 0;
        for (int i = startseg; i < nsegs; i++)
        {
            Eigen::Vector3d v0 = state.centerline.row((nverts + i - 1) % nverts).transpose();
            Eigen::Vector3d v1 = state.centerline.row(i).transpose();
            Eigen::Vector3d v2 = state.centerline.row((i + 1) % nverts).transpose();
            Eigen::Vector3d t01 = (v1 - v0) / (v1 - v0).norm();
            Eigen::Vector3d t12 = (v2 - v1) / (v2 - v1).norm();
            Eigen::Vector3d kb = 2.0*t01.cross(t12) / (1.0 + t01.dot(t12));
            Eigen::Vector3d db11 = state.directors.row((nsegs + i - 1) % nsegs);
            Eigen::Vector3d db21 = t01.cross(db11);
            Eigen::Vector3d db12 = state.directors.row(i);
            Eigen::Vector3d db22 = t12.cross(db12);
            double theta1 = state.thetas[(nsegs + i - 1) % nsegs];
            double theta2 = state.thetas[i];
            Eigen::Vector3d d11 = db11*cos(theta1) + db21*sin(theta1);
            Eigen::Vector3d d21 = -db11*sin(theta1) + db21*cos(theta1);
            Eigen::Vector3d d12 = db12*cos(theta2) + db22*sin(theta2);
            Eigen::Vector3d d22 = -db12*sin(theta2) + db22*cos(theta2);
            double k1 = 0.5*(d21 + d22).dot(kb);
            double k2 = 0.5*(d11 + d12).dot(kb);
            double len = 0.5*(rod.restlens[(nsegs + i - 1) % nsegs] + rod.restlens[i]);
            double width = 0.5*(rod.widths[(nsegs + i - 1) % nsegs] + rod.widths[i]);
            double factor1 = 0.5*rod.params.kbending*rod.params.thickness*rod.params.thickness*rod.params.thickness / len;
            double factor2 = 0.5*rod.params.kbending*rod.params.thickness / len;
            double vertr1 = sqrt(factor1*width)*k1;
            double vertr2 = sqrt(factor2*width*width*width)*k2;
            r[roffset + 2 * idx] = vertr1;
            r[roffset + 2 * idx + 1] = vertr2;
            if (Jr)
            {
                Eigen::Vector3d ttilde = (t01 + t12) / (1.0 + t01.dot(t12));
                Eigen::Vector3d d1tilde = (d11 + d12) / (1.0 + t01.dot(t12));
                Eigen::Vector3d d2tilde = (d21 + d22) / (1.0 + t01.dot(t12));
                Eigen::Vector3d dk11 = 1.0 / (v1 - v0).norm() * (-k1 * ttilde + t12.cross(d2tilde));
                Eigen::Vector3d dk12 = 1.0 / (v2 - v1).norm() * (-k1 * ttilde - t01.cross(d2tilde));
                for (int j = 0; j < 3; j++)
                {
                    J.push_back(Eigen::Triplet<double>(roffset + 2 * idx, dofoffset + 3 * ((nverts + i - 1) % nverts) + j, -sqrt(factor1*width) * dk11[j]));
                    J.push_back(Eigen::Triplet<double>(roffset + 2 * idx, dofoffset + 3 * i + j, sqrt(factor1*width) * dk11[j] - sqrt(factor1*width) * dk12[j]));
                    J.push_back(Eigen::Triplet<double>(roffset + 2 * idx, dofoffset + 3 * ((i + 1) % nverts) + j, sqrt(factor1*width) * dk12[j]));
                }
                Eigen::Vector3d dk21 = 1.0 / (v1 - v0).norm() * (-k2 * ttilde + t12.cross(d1tilde));
                Eigen::Vector3d dk22 = 1.0 / (v2 - v1).norm() * (-k2 * ttilde - t01.cross(d1tilde));
                for (int j = 0; j < 3; j++)
                {
                    J.push_back(Eigen::Triplet<double>(roffset + 2 * idx + 1, dofoffset + 3 * ((nverts + i - 1) % nverts) + j, -sqrt(factor2*width*width*width) * dk21[j]));
                    J.push_back(Eigen::Triplet<double>(roffset + 2 * idx + 1, dofoffset + 3 * i + j, sqrt(factor2*width*width*width) * dk21[j] - sqrt(factor2*width*width*width) * dk22[j]));
                    J.push_back(Eigen::Triplet<double>(roffset + 2 * idx + 1, dofoffset + 3 * ((i + 1) % nverts) + j, sqrt(factor2*width*width*width) * dk22[j]));
                }

                Eigen::Vector3d Dd11 = -db11*sin(theta1) + db21*cos(theta1);
                Eigen::Vector3d Dd21 = -db11*cos(theta1) - db21*sin(theta1);
                Eigen::Vector3d Dd12 = -db12*sin(theta2) + db22*cos(theta2);
                Eigen::Vector3d Dd22 = -db12*cos(theta2) - db22*sin(theta2);
                J.push_back(Eigen::Triplet<double>(roffset + 2 * idx, dofoffset + thetaoffset + (nsegs + i - 1) % nsegs, 0.5*sqrt(factor1*width) * Dd21.dot(kb)));
                J.push_back(Eigen::Triplet<double>(roffset + 2 * idx + 1, dofoffset + thetaoffset + (nsegs + i - 1) % nsegs, 0.5*sqrt(factor2*width*width*width) * Dd11.dot(kb)));
                J.push_back(Eigen::Triplet<double>(roffset + 2 * idx, dofoffset + thetaoffset + i, 0.5*sqrt(factor1*width) * Dd22.dot(kb)));
                J.push_back(Eigen::Triplet<double>(roffset + 2 * idx + 1, dofoffset + thetaoffset + i, 0.5*sqrt(factor2*width*width*width) * Dd12.dot(kb)));
            }
            idx++;
        }

        roffset += 2 * idx;

        idx = 0;
        for (int i = startseg; i < nsegs; i++)
        {
            Eigen::Vector3d v0 = state.centerline.row((nverts + i - 1) % nverts).transpose();
            Eigen::Vector3d v1 = state.centerline.row(i).transpose();
            Eigen::Vector3d v2 = state.centerline.row((i + 1) % nverts).transpose();
            Eigen::Vector3d t01 = (v1 - v0) / (v1 - v0).norm();
            Eigen::Vector3d t12 = (v2 - v1) / (v2 - v1).norm();
            Eigen::Vector3d kb = 2.0*t01.cross(t12) / (1.0 + t01.dot(t12));
            Eigen::Vector3d db11 = state.directors.row((nsegs + i - 1) % nsegs);
            Eigen::Vector3d db21 = t01.cross(db11);
            Eigen::Vector3d db12 = state.directors.row(i);
            Eigen::Vector3d db22 = t12.cross(db12);
            double theta1 = state.thetas[(nsegs + i - 1) % nsegs];
            double theta2 = state.thetas[i];
            Eigen::Vector3d d1 = db11*cos(theta1) + db21*sin(theta1);
            Eigen::Vector3d d2 = db12*cos(theta2) + db22*sin(theta2);
            Eigen::Vector3d d1t = parallelTransport(d1, v1 - v0, v2 - v1);
            double theta = angle(d1t, d2, t12);
            double len = 0.5*(rod.restlens[(nsegs + i - 1) % nsegs] + rod.restlens[i]);
            double width = 0.5*(rod.widths[(nsegs + i - 1) % nsegs] + rod.widths[i]);
            double factor = 0.5*rod.params.ktwist*rod.params.thickness*rod.params.thickness*rod.params.thickness / len;
            r[roffset + idx] = sqrt(factor*width)*theta;
            if (Jr)
            {
                Eigen::Vector3d dtheta1 = 0.5*kb / (v1 - v0).norm();
                Eigen::Vector3d dtheta2 = 0.5*kb / (v2 - v1).norm();
                for (int j = 0; j < 3; j++)
                {
                    J.push_back(Eigen::Triplet<double>(roffset + idx, dofoffset + 3 * ((nverts + i - 1) % nverts) + j, -sqrt(factor*width)*dtheta1[j]));
                    J.push_back(Eigen::Triplet<double>(roffset + idx, dofoffset + 3 * i + j, sqrt(factor*width)*dtheta1[j] - sqrt(factor*width)*dtheta2[j]));
                    J.push_back(Eigen::Triplet<double>(roffset + idx, dofoffset + 3 * ((i + 1) % nverts) + j, sqrt(factor*width)*dtheta2[j]));
                }

                J.push_back(Eigen::Triplet<double>(roffset + idx, dofoffset + thetaoffset + (nsegs + i - 1) % nsegs, -sqrt(factor*width)));
                J.push_back(Eigen::Triplet<double>(roffset + idx, dofoffset + thetaoffset + i, sqrt(factor*width)));
            }
            idx++;
        }

        roffset += idx;
        dofoffset += 3 * rod.numVertices() + rod.numSegments();
    }

    int nconstraints = (int)config.constraints.size();
    // constraint energy rod1 -> rod2
    for (int i = 0; i < nconstraints; i++)
    {
        const Constraint &c = config.constraints[i];
        const RodState &rs1 = config.rods[c.rod1]->curState;
        const RodState &rs2 = config.rods[c.rod2]->curState;

        int nverts1 = config.rods[c.rod1]->numVertices();
        Eigen::Vector3d p1 = config.rods[c.rod1]->curState.centerline.row(c.seg1);
        Eigen::Vector3d p2 = config.rods[c.rod1]->curState.centerline.row((c.seg1 + 1) % nverts1);
        int nverts2 = config.rods[c.rod2]->numVertices();
        Eigen::Vector3d q1 = config.rods[c.rod2]->curState.centerline.row(c.seg2);
        Eigen::Vector3d q2 = config.rods[c.rod2]->curState.centerline.row((c.seg2 + 1) % nverts2);
        Eigen::Vector3d t01 = (p2 - p1) / (p2 - p1).norm();
        Eigen::Vector3d t12 = (q2 - q1) / (q2 - q1).norm();

        Eigen::Vector3d db11 = rs1.directors.row(c.seg1);
        Eigen::Vector3d db21 = t01.cross(db11);
        Eigen::Vector3d db12 = rs2.directors.row(c.seg2);
        Eigen::Vector3d db22 = t12.cross(db12);
        double theta1 = config.rods[c.rod1]->curState.thetas[c.seg1];
        double theta2 = config.rods[c.rod2]->curState.thetas[c.seg2];
        Eigen::Vector3d d1 = db11*cos(theta1) + db21*sin(theta1);
        Eigen::Vector3d d2 = db12*cos(theta2) + db22*sin(theta2);

        Eigen::Vector3d pt1 = (1.0 - c.bary1)*p1 + c.bary1*p2;
        Eigen::Vector3d pt2 = (1.0 - c.bary2)*q1 + c.bary2*q2;
        for (int j = 0; j < 3; j++)
        {
            r[roffset + 3 * i + j] = sqrt(0.5 * c.stiffness) * (pt1 + c.assignment * d1 * config.rods[c.rod1]->params.thickness - pt2)[j];
        }

        if (Jr)
        {
            int rod1offset = 0;
            int rod2offset = 0;
            for (int rod = 0; rod < c.rod1; rod++)
            {
                rod1offset += 3 * config.rods[rod]->numVertices() + config.rods[rod]->numSegments();
            }
            for (int rod = 0; rod < c.rod2; rod++)
            {
                rod2offset += 3 * config.rods[rod]->numVertices() + config.rods[rod]->numSegments();
            }
            for (int j = 0; j < 3; j++)
            {
                J.push_back(Eigen::Triplet<double>(roffset + 3 * i + j, rod1offset + 3 * c.seg1 + j, sqrt(0.5*c.stiffness) * (1.0 - c.bary1)));                
                J.push_back(Eigen::Triplet<double>(roffset + 3 * i + j, rod1offset + 3 * ((c.seg1 + 1) % nverts1) + j, sqrt(0.5*c.stiffness) * c.bary1));
                J.push_back(Eigen::Triplet<double>(roffset + 3 * i + j, rod2offset + 3 * c.seg2 + j, -sqrt(0.5*c.stiffness) * (1.0 - c.bary2)));
                J.push_back(Eigen::Triplet<double>(roffset + 3 * i + j, rod2offset + 3 * ((c.seg2 + 1) % nverts2) + j, -sqrt(0.5*c.stiffness) * c.bary2));
                for (int k = 0; k < 3; k++)
                {
                    J.push_back(Eigen::Triplet<double>(roffset + 3 * i + j, rod1offset + 3 * c.seg1 + k, sqrt(0.5*c.stiffness) * config.rods[c.rod1]->params.thickness * c.assignment * t01[j]*d1[k]/ (p2 - p1).norm()));
                    J.push_back(Eigen::Triplet<double>(roffset + 3 * i + j, rod1offset + 3 * ((c.seg1 + 1) % nverts1) + k, -sqrt(0.5*c.stiffness) * config.rods[c.rod1]->params.thickness * c.assignment * t01[j]*d1[k]/ (p2 - p1).norm()));                    
                }
                Eigen::Vector3d Dd1 = -db11*sin(theta1) + db21*cos(theta1);
                Eigen::Vector3d Dd2 = -db12*sin(theta2) + db22*cos(theta2);
                J.push_back(Eigen::Triplet<double>(roffset + 3 * i + j, rod1offset + 3 * config.rods[c.rod1]->numVertices() + c.seg1, c.assignment * sqrt(0.5*c.stiffness) * config.rods[c.rod1]->params.thickness * Dd1[j]));
                if (allowSliding)
                {
                    J.push_back(Eigen::Triplet<double>(roffset + 3 * i + j, baryoffset + 2 * i, sqrt(0.5*c.stiffness) * (p2 - p1)[j]));
                    J.push_back(Eigen::Triplet<double>(roffset + 3 * i + j, baryoffset + 2 * i + 1, -sqrt(0.5*c.stiffness) * (q2 - q1)[j]));
                }
            }
        }
    }
    roffset += 3 * nconstraints;

    // constraint energy rod2 -> rod1
    for (int i = 0; i < nconstraints; i++)
    {
        const Constraint &c = config.constraints[i];
        const RodState &rs1 = config.rods[c.rod1]->curState;
        const RodState &rs2 = config.rods[c.rod2]->curState;

        int nverts1 = config.rods[c.rod1]->numVertices();
        Eigen::Vector3d p1 = config.rods[c.rod1]->curState.centerline.row(c.seg1);
        Eigen::Vector3d p2 = config.rods[c.rod1]->curState.centerline.row((c.seg1 + 1) % nverts1);
        int nverts2 = config.rods[c.rod2]->numVertices();
        Eigen::Vector3d q1 = config.rods[c.rod2]->curState.centerline.row(c.seg2);
        Eigen::Vector3d q2 = config.rods[c.rod2]->curState.centerline.row((c.seg2 + 1) % nverts2);
        Eigen::Vector3d t01 = (p2 - p1) / (p2 - p1).norm();
        Eigen::Vector3d t12 = (q2 - q1) / (q2 - q1).norm();

        Eigen::Vector3d db11 = rs1.directors.row(c.seg1);
        Eigen::Vector3d db21 = t01.cross(db11);
        Eigen::Vector3d db12 = rs2.directors.row(c.seg2);
        Eigen::Vector3d db22 = t12.cross(db12);
        double theta1 = config.rods[c.rod1]->curState.thetas[c.seg1];
        double theta2 = config.rods[c.rod2]->curState.thetas[c.seg2];
        Eigen::Vector3d d1 = db11*cos(theta1) + db21*sin(theta1);
        Eigen::Vector3d d2 = db12*cos(theta2) + db22*sin(theta2);

        Eigen::Vector3d pt1 = (1.0 - c.bary1)*p1 + c.bary1*p2;
        Eigen::Vector3d pt2 = (1.0 - c.bary2)*q1 + c.bary2*q2;
        for (int j = 0; j < 3; j++)
        {
            r[roffset + 3 * i + j] = sqrt(0.5 * c.stiffness) * (pt1 + c.assignment * d2 * config.rods[c.rod2]->params.thickness - pt2)[j];
        }
            

        if (Jr)
        {
            int rod1offset = 0;
            int rod2offset = 0;
            for (int rod = 0; rod < c.rod1; rod++)
            {
                rod1offset += 3 * config.rods[rod]->numVertices() + config.rods[rod]->numSegments();
            }
            for (int rod = 0; rod < c.rod2; rod++)
            {
                rod2offset += 3 * config.rods[rod]->numVertices() + config.rods[rod]->numSegments();
            }
            for (int j = 0; j < 3; j++)
            {
                J.push_back(Eigen::Triplet<double>(roffset + 3 * i + j, rod1offset + 3 * c.seg1 + j, sqrt(0.5*c.stiffness) * (1.0 - c.bary1)));
                J.push_back(Eigen::Triplet<double>(roffset + 3 * i + j, rod1offset + 3 * ((c.seg1 + 1) % nverts1) + j, sqrt(0.5*c.stiffness) * c.bary1));
                J.push_back(Eigen::Triplet<double>(roffset + 3 * i + j, rod2offset + 3 * c.seg2 + j, -sqrt(0.5*c.stiffness) * (1.0 - c.bary2)));
                J.push_back(Eigen::Triplet<double>(roffset + 3 * i + j, rod2offset + 3 * ((c.seg2 + 1) % nverts2) + j, -sqrt(0.5*c.stiffness) * c.bary2));
                for (int k = 0; k < 3; k++)
                {
                    J.push_back(Eigen::Triplet<double>(roffset + 3 * i + j, rod2offset + 3 * c.seg2 + k, sqrt(0.5*c.stiffness) * config.rods[c.rod2]->params.thickness * c.assignment * t12[j]*d2[k]/ (q2 - q1).norm()));
                    J.push_back(Eigen::Triplet<double>(roffset + 3 * i + j, rod2offset + 3 * ((c.seg2 + 1) % nverts2) + k, -sqrt(0.5*c.stiffness) * config.rods[c.rod2]->params.thickness * c.assignment * t12[j]*d2[k]/ (q2 - q1).norm()));                    
                }

                Eigen::Vector3d Dd1 = -db11*sin(theta1) + db21*cos(theta1);
                Eigen::Vector3d Dd2 = -db12*sin(theta2) + db22*cos(theta2);
                J.push_back(Eigen::Triplet<double>(roffset + 3 * i + j, rod2offset + 3 * config.rods[c.rod2]->numVertices() + c.seg2, c.assignment * sqrt(0.5*c.stiffness) * config.rods[c.rod2]->params.thickness * Dd2[j]));

                if (allowSliding)
                {
                    J.push_back(Eigen::Triplet<double>(roffset + 3 * i + j, baryoffset + 2 * i, sqrt(0.5*c.stiffness) * (p2 - p1)[j]));
                    J.push_back(Eigen::Triplet<double>(roffset + 3 * i + j, baryoffset + 2 * i + 1, -sqrt(0.5*c.stiffness) * (q2 - q1)[j]));
                }
            }
        }
    }
    roffset += 3 * nconstraints;

    for (int i = 0; i < nconstraints; i++)
    {
        const Constraint &c = config.constraints[i];
        const RodState &rs1 = config.rods[c.rod1]->curState;
        const RodState &rs2 = config.rods[c.rod2]->curState;

        Eigen::Vector3d v0 = rs1.centerline.row(c.seg1).transpose();
        Eigen::Vector3d v1 = rs1.centerline.row((c.seg1 + 1) % config.rods[c.rod1]->numVertices()).transpose();

        Eigen::Vector3d w0 = rs2.centerline.row(c.seg2).transpose();
        Eigen::Vector3d w1 = rs2.centerline.row((c.seg2 + 1) % config.rods[c.rod2]->numVertices()).transpose();

        Eigen::Vector3d t01 = (v1 - v0) / (v1 - v0).norm();
        Eigen::Vector3d t12 = (w1 - w0) / (w1 - w0).norm();
        Eigen::Vector3d db11 = rs1.directors.row(c.seg1);
        Eigen::Vector3d db21 = t01.cross(db11);
        Eigen::Vector3d db12 = rs2.directors.row(c.seg2);
        Eigen::Vector3d db22 = t12.cross(db12);
        double theta1 = rs1.thetas[c.seg1];
        double theta2 = rs2.thetas[c.seg2];
        Eigen::Vector3d d1 = db11 * cos(theta1) + db21 * sin(theta1);
        Eigen::Vector3d d2 = db12 * cos(theta2) + db22 * sin(theta2);
        Eigen::Vector3d axis = d1.cross(d2);

        if (axis.norm() == 0.0)
        {
            r[roffset + i] = 0;
            continue;
        }
        axis /= axis.norm();
        double theta = angle(d1, d2, axis);
        double factor = angleWeight * 0.5 * c.stiffness;
        r[roffset + i] = sqrt(factor)*theta;
        if (Jr)
        {
            int rod1offset = 0;
            int rod2offset = 0;
            for (int rod = 0; rod < c.rod1; rod++)
            {
                rod1offset += 3 * config.rods[rod]->numVertices() + config.rods[rod]->numSegments();
            }
            for (int rod = 0; rod < c.rod2; rod++)
            {
                rod2offset += 3 * config.rods[rod]->numVertices() + config.rods[rod]->numSegments();
            }

            Eigen::Vector3d dtheta1 = -axis.cross(d1);
            Eigen::Vector3d dtheta2 = axis.cross(d2);
            for (int j = 0; j < 3; j++)
            {
                J.push_back(Eigen::Triplet<double>(roffset + i, rod1offset + 3 * c.seg1 + j, sqrt(factor)*dtheta1.dot(t01)*d1[j] / (v1-v0).norm()));
                J.push_back(Eigen::Triplet<double>(roffset + i, rod1offset + 3 * ((c.seg1 + 1) % config.rods[c.rod1]->numVertices()) + j, -sqrt(factor)*dtheta1.dot(t01)*d1[j] / (v1-v0).norm()));
                J.push_back(Eigen::Triplet<double>(roffset + i, rod2offset + 3 * c.seg2 + j, sqrt(factor)*dtheta2.dot(t12)*d2[j] / (w1-w0).norm()));
                J.push_back(Eigen::Triplet<double>(roffset + i, rod2offset + 3 * ((c.seg2 + 1) % config.rods[c.rod2]->numVertices()) + j, -sqrt(factor)*dtheta2.dot(t12)*d2[j] / (w1-w0).norm()));
            }

            int rod1thetaoffset = rod1offset + 3 * config.rods[c.rod1]->numVertices();
            int rod2thetaoffset = rod2offset + 3 * config.rods[c.rod2]->numVertices();
            Eigen::Vector3d d1dtheta = -db11 * sin(theta1) + db21 * cos(theta1);
            Eigen::Vector3d d2dtheta = -db12 * sin(theta2) + db22 * cos(theta2);

            J.push_back(Eigen::Triplet<double>(roffset + i, rod1thetaoffset + c.seg1, sqrt(factor)*dtheta1.dot(d1dtheta)));
            J.push_back(Eigen::Triplet<double>(roffset + i, rod2thetaoffset + c.seg2, sqrt(factor)*dtheta2.dot(d2dtheta)));
        }

    }

    roffset += nconstraints;
    double penstiffness = 1e4;
    for (int i = 0; i < nconstraints; i++)
    {
        const Constraint &c = config.constraints[i];
        if (!allowSliding || c.rod1 == c.rod2)
        {
            for(int j=0; j<4; j++)
                r[roffset + 4 * i + j] = 0;            
        }
        else
        {
            Constraint &c = config.constraints[i];
            if (c.bary1 > 0)
            {
                r[roffset + 4 * i + 0] = 0;
            }
            else
            {
                r[roffset + 4 * i + 0] = penstiffness * c.bary1;
                J.push_back(Eigen::Triplet<double>(roffset + 4 * i + 0, baryoffset + 2 * i, penstiffness));
            }
            if (c.bary1 < 1.0)
            {
                r[roffset + 4 * i + 1] = 0;
            }            
            else
            {
                r[roffset + 4 * i + 1] = penstiffness * (c.bary1 - 1.0);
                J.push_back(Eigen::Triplet<double>(roffset + 4 * i + 1, baryoffset + 2 * i, penstiffness));
            }
            if (c.bary2 > 0)
            {
                r[roffset + 4 * i + 2] = 0;
            }
            else
            {
                r[roffset + 4 * i + 2] = penstiffness * c.bary2;
                J.push_back(Eigen::Triplet<double>(roffset + 4 * i + 2, baryoffset + 2 * i + 1, penstiffness));
            }
            if (c.bary2 < 1.0)
            {
                r[roffset + 4 * i + 3] = 0;
            }            
            else
            {
                r[roffset + 4 * i + 3] = penstiffness * (c.bary2 - 1.0);
                J.push_back(Eigen::Triplet<double>(roffset + 4 * i + 3, baryoffset + 2 * i + 1, penstiffness));
            }
        }
    }
    roffset += 4 * nconstraints;

    if (anchorPoints)
    {
        int idx = 0;
        int dofoffset = 0;
        for (int rodidx = 0; rodidx < config.numRods(); rodidx++)
        {
            Rod *rod = config.rods[rodidx];
            double normalStiffness = 1000.0;
            double tangentStiffness = 1.0;
            int nverts = rod->numVertices();
            for (int i = 0; i < nverts; i++)
            {
                Eigen::Vector3d disp = rod->curState.centerline.row(i).transpose() - anchorPoints->row(idx/3).transpose();
                Eigen::Vector3d n = anchorNormals->row(idx/3);                
                Eigen::Vector3d t1 = perpToVector(n);
                Eigen::Vector3d t2 = n.cross(t1);
                r[roffset + idx] = normalStiffness * disp.dot(n);
                r[roffset + idx + 1] = tangentStiffness * disp.dot(t1);
                r[roffset + idx + 2] = tangentStiffness * disp.dot(t2);
                if (Jr)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        J.push_back(Eigen::Triplet<double>(roffset + idx, dofoffset + 3 * i + j, normalStiffness*n[j]));
                        J.push_back(Eigen::Triplet<double>(roffset + idx + 1, dofoffset + 3 * i + j, tangentStiffness*t1[j]));
                        J.push_back(Eigen::Triplet<double>(roffset + idx + 2, dofoffset + 3 * i + j, tangentStiffness*t2[j]));
                    }
                }
                idx+=3;
            }
            dofoffset += 3 * rod->numVertices() + rod->numSegments();
        }
        roffset += idx;
    }

    if (roffset != nterms)
        exit(-1);

    if (Jr)
        Jr->setFromTriplets(J.begin(), J.end());
}

Eigen::Vector3d parallelTransport(const Eigen::Vector3d &v, const Eigen::Vector3d &e1, const Eigen::Vector3d &e2)
{
    Eigen::Vector3d t1 = e1 / e1.norm();
    Eigen::Vector3d t2 = e2 / e2.norm();
    Eigen::Vector3d n = t1.cross(t2);
    if (n.norm() < 1e-8)
        return v;
    n /= n.norm();
    Eigen::Vector3d p1 = n.cross(t1);
    Eigen::Vector3d p2 = n.cross(t2);
    return v.dot(n)*n + v.dot(t1)*t2 + v.dot(p1)*p2;
}
