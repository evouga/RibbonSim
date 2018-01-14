#include "RodEnergy.h"
#include <Eigen/Geometry>
#include <iostream>
#include <Eigen/Sparse>



double projectionEnergy(const Rod &rod, const RodState &state, Eigen::MatrixXd *dE)
{
    int nverts = state.centerline.rows();
    int nsegs = rod.isClosed() ? nverts: nverts - 1;
    double energy = 0;
    for (int i = 0; i < nsegs; i++)
    {
        Eigen::Vector3d v1 = state.centerline.row(i).transpose();
        Eigen::Vector3d v2 = state.centerline.row((i+1)%nverts).transpose();
        double len = (v1 - v2).norm();
        double restlen = rod.restlens[i];
        double factor = 0.5 * rod.params.kstretching * rod.widths[i] * rod.params.thickness / restlen;
        double segenergy = factor * (len - restlen)*(len - restlen);
        energy += segenergy;
        if (dE)
        {
            Eigen::Vector3d dlen = (v1 - v2) / len;
            dE->row(i) += 2.0*factor*(len - restlen)*dlen;
            dE->row((i + 1) % nverts) -= 2.0*factor*(len - restlen)*dlen;
        }
    }
    return energy;
}

double stretchingEnergy(const Rod &rod, const RodState &state, Eigen::MatrixXd *dE)
{
    int nverts = state.centerline.rows();
    int nsegs = rod.isClosed() ? nverts: nverts - 1;
    double energy = 0;
    for (int i = 0; i < nsegs; i++)
    {
        Eigen::Vector3d v1 = state.centerline.row(i).transpose();
        Eigen::Vector3d v2 = state.centerline.row((i+1)%nverts).transpose();
        double len = (v1 - v2).norm();
        double restlen = rod.restlens[i];
        double factor = 0.5 * rod.params.kstretching * rod.params.thickness / restlen;
        double segenergy = factor * rod.widths[i] * (len - restlen)*(len - restlen);
        energy += segenergy;
        if (dE)
        {
            Eigen::Vector3d dlen = (v1 - v2) / len;
            dE->segment<3>(3 * i) += 2.0*factor*rod.widths[i] * (len - restlen)*dlen;
            dE->segment<3>(3 * ((i + 1) % nverts)) -= 2.0*factor*rod.widths[i] * (len - restlen)*dlen;
        }
        if (dEdw)
        {
            Eigen::Vector3d dlen = (v1 - v2) / len;
            for (int j = 0; j < 3; j++)
            {
                dEdw->push_back(Eigen::Triplet<double>(3 * i + j, i, 2.0*factor* (len - restlen)*dlen[j]));
                dEdw->push_back(Eigen::Triplet<double>(3 * (i + 1) + j, i, -2.0*factor * (len - restlen)*dlen[j]));
            }
        }
    }
    return energy;
}

double bendingEnergy(const Rod &rod, const RodState &state, Eigen::VectorXd *dE, Eigen::VectorXd *dtheta, std::vector<Eigen::Triplet<double> > *dEdw, std::vector<Eigen::Triplet<double> > *dthetadw)
{
    double energy = 0;
    int nverts = state.centerline.rows();
    int nsegs = rod.isClosed() ? nverts : nverts - 1;
    int startseg = rod.isClosed() ? 0 : 1;
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
        double vertenergy = factor1*width*k1*k1 + factor2*width*width*width*k2*k2;
        energy += vertenergy;
        if (dE)
        {
            Eigen::Vector3d ttilde = (t01 + t12) / (1.0 + t01.dot(t12));
            Eigen::Vector3d d1tilde = (d11 + d12) / (1.0 + t01.dot(t12));
            Eigen::Vector3d d2tilde = (d21 + d22) / (1.0 + t01.dot(t12));
            Eigen::Vector3d dk11 = 1.0 / (v1-v0).norm() * (-k1 * ttilde + t12.cross(d2tilde));
            Eigen::Vector3d dk12 = 1.0 / (v2-v1).norm() * (-k1 * ttilde - t01.cross(d2tilde));
            dE->segment<3>(3 * ((nverts + i - 1) % nverts)) -= 2.0*factor1*width*k1*dk11.transpose();
            dE->segment<3>(3*i) += 2.0*factor1*width*k1*dk11.transpose() - 2.0*factor1*width*k1*dk12.transpose();
            dE->segment<3>(3 * ((i + 1) % nverts)) += 2.0*factor1*width*k1*dk12.transpose();
            Eigen::Vector3d dk21 = 1.0 / (v1-v0).norm() * (-k2 * ttilde + t12.cross(d1tilde));
            Eigen::Vector3d dk22 = 1.0 / (v2-v1).norm() * (-k2 * ttilde - t01.cross(d1tilde));
            dE->segment<3>(3 * ((nverts + i - 1) % nverts)) -= 2.0*factor2*width*width*width*k2*dk21.transpose();
            dE->segment<3>(3 * i) += 2.0*factor2*width*width*width*k2*dk21.transpose() - 2.0*factor2*width*width*width*k2*dk22.transpose();
            dE->segment<3>(3 * ((i + 1) % nverts)) += 2.0*factor2*width*width*width*k2*dk22.transpose();
        }
        if(dEdw)
        { 
            Eigen::Vector3d ttilde = (t01 + t12) / (1.0 + t01.dot(t12));
            Eigen::Vector3d d1tilde = (d11 + d12) / (1.0 + t01.dot(t12));
            Eigen::Vector3d d2tilde = (d21 + d22) / (1.0 + t01.dot(t12));
            Eigen::Vector3d dk11 = 1.0 / (v1-v0).norm() * (-k1 * ttilde + t12.cross(d2tilde));
            Eigen::Vector3d dk12 = 1.0 / (v2-v1).norm() * (-k1 * ttilde - t01.cross(d2tilde));
            Eigen::Vector3d dk21 = 1.0 / (v1-v0).norm() * (-k2 * ttilde + t12.cross(d1tilde));
            Eigen::Vector3d dk22 = 1.0 / (v2-v1).norm() * (-k2 * ttilde - t01.cross(d1tilde));

            for (int j = 0; j < 3; j++)
            {
                dEdw->push_back(Eigen::Triplet<double>(3 * ((nverts + i - 1) % nverts) + j, (nsegs + i - 1) % nsegs, -factor1*k1*dk11[j]));
                dEdw->push_back(Eigen::Triplet<double>(3 * ((nverts + i - 1) % nverts) + j, i, -factor1*k1*dk11[j]));
                dEdw->push_back(Eigen::Triplet<double>(3 * i + j, (nsegs + i - 1) % nsegs, factor1*k1*dk11[j] - factor1*k1*dk12[j]));
                dEdw->push_back(Eigen::Triplet<double>(3 * i + j, i, factor1*k1*dk11[j] - factor1*k1*dk12[j]));
                dEdw->push_back(Eigen::Triplet<double>(3 * ((i + 1) % nverts) + j, (nsegs + i - 1) % nsegs, factor1*k1*dk12[j]));
                dEdw->push_back(Eigen::Triplet<double>(3 * ((i + 1) % nverts) + j, i, factor1*k1*dk12[j]));
                dEdw->push_back(Eigen::Triplet<double>(3 * ((nverts + i - 1) % nverts) + j, (nsegs + i - 1) % nsegs, -3.0*factor2*width*width*k2*dk21[j]));
                dEdw->push_back(Eigen::Triplet<double>(3 * ((nverts + i - 1) % nverts) + j, i, -3.0*factor2*width*width*k2*dk21[j]));
                dEdw->push_back(Eigen::Triplet<double>(3 * i + j, (nsegs + i - 1) % nsegs, 3.0*factor2*width*width*k2*dk21[j] - 3.0*factor2*width*width*k2*dk22[j]));
                dEdw->push_back(Eigen::Triplet<double>(3 * i + j, i, 3.0*factor2*width*width*k2*dk21[j] - 3.0*factor2*width*width*k2*dk22[j]));
                dEdw->push_back(Eigen::Triplet<double>(3 * ((i + 1) % nverts) + j, (nsegs + i - 1) % nsegs, 3.0*factor2*width*width*k2*dk22[j]));
                dEdw->push_back(Eigen::Triplet<double>(3 * ((i + 1) % nverts) + j, i, 3.0*factor2*width*width*k2*dk22[j]));
            }
        }
        if (dtheta)
        {
            Eigen::Vector3d Dd11 = -db11*sin(theta1) + db21*cos(theta1);
            Eigen::Vector3d Dd21 = -db11*cos(theta1) - db21*sin(theta1);
            Eigen::Vector3d Dd12 = -db12*sin(theta2) + db22*cos(theta2);
            Eigen::Vector3d Dd22 = -db12*cos(theta2) - db22*sin(theta2);
            (*dtheta)[(nsegs + i - 1) % nsegs] += factor1*width*k1*Dd21.dot(kb);
            (*dtheta)[(nsegs + i - 1) % nsegs] += factor2*width*width*width*k2*Dd11.dot(kb);
            (*dtheta)[i] += factor1*width*k1*Dd22.dot(kb);
            (*dtheta)[i] += factor2*width*width*width*k2*Dd12.dot(kb);
        }
        if (dthetadw)
        {
            Eigen::Vector3d Dd11 = -db11*sin(theta1) + db21*cos(theta1);
            Eigen::Vector3d Dd21 = -db11*cos(theta1) - db21*sin(theta1);
            Eigen::Vector3d Dd12 = -db12*sin(theta2) + db22*cos(theta2);
            Eigen::Vector3d Dd22 = -db12*cos(theta2) - db22*sin(theta2);
            dthetadw->push_back(Eigen::Triplet<double>((nsegs + i - 1) % nsegs, (nsegs + i - 1) % nsegs, 0.5 * factor1*k1*Dd21.dot(kb)));
            dthetadw->push_back(Eigen::Triplet<double>((nsegs + i - 1) % nsegs, i, 0.5 * factor1*k1*Dd21.dot(kb)));
            dthetadw->push_back(Eigen::Triplet<double>((nsegs + i - 1) % nsegs, (nsegs + i - 1) % nsegs, 3.0 / 2.0 * factor2*width*width*k2*Dd11.dot(kb)));
            dthetadw->push_back(Eigen::Triplet<double>((nsegs + i - 1) % nsegs, i, 3.0 / 2.0 * factor2*width*width*k2*Dd11.dot(kb)));
            dthetadw->push_back(Eigen::Triplet<double>(i, (nsegs + i - 1) % nsegs, 0.5 * factor1*k1*Dd22.dot(kb)));
            dthetadw->push_back(Eigen::Triplet<double>(i, i, 0.5 * factor1*k1*Dd22.dot(kb)));
            dthetadw->push_back(Eigen::Triplet<double>(i, (nsegs + i - 1) % nsegs, 3.0 / 2.0 * factor2*width*width*k2*Dd12.dot(kb)));
            dthetadw->push_back(Eigen::Triplet<double>(i, i, 3.0 / 2.0 * factor2*width*width*k2*Dd12.dot(kb)));
        }
    }
    return energy;
}

double angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2, const Eigen::Vector3d &axis)
{
    double result =  2.0 * atan2((v1.cross(v2)).dot(axis), v1.norm()*v2.norm() + v1.dot(v2));
    return result;
}

double twistingEnergy(const Rod &rod, const RodState &state, Eigen::VectorXd *dE, Eigen::VectorXd *dtheta, std::vector<Eigen::Triplet<double> > *dEdw, std::vector<Eigen::Triplet<double> > *dthetadw)
{
    double energy = 0;
    int nverts = state.centerline.rows();
    int nsegs = rod.isClosed() ? nverts: nverts - 1;
    int startseg = rod.isClosed() ? 0 : 1;
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
        energy += factor*theta*theta;
        if (dE)
        {
            Eigen::Vector3d dtheta1 = 0.5*kb / (v1 - v0).norm();
            Eigen::Vector3d dtheta2 = 0.5*kb / (v2 - v1).norm();
            dE->segment<3>(3 * ((nverts + i - 1) % nverts)) += -2.0*factor*width*theta*dtheta1;
            dE->segment<3>(3 * i) += 2.0*factor*width*theta*dtheta1 - 2.0*factor*width*theta*dtheta2;
            dE->segment<3>(3 * ((i + 1) % nverts)) += 2.0*factor*width*theta*dtheta2;
        }
        if (dEdw)
        {
            Eigen::Vector3d dtheta1 = 0.5*kb / (v1 - v0).norm();
            Eigen::Vector3d dtheta2 = 0.5*kb / (v2 - v1).norm();
            for (int j = 0; j < 3; j++)
            {
                dEdw->push_back(Eigen::Triplet<double>(3 * ((nverts + i - 1) % nverts) + j, (nsegs + i - 1) % nsegs, -factor*theta*dtheta1[j]));
                dEdw->push_back(Eigen::Triplet<double>(3 * ((nverts + i - 1) % nverts) + j, i, -factor*theta*dtheta1[j]));
                dEdw->push_back(Eigen::Triplet<double>(3 * i + j, (nsegs + i - 1) % nsegs, factor*theta*dtheta1[j] - factor*theta*dtheta2[j]));
                dEdw->push_back(Eigen::Triplet<double>(3 * i + j, i, factor*theta*dtheta1[j] - factor*theta*dtheta2[j]));
                dEdw->push_back(Eigen::Triplet<double>(3 * ((i + 1) % nverts) + j, (nsegs + i - 1) % nsegs, factor*theta*dtheta2[j]));
                dEdw->push_back(Eigen::Triplet<double>(3 * ((i + 1) % nverts) + j, i, factor*theta*dtheta2[j]));
            }
        }
        if (dtheta)
        {
            (*dtheta)[(nsegs + i - 1) % nsegs] -= 2.0*factor*width*theta;
            (*dtheta)[i] += 2.0*factor*width*theta;
        }
        if (dthetadw)
        {
            dthetadw->push_back(Eigen::Triplet<double>((nsegs + i - 1) % nsegs, (nsegs + i - 1) % nsegs, -factor*theta));
            dthetadw->push_back(Eigen::Triplet<double>((nsegs + i - 1) % nsegs, i, -factor*theta));
            dthetadw->push_back(Eigen::Triplet<double>(i, (nsegs + i - 1) % nsegs, factor*theta));
            dthetadw->push_back(Eigen::Triplet<double>(i, i, factor*theta));
        }
    }
    return energy;
}

double rodEnergy(const Rod &rod, const RodState &state, Eigen::VectorXd *dE, Eigen::VectorXd *dtheta, std::vector<Eigen::Triplet<double> > *dEdw, std::vector<Eigen::Triplet<double> > *dthetadw)
{
    if (dE)
    {
        dE->resize(3 * state.centerline.rows());
        dE->setZero();
    }
    if (dtheta)
    {
        dtheta->resize(state.thetas.size());
        dtheta->setZero();
    }
    if (dEdw)
        dEdw->clear();
    if (dthetadw)
        dthetadw->clear();

    double senergy = stretchingEnergy(rod, state, dE, dEdw);
    double benergy = bendingEnergy(rod, state, dE, dtheta, dEdw, dthetadw);
    double tenergy = twistingEnergy(rod, state, dE, dtheta, dEdw, dthetadw);
    return senergy + benergy + tenergy;
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

double positionConstraintEnergy(RodConfig &config, std::vector<Eigen::VectorXd> *dEs)
{
    
    int nconstraints = (int)config.constraints.size();
    double totenergy = 0;
    for (int i = 0; i < nconstraints; i++)
    {
        const Constraint &c = config.constraints[i];
        int nverts1 = config.rods[c.rod1]->numVertices();
        Eigen::Vector3d p1 = config.rods[c.rod1]->curState.centerline.row(c.seg1);
        Eigen::Vector3d p2 = config.rods[c.rod1]->curState.centerline.row((c.seg1+1)%nverts1);
        int nverts2 = config.rods[c.rod2]->numVertices();
        Eigen::Vector3d q1 = config.rods[c.rod2]->curState.centerline.row(c.seg2);
        Eigen::Vector3d q2 = config.rods[c.rod2]->curState.centerline.row((c.seg2+1)%nverts2);

        Eigen::Vector3d pt1 = (1.0 - c.bary1)*p1 + c.bary1*p2;
        Eigen::Vector3d pt2 = (1.0 - c.bary2)*q1 + c.bary2*q2;
        double energy = 0.5*c.stiffness*(pt1 - pt2).dot(pt1 - pt2);
        totenergy += energy;

        if (dEs)
        {
            (*dEs)[c.rod1].segment<3>(3 * c.seg1) += c.stiffness*(1.0 - c.bary1)*(pt1 - pt2);
            (*dEs)[c.rod1].segment<3>(3 * ((c.seg1 + 1) % nverts1)) += c.stiffness*c.bary1*(pt1 - pt2);
            (*dEs)[c.rod2].segment<3>(3 * c.seg2) -= c.stiffness*(1.0 - c.bary2)*(pt1 - pt2);
            (*dEs)[c.rod2].segment<3>(3 * ((c.seg2 + 1) % nverts1)) -= c.stiffness*c.bary2*(pt1 - pt2);
        }
    }
    return totenergy;
}

double directorConstraintEnergy(RodConfig &config, std::vector<Eigen::VectorXd> *dEs, std::vector<Eigen::VectorXd> *dthetas)
{
    int nconstraints = (int)config.constraints.size();
    double totenergy = 0;
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
        Eigen::Vector3d kb = 2.0*t01.cross(t12) / (1.0 + t01.dot(t12));
        Eigen::Vector3d db11 = rs1.directors.row(c.seg1);
        Eigen::Vector3d db21 = t01.cross(db11);
        Eigen::Vector3d db12 = rs2.directors.row(c.seg2);
        Eigen::Vector3d db22 = t12.cross(db12);
        double theta1 = rs1.thetas[c.seg1];
        double theta2 = rs2.thetas[c.seg2];
        Eigen::Vector3d d1 = db11*cos(theta1) + db21*sin(theta1);
        Eigen::Vector3d d2 = db12*cos(theta2) + db22*sin(theta2);
        Eigen::Vector3d d1t = parallelTransport(d1, v1 - v0, w1 - w0);
        double theta = angle(d1t, d2, t12);
        double factor = 0.5 * c.stiffness;
        totenergy += factor*theta*theta;
        if (dEs)
        {
            Eigen::Vector3d dtheta1 = 0.5*kb / (v1 - v0).norm();
            Eigen::Vector3d dtheta2 = 0.5*kb / (w1 - w0).norm();
            (*dEs)[c.rod1].segment<3>(3 * c.seg1) += -2.0*factor*theta*dtheta1;
            (*dEs)[c.rod1].segment<3>(3 * ((c.seg1 + 1) % config.rods[c.rod1]->numVertices())) += 2.0*factor*theta*dtheta1;
            
            (*dEs)[c.rod2].segment<3>(3 * c.seg2) += -2.0*factor*theta*dtheta2;
            (*dEs)[c.rod2].segment<3>(3 * ((c.seg2 + 1) % config.rods[c.rod2]->numVertices())) += 2.0*factor*theta*dtheta2;
        }
        if (dthetas)
        {
            (*dthetas)[c.rod1][c.seg1] -= 2.0*factor*theta;
            (*dthetas)[c.rod2][c.seg2] += 2.0*factor*theta;
        }
    }
    return totenergy;
}

double constraintEnergy(RodConfig &config, std::vector<Eigen::VectorXd> *dEs, std::vector<Eigen::VectorXd> *dthetas)
{
    int nrods = config.numRods();
    double result = 0;
    if (dEs)
    {
        dEs->resize(nrods);
        for (int i = 0; i < nrods; i++)
        {
            int nverts = config.rods[i]->numVertices();
            (*dEs)[i].resize(3 * nverts);
            (*dEs)[i].setZero();
        }
    }
    if (dthetas)
    {
        dthetas->resize(nrods);
        for (int i = 0; i < nrods; i++)
        {
            int nsegs = config.rods[i]->numSegments();
            (*dthetas)[i].resize(nsegs);
            (*dthetas)[i].setZero();
        }
    }
    result += positionConstraintEnergy(config, dEs);
    result += directorConstraintEnergy(config, dEs, dthetas);

    return result;
}

