#include "RodEnergy.h"
#include <Eigen/Geometry>
#include <iostream>

Rod::Rod(const RodState &startState, const RodParams &params) : startState(startState), params(params)
{
    int nverts = startState.centerline.rows();
    assert(startState.ceterlineVel.rows() == nverts);
    assert(startState.directors.rows() == nverts - 1);
    assert(startState.directorAngVel.size() == nverts - 1);
    assert(startState.thetas.size() == nverts - 1);

    curState = startState;
    initializeRestQuantities();
}

void Rod::initializeRestQuantities()
{
    int nverts = startState.centerline.rows();
    restlens.resize(nverts-1);
    for (int i = 0; i < nverts - 1; i++)
    {
        Eigen::Vector3d v1 = startState.centerline.row(i).transpose();
        Eigen::Vector3d v2 = startState.centerline.row(i+1).transpose();
        double len = (v1 - v2).norm();
        restlens[i] = len;
    }

    masses.resize(nverts);
    masses.setZero();
    for (int i = 0; i < nverts - 1; i++)
    {
        double len = restlens[i];
        double totmass = params.width*params.thickness*len*params.rho;
        masses[i] += totmass / 2.0;
        masses[i + 1] += totmass / 2.0;
    }

    momInertia.resize(nverts - 1);
    for (int i = 0; i < nverts - 1; i++)
    {
        double len = restlens[i];
        double mass = params.width*params.thickness*len*params.rho;
        momInertia[i] = mass / 12.0 * (params.width*params.width + params.thickness*params.thickness);
    }
}

double stretchingEnergy(const Rod &rod, const RodState &state, Eigen::MatrixXd *dE)
{
    int nsegs = state.centerline.rows() - 1;
    double energy = 0;
    for (int i = 0; i < nsegs; i++)
    {
        Eigen::Vector3d v1 = state.centerline.row(i).transpose();
        Eigen::Vector3d v2 = state.centerline.row(i+1).transpose();
        double len = (v1 - v2).norm();
        double restlen = rod.restlens[i];
        double factor = 0.5 * rod.params.kstretching * rod.params.width * rod.params.thickness / restlen;
        double segenergy = factor * (len - restlen)*(len - restlen);
        energy += segenergy;
        if (dE)
        {
            Eigen::Vector3d dlen = (v1 - v2) / len;
            dE->row(i) += 2.0*factor*(len - restlen)*dlen;
            dE->row(i+1) -= 2.0*factor*(len - restlen)*dlen;
        }
    }
    return energy;
}

double bendingEnergy(const Rod &rod, const RodState &state, Eigen::MatrixXd *dE, Eigen::VectorXd *dtheta)
{
    double energy = 0;
    int nverts = state.centerline.rows();
    for (int i = 1; i < nverts - 1; i++)
    {
        Eigen::Vector3d v0 = state.centerline.row(i - 1).transpose();
        Eigen::Vector3d v1 = state.centerline.row(i).transpose();
        Eigen::Vector3d v2 = state.centerline.row(i + 1).transpose();
        Eigen::Vector3d t01 = (v1 - v0) / (v1 - v0).norm();
        Eigen::Vector3d t12 = (v2 - v1) / (v2 - v1).norm();
        Eigen::Vector3d kb = 2.0*t01.cross(t12) / (1.0 + t01.dot(t12));
        Eigen::Vector3d db11 = state.directors.row(i - 1);
        Eigen::Vector3d db21 = t01.cross(db11);
        Eigen::Vector3d db12 = state.directors.row(i);
        Eigen::Vector3d db22 = t12.cross(db12);
        double theta1 = state.thetas[i-1];
        double theta2 = state.thetas[i];
        Eigen::Vector3d d11 = db11*cos(theta1) + db21*sin(theta1);
        Eigen::Vector3d d21 = -db11*sin(theta1) + db21*cos(theta1);
        Eigen::Vector3d d12 = db12*cos(theta2) + db22*sin(theta2);
        Eigen::Vector3d d22 = -db12*sin(theta2) + db22*cos(theta2);
        double k1 = 0.5*(d21 + d22).dot(kb);
        double k2 = 0.5*(d11 + d12).dot(kb);
        double len = 0.5*(rod.restlens[i - 1] + rod.restlens[i]);
        double factor1 = 0.5*rod.params.kbending*rod.params.width*rod.params.thickness*rod.params.thickness*rod.params.thickness / len;
        double factor2 = 0.5*rod.params.kbending*rod.params.width*rod.params.width*rod.params.width*rod.params.thickness / len;
        double vertenergy = factor1*k1*k1 + factor2*k2*k2;
        energy += vertenergy;
        if (dE)
        {
            Eigen::Vector3d ttilde = (t01 + t12) / (1.0 + t01.dot(t12));
            Eigen::Vector3d d1tilde = (d11 + d12) / (1.0 + t01.dot(t12));
            Eigen::Vector3d d2tilde = (d21 + d22) / (1.0 + t01.dot(t12));
            Eigen::Vector3d dk11 = 1.0 / (v1-v0).norm() * (-k1 * ttilde + t12.cross(d2tilde));
            Eigen::Vector3d dk12 = 1.0 / (v2-v1).norm() * (-k1 * ttilde - t01.cross(d2tilde));
            dE->row(i - 1) -= 2.0*factor1*k1*dk11.transpose();
            dE->row(i) += 2.0*factor1*k1*dk11.transpose() - 2.0*factor1*k1*dk12.transpose();
            dE->row(i + 1) += 2.0*factor1*k1*dk12.transpose();
            Eigen::Vector3d dk21 = 1.0 / (v1-v0).norm() * (-k2 * ttilde + t12.cross(d1tilde));
            Eigen::Vector3d dk22 = 1.0 / (v2-v1).norm() * (-k2 * ttilde - t01.cross(d1tilde));
            dE->row(i - 1) -= 2.0*factor2*k2*dk21.transpose();
            dE->row(i) += 2.0*factor2*k2*dk21.transpose() - 2.0*factor2*k2*dk22.transpose();
            dE->row(i + 1) += 2.0*factor2*k2*dk22.transpose();
        }
        if (dtheta)
        {
            Eigen::Vector3d Dd11 = -db11*sin(theta1) + db21*cos(theta1);
            Eigen::Vector3d Dd21 = -db11*cos(theta1) - db21*sin(theta1);
            Eigen::Vector3d Dd12 = -db12*sin(theta2) + db22*cos(theta2);
            Eigen::Vector3d Dd22 = -db12*cos(theta2) - db22*sin(theta2);
            (*dtheta)[i - 1] += factor1*k1*Dd21.dot(kb);
            (*dtheta)[i - 1] += factor2*k2*Dd11.dot(kb);
            (*dtheta)[i] += factor1*k1*Dd22.dot(kb);
            (*dtheta)[i] += factor2*k2*Dd12.dot(kb);
        }
    }
    return energy;
}

double angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2, const Eigen::Vector3d &axis)
{
    return 2.0 * atan2((v1.cross(v2)).dot(axis), v1.norm()*v2.norm() + v1.dot(v2));
}

double twistingEnergy(const Rod &rod, const RodState &state, Eigen::MatrixXd *dE, Eigen::VectorXd *dtheta)
{
    double energy = 0;
    int nverts = state.centerline.rows();
    for (int i = 1; i < nverts - 1; i++)
    {
        Eigen::Vector3d v0 = state.centerline.row(i - 1).transpose();
        Eigen::Vector3d v1 = state.centerline.row(i).transpose();
        Eigen::Vector3d v2 = state.centerline.row(i + 1).transpose();
        Eigen::Vector3d t01 = (v1 - v0) / (v1 - v0).norm();
        Eigen::Vector3d t12 = (v2 - v1) / (v2 - v1).norm();
        Eigen::Vector3d kb = 2.0*t01.cross(t12) / (1.0 + t01.dot(t12));
        Eigen::Vector3d db11 = state.directors.row(i - 1);
        Eigen::Vector3d db21 = t01.cross(db11);
        Eigen::Vector3d db12 = state.directors.row(i);
        Eigen::Vector3d db22 = t12.cross(db12);
        double theta1 = state.thetas[i - 1];
        double theta2 = state.thetas[i];
        Eigen::Vector3d d1 = db11*cos(theta1) + db21*sin(theta1);
        Eigen::Vector3d d2 = db12*cos(theta2) + db22*sin(theta2);
        Eigen::Vector3d d1t = parallelTransport(d1, v1 - v0, v2 - v1);
        double theta = angle(d1t, d2, t12);
        double len = 0.5*(rod.restlens[i - 1] + rod.restlens[i]);
        double factor = 0.5*rod.params.ktwist*rod.params.width*rod.params.thickness*rod.params.thickness*rod.params.thickness / len;
        energy += factor*theta*theta;
        if (dE)
        {
            Eigen::Vector3d dtheta1 = 0.5*kb / (v1 - v0).norm();
            Eigen::Vector3d dtheta2 = 0.5*kb / (v2 - v1).norm();
            dE->row(i - 1) += -2.0*factor*theta*dtheta1;
            dE->row(i) += 2.0*factor*theta*dtheta1 - 2.0*factor*theta*dtheta2;
            dE->row(i + 1) += 2.0*factor*theta*dtheta2;
        }
        if (dtheta)
        {
            (*dtheta)[i - 1] -= 2.0*factor*theta;
            (*dtheta)[i] += 2.0*factor*theta;
        }
    }
    return energy;
}

double rodEnergy(const Rod &rod, const RodState &state, Eigen::MatrixXd *dE, Eigen::VectorXd *dtheta)
{
    if (dE)
    {
        dE->resize(state.centerline.rows(), 3);
        dE->setZero();
    }
    if (dtheta)
    {
        dtheta->resize(state.thetas.size());
        dtheta->setZero();
    }
    double energy = 0;
    energy += stretchingEnergy(rod, state, dE);
    energy += bendingEnergy(rod, state, dE, dtheta);
    energy += twistingEnergy(rod, state, dE, dtheta);
    return energy;
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