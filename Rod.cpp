#include "Rod.h"
#include <iostream>
#include <sstream>
#include <igl/writeOBJ.h>
#include <igl/read_triangle_mesh.h>
#include <igl/barycentric_coordinates.h>

#include <Eigen/Geometry>


Rod::Rod(const RodState &startState, RodParams &params, bool isClosed) : startState(startState), params(params), isClosed_(isClosed)
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
    initializeRestQuantities();
}


Eigen::Vector3d faceNormal(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, int faceidx)
{
    Eigen::Vector3d p0 = V.row(F(faceidx, 0));
    Eigen::Vector3d p1 = V.row(F(faceidx, 1));
    Eigen::Vector3d p2 = V.row(F(faceidx, 2));
    Eigen::Vector3d n = (p1 - p0).cross(p2 - p0);
    n /= n.norm();
    return n;
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


void Rod::updateProjectionVars(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
    Eigen::VectorXd sqrD;
    Eigen::VectorXi closestFaces;
    Eigen::MatrixXd c_point;
    tree->squared_distance(V, F, curState.centerline, sqrD, closestFaces, c_point);
    curState.closestFaceNormals   = Eigen::MatrixXd::Zero(closestFaces.rows(), 3);
    curState.closestFaceCentroids = Eigen::MatrixXd::Zero(closestFaces.rows(), 3);
    std::cout << curState.centerline.rows();
    for (int i = 0; i < closestFaces.rows(); i++)
    {
        int faceidx = closestFaces(i);
	curState.closestFaceNormals.row(i) = faceNormal( V, F, faceidx );
	Eigen::Vector3d p0 = V.row(F(faceidx, 0));
	Eigen::Vector3d p1 = V.row(F(faceidx, 1));
	Eigen::Vector3d p2 = V.row(F(faceidx, 2));
        curState.closestFaceCentroids.row(i) = ( p0 + p1 + p2 ) / 3.;
    }
}


void Rod::projectToMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
    Eigen::VectorXd sqrD;
    Eigen::VectorXi closestFaces;
    Eigen::MatrixXd c_point;
    tree->squared_distance(V, F, curState.centerline, sqrD, closestFaces, c_point);
    curState.closestFaceNormals   = Eigen::MatrixXd::Zero(closestFaces.rows(), 3);
    curState.closestFaceCentroids = Eigen::MatrixXd::Zero(closestFaces.rows(), 3);
    

    int nverts = curState.centerline.rows();
    int nsegs = isClosed() ? nverts : nverts - 1;
   
    curState.centerlineVel.setZero();
    curState.directorAngVel.setZero();
 //   curState.thetas.setZero();
    for (int i = 0; i < nsegs; i++)
    { 
	Eigen::Vector3d oldv1 = curState.centerline.row(i);
	Eigen::Vector3d oldv2 = curState.centerline.row((i + 1) % nverts);
	Eigen::Vector3d v1 = c_point.row(i);
	Eigen::Vector3d v2 = c_point.row((i + 1) % nverts);
        curState.directors.row(i) = parallelTransport(curState.directors.row(i), oldv2 - oldv1, v2 - v1); 
    }

    c_points = c_point; // for visualization
    curState.centerline = c_point;

    for (int i = 0; i < nsegs; i++)
    {
        Eigen::Vector3d v1 = curState.centerline.row(i).transpose();
        Eigen::Vector3d v2 = curState.centerline.row((i+1)%nverts).transpose();
        double len = (v1 - v2).norm();
        restlens[i] = len;
    }

    masses.resize(nverts);
    masses.setZero();
    for (int i = 0; i < nsegs; i++)
    {
        double len = restlens[i];
        double totmass = curState.widths[i]*params.thickness*len*params.rho;
        masses[i] += totmass / 2.0;
        masses[(i + 1)%nverts] += totmass / 2.0;
    }

    momInertia.resize(nsegs);
    for (int i = 0; i < nsegs; i++)
    {
        double len = restlens[i];
        double mass = curState.widths[i]*params.thickness*len*params.rho;
        momInertia[i] = mass / 12.0 * (curState.widths[i]*curState.widths[i] + params.thickness*params.thickness);
    }

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
        double totmass = startState.widths[i]*params.thickness*len*params.rho;
        masses[i] += totmass / 2.0;
        masses[(i + 1)%nverts] += totmass / 2.0;
    }

    momInertia.resize(nsegs);
    for (int i = 0; i < nsegs; i++)
    {
        double len = restlens[i];
        double mass = startState.widths[i]*params.thickness*len*params.rho;
        momInertia[i] = mass / 12.0 * (startState.widths[i]*startState.widths[i] + params.thickness*params.thickness);
    }
}



