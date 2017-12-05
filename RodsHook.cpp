#include "RodsHook.h"

void RodsHook::initGUI(igl::viewer::Viewer &viewer)
{
    params.kbending = 1.0;
    params.kstretching = 1.0;
    params.ktwist = 1.0 / 3.0;
    params.rho = 1.0;
    params.thickness = 1e-4;
    params.width = 0.01;
    dt = 1e-3;

    viewer.ngui->addGroup("Rod Parameters");
    viewer.ngui->addVariable("Thickness", params.thickness);
    viewer.ngui->addVariable("Width", params.width);
    viewer.ngui->addVariable("Stretching k", params.kstretching);
    viewer.ngui->addVariable("Bending k", params.kbending);
    viewer.ngui->addVariable("Twisting k", params.ktwist);

    viewer.ngui->addGroup("Sim Options");
    viewer.ngui->addVariable("Time Step", dt);
    
}

void RodsHook::createVisualizationMesh()
{
    int nverts = centerline.rows();
    Eigen::MatrixXd B(nverts, 3);
    for (int i = 0; i < nverts - 2; i++)
    {
        Eigen::Vector3d v0 = centerline.row(i).transpose();
        Eigen::Vector3d v1 = centerline.row(i+1).transpose();
        Eigen::Vector3d v2 = centerline.row(i+2).transpose();
        Eigen::Vector3d b = (v1 - v0).cross(v2 - v1);
        b /= b.norm();
        B.row(i + 1) = b.transpose();
    }
    B.row(0) = B.row(1);
    B.row(nverts - 1) = B.row(nverts - 2);
    int nedges = nverts - 1;
    Q.resize(8 * nedges, 3);
    F.resize(8 * nedges, 3);
    for (int i = 0; i < nedges; i++)
    {
        Eigen::Vector3d v0 = centerline.row(i).transpose();
        Eigen::Vector3d v1 = centerline.row(i+1).transpose();
        Eigen::Vector3d T = v1 - v0;
        T /= T.norm();
        Eigen::Vector3d B1 = B.row(i).transpose();
        Eigen::Vector3d B2 = B.row(i+1).transpose();
        Eigen::Vector3d N1 = B1.cross(T);
        Eigen::Vector3d N2 = B2.cross(T);
        Q.row(8 * i + 0) = (v0 + params.thickness / 2.0 * N1 - params.width / 2.0 * B1).transpose();
        Q.row(8 * i + 1) = (v0 + params.thickness / 2.0 * N1 + params.width / 2.0 * B1).transpose();
        Q.row(8 * i + 2) = (v0 - params.thickness / 2.0 * N1 + params.width / 2.0 * B1).transpose();
        Q.row(8 * i + 3) = (v0 - params.thickness / 2.0 * N1 - params.width / 2.0 * B1).transpose();
        Q.row(8 * i + 4) = (v1 + params.thickness / 2.0 * N2 - params.width / 2.0 * B2).transpose();
        Q.row(8 * i + 5) = (v1 + params.thickness / 2.0 * N2 + params.width / 2.0 * B2).transpose();
        Q.row(8 * i + 6) = (v1 - params.thickness / 2.0 * N2 + params.width / 2.0 * B2).transpose();
        Q.row(8 * i + 7) = (v1 - params.thickness / 2.0 * N2 - params.width / 2.0 * B2).transpose();
        for (int j = 0; j < 4; j++)
        {
            F(8 * i + 2 * j, 0) = 8 * i + j;
            F(8 * i + 2 * j, 2) = 8 * i + 4 + j;
            F(8 * i + 2 * j, 1) = 8 * i + 4 + ((j + 1) % 4);
            F(8 * i + 2 * j + 1, 0) = 8 * i + (j + 1) % 4;
            F(8 * i + 2 * j + 1, 2) = 8 * i + j;
            F(8 * i + 2 * j + 1, 1) = 8 * i + 4 + ((j + 1) % 4);
        }
    }
}