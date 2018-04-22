#include "RodsHook.h"
#include "RodParser.h"
#include "utility/simple_svg_1.0.0.hpp"
#include <igl/unproject_onto_mesh.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <Eigen/Dense>
#include <igl/point_mesh_squared_distance.h>

RodsHook::RodsHook() : PhysicsHook(), iter(0), forceResidual(0.0), angleWeight(1e3), newWidth(0.01), dirty(true), config(NULL) 
{
    savePrefix = "rod_";
    loadName = "../configs/torus7.rod";
    targetMeshName = "../meshes/torus.obj";
    visualizeConstraints = true;
    allowSliding = false;
    stickToMesh = false;

    Q.resize(0, 3);
    F.resize(0, 3);
    renderQ.resize(0, 3);
    renderF.resize(0, 3);
}

void RodsHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
{
    if (ImGui::CollapsingHeader("Configuration", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputText("Config File", loadName);
    }
    if (ImGui::CollapsingHeader("Rod Controls", ImGuiTreeNodeFlags_DefaultOpen))
    {
        if (ImGui::Button("Save Geometry", ImVec2(-1, 0)))
        {
            saveRods();
        }
        ImGui::InputText("Save Prefix", savePrefix);
        if (ImGui::Button("Subdivide", ImVec2(-1, 0)))
        {
            linearSubdivision();
        }
        if (ImGui::Button("Export Weave", ImVec2(-1, 0)))
        {
            exportWeave();
        }
        ImGui::Checkbox("Show Constraints", &visualizeConstraints);
        bool oldVisuaizeTargetMesh = visualizeTargetMesh;
        ImGui::Checkbox("Show Target Mesh", &visualizeTargetMesh);
        if (oldVisuaizeTargetMesh != visualizeTargetMesh)
            updateRenderGeometry();

        if (ImGui::Button("Set Widths", ImVec2(-1, 0)))
        {
            setWidths();
        }
        ImGui::InputFloat("New Widths", &newWidth);
    }

    if (ImGui::CollapsingHeader("Sim Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputFloat("Orientation Weight", &angleWeight);
        ImGui::Checkbox("Allow Sliding", &allowSliding);
        ImGui::InputText("Target Mesh", targetMeshName);
        if (ImGui::Button("Load Mesh", ImVec2(-1, 0)))
        {
            loadTargetMesh();
        }
        ImGui::Checkbox("Stick to Target Mesh", &stickToMesh);
    }

    if (ImGui::CollapsingHeader("Sim Status", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Text("Iteration %d", iter);
        ImGui::Text("Force Residual %f", forceResidual);
    }
}

bool RodsHook::mouseClicked(igl::opengl::glfw::Viewer &viewer, int button)
{
    int fid;
    Eigen::Vector3f bc;
    // Cast a ray in the view direction starting from the mouse position
    double x = viewer.current_mouse_x;
    double y = viewer.core.viewport(3) - viewer.current_mouse_y;
    if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.core.model,
        viewer.core.proj, viewer.core.viewport, this->Q, this->F, fid, bc))
    {
        std::cout << fid << " - clicked on face #\n";
        int prevId = 0;
        int nextId = 0;
        for (int i = 0; i < config->numRods(); i++)
        {
            nextId += config->rods[i]->numVertices() * 8;
            if (fid < nextId && fid > prevId)
            {
                config->rods[i]->visible = !config->rods[i]->visible;
                break;
            }
            prevId = nextId;
        }

        return true;
    }
    return false;
}

void RodsHook::initSimulation()
{
    iter = 0;
    forceResidual = 0;

    if (config)
        delete config;

    config = readRod(loadName.c_str());
    if (!config)
        exit(-1);
    centerScene();
   
    createVisualizationMesh();
    dirty = true;
}

void RodsHook::createVisualizationMesh()
{
    config->createVisualizationMesh(Q, F);
    exportWeave();
}

void RodsHook::loadTargetMesh()
{
    if (!igl::readOBJ(targetMeshName, targetV, targetF))
    {
        std::cerr << "Warning: couldn't load target mesh " << targetMeshName << std::endl;
        targetV.resize(0, 3);
        targetF.resize(0, 3);
    }
    updateRenderGeometry();
}

double lineSearch(RodConfig &config, 
    const Eigen::VectorXd &update, 
    double angleWeight, 
    bool allowSliding,
    Eigen::MatrixXd *anchorPoints,
    Eigen::MatrixXd *anchorNormals)
{
    double t = 1.0;
    double c1 = 0.1;
    double c2 = 0.9;
    double alpha = 0;
    double infinity = 1e6;
    double beta = infinity;

    Eigen::VectorXd r;
    Eigen::SparseMatrix<double> J;
    rAndJ(config, r, &J, angleWeight, allowSliding, anchorPoints, anchorNormals);  

    Eigen::VectorXd dE;
    Eigen::VectorXd newdE;
    std::vector<RodState> start;
    std::vector<Constraint> startC;
    for (int i = 0; i < config.numRods(); i++)
    {
        start.push_back(config.rods[i]->curState);
    }
    startC = config.constraints;
    double orig = 0.5 * r.transpose() * r;
    dE = J.transpose() * r;
    double deriv = -dE.dot(update);
    assert(deriv < 0);

    std::cout << "Starting line search, original energy " << orig << ", descent magnitude " << deriv << std::endl;

    while (true)
    {
        int nrods = config.numRods();
        int nconstraints = config.numConstraints();
        int dofoffset = 0;
        for (int rod = 0; rod < nrods; rod++)
        {
            int nverts = config.rods[rod]->numVertices();
            int nsegs = config.rods[rod]->numSegments();
            for (int i = 0; i < nsegs; i++)
            {
                Eigen::Vector3d oldv1 = start[rod].centerline.row(i);
                Eigen::Vector3d oldv2 = start[rod].centerline.row((i + 1) % nverts);
                Eigen::Vector3d v1 = oldv1 - t * update.segment<3>(dofoffset + 3 * i);
                Eigen::Vector3d v2 = oldv2 - t * update.segment<3>(dofoffset + 3 * ((i + 1) % nverts));

                config.rods[rod]->curState.directors.row(i) = parallelTransport(start[rod].directors.row(i), oldv2 - oldv1, v2 - v1);
            }
            for (int i = 0; i < nverts; i++)
                config.rods[rod]->curState.centerline.row(i) = start[rod].centerline.row(i) - t * update.segment<3>(dofoffset + 3 * i).transpose();
            for (int i = 0; i < nsegs; i++)
            {
                config.rods[rod]->curState.thetas[i] = start[rod].thetas[i] - t * update[dofoffset + 3 * nverts + i];
            }
            dofoffset += 3 * nverts + nsegs;
        }
        for (int i = 0; i < nconstraints; i++)
        {
            config.constraints[i].bary1 = startC[i].bary1 - t * update[dofoffset + 2 * i];           
            config.constraints[i].bary2 = startC[i].bary2 - t * update[dofoffset + 2 * i + 1];            
        }

        rAndJ(config, r, &J, angleWeight, allowSliding, anchorPoints, anchorNormals);

        double newenergy = 0.5 * r.transpose() * r;
        newdE = J.transpose() * r;

        std::cout << "Trying t = " << t << ", energy now " << newenergy << std::endl;

        if (std::isnan(newenergy) || newenergy > orig + t*deriv*c1)
        {
            beta = t;
            t = 0.5*(alpha + beta);
        }
        else if (-newdE.dot(update) < c2*deriv)
        {
            alpha = t;
            if (beta == infinity)
            {
                t = 2 * alpha;
            }
            else
            {
                t = 0.5*(alpha + beta);
            }

            if (beta - alpha < 1e-8)
            {
                return newdE.squaredNorm();
            }
        }
        else
        {
            return newdE.squaredNorm();
        }
    }
}

void RodsHook::centerScene()
{
    Eigen::Vector3d centroid(0, 0, 0);
    Eigen::Vector3d origcentroid(0, 0, 0);
    int nverts = 0;
    for (int i = 0; i < config->numRods(); i++)
    {
        for (int j = 0; j < config->rods[i]->numVertices(); j++)
        {
            centroid += config->rods[i]->curState.centerline.row(j).transpose();
            origcentroid += config->rods[i]->startState.centerline.row(j).transpose();
            nverts++;
        }
    }
    centroid /= nverts;
    origcentroid /= nverts;
    Eigen::MatrixXd A(3, nverts);
    Eigen::MatrixXd B(3, nverts);
    int idx = 0;
    for (int i = 0; i < config->numRods(); i++)
    {
        for (int j = 0; j < config->rods[i]->numVertices(); j++)
        {
            config->rods[i]->curState.centerline.row(j) -= centroid;
            A.col(idx) = config->rods[i]->curState.centerline.row(j).transpose();
            B.col(idx) = config->rods[i]->startState.centerline.row(j).transpose() - origcentroid;
            idx++;
        }
    }
    Eigen::Matrix3d M = B*A.transpose();
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d R = svd.matrixU() * svd.matrixV().transpose();
    if (R.determinant() < 0)
    {
        Eigen::Matrix3d sigma;
        sigma.setIdentity();
        sigma(2, 2) = -1;
        R = svd.matrixU() * sigma * svd.matrixV().transpose();
    }
    for (int i = 0; i < config->numRods(); i++)
    {
        for (int j = 0; j < config->rods[i]->numVertices(); j++)
        {
            config->rods[i]->curState.centerline.row(j) = config->rods[i]->curState.centerline.row(j) * R.transpose();            
        }
        for (int j = 0; j < config->rods[i]->numSegments(); j++)
        {
            config->rods[i]->curState.directors.row(j) = config->rods[i]->curState.directors.row(j) * R.transpose();
        }
    }
}

void RodsHook::slideConstraints()
{
    int nconstraints = config->numConstraints();
    std::vector<int> todelete;
    for (int i = 0; i < nconstraints; i++)
    {
        Constraint &c = config->constraints[i];
        if (c.bary1 < 0)
        {
            if (config->rods[c.rod1]->isClosed() || c.seg1 > 0)
            {
                c.seg1 = (c.seg1 + config->rods[c.rod1]->numSegments() - 1) % config->rods[c.rod1]->numSegments();
                c.bary1 = 1.0;
            }
            else
            {
                todelete.push_back(i);
                continue;
            }
        }
        else if (c.bary1 > 1.0)
        {
            if (config->rods[c.rod1]->isClosed() || c.seg1 < config->rods[c.rod1]->numSegments()-1)
            {
                c.seg1 = (c.seg1 + 1) % config->rods[c.rod1]->numSegments();
                c.bary1 = 0.0;
            }
            else
            {
                todelete.push_back(i);
                continue;
            }
        }
        if (c.bary2 < 0)
        {
            if (config->rods[c.rod2]->isClosed() || c.seg2 > 0)
            {
                c.seg2 = (c.seg2 + config->rods[c.rod2]->numSegments() - 1) % config->rods[c.rod2]->numSegments();
                c.bary2 = 1.0;
            }
            else
            {
                todelete.push_back(i);
                continue;
            }
        }
        else if (c.bary2 > 1.0)
        {
            if (config->rods[c.rod2]->isClosed() || c.seg2 < config->rods[c.rod2]->numSegments()-1)
            {
                c.seg2 = (c.seg2 + 1) % config->rods[c.rod2]->numSegments();
                c.bary2 = 0.0;
            }
            else
            {
                todelete.push_back(i);
                continue;
            }
        }
    }

    for (std::vector<int>::reverse_iterator it = todelete.rbegin(); it != todelete.rend(); ++it)
    {
        config->constraints.erase(config->constraints.begin() + *it);
    }
}

void RodsHook::findAnchorPoints(Eigen::MatrixXd &anchorPoints, Eigen::MatrixXd &anchorNormals)
{
    int nrods = config->numRods();
    int nverts = 0;
    for (int i = 0; i < nrods; i++)
        nverts += config->rods[i]->numVertices();
    anchorPoints.resize(nverts, 3);
    anchorNormals.resize(nverts, 3);
    int idx = 0;
    for (int i = 0; i < nrods; i++)
    {
        Eigen::VectorXd sqrD;
        Eigen::VectorXi I;
        Eigen::MatrixXd C;
        igl::point_mesh_squared_distance(
            config->rods[i]->curState.centerline,
            targetV,
            targetF,
            sqrD,
            I,
            C
        );
        int nverts = config->rods[i]->numVertices();
        for (int j = 0; j < nverts; j++)
        {
            anchorPoints.row(idx) = C.row(j);
            int face = I[j];
            Eigen::Vector3d e1 = targetV.row(targetF(face, 1)).transpose() - targetV.row(targetF(face, 0)).transpose();
            Eigen::Vector3d e2 = targetV.row(targetF(face, 2)).transpose() - targetV.row(targetF(face, 0)).transpose();
            Eigen::Vector3d n = e1.cross(e2);
            anchorNormals.row(idx) = n / n.norm();
            idx++;
        }
    }
}

void RodsHook::testFiniteDifferences()
{
    Eigen::VectorXd r;
    Eigen::SparseMatrix<double> Jr;
    rAndJ(*config,
        r,
        &Jr,
        angleWeight,
        allowSliding,
        NULL,
        NULL);
    
    Eigen::VectorXd newr;
    double eps = 1e-6;
    int ndofs = Jr.cols();
    Eigen::VectorXd pert(ndofs);
    pert.setRandom();

    int idx2 = 0;
    for (int i = 0; i < config->numRods(); i++)
    {
        idx2 += 3*config->rods[i]->numVertices();
        idx2 += config->rods[i]->numSegments();
    }
    for (int i = 0; i < config->numConstraints(); i++)
    {
        if (allowSliding)
        {
            idx2+=2;
        }
        else
        {
            pert[idx2++] = 0;
            pert[idx2++] = 0;
        }
    }

    Eigen::VectorXd exact = Jr * pert;
    int idx = 0;
    for (int i = 0; i < config->numRods(); i++)
    {
        int nverts = config->rods[i]->numVertices();
        Eigen::MatrixXd newv = config->rods[i]->curState.centerline;

        for (int j = 0; j < config->rods[i]->numVertices(); j++)
        {
            for (int k = 0; k < 3; k++)
            {
                newv(j, k) += eps*pert[idx++];
            }
        }
        for (int j = 0; j < config->rods[i]->numSegments(); j++)
        {
            Eigen::Vector3d newseg = newv.row((j + 1) % config->rods[i]->numVertices()) - newv.row(j);
            Eigen::Vector3d oldseg = config->rods[i]->curState.centerline.row((j + 1) % config->rods[i]->numVertices()) - config->rods[i]->curState.centerline.row(j);
            config->rods[i]->curState.directors.row(j) = parallelTransport(config->rods[i]->curState.directors.row(j), oldseg, newseg).transpose();
        }
        config->rods[i]->curState.centerline = newv;

        for (int j = 0; j < config->rods[i]->numSegments(); j++)
        {
            config->rods[i]->curState.thetas[j] += eps*pert[idx++];
        }
    }
    for (int i = 0; i < config->numConstraints(); i++)
    {
        config->constraints[i].bary1 += eps * pert[idx++];
        config->constraints[i].bary2 += eps * pert[idx++];
    }
    rAndJ(*config,
        newr,
        NULL,
        angleWeight,
        allowSliding,
        NULL,
        NULL);
    std::ofstream ofs("findiff.txt");
    ofs << "exact:" << std::endl;
    ofs << exact << std::endl << std::endl;
    ofs << "findiff: " << std::endl;
    ofs << (newr - r) / eps << std::endl;
    exit(0);
}

bool RodsHook::simulateOneStep()
{
    //testFiniteDifferences();
    Eigen::VectorXd r;
    Eigen::SparseMatrix<double> Jr;
    Eigen::MatrixXd anchorPoints;
    Eigen::MatrixXd anchorNormals;
    bool useanchors = false;
    if (stickToMesh && targetV.rows() > 0)
    {
        findAnchorPoints(anchorPoints, anchorNormals);
        useanchors = true;
    }
    rAndJ(*config, 
        r, 
        &Jr, 
        angleWeight, 
        allowSliding, 
        useanchors ? &anchorPoints : NULL, 
        useanchors ? &anchorNormals : NULL);
    std::cout << "Orig energy: " << 0.5 * r.transpose() * r << std::endl;
    Eigen::SparseMatrix<double> mat = Jr.transpose() * Jr;
    Eigen::SparseMatrix<double> I(mat.rows(), mat.cols());
    I.setIdentity();    
    double Treg = 1e-6;
    mat += Treg*I;

    Eigen::VectorXd rhs = Jr.transpose() * r;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver(mat);
    Eigen::VectorXd delta = solver.solve(rhs);
    if (solver.info() != Eigen::Success)
        exit(-1);
    std::cout << "Solver residual: " << (mat*delta - rhs).norm() << std::endl;

    forceResidual = lineSearch(*config, 
        delta, 
        angleWeight, 
        allowSliding, 
        useanchors ? &anchorPoints : NULL, 
        useanchors ? &anchorNormals : NULL);

    slideConstraints();
    centerScene();

    iter++;

    createVisualizationMesh();    

    return false;
}

void RodsHook::saveRods()
{
    config->saveRodGeometry(savePrefix);
}


using namespace svg;
void RodsHook::exportWeave()
{
    double strip_width = 45.;
    double strip_space = 10.;
    double strip_stretch = 600.;
    double label_spacing = 90.;

    enum Defaults { Transparent = -1, Aqua, Black, Blue, Brown, Cyan, Fuchsia,
        Green, Lime, Magenta, Orange, Purple, Red, Silver, White, Yellow };

    Color clist[] = {Color::Blue, Color::Red, Color::Yellow, Color::Lime, Color::Orange, Color::Purple, Color::Cyan, Color::Black};
    int colorlen = 7;

    double maxlen = 0;
    double maxseg = 0;
    for (int i = 0; i < config->numRods(); i++) 
    {  
        if (config->rods[i]->numSegments() > maxlen) 
            maxlen = config->rods[i]->numSegments();

        for (int j = 0; j < config->rods[i]->numSegments() - 1; j++)
        {
            Eigen::Vector3d seg = config->rods[i]->curState.centerline.row(j) - config->rods[i]->curState.centerline.row(j + 1);
            if ( seg.norm() > maxseg )
                maxseg = seg.norm(); // Can opt w/ norm squared if slow...
        }
    } 

    Dimensions dimensions(strip_stretch * maxlen * maxseg, (config->numRods() + 1) * (strip_width + strip_space));
    Document doc("my_svg.svg", Layout(dimensions, Layout::BottomLeft));

    Eigen::MatrixXi collisions = Eigen::MatrixXi::Constant(config->numRods(), maxlen, 0);
    Eigen::MatrixXi collisions_strip_match = Eigen::MatrixXi::Constant(config->numRods(), maxlen, 0);
    Eigen::MatrixXi collisions_circ_match = Eigen::MatrixXi::Constant(config->numRods(), maxlen, 0);
    Eigen::MatrixXi self_intersect = Eigen::MatrixXi::Constant(config->numRods(), maxlen, 0);
    Eigen::MatrixXi angles = Eigen::MatrixXi::Constant(config->numRods(), maxlen, 0.);
    std::vector<Eigen::Matrix3d> rotations;
    rotations.push_back( Eigen::MatrixXd::Identity(3,3) );
    int match_iter = 0;
    for (int i = 0; i < config->numConstraints(); i++)
    {
        Constraint c = config->constraints[i];
        collisions(c.rod1, c.seg1) = (1 + c.rod2) * c.assignment;
        collisions(c.rod2, c.seg2) = (1 + c.rod1) * c.assignment * -1;
        collisions_strip_match(c.rod1, c.seg1) = c.rod2 % colorlen;//match_iter;
        collisions_strip_match(c.rod2, c.seg2) = c.rod1 % colorlen;//match_iter;
        collisions_circ_match(c.rod1, c.seg1) = match_iter;
        collisions_circ_match(c.rod2, c.seg2) = match_iter;
        config->constraints[i].color = match_iter;
        if (c.rod1 == c.rod2)
        {
            self_intersect(c.rod1, c.seg1) = 1;
            self_intersect(c.rod2, c.seg2) = 1;
        } 

        if ( !config->rods[c.rod1]->visible || !config->rods[c.rod2]->visible )
        {
            self_intersect(c.rod1, c.seg1) = -1;
            self_intersect(c.rod2, c.seg2) = -1;
        }
        match_iter = (match_iter + 3) % (colorlen * 3);
        Eigen::Vector3d r1 = config->rods[c.rod1]->curState.centerline.row(c.seg1) - 
                             config->rods[c.rod1]->curState.centerline.row(c.seg1 + 1); 
        Eigen::Vector3d r2 = config->rods[c.rod2]->curState.centerline.row(c.seg2) - 
                             config->rods[c.rod2]->curState.centerline.row(c.seg2 + 1);
        r1.normalize();
        r2.normalize();


        Eigen::Vector3d d1 = config->rods[c.rod1]->curState.directors.row(c.seg1);
        Eigen::Vector3d d2 = r1.cross(d1);
        double theta = config->rods[c.rod1]->curState.thetas[c.seg1];
        Eigen::Vector3d n_r1 = d1*cos(theta) + d2*sin(theta);

        Eigen::Vector3d n = r1.cross(r2).normalized();
        if ( n.dot(n_r1) < .01 )
            n *= -1;

        Eigen::Matrix3d f1;
        f1.col(0) = r1;
        f1.col(1) = n.cross(r1);
        f1.col(2) = n;
        
        Eigen::Matrix3d f2;
        f2.col(0) = r2;
        f2.col(1) = n.cross(r2);
        f2.col(2) = n;

        // Eigen::Matrix3d f3 = Eigen::Matrix3d::Identity();
        // f3.col(2) = n;
        // Eigen::Matrix3d toPlane = f3.inverse();
        // f1 = toPlane * f1;
        // f2 = toPlane * f2;

  //      std::cout << f1 << "\n";

        if ( c.assignment > 0 )
        {
            angles(c.rod1, c.seg1) = 0;
            angles(c.rod2, c.seg2) = rotations.size();
            rotations.push_back(  f1.inverse() * f2 );
        }
        else 
        {
            angles(c.rod2, c.seg2) = 0;
            angles(c.rod1, c.seg1) = rotations.size();
            rotations.push_back(  f2.inverse() * f1 );
        }
    } 


    // Draw lines
    for (int i = 0; i < config->numRods(); i++) 
    {
        Rod *r = config->rods[i];
    //    Polyline pl_center(Fill(Color::Transparent), Stroke(strip_width - 6., clist[i % colorlen]));
        svg::Polyline pl_l(Fill(Color::Transparent), Stroke(1., Color::Black));
        svg::Polyline pl_r(Stroke(.5, Color::Black));
        double startpoint = 0.;
        double endpoint;
        double x_shift = i * (strip_width + strip_space);
        pl_l << Point(startpoint, x_shift);
  //      pl_center << Point(startpoint, x_shift);
        for (int j = 0; j < r->numSegments(); j++)
        { 
            Eigen::Vector3d seg = r->curState.centerline.row(j) - r->curState.centerline.row(j + 1);
            endpoint = seg.norm() * strip_stretch + startpoint;
            pl_l << Point(endpoint, x_shift);
            svg::Color c( r->colors(j, 0 ) * 255, r->colors(j, 1 ) * 255, r->colors(j, 2 ) * 255);
            doc << Line( Point(startpoint, x_shift + strip_width / 2.), Point(endpoint, x_shift + strip_width / 2.), Stroke(strip_width - 6., c) );
   //         pl_center << Point(endpoint, x_shift);
            startpoint = endpoint; 
        }

        pl_r = pl_l;
        pl_r.offset( Point (0, strip_width) );
   //     pl_center.offset( Point (0, strip_width / 2.) );
        doc << pl_r << pl_l;
    }

    // Draw crossings and labels 
    for (int i = 0; i < config->numRods(); i++) 
    {
        Rod *r = config->rods[i];
        double startpoint = 0.;
        double endpoint;
        double x_shift = i * (strip_width + strip_space) + strip_space;
        double lasttext = -200.;

        for (int j = 0; j < r->numSegments() - 1; j++)
        { 
            Eigen::Vector3d seg = r->curState.centerline.row(j) - r->curState.centerline.row(j + 1);
            endpoint = seg.norm() * strip_stretch + startpoint;
         //   pl_l << Point(endpoint, x_shift);
            startpoint = endpoint;
            if ( collisions(i,j) != 0) 
            { 
                double y_center = i * (strip_width + strip_space) + .5 * strip_width;
                if ( self_intersect(i,j) == -1 ) {}
                else if ( self_intersect(i,j) > 0 )
                {
                    doc << Circle( Point(endpoint, y_center), strip_width / 2., Fill(Color::Fuchsia), Stroke(3., Color::Fuchsia));
                }
                else 
                {
                    Color circ_col = clist[collisions_circ_match(i,j) % colorlen];
                    if ( (collisions_circ_match(i,j) % colorlen) == (i % colorlen) ) 
                        Color circ_col = Color::White;
                    doc << Circle( Point(endpoint, y_center), strip_width / 3., Fill(circ_col), Stroke(3., circ_col));
                }
                
                svg::Polyline mark_crossing(Fill(Color::Transparent), Stroke(3., Color::Black));

                int mark_cidx = collisions_strip_match(i,j);
                if (mark_cidx == (i % colorlen))
                    mark_cidx = colorlen;
                if (collisions(i,j) < 0)
                {
                    mark_crossing = svg::Polyline(Fill(Color::Transparent), Stroke(7., clist[mark_cidx]));
                }
                else 
                {
                    mark_crossing = svg::Polyline(Fill(Color::Transparent), Stroke(2., clist[mark_cidx]));
                }

                Eigen::Matrix3d orient = rotations[angles(i,j)] * Eigen::MatrixXd::Identity(3,3) * strip_width / 2.;

                Eigen::Vector2d pos_x = Eigen::Vector2d(orient(0, 0), orient(0, 1));
                Eigen::Vector2d pos_y = Eigen::Vector2d(orient(1, 0), orient(1, 1)); // 1. / cos(angles(i, j)) *

        //        std::cout << orient.col(2) << "\n";

                mark_crossing << Point( endpoint + pos_x(0) + pos_y(0) * (1.),   pos_x(1) + pos_y(1) * (1.) + y_center )
                              << Point( endpoint + pos_x(0) - pos_y(0) * (1.),   pos_x(1) - pos_y(1) * (1.) + y_center )
                              << Point( endpoint - pos_x(0) - pos_y(0) * (1.), - pos_x(1) - pos_y(1) * (1.) + y_center )
                              << Point( endpoint - pos_x(0) + pos_y(0) * (1.), - pos_x(1) + pos_y(1) * (1.) + y_center );


                // mark_crossing << Point( endpoint + pos_x(0),   pos_x(1) + y_center )
                //               << Point( endpoint - pos_y(0) - pos_x(0), - pos_y(1)- pos_x(1) + y_center )
                //               << Point( endpoint - pos_x(0), - pos_x(1) + y_center );



/*                mark_crossing << Point( endpoint + pos_x(0) + pos_y(0) * (1 + (tan(angles(i, j)) * strip_width / 2.)),   pos_x(1) + pos_y(1) * (1 + (tan(angles(i, j)) * strip_width / 2.)) + y_center )
                              << Point( endpoint + pos_x(0) - pos_y(0) * (1 - (tan(angles(i, j)) * strip_width / 2.)),   pos_x(1) - pos_y(1) * (1 - (tan(angles(i, j)) * strip_width / 2.)) + y_center )
                              << Point( endpoint - pos_x(0) - pos_y(0) * (1 + (tan(angles(i, j)) * strip_width / 2.)), - pos_x(1) - pos_y(1) * (1 + (tan(angles(i, j)) * strip_width / 2.)) + y_center )
                              << Point( endpoint - pos_x(0) + pos_y(0) * (1 - (tan(angles(i, j)) * strip_width / 2.)), - pos_x(1) + pos_y(1) * (1 - (tan(angles(i, j)) * strip_width / 2.)) + y_center );*/
                if ( self_intersect(i,j) == -1 ) {}
                else {
                    doc << mark_crossing;

                    char shift = collisions_circ_match(i,j) + 'a';
                    doc << svg::Text(Point(endpoint - 10, x_shift  - strip_space / 2), std::to_string(abs(collisions(i,j))) + shift, Color::Black, Font(10, "Verdana"));
                    doc << svg::Text(Point(endpoint + 10, x_shift - strip_space / 2 + strip_width / 2.), std::to_string(abs(collisions(i,j))) + shift, Color::White, Font(10, "Verdana"));
                    lasttext = endpoint;
                }
            } 
            if (lasttext < endpoint - label_spacing) 
            {
     //           std::cout << endpoint << "\n";
                doc << svg::Text(Point(endpoint, x_shift - strip_space / 2 + strip_width / 2. - 6), std::to_string(i + 1), Color::White, Font(14, "Verdana"));
                lasttext = endpoint;
            }
        }


    } 
    
    doc.save();

}

void RodsHook::showConstraints()
{
    int nconstraints = config->constraints.size();
    constraintPoints.resize(2 *nconstraints,3);
    constraintEdges.resize(nconstraints,2);
    constraintColors.resize(nconstraints, 3);
    constraintColors.col(0).setConstant(1);
    constraintColors.col(1).setConstant(1);
    constraintColors.col(2).setConstant(0);
    for (int i = 0; i < nconstraints; i++)
    {
        Constraint &c = config->constraints[i];
        Eigen::Vector3d v1 = config->rods[c.rod1]->curState.centerline.row(c.seg1);
        Eigen::Vector3d v2 = config->rods[c.rod1]->curState.centerline.row((c.seg1 + 1) % config->rods[c.rod1]->numVertices());
        Eigen::Vector3d w1 = config->rods[c.rod2]->curState.centerline.row(c.seg2);
        Eigen::Vector3d w2 = config->rods[c.rod2]->curState.centerline.row((c.seg2 + 1) % config->rods[c.rod2]->numVertices());
        constraintPoints.row(2 * i) = (1.0 - c.bary1)*v1 + c.bary1 * v2;
        constraintPoints.row(2 * i + 1) = (1.0 - c.bary2)*w1 + c.bary2 * w2;
        constraintEdges(i, 0) = 2 * i;
        constraintEdges(i, 1) = 2 * i + 1;
    }
}

void RodsHook::linearSubdivision()
{
    int nrods = config->numRods();
    for(int i=0; i<nrods; i++)
    {
        int nverts = config->rods[i]->numVertices();        
        RodState newstate;
        newstate.centerline.resize(2*nverts - 1,3);
        for(int j=0; j<nverts; j++)
        {
            newstate.centerline.row(2*j) = config->rods[i]->startState.centerline.row(j);
            if(j != nverts-1)
            {
                newstate.centerline.row(2*j+1) = 0.5 * (config->rods[i]->startState.centerline.row(j) + config->rods[i]->startState.centerline.row(j+1) );
            }
        }
        newstate.centerlineVel.resize(2*nverts - 1,3);
        newstate.centerlineVel.setZero();
        int nsegs = config->rods[i]->numSegments();
        newstate.directors.resize(2*nsegs, 3);
        for(int j=0; j<nsegs; j++)
        {
            newstate.directors.row(2*j) = config->rods[i]->startState.directors.row(j);
            newstate.directors.row(2*j+1) = config->rods[i]->startState.directors.row(j);
        }
        newstate.thetas.resize(2*nsegs);
        newstate.thetas.setZero();
        newstate.directorAngVel.resize(2*nsegs);
        newstate.directorAngVel.setZero();
        config->rods[i]->startState = newstate;
        config->rods[i]->curState = newstate;
        Eigen::VectorXd newwidths(2*nsegs);
        for(int j=0; j<nsegs; j++)
        {
            newwidths[2*j] = config->rods[i]->widths[j];
            newwidths[2*j+1] = config->rods[i]->widths[j];
        }
        config->rods[i]->widths = newwidths;
        config->rods[i]->initializeRestQuantities();
    }

    std::vector<Constraint> newconstraints;
    int nconstraints = config->constraints.size();
    for(int i=0; i<nconstraints; i++)
    {
        Constraint c = config->constraints[i];
        Constraint newc;
        newc.rod1 = c.rod1;
        newc.rod2 = c.rod2;
        newc.stiffness = c.stiffness;
        newc.assignment = c.assignment;
        newc.color = 0;
        if(c.bary1 < 0.5)
        {
            newc.seg1 = 2*c.seg1;
            newc.bary1 = c.bary1 * 2.0;
        }
        else
        {
            newc.seg1 = 2*c.seg1 + 1;
            newc.bary1 = c.bary1 * 2.0 - 1.0;
        }
        if(c.bary2 < 0.5)
        {
            newc.seg2 = 2*c.seg2;
            newc.bary2 = c.bary2 * 2.0;
        }
        else
        {
            newc.seg2 = 2*c.seg2 + 1;
            newc.bary2 = c.bary2 * 2.0 - 1.0;
        }
        newconstraints.push_back(newc);
    }
    config->constraints = newconstraints;
    config->initWeave();
    centerScene();
    createVisualizationMesh();
    updateRenderGeometry();
}

void RodsHook::setWidths()
{
    for(int i=0; i<config->numRods(); i++)
    {
        for(int j=0; j<config->rods[i]->numSegments(); j++)
        {
            config->rods[i]->widths[j] = newWidth;
        }
    }
    createVisualizationMesh();
    updateRenderGeometry();  
}
