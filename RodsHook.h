#include "PhysicsHook.h"
#include "RodEnergy.h"

static double dot_colors[7][3] = { { .001, .001, .95 },
                                          { 0.95, 0.001, 0.001 },
                                          { 0.95, 0.95, 0.001 },
                                          { 0.01, 0.95, 0.01 },
                                          { 0.95, 0.7, 0.01 },
                                          { 0.5, 0.01, 0.5 },
                                          { 0.001, 0.95, 0.95 } };

class RodsHook : public PhysicsHook
{
public:
    RodsHook();

    virtual void drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &viewer);

    virtual bool mouseClicked(igl::opengl::glfw::Viewer &viewer, int button);

    virtual void initSimulation();

    virtual void tick() {}
    
    virtual void updateRenderGeometry()
    {
        int numRenderVerts = Q.rows();
        if (visualizeTargetMesh)
            numRenderVerts += targetV.rows();
        int numRenderFaces = F.rows();
        if (visualizeTargetMesh)
            numRenderFaces += targetF.rows();
        if (numRenderVerts != renderQ.rows() || numRenderFaces != renderF.rows())
            dirty = true;

        renderQ.resize(numRenderVerts, 3);
        for (int i = 0; i < Q.rows(); i++)
        {
            renderQ.row(i) = Q.row(i);
        }
        if (visualizeTargetMesh)
        {
            int offset = Q.rows();
            for (int i = 0; i < targetV.rows(); i++)
            {
                renderQ.row(offset + i) = targetV.row(i);
            }
        }
        
        renderF.resize(numRenderFaces, 3);
        for (int i = 0; i < F.rows(); i++)
        {
            renderF.row(i) = F.row(i);
        }
        if (visualizeTargetMesh)
        {
            int offset = F.rows();
            int voffset = Q.rows();
            for (int i = 0; i < targetF.rows(); i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    renderF(offset + i, j) = targetF(i, j) + voffset;
                }
            }
        }
        
        showConstraints();
    }

    virtual bool simulateOneStep();
    
    virtual void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer)
    {
        if (dirty)
        {
            viewer.data().clear();
            dirty = false;
        }
        viewer.data().set_mesh(renderQ, renderF);

        int faces = renderF.rows();
        faceColors.resize(faces, 4);
        faceColors.col(0).setConstant(0.7);
        faceColors.col(1).setConstant(0.7);
        faceColors.col(2).setConstant(0.7);
        faceColors.col(3).setConstant(1.0);
        int pos = 0;

        double transp = 1.;

        for (int i = 0; i < config->numRods(); i++)
        {
            if (config->rods[i]->visible)
            {
                transp = 1.;
            }
            else
            {
                transp = -10.;
            }

            for (int j = 0; j < config->rods[i]->numSegments(); j++)
            {
                Eigen::Vector3d col = config->rods[i]->colors.row(j);
                for ( int f = 0; f < 8; f++)
                {
                    faceColors.row(pos) = Eigen::Vector4d(col(0), col(1), col(2), transp);
                    pos++;
                }
            }
        }


        int numConst = config->numConstraints();
        Eigen::MatrixXd P = Eigen::MatrixXd::Zero(numConst, 3);
        Eigen::MatrixXd C = Eigen::MatrixXd::Zero(numConst, 3);
        viewer.data().set_points(P, C);

        int num_colors = 7;
        bool show = false;

        for (int i = 0; i < numConst; i++)
        {
            Constraint c = config->constraints[i];
            double* col = dot_colors[(c.color % num_colors)];

            if (config->rods[c.rod1]->visible && config->rods[c.rod2]->visible)
            {
                show = true;
                P.row(i) = constraintPoints.row(2 * i);
                C.row(i) = Eigen::Vector3d(col[0], col[1], col[2]);
            }
        }
        if (visualizeConstraints)
            viewer.data().set_points(P, C);

        viewer.core.lighting_factor = 0.;
        viewer.data().set_colors(faceColors);

        if(constraintEdges.rows() > 0)
            viewer.data().set_edges(constraintPoints, constraintEdges, constraintColors);
    }

    void saveRods();
    void setWidths();

private:    
    void createVisualizationMesh();
    void showConstraints();
    void linearSubdivision();
    void exportWeave();    
    void centerScene();
    void slideConstraints();
    void loadTargetMesh();
    void findAnchorPoints(Eigen::MatrixXd &anchorPoints, Eigen::MatrixXd &anchorNormals);
    void testFiniteDifferences();

    std::string loadName;

    int iter;
    double forceResidual;
    float angleWeight;
    bool allowSliding;
    float newWidth;
    float expLenScale;

    RodConfig *config;

    std::string savePrefix;
    bool visualizeConstraints;
    bool visualizeTargetMesh;

    std::string targetMeshName;
    bool stickToMesh;

    Eigen::MatrixXd targetV;
    Eigen::MatrixXi targetF;

    // for visualization
    Eigen::MatrixXd Q;
    Eigen::MatrixXi F;
    Eigen::MatrixXd constraintPoints;
    Eigen::MatrixXi constraintEdges;
    Eigen::MatrixXd constraintColors;

    Eigen::MatrixXd faceColors;

    int segsPerEdge;
    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
    bool dirty;
};
