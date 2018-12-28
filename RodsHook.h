#include "PhysicsHook.h"
#include "RodEnergy.h"

struct SceneStats
{
    int numRods;
    int numCrossings;
    double totalLength;
    double dimensions[3];

    double meanWidth;
    double meanThickness;
    double meanModulus;
    double meanDensity;
};

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
        recomputeStats();
    }

    virtual bool simulateOneStep();
    
    virtual void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer);

    void saveRods();
    void setWidths();

private:    
    void createVisualizationMesh();
    void showConstraints();
    void deleteInvisibleRods();
    void linearSubdivision();
    void trimLooseEnds();
    void exportSomeRods(const char*filename, int firstRod, int numRods);
    void exportWeave();    
    void centerScene();
    void slideConstraints();
    void loadTargetMesh();
    void findAnchorPoints(Eigen::MatrixXd &anchorPoints, Eigen::MatrixXd &anchorNormals);
    void testFiniteDifferences();
    void saveConfig();
    void hideLongRods();
    void fitFloorHeight();
    void rescaleRods(double factor);
    void recomputeStats();

    std::string loadName;

    int iter;
    double forceResidual;
    float constraintWeight;
    bool enableGravity;
    Eigen::Vector3d gravityDir;
    double floorHeight;
    double floorWeight;
    bool allowSliding;
    float newWidth;
    float expLenScale;
    bool limitRenderLen;
    float maxRenderLen;
    double rescaleFactor;
    SceneStats stats;

    RodConfig *config;

    std::string savePrefix;
    bool visualizeConstraints;
    bool visualizeTargetMesh;

    int rodsPerSVG;
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
    Eigen::MatrixXd floorQ;
    Eigen::MatrixXi floorF;
    Eigen::MatrixXd floorColors;
    bool dirty;
};
