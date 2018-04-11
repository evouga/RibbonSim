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
    RodsHook() : PhysicsHook(), iter(0), forceResidual(0.0), angleWeight(1e3), newWidth(0.01), dirty(true), config(NULL) {
    }

    virtual void initGUI(igl::viewer::Viewer &viewer);

    virtual void initSimulation();
    
    virtual void updateRenderGeometry()
    {
        if (Q.rows() != renderQ.rows() || F.rows() != renderF.rows())
            dirty = true;
        renderQ = Q;
        renderF = F;        
        showConstraints();
    }

    virtual bool simulateOneStep();
    
    virtual void renderRenderGeometry(igl::viewer::Viewer &viewer)
    {
        if (dirty)
        {
            viewer.data.clear();
            dirty = false;
        }
        viewer.data.set_mesh(renderQ, renderF);

        int faces = renderF.rows();
        faceColors.resize(faces, 4);
        int pos = 0;

        faceColors.setConstant(0.7);

        for (int i = 0; i < config->numRods(); i++)
        {
            double alpha = 1.;
            for (int j = 0; j < config->rods[i]->numSegments(); j++)
            {
                Eigen::Vector3d col = config->rods[i]->colors.row(j);
                for ( int f = 0; f < 8; f++)
                {
                    faceColors.row(pos) = Eigen::Vector4d(col(0), col(1), col(2), alpha);
                    pos++;
                }
            }
        }


        int numConst = config->numConstraints();
        Eigen::MatrixXd P(numConst, 3);
        Eigen::MatrixXd C(numConst, 3);

        int acc_hack = 0;
        int num_colors = 7;

        for (int i = 0; i < numConst; i++)
        {
            Constraint c = config->constraints[i];
            double* col = dot_colors[(c.color % num_colors)];
            acc_hack += c.color;
            P.row(i) = config->rods[c.rod1]->curState.centerline.row(c.seg1);
            C.row(i) = Eigen::Vector3d(col[0], col[1], col[2]);
        }
        if (acc_hack > 3)
            viewer.data.set_points(P, C);

        viewer.core.lighting_factor = 0.;
        viewer.data.set_colors(faceColors);




        if(forceEdges.rows() > 0)
            viewer.data.set_edges(forcePoints, forceEdges, forceColors);
    }

    void saveRods();
    void setWidths();

private:    
    void createVisualizationMesh();
    void showConstraints();
    void linearSubdivision();
    void exportWeave();

    std::string loadName;

    double iter;
    double forceResidual;
    double angleWeight;
    double newWidth;
    
    RodConfig *config;

    std::string savePrefix;

    // for visualization
    Eigen::MatrixXd Q;
    Eigen::MatrixXi F;
    Eigen::MatrixXd forcePoints;
    Eigen::MatrixXi forceEdges;
    Eigen::MatrixXd forceColors;

    Eigen::MatrixXd faceColors;

    int segsPerEdge;
    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
    bool dirty;
};
