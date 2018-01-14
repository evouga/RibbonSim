#include "PhysicsHook.h"
#include "RodEnergy.h"

class RodsHook : public PhysicsHook
{
public:
    RodsHook() : PhysicsHook(), iter(0), forceResidual(0.0), dirty(true), config(NULL) {
    }

    virtual void initGUI(igl::viewer::Viewer &viewer);

    virtual void initSimulation();
    
    virtual void updateRenderGeometry()
    {
        if (Q.rows() != renderQ.rows() || F.rows() != renderF.rows())
            dirty = true;
        renderQ = Q;
        renderF = F;        
    }

    virtual bool simulateOneStep();
    
    virtual void renderRenderGeometry(igl::viewer::Viewer &viewer)
    {
        if (dirty)
        {
            viewer.data.clear();
            dirty = false;
        }
        
        if (showTarget) 
	{
            viewer.data.clear();	    
	    Eigen::MatrixXd V(renderQ.rows() + config->V_mesh.rows(), renderQ.cols());
	    V<<renderQ, config->V_mesh;
	    Eigen::MatrixXi F(renderF.rows() + config->F_mesh.rows(), renderF.cols());
	    F<<renderF,(config->F_mesh.array()+renderF.rows());

	    viewer.data.set_mesh(V, F);
	    Eigen::MatrixXd C(F.rows(),3);
	    C << Eigen::RowVector3d(0.2,0.3,0.8).replicate(renderF.rows(),1),
		 Eigen::RowVector3d(1.0,0.7,0.2).replicate(config->F_mesh.rows(),1);  
	    viewer.data.set_colors(C);
	}
        else 
	{
	    viewer.data.clear();
            viewer.data.set_mesh(renderQ, renderF);
	}


	if(forceEdges.rows() > 0)
            viewer.data.set_edges(forcePoints, forceEdges, forceColors);
    }

    void saveRods();
    void refreshLoadBuffers(); 
 
private:    
    void createVisualizationMesh();
    void showForces(int rod, const Eigen::VectorXd &dE);


    char loadBuffer [100];
    char meshBuffer [100];
    std::string loadName;

    bool showTarget = true;

    double dt;
    double damp;

    double iter;
    double forceResidual;
    
    RodConfig *config;

    std::string savePrefix;

    // for visualization
    Eigen::MatrixXd Q;
    Eigen::MatrixXi F;
    Eigen::MatrixXd forcePoints;
    Eigen::MatrixXi forceEdges;
    Eigen::MatrixXd forceColors;

    int segsPerEdge;
    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
    bool dirty;
};
