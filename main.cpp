#include <igl/viewer/Viewer.h>
#include <thread>
#include "PhysicsHook.h"
#include "RodsHook.h"

static PhysicsHook *hook = NULL;

void toggleSimulation()
{
    if (!hook)
        return;

    if (hook->isPaused())
        hook->run();
    else
        hook->pause();
}

void resetSimulation()
{
    if (!hook)
        return;

    hook->reset();
}

bool drawCallback(igl::viewer::Viewer &viewer)
{
    if (!hook)
        return false;

    hook->render(viewer);
    return false;
}

bool keyCallback(igl::viewer::Viewer& viewer, unsigned int key, int modifiers)
{
    if (key == ' ')
    {
        toggleSimulation();
        return true;
    }
    return false;
}

bool initGUI(igl::viewer::Viewer &viewer)
{
    viewer.ngui->window()->setVisible(false);
    viewer.ngui->addWindow(Eigen::Vector2i(10, 10), "Sim Settings");
    viewer.ngui->addButton("Run/Pause Sim", toggleSimulation);
    viewer.ngui->addButton("Reset Sim", resetSimulation);
    hook->initGUI(viewer);    
    viewer.screen->performLayout();
    hook->reset();
    return false;
}

int main(int argc, char *argv[])
{
  igl::viewer::Viewer viewer;

  hook = new RodsHook();

  viewer.data.set_face_based(true);
  viewer.core.is_animating = true;
  viewer.callback_key_pressed = keyCallback;
  viewer.callback_pre_draw = drawCallback;
  viewer.callback_init = initGUI;
  viewer.launch();
}
