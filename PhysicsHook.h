#ifndef PHYSICSHOOK_H
#define PHYSICSHOOK_H

#include <mutex>
#include <thread>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>

class PhysicsHook
{
public:
    PhysicsHook() : sim_thread(NULL), please_pause(false), please_die(false), running(false) 
    {
    }

    virtual ~PhysicsHook()
    {
        killSimThread();        
    }

    /*
    * Runs when the user redraws/interacts with the GUI.
    */
    virtual void drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu) = 0;

    /*
     * Runs once when the simulation is initialized, and again each time the user resets the simulation.
     */
    virtual void initSimulation() = 0;

    /*
    * Called every once in a while (and always before every simulation step) even if the simulation is paused.
    * Use this to update the visualization in response to user input etc.
    */
    virtual void tick() {};

    /*
     * Takes one simulation "step." You can do whatever you want here, but the granularity of your computation should 
     * be small enough that the user can view/pause/kill the simulation at interactive rates.
     * This method *must* be thread-safe with respect to renderRenderGeometry() (easiest is to not touch any rendering
     * data structures at all).
     */
    virtual bool simulateOneStep() = 0;

    /*
     * Update the rendering data structures here. This method will be called in alternation with simulateOneStep().
     * This method blocks rendering in the viewer, so do *not* do extensive computation here (leave it to 
     * simulateOneStep()).
     */
    virtual void updateRenderGeometry() = 0;

    /*
     * Perform any actual rendering here. This method *must* be thread-safe with respect to simulateOneStep().
     * This method runs in the same thread as the viewer and blocks user IO, so there really should not be any
     * extensive computation here or the UI will lag/become unresponsive (the whole reason the simulation itself
     * is in its own thread.)
     */
    virtual void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer) = 0;

    /*
    * Called when the user clicks on the simulation panel.
    * This method is called in the *rendering thread*. If you need to make changes to the simulation state, you
    * should stash the mouse click in a message queue and deal with it in the simulation thread.
    */
    virtual bool mouseClicked(igl::opengl::glfw::Viewer &viewer, int button) { return false; }

    /*
    * Called when the user unclicks the mouse on the simulation panel.
    * This method is called in the *rendering thread*. If you need to make changes to the simulation state, you
    * should stash the mouse click in a message queue and deal with it in the simulation thread.
    */
    virtual bool mouseReleased(igl::opengl::glfw::Viewer &viewer, int button) { return false; }

    /*
    * Called when the user drags the mouse on the simulation panel.
    * This method is called in the *rendering thread*. If you need to make changes to the simulation state, you
    * should stash the mouse click in a message queue and deal with it in the simulation thread.
    */
    virtual bool mouseMoved(igl::opengl::glfw::Viewer &viewer, int button) { return false; }

    /*
     * Runs the simulation, if it has been paused (or never started).
     */
    void run()
    {
        status_mutex.lock();
        please_pause = false;
        status_mutex.unlock();
    }

    /*
     * Resets the simulation (and leaves it in a paused state; call run() to start it).
     */
    void reset()
    {
        killSimThread();
        please_die = running = false;
        please_pause = true;
        initSimulation();
        updateRenderGeometry();
        sim_thread = new std::thread(&PhysicsHook::runSimThread, this);
    }

    /*
     * Pause a running simulation. The simulation will pause at the end of its current "step"; this method will not
     * interrupt simulateOneStep mid-processing.
     */
    void pause()
    {
        status_mutex.lock();
        please_pause = true;
        status_mutex.unlock();
    }

    bool isPaused()
    {
        bool ret = false;
        status_mutex.lock();
        if(running && please_pause)
            ret = true;
        status_mutex.unlock();
        return ret;
    }

    void render(igl::opengl::glfw::Viewer &viewer)
    {
        render_mutex.lock();
        renderRenderGeometry(viewer);
        render_mutex.unlock();
    }
    
protected:
    void runSimThread()
    {
        status_mutex.lock();
        running = true;
        status_mutex.unlock();

        bool done = false;
        while (!done)
        {
            tick();

            status_mutex.lock();
            bool pausenow = please_pause;
            status_mutex.unlock();
            if (pausenow)
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }
            else
            {
                done = simulateOneStep();
                render_mutex.lock();
                updateRenderGeometry();
                render_mutex.unlock();
            }
            status_mutex.lock();
            if (please_die)
                done = true;
            status_mutex.unlock();
        }

        status_mutex.lock();
        running = false;
        status_mutex.unlock();
    }

    void killSimThread()
    {
        if (sim_thread)
        {
            status_mutex.lock();
            please_die = true;
            status_mutex.unlock();
            sim_thread->join();
            delete sim_thread;
            sim_thread = NULL;
        }
    }

    std::thread *sim_thread;
    bool please_pause;
    bool please_die;
    bool running;
    std::mutex render_mutex;
    std::mutex status_mutex;
};

#endif