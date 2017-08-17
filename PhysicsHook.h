#ifndef PHYSICSHOOK_H
#define PHYSICSHOOK_H

#include <mutex>
#include <thread>
#include <igl/viewer/Viewer.h>

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

    virtual void initSimulation() = 0;
    virtual bool simulateOneStep() = 0;
    virtual void updateRenderGeometry() = 0;
    virtual void renderRenderGeometry(igl::viewer::Viewer &viewer) = 0;

    void run()
    {
        status_mutex.lock();
        please_pause = false;
        status_mutex.unlock();
    }

    void reset()
    {
        killSimThread();
        please_die = running = false;
        please_pause = true;
        initSimulation();
        updateRenderGeometry();
        sim_thread = new std::thread(&PhysicsHook::runSimThread, this);
    }

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

    void render(igl::viewer::Viewer &viewer)
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