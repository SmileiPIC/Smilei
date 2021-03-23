#ifndef TIMER_H
#define TIMER_H

#include <string>
#include <vector>

#include "SmileiMPI.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class Timer
//  --------------------------------------------------------------------------------------------------------------------
class Timer
{
    friend class DiagnosticPerformances;
public:
    //! Timer creator, (default no MPI set)
    Timer( std::string name );
    //! Timer destructor
    ~Timer();
    //! Init t0 of the timer, synchronized through MPI
    void init( SmileiMPI *smpi );
    //! Accumulate time couting from last init/restart
    void update( bool store = false );
    
    //! Accumulate time couting from last init/restart without omp master for tasking
    void updateInTask( bool store = false );

    
#ifdef __DETAILED_TIMERS
    //! Accumulate time couting from last init/restart using patch detailed timers
    void update( VectorPatch &vecPatches, bool store = false );
    
    //! Accumulate time couting from last init/restart using patch detailed timers spreaded between threads
    void updateThreaded( VectorPatch &vecPatches, bool store = false );
#endif
    
    //! Start a new cumulative period
    void restart();
    //! Start a new cumulative period without omp master for tasking
    void restartInTask();

    
    //! Start a new cumulative period
    void reboot();
    //! Return accumulated time
    double getTime()
    {
        return time_acc_;
    }
    //! Print accumulated time in stdout
    void print( double tot );
    //! name of the timer
    inline std::string name()
    {
        return name_;
    }
    
    //! Timer name
    std::string name_;
    
    //! Accumulated time in current timer
    double time_acc_;
    
    std::vector<double> register_timers;
    
#ifdef __DETAILED_TIMERS
    //! Id of the associated timer in the patch timer array
    unsigned int patch_timer_id;
#endif
    
private:
    //! Last timer start
    double last_start_;
    //! MPI process timer synchronized through MPI
    SmileiMPI *smpi_;
    
};


#endif
