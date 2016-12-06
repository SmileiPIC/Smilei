#ifndef TIMER_H
#define TIMER_H

#include <string>
#include <vector>

class SmileiMPI;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Timer
//  --------------------------------------------------------------------------------------------------------------------
class Timer {
public:
    //! Timer creator, (default no MPI set)
    Timer();
    //! Timer destructor
    ~Timer();
    //! Init t0 of the timer, no MPI
    void init(std::string name);
    //! Init t0 of the timer, synchronized through MPI
    void init(SmileiMPI *smpi, std::string name);
    //! Accumulate time couting from last init/restart
    void update(bool store = false);
    //! Start a new cumulative period
    void restart();
    //! Start a new cumulative period
    void reboot();
    //! Return accumulated time
    double getTime(){return time_acc_;}
    //! Print accumulated time in stdout
    void print(double tot);
    //! name of the timer
    inline std::string name() {return name_;}

    //! Timer name 
    std::string name_;
    std::vector<double> register_timers;

    static std::vector<Timer> consolidate_timers( std::vector<Timer> timers );


private:
    //! Accumulated time in current timer
    double time_acc_;
    //! Last timer start
    double last_start_;
    //! MPI process timer synchronized through MPI
    SmileiMPI* smpi_;

};


#endif

