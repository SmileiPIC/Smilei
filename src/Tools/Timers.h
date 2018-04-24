#ifndef TIMERS_H
#define TIMERS_H

#include <string>
#include <vector>

#include "Timer.h"

class SmileiMPI;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Timers
//  --------------------------------------------------------------------------------------------------------------------
class Timers {
public:
    //! Constructor
    Timers(SmileiMPI *smpi );
    //! Destructor
    ~Timers();

    Timer global    ;
    Timer particles ;
    Timer maxwell   ;
    Timer diags     ;
    Timer densities ;
    Timer collisions;
    Timer movWindow ;
    Timer loadBal   ;
    Timer syncPart  ;
    Timer syncField ;
    Timer syncDens  ;
    Timer diagsNEW  ;
    Timer reconfiguration  ;

    void profile(SmileiMPI * smpi);
    std::vector<Timer*> consolidate(SmileiMPI * smpi);
    void reboot();

private:
    std::vector<Timer*> timers;

};


#endif
