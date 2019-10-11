#ifndef TIMERS_H
#define TIMERS_H

#include <string>
#include <vector>

#include "Timer.h"

class SmileiMPI;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Timers
//  --------------------------------------------------------------------------------------------------------------------
class Timers
{
public:
    //! Constructor
    Timers( SmileiMPI *smpi );
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
    Timer particleMerging;
    Timer particleInjection;
    Timer diagsNEW  ;
    Timer reconfiguration  ;
    Timer envelope  ;
    Timer susceptibility ;
    Timer grids ;
#ifdef __DETAILED_TIMERS
    Timer interpolator  ;
    Timer pusher  ;
    Timer projector  ;
    Timer cell_keys  ;
    Timer ionization  ;
    Timer radiation  ;
    Timer multiphoton_Breit_Wheeler_timer  ;
    
    Timer interp_fields_env  ;
    Timer proj_susceptibility  ;
    Timer push_mom ;
    Timer interp_env_old  ;
    Timer proj_currents  ;
    Timer push_pos ;
    
    Timer sorting ;
    
#endif
    
    // Where the patch timers start in the timer vector
    unsigned int patch_timer_id_start ;
    
    //! Output the timer profile
    void profile( SmileiMPI *smpi );
    
    //! Perform the required processing on the timers for output
    std::vector<Timer *> consolidate( SmileiMPI *smpi, bool final_profile = false );
    
    void reboot();
    
private:
    std::vector<Timer *> timers;
    
};


#endif
