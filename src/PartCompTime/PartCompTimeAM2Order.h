#ifndef PARTCOMPTIMEAM2ORDER_H
#define PARTCOMPTIMEAM2ORDER_H

#include "PartCompTime.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PartCompTimeAM2Order
//! Evaluation of the particle time in AM at Order 2
//  --------------------------------------------------------------------------------------------------------------------
class PartCompTimeAM2Order final : public PartCompTime
{
public:
    PartCompTimeAM2Order();
    ~PartCompTimeAM2Order() override final {};
    
    // -------------------------------------------------------------------------
    //! Evaluate the time (simple precision) to compute all particles
    //! in the current patch with vectorized operators
    //! @param count the numer of particles per cell
    //! @aram vecto_time time in vector mode
    //! @aram scalar_time time in scalar mode
    // -------------------------------------------------------------------------
    virtual void operator() (   const std::vector<int> &count,
                        float &vecto_time,
                        float &scalar_time ) override final;
    
    inline float __attribute__((always_inline)) getParticleComputationTimeVecto( const float log_particle_number );
    inline float __attribute__((always_inline)) getParticleComputationTimeScalar( const float log_particle_number );
    
private:

};//END class PartCompTimeAM2Order

#endif
