#ifndef PARTCOMPTIME2D4ORDER_H
#define PARTCOMPTIME2D4ORDER_H

#include "PartCompTime.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PartCompTime2D4Order
//! Evaluation of the particle time in 2D at Order 4
//  --------------------------------------------------------------------------------------------------------------------
class PartCompTime2D4Order final : public PartCompTime
{
public:
    PartCompTime2D4Order();
    ~PartCompTime2D4Order() override final {};
    
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

};//END class PartCompTime2D4Order

#endif
