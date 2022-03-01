#ifndef PARTCOMPTIME_H
#define PARTCOMPTIME_H

#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

//  ----------------------------------------------------------------------------
//! Class PartCompTime
//  ----------------------------------------------------------------------------
class PartCompTime
{
public:
    PartCompTime( );
    virtual ~PartCompTime() {};
    
    // -------------------------------------------------------------------------
    //! Evaluate the time (simple precision) to compute all particles
    //! in the current patch with vectorized operators
    //! @param count the numer of particles per cell
    //! @aram vecto_time time in vector mode
    //! @aram scalar_time time in scalar mode
    // -------------------------------------------------------------------------
    virtual void operator() (   const std::vector<int> &count,
                        float &vecto_time,
                        float &scalar_time );
    
    // -------------------------------------------------------------------------
    //! Evaluate the time necessary to compute `log_particle_number` particles
    //! using vectorized operators
    //! @aram log_particle_number number of particles to compute
    //! @return normalized time
    // -------------------------------------------------------------------------
    inline float __attribute__((always_inline)) getParticleComputationTimeVecto( const float log_particle_number );
    
    // -------------------------------------------------------------------------
    //! Evaluate the time necessary to compute `log_particle_number` particles
    //! using scalar operators
    //! @aram log_particle_number number of particles to compute
    //! @return normalized time
    // -------------------------------------------------------------------------
    inline float __attribute__((always_inline)) getParticleComputationTimeScalar( const float log_particle_number );
    
private:

};//END class PartCompTime

#endif
