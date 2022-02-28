#ifndef PARTCOMPTIME3D4ORDER_H
#define PARTCOMPTIME3D4ORDER_H

#include "PartCompTime.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PartCompTime3D4Order
//  --------------------------------------------------------------------------------------------------------------------
class PartCompTime3D4Order final : public PartCompTime
{
public:
    PartCompTime3D4Order();
    ~PartCompTime3D4Order() override final {};
    
    inline float __attribute__((always_inline)) getParticleComputationTimeVecto( const float log_particle_number );
    inline float __attribute__((always_inline)) getParticleComputationTimeScalar( const float log_particle_number );
    
private:

};//END class PartCompTime3D2Order

#endif
