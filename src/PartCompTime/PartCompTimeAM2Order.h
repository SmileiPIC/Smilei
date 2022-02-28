#ifndef PARTCOMPTIMEAM2ORDER_H
#define PARTCOMPTIMEAM2ORDER_H

#include "PartCompTime.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PartCompTimeAM2Order
//! Evaluation du temps de calcul des macro-particules en mode AM Order 2
//  --------------------------------------------------------------------------------------------------------------------
class PartCompTimeAM2Order final : public PartCompTime
{
public:
    PartCompTimeAM2Order();
    ~PartCompTimeAM2Order() override final {};
    
    inline float __attribute__((always_inline)) getParticleComputationTimeVecto( const float log_particle_number );
    inline float __attribute__((always_inline)) getParticleComputationTimeScalar( const float log_particle_number );
    
private:

};//END class PartCompTimeAM2Order

#endif
