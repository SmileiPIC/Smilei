#ifndef BIRTHRECORDS_H
#define BIRTHRECORDS_H

#include <vector>

#include "Particles.h"
#include "Ionization.h"

//! Records some information when the particles are created by ionization and other mechanisms
//! Required for DiagNewParticles
struct BirthRecords
{
    BirthRecords( Particles &source_particles ) {
        p_.initialize( 0, source_particles );
        p_.double_prop_.push_back( &birth_time_ ); // the birth time is the last property
    };
    void clear() {
        birth_time_.resize( 0 );
        p_.resize( 0 );
    };
    void update( Particles &source_particles, size_t npart, double time_dual, Ionization *I ) {
        size_t prev_size = birth_time_.size();
        birth_time_.resize( prev_size + npart, time_dual );
        source_particles.copyParticles( 0, npart, p_, prev_size );
        // If electrons come from ionization, then the charge is replaced by that of the ionized ions
        if( I && I->save_ion_charge_ ) {
            copy( I->ion_charge_.begin(), I->ion_charge_.end(), &p_.Charge[prev_size] );
            I->ion_charge_.clear();
        }
    }
    
    std::vector<double> birth_time_;
    Particles p_;
};

#endif
