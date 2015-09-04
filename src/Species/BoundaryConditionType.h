/*! @file BoundaryConditionType.h
 
 @brief BoundaryConditionType.h for particle boundary conditions

 */

#ifndef BOUNDARYCONDITIONTYPE_H
#define BOUNDARYCONDITIONTYPE_H

#include <cmath>
#include <cstdlib>

#include "Particles.h"
#include "PicParams.h"
//#include "tabulatedFunctions.h"
#include "userFunctions.h"

//!
//! int function( Particles &particles, int ipart, int direction, double limit_pos )
//!     returns :
//!         0 if particle ipart have to be deleted from current process (MPI or BC)
//!         1 otherwise
//!

inline int refl_particle( Particles &particles, int ipart, int direction, double limit_pos, SpeciesStructure &params, double &nrj_iPart) {
    nrj_iPart = 0.;     // no energy loss during reflection
    particles.position(direction, ipart) = limit_pos - particles.position(direction, ipart);
    particles.momentum(direction, ipart) = -particles.momentum(direction, ipart);
    return 1;
}

inline int supp_particle( Particles &particles, int ipart, int direction, double limit_pos, SpeciesStructure &params, double &nrj_iPart) {
    nrj_iPart = particles.weight(ipart)*(particles.lor_fac(ipart)-1.0); // energy lost
    particles.position(direction, ipart) = particles.position_old(direction, ipart);
    particles.charge(ipart) = 0;
    return 0;
}

inline int stop_particle( Particles &particles, int ipart, int direction, double limit_pos, SpeciesStructure &params, double &nrj_iPart) {
    nrj_iPart = particles.weight(ipart)*(particles.lor_fac(ipart)-1.0); // energy lost
    particles.position(direction, ipart) = particles.position_old(direction, ipart);
    particles.momentum(0, ipart) = 0.;
    particles.momentum(1, ipart) = 0.;
    particles.momentum(2, ipart) = 0.;
    return 1;

}

//!\todo (MG) at the moment the particle is thermalize whether or not there is a plasma initially at the boundary.
// ATTENTION: here the thermalization assumes a Maxwellian distribution, maybe we should add some checks on thermT (MG)!
inline int thermalize_particle( Particles &particles, int ipart, int direction, double limit_pos, SpeciesStructure &params, double &nrj_iPart) {
    
    // checking the particle's velocity compared to the thermal one
    double p2 = 0.;
    for (unsigned int i=0; i<3; i++) {
        p2 += pow( particles.momentum(i,ipart) , 2);
    }
    double v = sqrt(p2)/particles.lor_fac(ipart);
    
    // energy before thermalization
    nrj_iPart = particles.weight(ipart)*(particles.lor_fac(ipart)-1.0);
    
    // Apply bcs depending on the particle velocity
    // --------------------------------------------
    if ( v>3.0*params.thermalVelocity[0] ) {    //IF VELOCITY > 3*THERMAL VELOCITY THEN THERMALIZE IT

        // velocity of the particle after reflection (unchanged in the directions that are not resolved in the simulations)
        for (int i=0; i<params.nDim_fields; i++) {
            
            if (i==direction) {
                // change of velocity in the direction normal to the reflection plane
                double sign_vel = -(particles.position(direction, ipart)-0.5*limit_pos)
                /          std::abs(particles.position(direction, ipart)-0.5*limit_pos);
                particles.momentum(i,ipart) = sign_vel * params.thermalMomentum[i]
                *                             std::sqrt( -std::log(1.0-((double)rand() / RAND_MAX)) );
                
            } else {
                // change of momentum in the direction(s) along the reflection plane
                double sign_rnd = (double)rand() / RAND_MAX - 0.5; sign_rnd = (sign_rnd)/std::abs(sign_rnd);
                particles.momentum(i,ipart) = sign_rnd * params.thermalMomentum[i]
                *                               userFunctions::erfinv( (double)rand() / RAND_MAX );
                //*                             erfinv::instance().call( (double)rand() / RAND_MAX );
            }//if
            
        }//i

    } else {                                    // IF VELOCITY < 3*THERMAL SIMPLY REFLECT IT
        particles.momentum(direction, ipart) = -particles.momentum(direction, ipart);
        
    }// endif on v vs. thermalVelocity
    
    // position of the particle after reflection
    particles.position(direction, ipart) = limit_pos - particles.position(direction, ipart);
    
    // energy lost during thermalization
    nrj_iPart -= particles.weight(ipart)*(particles.lor_fac(ipart)-1.0);
    
    
    /* HERE IS AN ATTEMPT TO INTRODUCE A SPACE DEPENDENCE ON THE BCs
    // double val_min(params.dens_profile.vacuum_length[1]), val_max(params.dens_profile.vacuum_length[1]+params.dens_profile.length_params_y[0]);

    if ( ( particles.position(1,ipart) >= val_min ) && ( particles.position(1,ipart) <= val_max ) ) {
	// nrj computed during diagnostics
	particles.position(direction, ipart) = limit_pos - particles.position(direction, ipart);
        particles.momentum(direction, ipart) = sqrt(params.thermalVelocity[direction]) * tabFcts.erfinv( (double)rand() / RAND_MAX );
    }
    else {
	stop_particle( particles, ipart, direction, limit_pos, params, nrj_iPart );
    }
     */
    
    return 1;

}

#endif
