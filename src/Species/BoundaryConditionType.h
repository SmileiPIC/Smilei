/*! @file BoundaryConditionType.h
 
 @brief BoundaryConditionType.h for particle boundary conditions

 */

#ifndef BOUNDARYCONDITIONTYPE_H
#define BOUNDARYCONDITIONTYPE_H

#include <cmath>
#include <cstdlib>

#include "Particles.h"
#include "PicParams.h"
#include "tabulatedFunctions.h"
//!
//! int function( Particles &particles, int ipart, int direction, double limit_pos )
//!     returns :
//!         0 if particle ipart have to be deleted from current process (MPI or BC)
//!         1 otherwise
//!

inline int refl_particle( Particles &particles, int ipart, int direction, double limit_pos, SpeciesStructure &params, double &nrj_iPart, tabulatedFunctions &tabFcts ) {
    nrj_iPart = 0.;     // no energy loss during reflection
    particles.position(direction, ipart) = limit_pos - particles.position(direction, ipart);
    particles.momentum(direction, ipart) = -particles.momentum(direction, ipart);
    return 1;
}

inline int supp_particle( Particles &particles, int ipart, int direction, double limit_pos, SpeciesStructure &params, double &nrj_iPart, tabulatedFunctions &tabFcts ) {
    nrj_iPart = particles.weight(ipart)*(particles.lor_fac(ipart)-1.0); // energy lost
    particles.position(direction, ipart) = particles.position_old(direction, ipart);
    particles.charge(ipart) = 0;
    return 0;
}

inline int stop_particle( Particles &particles, int ipart, int direction, double limit_pos, SpeciesStructure &params, double &nrj_iPart, tabulatedFunctions &tabFcts ) {
    nrj_iPart = particles.weight(ipart)*(particles.lor_fac(ipart)-1.0); // energy lost
    particles.position(direction, ipart) = particles.position_old(direction, ipart);
    particles.momentum(0, ipart) = 0.;
    particles.momentum(1, ipart) = 0.;
    particles.momentum(2, ipart) = 0.;
    return 1;

}

//!\todo (MG) at the moment the particle is thermalize whether or not there is a plasma initially at the boundary
inline int thermalize_particle( Particles &particles, int ipart, int direction, double limit_pos, SpeciesStructure &params, double &nrj_iPart, tabulatedFunctions &tabFcts ) {
    nrj_iPart = particles.weight(ipart)*(particles.lor_fac(ipart)-1.0); // energy before thermalization
    particles.position(direction, ipart) = limit_pos - particles.position(direction, ipart);
    for (int i=0; i<params.nDim_fields; i++) {
        if (i==direction) {
            particles.momentum(i,ipart) = params.thermalVelocity[direction]
            *                             sqrt(-2.0 * log(1.0-((double)rand() / RAND_MAX)) );
        } else {
            particles.momentum(i,ipart) = params.thermalVelocity[direction]
            *                             sqrt( 2.0 * tabFcts.erfinv( (double)rand() / RAND_MAX ) );
        }//if
    }//i
    nrj_iPart -= particles.weight(ipart)*(particles.lor_fac(ipart)-1.0); // energy lost
    
    /* HERE IS AN ATTEMPT TO INTRODUCE A SPACE DEPENDENCE ON THE BCs
    // double val_min(params.dens_profile.vacuum_length[1]), val_max(params.dens_profile.vacuum_length[1]+params.dens_profile.length_params_y[0]);

    if ( ( particles.position(1,ipart) >= val_min ) && ( particles.position(1,ipart) <= val_max ) ) {
	// nrj computed during diagnostics
	particles.position(direction, ipart) = limit_pos - particles.position(direction, ipart);
        particles.momentum(direction, ipart) = sqrt(params.thermalVelocity[direction]) * tabFcts.erfinv( (double)rand() / RAND_MAX );
    }
    else {
	stop_particle( particles, ipart, direction, limit_pos, params, nrj_iPart, tabFcts );
    }
     */
    
    return 1;

}

#endif
