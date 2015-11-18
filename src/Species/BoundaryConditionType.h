/*! @file BoundaryConditionType.h
 
 @brief BoundaryConditionType.h for particle boundary conditions

 */

#ifndef BOUNDARYCONDITIONTYPE_H
#define BOUNDARYCONDITIONTYPE_H

#include "Particles.h"
#include "PicParams.h"

//!
//! int function( Particles &particles, int ipart, int direction, double limit_pos )
//!     returns :
//!         0 if particle ipart have to be deleted from current process (MPI or BC)
//!         1 otherwise
//!

inline int refl_particle( Particles &particles, int ipart, int direction, double limit_pos, SpeciesStructure &params, double &nrj_iPart ) {
    // nrj computed during diagnostics
    nrj_iPart = 0.;//particles.weight(ipart)*(particles.lor_fac(ipart)-1.0);
    particles.position(direction, ipart) = limit_pos - particles.position(direction, ipart);
    particles.momentum(direction, ipart) = -particles.momentum(direction, ipart);
    return 1;

}

inline int supp_particle( Particles &particles, int ipart, int direction, double limit_pos, SpeciesStructure &params, double &nrj_iPart ) {
    //nrj_iPart = particles.weight(ipart)*(particles.lor_fac(ipart)-1.0);
    //particles.position(direction, ipart) = particles.position_old(direction, ipart);
    //particles.charge(ipart) = 0;
    return 0;
}

inline int stop_particle( Particles &particles, int ipart, int direction, double limit_pos, SpeciesStructure &params, double &nrj_iPart ) {
    nrj_iPart = particles.weight(ipart)*(particles.lor_fac(ipart)-1.0);
    particles.position(direction, ipart) = particles.position_old(direction, ipart);
    particles.momentum(0, ipart) = 0.;
    particles.momentum(1, ipart) = 0.;
    particles.momentum(2, ipart) = 0.;
    return 1;

}

inline int adrien_particle( Particles &particles, int ipart, int direction, double limit_pos, SpeciesStructure &params, double &nrj_iPart ) {
    double val_min(params.vacuum_length[1]), val_max(params.vacuum_length[1]+params.dens_length_y[0]);
    if ( ( particles.position(1,ipart) >= val_min ) && ( particles.position(1,ipart) <= val_max ) ) {
	// nrj computed during diagnostics
	particles.position(direction, ipart) = limit_pos - particles.position(direction, ipart);
	//particles.momentum(X, ipart) = random,params.temperature[]
    }
    else {
	stop_particle( particles, ipart, direction, limit_pos, params, nrj_iPart );
    }

    return 1;

}

#endif
