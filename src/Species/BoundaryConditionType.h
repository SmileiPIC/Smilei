/*! @file BoundaryConditionType.h

 @brief BoundaryConditionType.h for particle boundary conditions

 */

#ifndef BOUNDARYCONDITIONTYPE_H
#define BOUNDARYCONDITIONTYPE_H

#include <cmath>
#include <cstdlib>

#include "Particles.h"
#include "Species.h"
#include "Params.h"
#include "tabulatedFunctions.h"
#include "userFunctions.h"

//!
//! int function( Particles &particles, int ipart, int direction, double limit_pos )
//!     returns :
//!         0 if particle ipart have to be deleted from current process (MPI or BC)
//!         1 otherwise
//!

void internal_inf( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_inf, double dt, Species *species,
                   int ithread, double &nrj_iPart );

void internal_sup( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_sup, double dt, Species *species,
                   int ithread, double &nrj_iPart );

void reflect_particle_inf( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_inf, double dt, Species *species,
                           int ithread, double &nrj_iPart );

void reflect_particle_sup( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_sup, double dt, Species *species,
                           int ithread, double &nrj_iPart );

void reflect_particle_wall( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_sup, double dt, Species *species,
                           int ithread, double &nrj_iPart );

// direction not used below, direction is "r"
void refl_particle_AM( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_sup, double dt, Species *species,
                      int ithread, double &nrj_iPart );

void remove_particle_inf( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_inf, double dt, Species *species,
                         int ithread, double &nrj_iPart );

void remove_particle_sup( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_sup, double dt, Species *species,
                         int ithread, double &nrj_iPart );

void remove_particle_wall( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_sup, double dt, Species *species,
                         int ithread, double &nrj_iPart );

//! Delete photon (mass_==0) at the boundary and keep the energy for diagnostics
void remove_photon_inf( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_inf, double dt, Species *species,
                       int ithread, double &nrj_iPart );

void remove_photon_sup( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_sup, double dt, Species *species,
                       int ithread, double &nrj_iPart );

void stop_particle_inf( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_inf, double dt, Species *species,
                       int ithread, double &nrj_iPart );

void stop_particle_sup( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_sup, double dt, Species *species,
                       int ithread, double &nrj_iPart );

void stop_particle_wall( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_sup, double dt, Species *species,
                       int ithread, double &nrj_iPart );

void stop_particle_AM( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_pos, double dt, Species *species,
                      int ithread, double &nrj_iPart );

//!\todo (MG) at the moment the particle is thermalize whether or not there is a plasma initially at the boundary.
// ATTENTION: here the thermalization assumes a Maxwellian distribution, maybe we should add some checks on thermal_boundary_temperature (MG)!
void thermalize_particle_inf( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_pos,
                          double dt, Species *species, int ithread, double &nrj_iPart );

void thermalize_particle_sup( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_pos,
                          double dt, Species *species, int ithread, double &nrj_iPart );

void thermalize_particle_wall( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_pos,
                          double dt, Species *species, int ithread, double &nrj_iPart );


#endif
