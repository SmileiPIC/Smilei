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

void internal_inf( Particles &particles, int imin, int imax, int direction, double limit_inf, Species *species,
                   double &nrj_iPart );

void internal_sup( Particles &particles, int imin, int imax, int direction, double limit_sup, Species *species,
                   double &nrj_iPart );

void reflect_particle_inf( Particles &particles, int imin, int imax, int direction, double limit_inf, Species *species,
                           double &nrj_iPart );

void reflect_particle_sup( Particles &particles, int imin, int imax, int direction, double limit_sup, Species *species,
                           double &nrj_iPart );

void reflect_particle_wall( Particles &particles, int imin, int imax, int direction, double limit_sup, Species *species,
                           double &nrj_iPart );

// direction not used below, direction is "r"
void refl_particle_AM( Particles &particles, int imin, int imax, int direction, double limit_sup, Species *species,
                      double &nrj_iPart );

void remove_particle_inf( Particles &particles, int imin, int imax, int direction, double limit_inf, Species *species,
                         double &nrj_iPart );

void remove_particle_sup( Particles &particles, int imin, int imax, int direction, double limit_sup, Species *species,
                         double &nrj_iPart );

void remove_particle_wall( Particles &particles, int imin, int imax, int direction, double limit_sup, Species *species,
                         double &nrj_iPart );

//! Delete photon (mass_==0) at the boundary and keep the energy for diagnostics
void remove_photon_inf( Particles &particles, int imin, int imax, int direction, double limit_inf, Species *species,
                       double &nrj_iPart );

void remove_photon_sup( Particles &particles, int imin, int imax, int direction, double limit_sup, Species *species,
                       double &nrj_iPart );

void stop_particle_inf( Particles &particles, int imin, int imax, int direction, double limit_inf, Species *species,
                       double &nrj_iPart );

void stop_particle_sup( Particles &particles, int imin, int imax, int direction, double limit_sup, Species *species,
                       double &nrj_iPart );

void stop_particle_wall( Particles &particles, int imin, int imax, int direction, double limit_sup, Species *species,
                       double &nrj_iPart );

void stop_particle_AM( Particles &particles, int imin, int imax, int direction, double limit_pos, Species *species,
                      double &nrj_iPart );

//!\todo (MG) at the moment the particle is thermalize whether or not there is a plasma initially at the boundary.
// ATTENTION: here the thermalization assumes a Maxwellian distribution, maybe we should add some checks on thermal_boundary_temperature (MG)!
inline void thermalize_particle( Particles &particles, int imin, int imax, int direction, double limit_pos,
                                 Species *species, double &nrj_iPart )
{
    int ipart = imin;
    // checking the particle's velocity compared to the thermal one
    double p2 = 0.;
    for( unsigned int i=0; i<3; i++ ) {
        p2 += particles.momentum( i, ipart )*particles.momentum( i, ipart );
    }
    double v = sqrt( p2 )/particles.LorentzFactor( ipart );
    
    // energy before thermalization
    nrj_iPart = particles.weight( ipart )*( particles.LorentzFactor( ipart )-1.0 );
    
    // Apply bcs depending on the particle velocity
    // --------------------------------------------
    if( v>3.0*species->thermal_velocity_[0] ) {     //IF VELOCITY > 3*THERMAL VELOCITY THEN THERMALIZE IT
    
        // velocity of the particle after thermalization/reflection
        //for (int i=0; i<species->nDim_fields; i++) {
        for( int i=0; i<3; i++ ) {
        
            if( i==direction ) {
                // change of velocity in the direction normal to the reflection plane
                double sign_vel = -particles.momentum( i, ipart )/std::abs( particles.momentum( i, ipart ) );
                particles.momentum( i, ipart ) = sign_vel * species->thermal_momentum_[i]
                                                 *                             std::sqrt( -std::log( 1.0-Rand::uniform1() ) );
                                                 
            } else {
                // change of momentum in the direction(s) along the reflection plane
                double sign_rnd = Rand::uniform() - 0.5;
                sign_rnd = ( sign_rnd )/std::abs( sign_rnd );
                particles.momentum( i, ipart ) = sign_rnd * species->thermal_momentum_[i]
                                                 *                             userFunctions::erfinv( Rand::uniform1() );
            }//if
            
        }//i
        
        // Adding the mean velocity (using relativistic composition)
        double vx, vy, vz, v2, g, gm1, Lxx, Lyy, Lzz, Lxy, Lxz, Lyz, gp, px, py, pz;
        // mean-velocity
        vx  = -species->thermal_boundary_velocity_[0];
        vy  = -species->thermal_boundary_velocity_[1];
        vz  = -species->thermal_boundary_velocity_[2];
        v2  = vx*vx + vy*vy + vz*vz;
        if( v2>0. ) {
        
            g   = 1.0/sqrt( 1.0-v2 );
            gm1 = g - 1.0;
            
            // compute the different component of the Matrix block of the Lorentz transformation
            Lxx = 1.0 + gm1 * vx*vx/v2;
            Lyy = 1.0 + gm1 * vy*vy/v2;
            Lzz = 1.0 + gm1 * vz*vz/v2;
            Lxy = gm1 * vx*vy/v2;
            Lxz = gm1 * vx*vz/v2;
            Lyz = gm1 * vy*vz/v2;
            
            // Lorentz transformation of the momentum
            gp = sqrt( 1.0 + pow( particles.momentum( 0, ipart ), 2 ) + pow( particles.momentum( 1, ipart ), 2 )
                       + pow( particles.momentum( 2, ipart ), 2 ) );
            px = -gp*g*vx + Lxx * particles.momentum( 0, ipart ) + Lxy * particles.momentum( 1, ipart ) + Lxz * particles.momentum( 2, ipart );
            py = -gp*g*vy + Lxy * particles.momentum( 0, ipart ) + Lyy * particles.momentum( 1, ipart ) + Lyz * particles.momentum( 2, ipart );
            pz = -gp*g*vz + Lxz * particles.momentum( 0, ipart ) + Lyz * particles.momentum( 1, ipart ) + Lzz * particles.momentum( 2, ipart );
            particles.momentum( 0, ipart ) = px;
            particles.momentum( 1, ipart ) = py;
            particles.momentum( 2, ipart ) = pz;
            
        }//ENDif vel != 0
        
    } else {                                    // IF VELOCITY < 3*THERMAL SIMPLY REFLECT IT
        particles.momentum( direction, ipart ) = -particles.momentum( direction, ipart );
        
    }// endif on v vs. thermal_velocity_
    
    // position of the particle after reflection
    particles.position( direction, ipart ) = limit_pos - particles.position( direction, ipart );
    
    // energy lost during thermalization
    nrj_iPart -= particles.weight( ipart )*( particles.LorentzFactor( ipart )-1.0 );
    
    
    /* HERE IS AN ATTEMPT TO INTRODUCE A SPACE DEPENDENCE ON THE BCs
    // double val_min(params.dens_profile.vacuum_length[1]), val_max(params.dens_profile.vacuum_length[1]+params.dens_profile.length_params_y[0]);
    
    if ( ( particles.position(1,ipart) >= val_min ) && ( particles.position(1,ipart) <= val_max ) ) {
        // nrj computed during diagnostics
        particles.position(direction, ipart) = limit_pos - particles.position(direction, ipart);
        particles.momentum(direction, ipart) = sqrt(params.thermal_velocity_[direction]) * tabFcts.erfinv( Rand::uniform() );
    }
    else {
        stop_particle( particles, ipart, direction, limit_pos, params, nrj_iPart );
    }
     */
    
}

#endif
