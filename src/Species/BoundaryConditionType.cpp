/*! @file BoundaryConditionType.cpp

 @brief BoundaryConditionType.cpp for particle boundary conditions

 */

#include <cmath>
#include <cstdlib>

#include "BoundaryConditionType.h"
#include "Params.h"
#include "tabulatedFunctions.h"
#include "userFunctions.h"

//!
//! int function( Particles &particles, int ipart, int direction, double limit_pos )
//!     returns :
//!         0 if particle ipart have to be deleted from current process (MPI or BC)
//!         1 otherwise
//!

void internal_inf( Particles &particles, SmileiMPI* smpi, int imin, int imax,
                          int direction, double limit_inf,
                          double dt, Species *species, int ithread, double &nrj_iPart )
{
    nrj_iPart = 0.;     // no energy loss during exchange
    double* position = &(particles.position(direction,0));
    int* cell_keys = &(particles.cell_keys[0]);
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] < limit_inf) {
            cell_keys[ ipart ] = -1;
        }
    }
}

void internal_sup( Particles &particles, SmileiMPI* smpi, int imin, int imax,
                          int direction, double limit_sup,
                          double dt, Species *species, int ithread, double &nrj_iPart )
{
    nrj_iPart = 0.;     // no energy loss during exchange
    double* position = &(particles.position(direction,0));
    int* cell_keys = &(particles.cell_keys[0]);
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] >= limit_sup) {
            cell_keys[ ipart ] = -1;
        }
    }
}


void reflect_particle_inf( Particles &particles, SmileiMPI* smpi, int imin, int imax,
                                  int direction, double limit_inf,
                                  double dt, Species *species, int ithread, double &nrj_iPart )
{
    nrj_iPart = 0.;     // no energy loss during reflection
    double* position = &(particles.position(direction,0));
    double* momentum = &(particles.momentum(direction,0));
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] < limit_inf ) {
            position[ ipart ] = 2.*limit_inf - position[ ipart ];
            momentum[ ipart ] = -momentum[ ipart ];
        }
    }
}

void reflect_particle_sup( Particles &particles, SmileiMPI* smpi, int imin, int imax,
                                  int direction, double limit_sup,
                                  double dt, Species *species, int ithread, double &nrj_iPart )
{
    nrj_iPart = 0.;     // no energy loss during reflection
    double* position = &(particles.position(direction,0));
    double* momentum = &(particles.momentum(direction,0));
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] >= limit_sup) {
            position[ ipart ] = 2*limit_sup - position[ ipart ];
            momentum[ ipart ] = -momentum[ ipart ];
        }
    }
}

void reflect_particle_wall( Particles &particles, SmileiMPI* smpi, int imin, int imax,
                                  int direction, double wall_position,
                                  double dt, Species *species, int ithread, double &nrj_iPart )
{
    nrj_iPart = 0.;     // no energy loss during reflection
    double* position = &(particles.position(direction,0));
    double* momentum = &(particles.momentum(direction,0));
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        double particle_position     = position[ipart];
        double particle_position_old = particle_position - dt*smpi->dynamics_invgf[ithread][ipart]*momentum[ipart]; 
        if ( ( wall_position-particle_position_old )*( wall_position-particle_position )<0) {
            position[ ipart ] = 2*wall_position - position[ ipart ];
            momentum[ ipart ] = -momentum[ ipart ];
        }
    }
}

// direction not used below, direction is "r"
void refl_particle_AM( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_sup, double dt, Species *species,
                             int ithread, double &nrj_iPart )
{
    nrj_iPart = 0.;     // no energy loss during reflection
    
    double* position_y = &(particles.position(1,0));
    double* position_z = &(particles.position(2,0));
    double* momentum_y = &(particles.momentum(1,0));
    double* momentum_z = &(particles.momentum(2,0));
    //limite_sup = 2*Rmax.
    //We look for the coordiunate of the point at which the particle crossed the boundary
    //We need to fine the parameter t at which (y + py*t)+(z+pz*t) = Rmax^2
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( particles.distance2ToAxis( ipart ) >= limit_sup*limit_sup ) {
            limit_sup *= 2.;
            double b = 2*( position_y[ ipart ]*momentum_y[ ipart ] + position_z[ ipart ]*momentum_z[ ipart ] );
            double r2 = ( position_y[ ipart ]*position_y[ ipart ] + position_z[ ipart ]*position_z[ ipart ] );
            double pr2 = ( momentum_y[ ipart ]*momentum_y[ ipart ] + momentum_z[ ipart ]*momentum_z[ ipart ] );
            double delta = b*b - pr2*( 4*r2 - limit_sup*limit_sup );
            
            //b and delta are neceseraliy >=0 otherwise there are no solution which means that something unsual happened
            if( b < 0 || delta < 0 ) {
                ERROR( "There are no solution to reflexion. This should never happen" );
            }
    
            double t = ( -b + sqrt( delta ) )/pr2*0.5;
    
            double y0, z0; //Coordinates of the crossing point 0
            y0 =  position_y[ ipart ] + momentum_y[ ipart ]*t ;
            z0 =  position_z[ ipart ] + momentum_z[ ipart ]*t ;
    
            //Update new particle position as a reflexion to the plane tangent to the circle at 0.
    
            position_y[ ipart ] -= 4*( position_y[ ipart ]-y0 )*y0/limit_sup ;
            position_z[ ipart ] -= 4*( position_z[ ipart ]-z0 )*z0/limit_sup ;
            
            momentum_y[ ipart ] *= 1 - 2*y0/limit_sup ;
            momentum_z[ ipart ] *= 1 - 2*z0/limit_sup ;
        }
    }    
}

void remove_particle_inf( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_inf, double dt, Species *species,
                                int ithread, double &nrj_iPart )
{
    nrj_iPart = 0.;
    double* position = &(particles.position(direction,0));
    short* charge = &(particles.charge(0));
    int* cell_keys = &(particles.cell_keys[0]);
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] < limit_inf) {
            nrj_iPart += particles.weight( ipart )*( particles.LorentzFactor( ipart )-1.0 ); // energy lost REDUCTION
            charge[ ipart ] = 0;
            cell_keys[ipart] = -1;
        }
    }
}

void remove_particle_sup( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_sup, double dt, Species *species,
                                int ithread, double &nrj_iPart )
{
    nrj_iPart = 0.;
    double* position = &(particles.position(direction,0));
    short* charge = &(particles.charge(0));
    int* cell_keys = &(particles.cell_keys[0]);
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] >= limit_sup) {
            nrj_iPart += particles.weight( ipart )*( particles.LorentzFactor( ipart )-1.0 ); // energy lost REDUCTION
            charge[ ipart ] = 0;
            cell_keys[ipart] = -1;
        }
    }
}

void remove_particle_wall( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double wall_position, double dt, Species *species,
                                int ithread, double &nrj_iPart )
{
    nrj_iPart = 0.;
    double* position = &(particles.position(direction,0));
    double* momentum = &(particles.momentum(direction,0));
    short* charge = &(particles.charge(0));
    int* cell_keys = &(particles.cell_keys[0]);
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        double particle_position     = position[ipart];
        double particle_position_old = particle_position - dt*smpi->dynamics_invgf[ithread][ipart]*momentum[ipart]; 
        if ( ( wall_position-particle_position_old )*( wall_position-particle_position )<0) {
            nrj_iPart += particles.weight( ipart )*( particles.LorentzFactor( ipart )-1.0 ); // energy lost REDUCTION
            charge[ ipart ] = 0;
            cell_keys[ipart] = -1;
        }
    }
}

//! Delete photon (mass_==0) at the boundary and keep the energy for diagnostics
void remove_photon_inf( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_inf, double dt, Species *species,
                              int ithread, double &nrj_iPart )
{
    nrj_iPart = 0.;
    double* position = &(particles.position(direction,0));
    short* charge = &(particles.charge(0));
    int* cell_keys = &(particles.cell_keys[0]);
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] < limit_inf) {
            nrj_iPart += particles.weight( ipart )*( particles.momentumNorm( ipart ) ); // energy lost REDUCTION
            charge[ ipart ] = 0;
            cell_keys[ipart] = -1;
        }
    }
}

void remove_photon_sup( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_sup, double dt, Species *species,
                              int ithread, double &nrj_iPart )
{
    nrj_iPart = 0.;
    double* position = &(particles.position(direction,0));
    short* charge = &(particles.charge(0));
    int* cell_keys = &(particles.cell_keys[0]);
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] >= limit_sup) {
            nrj_iPart += particles.weight( ipart )*( particles.momentumNorm( ipart ) ); // energy lost REDUCTION
            charge[ ipart ] = 0;
            cell_keys[ipart] = -1;
        }
    }
}

void stop_particle_inf( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_inf, double dt, Species *species,
                              int ithread, double &nrj_iPart )
{
    nrj_iPart = 0;
    double* position = &(particles.position(direction,0));
    double* momentum_x = &(particles.momentum(0,0));
    double* momentum_y = &(particles.momentum(1,0));
    double* momentum_z = &(particles.momentum(2,0));
    double* weight = &(particles.weight(0));
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] < limit_inf) {
            nrj_iPart = particles.weight( ipart )*( particles.LorentzFactor( ipart )-1.0 ); // energy lost REDUCTION
            position[ ipart ] = 2.*limit_inf - position[ ipart ];
            momentum_x[ ipart ] = 0.;
            momentum_y[ ipart ] = 0.;
            momentum_z[ ipart ] = 0.;
        }
    }
}

void stop_particle_sup( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_sup, double dt, Species *species,
                              int ithread, double &nrj_iPart )
{
    nrj_iPart = 0;
    double* position = &(particles.position(direction,0));
    double* momentum_x = &(particles.momentum(0,0));
    double* momentum_y = &(particles.momentum(1,0));
    double* momentum_z = &(particles.momentum(2,0));
    double* weight = &(particles.weight(0));
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] >= limit_sup) {
            nrj_iPart += particles.weight( ipart )*( particles.LorentzFactor( ipart )-1.0 ); // energy lost REDUCTION
            position[ ipart ] = 2.*limit_sup - position[ ipart ];
            momentum_x[ ipart ] = 0.;
            momentum_y[ ipart ] = 0.;
            momentum_z[ ipart ] = 0.;
        }
    }
}

void stop_particle_wall( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double wall_position, double dt, Species *species,
                              int ithread, double &nrj_iPart )
{
    nrj_iPart = 0;
    double* position = &(particles.position(direction,0));
    double* momentum_x = &(particles.momentum(0,0));
    double* momentum_y = &(particles.momentum(1,0));
    double* momentum_z = &(particles.momentum(2,0));
    double* weight = &(particles.weight(0));
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        double particle_position     = position[ipart];
        double particle_position_old = particle_position - dt*smpi->dynamics_invgf[ithread][ipart]*particles.Momentum[direction][ipart]; 
        if ( ( wall_position-particle_position_old )*( wall_position-particle_position )<0 ) {
            nrj_iPart += particles.weight( ipart )*( particles.LorentzFactor( ipart )-1.0 ); // energy lost REDUCTION
            position[ ipart ] = 2.*wall_position - position[ ipart ];
            momentum_x[ ipart ] = 0.;
            momentum_y[ ipart ] = 0.;
            momentum_z[ ipart ] = 0.;
        }
    }
}

void stop_particle_AM( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_pos, double dt, Species *species,
                             int ithread, double &nrj_iPart )
{
    double* position_y = &(particles.position(1,0));
    double* position_z = &(particles.position(2,0));
    double* momentum_y = &(particles.momentum(1,0));
    double* momentum_z = &(particles.momentum(2,0));


    nrj_iPart = 0;
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        nrj_iPart += particles.weight( ipart )*( particles.LorentzFactor( ipart )-1.0 ); // energy lost REDUCTION
        double distance_to_axis = sqrt( particles.distance2ToAxis( ipart ) );
        // limit_pos = 2*limit_pos
        double new_dist_to_axis = limit_pos - distance_to_axis;
    
        double delta = distance_to_axis - new_dist_to_axis;
        double cos = position_y[ ipart ] / distance_to_axis;
        double sin = position_z[ ipart ] / distance_to_axis;
    
        position_y[ ipart ] -= delta * cos ;
        position_z[ ipart ] -= delta * sin ;
    
        particles.momentum( 0, ipart ) = 0.;
        momentum_y[ ipart ] = 0.;
        momentum_z[ ipart ] = 0.;

    }
    
}
