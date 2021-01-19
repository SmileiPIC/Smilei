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
    double* position = particles.getPtrPosition(direction);
    int* cell_keys = particles.getPtrCellKeys();
#ifdef _GPU
    #pragma acc parallel deviceptr(position,cell_keys)
    #pragma acc loop gang worker vector
#endif
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
    double* position = particles.getPtrPosition(direction);
    int* cell_keys = particles.getPtrCellKeys();
#ifdef _GPU
    #pragma acc parallel deviceptr(position,cell_keys)
    #pragma acc loop gang worker vector
#endif
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] >= limit_sup) {
            cell_keys[ ipart ] = -1;
        }
    }
}

void internal_inf_AM( Particles &particles, SmileiMPI* smpi, int imin, int imax,
                          int direction, double limit_inf,
                          double dt, Species *species, int ithread, double &nrj_iPart )
{
    nrj_iPart = 0.;     // no energy loss during exchange
    double* position_y = particles.getPtrPosition(1);
    double* position_z = particles.getPtrPosition(2);
    int* cell_keys = particles.getPtrCellKeys();
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        double distance2ToAxis = position_y[ipart]*position_y[ipart]+position_z[ipart]*position_z[ipart];
        if ( distance2ToAxis < limit_inf*limit_inf ) {
            cell_keys[ ipart ] = -1;
        }
    }
}

void internal_sup_AM( Particles &particles, SmileiMPI* smpi, int imin, int imax,
                          int direction, double limit_sup,
                          double dt, Species *species, int ithread, double &nrj_iPart )
{
    nrj_iPart = 0.;     // no energy loss during exchange
    double* position_y = particles.getPtrPosition(1);
    double* position_z = particles.getPtrPosition(2);
    int* cell_keys = particles.getPtrCellKeys();
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        double distance2ToAxis = position_y[ipart]*position_y[ipart]+position_z[ipart]*position_z[ipart];
        if ( distance2ToAxis >= limit_sup*limit_sup ) {
            cell_keys[ ipart ] = -1;
        }
    }
}

void reflect_particle_inf( Particles &particles, SmileiMPI* smpi, int imin, int imax,
                                  int direction, double limit_inf,
                                  double dt, Species *species, int ithread, double &nrj_iPart )
{
    nrj_iPart = 0.;     // no energy loss during reflection
    double* position = particles.getPtrPosition(direction);
    double* momentum = particles.getPtrMomentum(direction);
#ifdef _GPU
    #pragma acc parallel deviceptr(position,momentum)
    #pragma acc loop gang worker vector
#endif
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
    double* position = particles.getPtrPosition(direction);
    double* momentum = particles.getPtrMomentum(direction);
#ifdef _GPU
    #pragma acc parallel deviceptr(position,momentum)
    #pragma acc loop gang worker vector
#endif
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
    double* position = particles.getPtrPosition(direction);
    double* momentum = particles.getPtrMomentum(direction);
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
    
    double* position_y = particles.getPtrPosition(1);
    double* position_z = particles.getPtrPosition(2);
    double* momentum_y = particles.getPtrMomentum(1);
    double* momentum_z = particles.getPtrMomentum(2);
    //limite_sup = 2*Rmax.
    //We look for the coordiunate of the point at which the particle crossed the boundary
    //We need to fine the parameter t at which (y + py*t)+(z+pz*t) = Rmax^2
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        double distance2ToAxis = position_y[ipart]*position_y[ipart]+position_z[ipart]*position_z[ipart];
        if ( distance2ToAxis >= limit_sup*limit_sup ) {
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
    double* position = particles.getPtrPosition(direction);
    double* momentum_x = particles.getPtrMomentum(0);
    double* momentum_y = particles.getPtrMomentum(1);
    double* momentum_z = particles.getPtrMomentum(2);
    short* charge    = particles.getPtrCharge();
    double* weight   = particles.getPtrWeight();
    int* cell_keys   = particles.getPtrCellKeys();
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] < limit_inf) {
            double LorentzFactor = sqrt( 1.+pow( momentum_x[ipart], 2 )+pow( momentum_y[ipart], 2 )+pow( momentum_z[ipart], 2 ) );
            nrj_iPart += weight[ ipart ]*( LorentzFactor-1.0 ); // energy lost REDUCTION
            charge[ ipart ] = 0;
            cell_keys[ipart] = -1;
        }
    }
}

void remove_particle_sup( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_sup, double dt, Species *species,
                                int ithread, double &nrj_iPart )
{
    nrj_iPart = 0.;
    double* position = particles.getPtrPosition(direction);
    double* momentum_x = particles.getPtrMomentum(0);
    double* momentum_y = particles.getPtrMomentum(1);
    double* momentum_z = particles.getPtrMomentum(2);
    short* charge    = particles.getPtrCharge();
    double* weight   = particles.getPtrWeight();
    int* cell_keys   = particles.getPtrCellKeys();
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] >= limit_sup) {
            double LorentzFactor = sqrt( 1.+pow( momentum_x[ipart], 2 )+pow( momentum_y[ipart], 2 )+pow( momentum_z[ipart], 2 ) );
            nrj_iPart += weight[ ipart ]*( LorentzFactor-1.0 ); // energy lost REDUCTION
            charge[ ipart ] = 0;
            cell_keys[ipart] = -1;
        }
    }
}

void remove_particle_wall( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double wall_position, double dt, Species *species,
                                int ithread, double &nrj_iPart )
{
    nrj_iPart = 0.;
    double* position = particles.getPtrPosition(direction);
    double* momentum = particles.getPtrMomentum(direction);
    double* momentum_x = particles.getPtrMomentum(0);
    double* momentum_y = particles.getPtrMomentum(1);
    double* momentum_z = particles.getPtrMomentum(2);
    short* charge    = particles.getPtrCharge();
    double* weight   = particles.getPtrWeight();
    int* cell_keys   = particles.getPtrCellKeys();
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        double particle_position     = position[ipart];
        double particle_position_old = particle_position - dt*smpi->dynamics_invgf[ithread][ipart]*momentum[ipart]; 
        if ( ( wall_position-particle_position_old )*( wall_position-particle_position )<0) {
            double LorentzFactor = sqrt( 1.+pow( momentum_x[ipart], 2 )+pow( momentum_y[ipart], 2 )+pow( momentum_z[ipart], 2 ) );
            nrj_iPart += weight[ ipart ]*( LorentzFactor-1.0 ); // energy lost REDUCTION
            charge[ ipart ] = 0;
            cell_keys[ipart] = -1;
        }
    }
}

void remove_particle_AM( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_sup, double dt, Species *species,
                                int ithread, double &nrj_iPart )
{
    nrj_iPart = 0.;
    double* position_y = particles.getPtrPosition(1);
    double* position_z = particles.getPtrPosition(2);
    double* momentum_x = particles.getPtrMomentum(0);
    double* momentum_y = particles.getPtrMomentum(1);
    double* momentum_z = particles.getPtrMomentum(2);
    short* charge    = particles.getPtrCharge();
    double* weight   = particles.getPtrWeight();
    int* cell_keys   = particles.getPtrCellKeys();
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        double distance2ToAxis = position_y[ipart]*position_y[ipart]+position_z[ipart]*position_z[ipart];
        if ( distance2ToAxis >= limit_sup*limit_sup ) {
            double LorentzFactor = sqrt( 1.+pow( momentum_x[ipart], 2 )+pow( momentum_y[ipart], 2 )+pow( momentum_z[ipart], 2 ) );
            nrj_iPart += weight[ ipart ]*( LorentzFactor-1.0 ); // energy lost REDUCTION
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
    double* position = particles.getPtrPosition(direction);
    double* momentum_x = particles.getPtrMomentum(0);
    double* momentum_y = particles.getPtrMomentum(1);
    double* momentum_z = particles.getPtrMomentum(2);
    double* weight     = particles.getPtrWeight();
    short* charge    = particles.getPtrCharge();
    int* cell_keys   = particles.getPtrCellKeys();
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] < limit_inf) {
            double momentumNorm = sqrt( pow( momentum_x[ipart], 2 )+pow( momentum_y[ipart], 2 )+pow( momentum_z[ipart], 2 ) );
            nrj_iPart += weight[ ipart ]*( momentumNorm ); // energy lost REDUCTION
            charge[ ipart ] = 0;
            cell_keys[ipart] = -1;
        }
    }
}

void remove_photon_sup( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_sup, double dt, Species *species,
                              int ithread, double &nrj_iPart )
{
    nrj_iPart = 0.;
    double* position = particles.getPtrPosition(direction);
    double* momentum_x = particles.getPtrMomentum(0);
    double* momentum_y = particles.getPtrMomentum(1);
    double* momentum_z = particles.getPtrMomentum(2);
    double* weight     = particles.getPtrWeight();
    short* charge    = particles.getPtrCharge();
    int* cell_keys   = particles.getPtrCellKeys();
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] >= limit_sup) {
            double momentumNorm = sqrt( pow( momentum_x[ipart], 2 )+pow( momentum_y[ipart], 2 )+pow( momentum_z[ipart], 2 ) );
            nrj_iPart += weight[ ipart ]*( momentumNorm ); // energy lost REDUCTION
            charge[ ipart ] = 0;
            cell_keys[ipart] = -1;
        }
    }
}

void stop_particle_inf( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_inf, double dt, Species *species,
                              int ithread, double &nrj_iPart )
{
    nrj_iPart = 0;
    double* position = particles.getPtrPosition(direction);
    double* momentum_x = particles.getPtrMomentum(0);
    double* momentum_y = particles.getPtrMomentum(1);
    double* momentum_z = particles.getPtrMomentum(2);
    double* weight     = particles.getPtrWeight();
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] < limit_inf) {
            double LorentzFactor = sqrt( 1.+pow( momentum_x[ipart], 2 )+pow( momentum_y[ipart], 2 )+pow( momentum_z[ipart], 2 ) );
            nrj_iPart = weight[ ipart ]*( LorentzFactor-1.0 ); // energy lost REDUCTION
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
    double* position = particles.getPtrPosition(direction);
    double* momentum_x = particles.getPtrMomentum(0);
    double* momentum_y = particles.getPtrMomentum(1);
    double* momentum_z = particles.getPtrMomentum(2);
    double* weight     = particles.getPtrWeight();
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] >= limit_sup) {
            double LorentzFactor = sqrt( 1.+pow( momentum_x[ipart], 2 )+pow( momentum_y[ipart], 2 )+pow( momentum_z[ipart], 2 ) );
            nrj_iPart = weight[ ipart ]*( LorentzFactor-1.0 ); // energy lost REDUCTION
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
    double* position = particles.getPtrPosition(direction);
    double* momentum = particles.getPtrMomentum(direction);
    double* momentum_x = particles.getPtrMomentum(0);
    double* momentum_y = particles.getPtrMomentum(1);
    double* momentum_z = particles.getPtrMomentum(2);
    double* weight     = particles.getPtrWeight();
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        double particle_position     = position[ipart];
        double particle_position_old = particle_position - dt*smpi->dynamics_invgf[ithread][ipart]*momentum[ipart];
        if ( ( wall_position-particle_position_old )*( wall_position-particle_position )<0 ) {
            double LorentzFactor = sqrt( 1.+pow( momentum_x[ipart], 2 )+pow( momentum_y[ipart], 2 )+pow( momentum_z[ipart], 2 ) );
            nrj_iPart += weight[ ipart ]*( LorentzFactor-1.0 ); // energy lost REDUCTION
            position[ ipart ] = 2.*wall_position - position[ ipart ];
            momentum_x[ ipart ] = 0.;
            momentum_y[ ipart ] = 0.;
            momentum_z[ ipart ] = 0.;
        }
    }
}

void stop_particle_AM( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_sup, double dt, Species *species,
                             int ithread, double &nrj_iPart )
{
    double* position_y = particles.getPtrPosition(1);
    double* position_z = particles.getPtrPosition(2);
    double* momentum_x = particles.getPtrMomentum(0);
    double* momentum_y = particles.getPtrMomentum(1);
    double* momentum_z = particles.getPtrMomentum(2);
    double* weight     = particles.getPtrWeight();

    nrj_iPart = 0;
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        double distance2ToAxis = position_y[ipart]*position_y[ipart]+position_z[ipart]*position_z[ipart];
        if ( distance2ToAxis >= limit_sup*limit_sup ) {
            double LorentzFactor = sqrt( 1.+pow( momentum_x[ipart], 2 )+pow( momentum_y[ipart], 2 )+pow( momentum_z[ipart], 2 ) );
            nrj_iPart += weight[ ipart ]*( LorentzFactor-1.0 ); // energy lost REDUCTION
            double distance_to_axis = sqrt( distance2ToAxis );
            // limit_pos = 2*limit_pos
            double new_dist_to_axis = limit_sup - distance_to_axis;
    
            double delta = distance_to_axis - new_dist_to_axis;
            double cos = position_y[ ipart ] / distance_to_axis;
            double sin = position_z[ ipart ] / distance_to_axis;
    
            position_y[ ipart ] -= delta * cos ;
            position_z[ ipart ] -= delta * sin ;
    
            momentum_x[ ipart ] = 0.;
            momentum_y[ ipart ] = 0.;
            momentum_z[ ipart ] = 0.;
        }
    }
    
}

void thermalize_particle_inf( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_inf,
                              double dt, Species *species, int ithread, double &nrj_iPart )
{
    int nDim = species->nDim_particle;
    double* position = particles.getPtrPosition(direction);
    double* momentum = particles.getPtrMomentum(direction);
    double* momentumRefl_2D = particles.getPtrMomentum((direction+1)%nDim);
    double* momentumRefl_3D = particles.getPtrMomentum((direction+2)%nDim);
    double* momentum_x = particles.getPtrMomentum(0);
    double* momentum_y = particles.getPtrMomentum(1);
    double* momentum_z = particles.getPtrMomentum(2);
    double* weight     = particles.getPtrWeight();

    nrj_iPart = 0;
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] < limit_inf) {
            // checking the particle's velocity compared to the thermal one
            double p2 = pow( momentum_x[ipart], 2 )+pow( momentum_y[ipart], 2 )+pow( momentum_z[ipart], 2 );
            double LorentzFactor = sqrt( 1.+p2 );
            double v = sqrt( p2 )/LorentzFactor;

            // energy before thermalization
            double initial_energy = LorentzFactor-1.0;

            // Apply bcs depending on the particle velocity
            // --------------------------------------------
            if( v>3.0*species->thermal_velocity_[0] ) {     //IF VELOCITY > 3*THERMAL VELOCITY THEN THERMALIZE IT

                // velocity of the particle after thermalization/reflection
                //for (int i=0; i<species->nDim_fields; i++) {

                // change of velocity in the direction normal to the reflection plane
                double sign_vel = -momentum[ ipart ]/std::abs( momentum[ ipart ] );
                momentum[ ipart ] = sign_vel * species->thermal_momentum_[direction]
                    *                             std::sqrt( -std::log( 1.0-Rand::uniform1() ) );

                // change of momentum in the direction(s) along the reflection plane
                if (nDim>1) {
                    double sign_rnd = Rand::uniform() - 0.5;
                    sign_rnd = ( sign_rnd )/std::abs( sign_rnd );
                    momentumRefl_2D[ ipart ] = sign_rnd * species->thermal_momentum_[(direction+1)%nDim]
                        *                             userFunctions::erfinv( Rand::uniform1() );
                    if (nDim>2) {
                        double sign_rnd = Rand::uniform() - 0.5;
                        sign_rnd = ( sign_rnd )/std::abs( sign_rnd );
                        momentumRefl_3D[ ipart ] = sign_rnd * species->thermal_momentum_[(direction+2)%nDim]
                            *                             userFunctions::erfinv( Rand::uniform1() );
                    }
                }

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
                    gp = sqrt( 1.0 + pow( momentum_x[ipart], 2 )+pow( momentum_y[ipart], 2 )+pow( momentum_z[ipart], 2 ) );
                    px = -gp*g*vx + Lxx * momentum_x[ ipart ] + Lxy * momentum_y[ ipart ] + Lxz * momentum_z[ ipart ];
                    py = -gp*g*vy + Lxy * momentum_x[ ipart ] + Lyy * momentum_y[ ipart ] + Lyz * momentum_z[ ipart ];
                    pz = -gp*g*vz + Lxz * momentum_x[ ipart ] + Lyz * momentum_y[ ipart ] + Lzz * momentum_z[ ipart ];
                    momentum_x[ ipart ] = px;
                    momentum_y[ ipart ] = py;
                    momentum_z[ ipart ] = pz;

                }//ENDif vel != 0

            } else {                                    // IF VELOCITY < 3*THERMAL SIMPLY REFLECT IT
                momentum[ ipart ] = -momentum[ ipart ];

            }// endif on v vs. thermal_velocity_

            // position of the particle after reflection
            position[ ipart ] = 2.*limit_inf - position[ ipart ];

            // energy lost during thermalization
            LorentzFactor = sqrt( 1.+pow( momentum_x[ipart], 2 )+pow( momentum_y[ipart], 2 )+pow( momentum_z[ipart], 2 ) );
            nrj_iPart += weight[ ipart ]*( initial_energy - LorentzFactor+1.0 );


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
    }
}

void thermalize_particle_sup( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double limit_sup,
                              double dt, Species *species, int ithread, double &nrj_iPart )
{
    int nDim = species->nDim_particle;
    double* position = particles.getPtrPosition(direction);
    double* momentum = particles.getPtrMomentum(direction);
    double* momentumRefl_2D = particles.getPtrMomentum((direction+1)%nDim);
    double* momentumRefl_3D = particles.getPtrMomentum((direction+2)%nDim);
    double* momentum_x = particles.getPtrMomentum(0);
    double* momentum_y = particles.getPtrMomentum(1);
    double* momentum_z = particles.getPtrMomentum(2);
    double* weight     = particles.getPtrWeight();

    nrj_iPart = 0;
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] >= limit_sup) {
            // checking the particle's velocity compared to the thermal one
            double p2 = pow( momentum_x[ipart], 2 )+pow( momentum_y[ipart], 2 )+pow( momentum_z[ipart], 2 );
            double LorentzFactor = sqrt( 1.+p2 );
            double v = sqrt( p2 )/LorentzFactor;

            // energy before thermalization
            double initial_energy = LorentzFactor-1.0;

            // Apply bcs depending on the particle velocity
            // --------------------------------------------
            if( v>3.0*species->thermal_velocity_[0] ) {     //IF VELOCITY > 3*THERMAL VELOCITY THEN THERMALIZE IT

                // velocity of the particle after thermalization/reflection
                //for (int i=0; i<species->nDim_fields; i++) {

                // change of velocity in the direction normal to the reflection plane
                double sign_vel = -momentum[ ipart ]/std::abs( momentum[ ipart ] );
                momentum[ ipart ] = sign_vel * species->thermal_momentum_[direction]
                    *                             std::sqrt( -std::log( 1.0-Rand::uniform1() ) );

                // change of momentum in the direction(s) along the reflection plane
                if (nDim>1) {
                    double sign_rnd = Rand::uniform() - 0.5;
                    sign_rnd = ( sign_rnd )/std::abs( sign_rnd );
                    momentumRefl_2D[ ipart ] = sign_rnd * species->thermal_momentum_[(direction+1)%nDim]
                        *                             userFunctions::erfinv( Rand::uniform1() );
                    if (nDim>2) {
                        double sign_rnd = Rand::uniform() - 0.5;
                        sign_rnd = ( sign_rnd )/std::abs( sign_rnd );
                        momentumRefl_3D[ ipart ] = sign_rnd * species->thermal_momentum_[(direction+2)%nDim]
                            *                             userFunctions::erfinv( Rand::uniform1() );
                    }
                }

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
                    gp = sqrt( 1.0 + pow( momentum_x[ipart], 2 )+pow( momentum_y[ipart], 2 )+pow( momentum_z[ipart], 2 ) );
                    px = -gp*g*vx + Lxx * momentum_x[ ipart ] + Lxy * momentum_y[ ipart ] + Lxz * momentum_z[ ipart ];
                    py = -gp*g*vy + Lxy * momentum_x[ ipart ] + Lyy * momentum_y[ ipart ] + Lyz * momentum_z[ ipart ];
                    pz = -gp*g*vz + Lxz * momentum_x[ ipart ] + Lyz * momentum_y[ ipart ] + Lzz * momentum_z[ ipart ];
                    momentum_x[ ipart ] = px;
                    momentum_y[ ipart ] = py;
                    momentum_z[ ipart ] = pz;

                }//ENDif vel != 0

            } else {                                    // IF VELOCITY < 3*THERMAL SIMPLY REFLECT IT
                momentum[ ipart ] = -momentum[ ipart ];

            }// endif on v vs. thermal_velocity_

            // position of the particle after reflection
            position[ ipart ] = 2.*limit_sup - position[ ipart ];

            // energy lost during thermalization
            LorentzFactor = sqrt( 1.+pow( momentum_x[ipart], 2 )+pow( momentum_y[ipart], 2 )+pow( momentum_z[ipart], 2 ) );
            nrj_iPart += weight[ ipart ]*( initial_energy - LorentzFactor+1.0 );


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
    }
}


void thermalize_particle_wall( Particles &particles, SmileiMPI* smpi, int imin, int imax, int direction, double wall_position,
                              double dt, Species *species, int ithread, double &nrj_iPart )
{
    int nDim = species->nDim_particle;
    double* position = particles.getPtrPosition(direction);
    double* momentum = particles.getPtrMomentum(direction);
    double* momentumRefl_2D = particles.getPtrMomentum((direction+1)%nDim);
    double* momentumRefl_3D = particles.getPtrMomentum((direction+2)%nDim);
    double* momentum_x = particles.getPtrMomentum(0);
    double* momentum_y = particles.getPtrMomentum(1);
    double* momentum_z = particles.getPtrMomentum(2);
    double* weight     = particles.getPtrWeight();

    nrj_iPart = 0;
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        double particle_position     = position[ipart];
        double particle_position_old = particle_position - dt*smpi->dynamics_invgf[ithread][ipart]*particles.Momentum[direction][ipart];
        if ( ( wall_position-particle_position_old )*( wall_position-particle_position )<0 ) {
            // checking the particle's velocity compared to the thermal one
            double p2 = pow( momentum_x[ipart], 2 )+pow( momentum_y[ipart], 2 )+pow( momentum_z[ipart], 2 );
            double LorentzFactor = sqrt( 1.+p2 );
            double v = sqrt( p2 )/LorentzFactor;

            // energy before thermalization
            double initial_energy = LorentzFactor-1.0;

            // Apply bcs depending on the particle velocity
            // --------------------------------------------
            if( v>3.0*species->thermal_velocity_[0] ) {     //IF VELOCITY > 3*THERMAL VELOCITY THEN THERMALIZE IT

                // velocity of the particle after thermalization/reflection
                //for (int i=0; i<species->nDim_fields; i++) {

                // change of velocity in the direction normal to the reflection plane
                double sign_vel = -momentum[ ipart ]/std::abs( momentum[ ipart ] );
                momentum[ ipart ] = sign_vel * species->thermal_momentum_[direction]
                    *                             std::sqrt( -std::log( 1.0-Rand::uniform1() ) );

                // change of momentum in the direction(s) along the reflection plane
                if (nDim>1) {
                    double sign_rnd = Rand::uniform() - 0.5;
                    sign_rnd = ( sign_rnd )/std::abs( sign_rnd );
                    momentumRefl_2D[ ipart ] = sign_rnd * species->thermal_momentum_[(direction+1)%nDim]
                        *                             userFunctions::erfinv( Rand::uniform1() );
                    if (nDim>2) {
                        double sign_rnd = Rand::uniform() - 0.5;
                        sign_rnd = ( sign_rnd )/std::abs( sign_rnd );
                        momentumRefl_3D[ ipart ] = sign_rnd * species->thermal_momentum_[(direction+2)%nDim]
                            *                             userFunctions::erfinv( Rand::uniform1() );
                    }
                }

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
                    gp = sqrt( 1.0 + pow( momentum_x[ipart], 2 )+pow( momentum_y[ipart], 2 )+pow( momentum_z[ipart], 2 ) );
                    px = -gp*g*vx + Lxx * momentum_x[ ipart ] + Lxy * momentum_y[ ipart ] + Lxz * momentum_z[ ipart ];
                    py = -gp*g*vy + Lxy * momentum_x[ ipart ] + Lyy * momentum_y[ ipart ] + Lyz * momentum_z[ ipart ];
                    pz = -gp*g*vz + Lxz * momentum_x[ ipart ] + Lyz * momentum_y[ ipart ] + Lzz * momentum_z[ ipart ];
                    momentum_x[ ipart ] = px;
                    momentum_y[ ipart ] = py;
                    momentum_z[ ipart ] = pz;

                }//ENDif vel != 0

            } else {                                    // IF VELOCITY < 3*THERMAL SIMPLY REFLECT IT
                momentum[ ipart ] = -momentum[ ipart ];

            }// endif on v vs. thermal_velocity_

            // position of the particle after reflection
            position[ ipart ] = 2.*wall_position - position[ ipart ];

            // energy lost during thermalization
            LorentzFactor = sqrt( 1.+pow( momentum_x[ipart], 2 )+pow( momentum_y[ipart], 2 )+pow( momentum_z[ipart], 2 ) );
            nrj_iPart += weight[ ipart ]*( initial_energy - LorentzFactor+1.0 );


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
    }
}
