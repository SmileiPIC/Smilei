/*! @file BoundaryConditionType.cpp

 @brief BoundaryConditionType.cpp for particle boundary conditions

 */

#include <cmath>
#include <cstdlib>

#include <string>
#include <algorithm>
#include <cctype>

#include "BoundaryConditionType.h"
#include "Params.h"
#include "userFunctions.h"


void internal_inf( Species *species, int imin, int imax, int direction, double limit_inf, double /*dt*/, std::vector<double> &/*invgf*/, Random * /*rand*/, double &energy_change )
{
    energy_change = 0.;     // no energy loss during exchange
    const double* const position  = species->particles->getPtrPosition( direction );
    int* const          cell_keys = species->particles->getPtrCellKeys();
#if defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel deviceptr(position,cell_keys)
    #pragma acc loop gang worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target is_device_ptr( /* to: */                           \
                                      position /* [imin:imax - imin] */ ) \
        is_device_ptr( /* tofrom: */                                      \
                       cell_keys /* [imin:imax - imin] */ )
    #pragma omp teams distribute parallel for
#endif
    for( int ipart=imin ; ipart<imax ; ipart++ ) {
        if( cell_keys[ ipart ] >= 0 && position[ ipart ] < limit_inf ) {
            cell_keys[ ipart ] = -2 - 2 * direction;
        }
    }
}

void internal_sup( Species *species, int imin, int imax, int direction, double limit_sup, double /*dt*/, std::vector<double> &/*invgf*/, Random * /*rand*/, double &energy_change )
{
    energy_change = 0.;     // no energy loss during exchange
    const double* const position  = species->particles->getPtrPosition( direction );
    int* const          cell_keys = species->particles->getPtrCellKeys();
#if defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel deviceptr(position,cell_keys)
    #pragma acc loop gang worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target is_device_ptr( /* to: */                           \
                                      position /* [imin:imax - imin] */ ) \
        is_device_ptr( /* tofrom: */                                      \
                       cell_keys /* [imin:imax - imin] */ )
    #pragma omp teams distribute parallel for
#endif
    for( int ipart=imin ; ipart<imax ; ipart++ ) {
        if( cell_keys[ ipart ] >= 0 && position[ ipart ] >= limit_sup ) {
            cell_keys[ ipart ] = -3 - 2 * direction;
        }
    }
}

void internal_inf_AM( Species *species, int imin, int imax, int /*direction*/, double limit_inf, double /*dt*/, std::vector<double> &/*invgf*/, Random * /*rand*/, double &energy_change )
{
    energy_change = 0.;     // no energy loss during exchange
    double* position_y = species->particles->getPtrPosition(1);
    double* position_z = species->particles->getPtrPosition(2);
    int* cell_keys = species->particles->getPtrCellKeys();
    double limit_inf2 = limit_inf*limit_inf;
    for( int ipart=imin ; ipart<imax ; ipart++ ) {
        double distance2ToAxis = position_y[ipart]*position_y[ipart]+position_z[ipart]*position_z[ipart];
        if( cell_keys[ ipart ] >= 0 && distance2ToAxis < limit_inf2 ) {
            cell_keys[ ipart ] = -4;
        }
    }
}

void internal_sup_AM( Species *species, int imin, int imax, int /*direction*/, double limit_sup, double /*dt*/, std::vector<double> &/*invgf*/, Random * /*rand*/, double &energy_change )
{
    energy_change = 0.;     // no energy loss during exchange
    double* position_y = species->particles->getPtrPosition(1);
    double* position_z = species->particles->getPtrPosition(2);
    int* cell_keys = species->particles->getPtrCellKeys();
    double limit_sup2 = limit_sup*limit_sup;
    for( int ipart=imin ; ipart<imax ; ipart++ ) {
        double distance2ToAxis = position_y[ipart]*position_y[ipart]+position_z[ipart]*position_z[ipart];
        if( cell_keys[ ipart ] >= 0 && distance2ToAxis >= limit_sup2 ) {
            cell_keys[ ipart ] = -5;
        }
    }
}

void reflect_particle_inf( Species *species, int imin, int imax, int direction, double limit_inf, double /*dt*/, std::vector<double> &/*invgf*/, Random * /*rand*/, double &energy_change )
{
    energy_change = 0.;     // no energy loss during reflection
    double* position = species->particles->getPtrPosition(direction);
    double* momentum = species->particles->getPtrMomentum(direction);
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc parallel deviceptr(position,momentum)
    #pragma acc loop gang worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target is_device_ptr( position, momentum )
    #pragma omp teams distribute parallel for
#endif
    for( int ipart=imin ; ipart<imax ; ipart++ ) {
        if( position[ ipart ] < limit_inf ) {
            position[ ipart ] = 2.*limit_inf - position[ ipart ];
            momentum[ ipart ] = -momentum[ ipart ];
        }
    }
}

void reflect_particle_sup( Species *species, int imin, int imax, int direction, double limit_sup, double /*dt*/, std::vector<double> &/*invgf*/, Random * /*rand*/, double &energy_change )
{
    energy_change = 0.;     // no energy loss during reflection
    double* position = species->particles->getPtrPosition(direction);
    double* momentum = species->particles->getPtrMomentum(direction);
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc parallel deviceptr(position,momentum)
    #pragma acc loop gang worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target is_device_ptr( position, momentum )
    #pragma omp teams distribute parallel for
#endif
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] >= limit_sup) {
            // The upper boundary does not belong to the domain so the reflection is done just before it.
            position[ ipart ] = 2*std::nextafter(limit_sup, 0) - position[ ipart ];
            momentum[ ipart ] = -momentum[ ipart ];
        }
    }
}

void reflect_particle_wall( Species *species, int imin, int imax, int direction, double wall_position, double dt, std::vector<double> &invgf, Random * /*rand*/, double &energy_change )
{
    energy_change = 0.;     // no energy loss during reflection
    double* position = species->particles->getPtrPosition(direction);
    double* momentum = species->particles->getPtrMomentum(direction);
    double* invgf_p  = invgf.data();
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc parallel deviceptr(position,momentum,invgf_p)
    #pragma acc loop gang worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target is_device_ptr(position,momentum,invgf_p)
    #pragma omp teams distribute parallel for
#endif
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        double particle_position     = position[ipart];
        double particle_position_old = particle_position - dt*invgf_p[ipart]*momentum[ipart]; 
        if ( ( wall_position - particle_position_old ) * ( wall_position - particle_position )<0) {
            position[ ipart ] = 2 * wall_position - position[ ipart ];
            momentum[ ipart ] = - momentum[ ipart ];
        }
    }
}

// direction not used below, direction is "r"
void refl_particle_AM( Species *species, int imin, int imax, int /*direction*/, double limit_sup, double /*dt*/, std::vector<double> &invgf, Random * /*rand*/, double &energy_change )
{
    energy_change = 0.;     // no energy loss during reflection
    
    double* position_y = species->particles->getPtrPosition(1);
    double* position_z = species->particles->getPtrPosition(2);
    double* momentum_y = species->particles->getPtrMomentum(1);
    double* momentum_z = species->particles->getPtrMomentum(2);
    //We look for the coordinate of the point at which the particle crossed the boundary
    //We need to fine the parameter t at which (y - vy*t)^2+(z-vz*t)^2 = Rmax^2
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        double r2 = position_y[ipart]*position_y[ipart]+position_z[ipart]*position_z[ipart];
        if ( r2 >= limit_sup*limit_sup ) {
            //solving v2*t2 + b*t + r2-rmax^2 = 0
            double b = -2.*( position_y[ ipart ]*momentum_y[ ipart ] + position_z[ ipart ]*momentum_z[ ipart ] )*invgf[ipart];
            double p2 = ( momentum_y[ ipart ]*momentum_y[ ipart ] + momentum_z[ ipart ]*momentum_z[ ipart ] );
            double v2 = p2 * invgf[ipart]*invgf[ipart];
            double delta = b*b - 4.*v2*( r2 - limit_sup*limit_sup );
            
            //delta is neceseraliy >=0 otherwise there is no solution which means that something unsual happened
            if( delta < 0 ) {
                ERROR( "There are no solution to reflexion. This should never happen" );
            }
    
            double t = ( -b - sqrt( delta ) )/v2*0.5;
    
            double y0, z0; //Coordinates of the crossing point 0
            y0 =  position_y[ ipart ] - momentum_y[ ipart ]*t*invgf[ipart] ;
            z0 =  position_z[ ipart ] - momentum_z[ ipart ]*t*invgf[ipart] ;
   
            //Update momentum after reflexion at crossing point
            //pr = p.er
            //pt = p.et
            //er = (y0/limit_sup, z0/limit_sup)
            //et = (-z0/limit_sup, y0/limit_sup)
            double pr = ( momentum_y[ ipart ]*y0 + momentum_z[ ipart ]*z0) / limit_sup;
            double pt = (-momentum_y[ ipart ]*z0 + momentum_z[ ipart ]*y0) / limit_sup;
           
            // New p = -pr.er + pt.et 
            momentum_y[ipart] = (-pr*y0 - pt*z0)/limit_sup;
            momentum_z[ipart] = (-pr*z0 + pt*y0)/limit_sup;
 
            //Update new particle position as a movement from (y0,z0) with new momentum during time t:
            // The upper boundary point (y0,z0) does not belong to the simulation domain so the reflection is done on the point just before towards the center of the domain.
            position_y[ ipart ] = std::nextafter(y0, 0.) + momentum_y[ipart]*invgf[ipart]*t ;
            position_z[ ipart ] = std::nextafter(z0, 0.) + momentum_z[ipart]*invgf[ipart]*t ;
        }
    }    
}

void remove_particle_inf( Species* species, 
                          int imin, int imax, 
                          int direction, 
                          double limit_inf, 
                          double /*dt*/, 
                          std::vector<double>& /*invgf*/, 
                          Random* /*rand*/, 
                          double& energy_change )
{

    double change_in_energy = 0.0;

    const double *const position   = species->particles->getPtrPosition( direction );
    const double *const momentum_x = species->particles->getPtrMomentum( 0 );
    const double *const momentum_y = species->particles->getPtrMomentum( 1 );
    const double *const momentum_z = species->particles->getPtrMomentum( 2 );
    short  *const charge     = species->particles->getPtrCharge();
    const double *const weight     = species->particles->getPtrWeight();
    int    *const cell_keys  = species->particles->getPtrCellKeys();

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target           is_device_ptr( position, momentum_x, momentum_y, momentum_z, charge, weight, cell_keys ) map( tofrom \
                                                                                                                               : change_in_energy )
    #pragma omp teams distribute parallel for reduction( + \
                                                         : change_in_energy )
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel deviceptr(position,momentum_x,momentum_y,momentum_z,weight,charge,cell_keys)
    #pragma acc loop gang worker vector reduction(+ : change_in_energy)
#else
    #pragma omp simd reduction(+ : change_in_energy)
#endif
    for( int ipart = imin; ipart < imax; ++ipart ) {
        if( position[ipart] < limit_inf ) {
            const double LorentzFactor = std::sqrt( 1.0 +
                                                    momentum_x[ipart] * momentum_x[ipart] +
                                                    momentum_y[ipart] * momentum_y[ipart] +
                                                    momentum_z[ipart] * momentum_z[ipart] );
            change_in_energy += weight[ipart] * ( LorentzFactor - 1.0 ); // energy lost REDUCTION
            charge[ipart]    = 0;
            cell_keys[ipart] = -1;
        }
    }

    energy_change = change_in_energy;
}

void remove_particle_sup( Species* species, 
                          int imin, int imax, 
                          int direction, 
                          double limit_sup, 
                          double /*dt*/, 
                          std::vector<double>& /*invgf*/, 
                          Random* /*rand*/, 
                          double& energy_change )
{

    double change_in_energy = 0.0;

    const double *const position   = species->particles->getPtrPosition( direction );
    const double *const momentum_x = species->particles->getPtrMomentum( 0 );
    const double *const momentum_y = species->particles->getPtrMomentum( 1 );
    const double *const momentum_z = species->particles->getPtrMomentum( 2 );
    short  *const charge     = species->particles->getPtrCharge();
    const double *const weight     = species->particles->getPtrWeight();
    int    *const cell_keys  = species->particles->getPtrCellKeys();

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target           is_device_ptr( position, momentum_x, momentum_y, momentum_z, charge, weight, cell_keys ) map( tofrom \
                                                                                                                               : change_in_energy )
    #pragma omp teams distribute parallel for reduction( + \
                                                         : change_in_energy )
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel deviceptr(position,momentum_x,momentum_y,momentum_z,weight,charge,cell_keys)
    #pragma acc loop gang worker vector reduction(+ : change_in_energy)
#else
    #pragma omp simd reduction(+ : change_in_energy)
#endif
    for( int ipart = imin; ipart < imax; ++ipart ) {
        if( position[ipart] >= limit_sup ) {
            const double LorentzFactor = std::sqrt( 1.0 +
                                                    momentum_x[ipart] * momentum_x[ipart] +
                                                    momentum_y[ipart] * momentum_y[ipart] +
                                                    momentum_z[ipart] * momentum_z[ipart] );
            change_in_energy += weight[ipart] * ( LorentzFactor - 1.0 ); // energy lost REDUCTION
            charge[ipart]    = 0;
            cell_keys[ipart] = -1;
        }
    }

    energy_change = change_in_energy;
}

void remove_particle_wall( Species *species, int imin, int imax, int direction, double wall_position, double dt, std::vector<double> &invgf, Random * /*rand*/, double &energy_change )
{
    energy_change = 0.;
    double* position = species->particles->getPtrPosition(direction);
    double* momentum = species->particles->getPtrMomentum(direction);
    double* momentum_x = species->particles->getPtrMomentum(0);
    double* momentum_y = species->particles->getPtrMomentum(1);
    double* momentum_z = species->particles->getPtrMomentum(2);
    short* charge    = species->particles->getPtrCharge();
    double* weight   = species->particles->getPtrWeight();
    int* cell_keys   = species->particles->getPtrCellKeys();
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        double particle_position     = position[ipart];
        double particle_position_old = particle_position - dt*invgf[ipart]*momentum[ipart]; 
        if ( ( wall_position-particle_position_old )*( wall_position-particle_position )<0) {
            double LorentzFactor = sqrt( 1. + momentum_x[ipart] * momentum_x[ipart] + momentum_y[ipart] * momentum_y[ipart] + momentum_z[ipart] * momentum_z[ipart] );
            energy_change += weight[ ipart ]*( LorentzFactor-1.0 ); // energy lost REDUCTION
            charge[ ipart ] = 0;
            cell_keys[ipart] = -1;
        }
    }
}

void remove_particle_AM( Species *species, int imin, int imax, int /*direction*/, double limit_sup, double /*dt*/, std::vector<double> &/*invgf*/, Random * /*rand*/, double &energy_change )
{
    energy_change = 0.;
    double* position_y = species->particles->getPtrPosition(1);
    double* position_z = species->particles->getPtrPosition(2);
    double* momentum_x = species->particles->getPtrMomentum(0);
    double* momentum_y = species->particles->getPtrMomentum(1);
    double* momentum_z = species->particles->getPtrMomentum(2);
    short* charge    = species->particles->getPtrCharge();
    double* weight   = species->particles->getPtrWeight();
    int* cell_keys   = species->particles->getPtrCellKeys();
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        double distance2ToAxis = position_y[ipart]*position_y[ipart]+position_z[ipart]*position_z[ipart];
        if ( distance2ToAxis >= limit_sup*limit_sup ) {
            double LorentzFactor = sqrt( 1. + momentum_x[ipart] * momentum_x[ipart] + momentum_y[ipart] * momentum_y[ipart] + momentum_z[ipart] * momentum_z[ipart] );
            energy_change += weight[ ipart ]*( LorentzFactor-1.0 ); // energy lost REDUCTION
            charge[ ipart ] = 0;
            cell_keys[ipart] = -1;
        }
    }
}

//! Delete photon (mass_==0) at the boundary and keep the energy for diagnostics
void remove_photon_inf( Species *species, int imin, int imax, int direction, double limit_inf, double /*dt*/, std::vector<double> &/*invgf*/, Random * /*rand*/, double &energy_change )
{
    energy_change = 0.;
    double* position = species->particles->getPtrPosition(direction);
    double* momentum_x = species->particles->getPtrMomentum(0);
    double* momentum_y = species->particles->getPtrMomentum(1);
    double* momentum_z = species->particles->getPtrMomentum(2);
    double* weight     = species->particles->getPtrWeight();
    short* charge    = species->particles->getPtrCharge();
    int* cell_keys   = species->particles->getPtrCellKeys();
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] < limit_inf) {
            double momentumNorm = sqrt( momentum_x[ipart] * momentum_x[ipart] + momentum_y[ipart] * momentum_y[ipart] + momentum_z[ipart] * momentum_z[ipart] );
            energy_change += weight[ ipart ]*( momentumNorm ); // energy lost REDUCTION
            charge[ ipart ] = 0;
            cell_keys[ipart] = -1;
        }
    }
}

void remove_photon_sup( Species *species, int imin, int imax, int direction, double limit_sup, double /*dt*/, std::vector<double> &/*invgf*/, Random * /*rand*/, double &energy_change )
{
    energy_change = 0.;
    double* position = species->particles->getPtrPosition(direction);
    double* momentum_x = species->particles->getPtrMomentum(0);
    double* momentum_y = species->particles->getPtrMomentum(1);
    double* momentum_z = species->particles->getPtrMomentum(2);
    double* weight     = species->particles->getPtrWeight();
    short* charge    = species->particles->getPtrCharge();
    int* cell_keys   = species->particles->getPtrCellKeys();
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] >= limit_sup) {
            double momentumNorm = sqrt( momentum_x[ipart] * momentum_x[ipart] + momentum_y[ipart] * momentum_y[ipart] + momentum_z[ipart] * momentum_z[ipart] );
            energy_change += weight[ ipart ]*( momentumNorm ); // energy lost REDUCTION
            charge[ ipart ] = 0;
            cell_keys[ipart] = -1;
        }
    }
}

void stop_particle_inf( Species *species, int imin, int imax, int direction, double limit_inf, double /*dt*/, std::vector<double> &/*invgf*/, Random * /*rand*/, double &energy_change )
{
    double change_in_energy = 0.0;
    double* position = species->particles->getPtrPosition(direction);
    double* momentum_x = species->particles->getPtrMomentum(0);
    double* momentum_y = species->particles->getPtrMomentum(1);
    double* momentum_z = species->particles->getPtrMomentum(2);
    double* weight     = species->particles->getPtrWeight();
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target is_device_ptr( position, momentum_x, momentum_y, momentum_z, weight ) map( tofrom : change_in_energy )
    #pragma omp teams distribute parallel for reduction( + : change_in_energy )
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel deviceptr(position,momentum_x,momentum_y,momentum_z,weight)
    #pragma acc loop gang worker vector reduction(+ : change_in_energy)
#else
    #pragma omp simd reduction(+ : change_in_energy)
#endif
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] < limit_inf) {
            double LorentzFactor = sqrt( 1. + momentum_x[ipart] * momentum_x[ipart] + momentum_y[ipart] * momentum_y[ipart] + momentum_z[ipart] * momentum_z[ipart] );
            change_in_energy += weight[ ipart ]*( LorentzFactor-1.0 ); // energy lost REDUCTION
            position[ ipart ] = 2.*limit_inf - position[ ipart ];
            momentum_x[ ipart ] = 0.;
            momentum_y[ ipart ] = 0.;
            momentum_z[ ipart ] = 0.;
        }
    }
    energy_change = change_in_energy;
}

void stop_particle_sup( Species *species, int imin, int imax, int direction, double limit_sup, double /*dt*/, std::vector<double> &/*invgf*/, Random * /*rand*/, double &energy_change )
{
    double change_in_energy = 0.0;
    double* position = species->particles->getPtrPosition(direction);
    double* momentum_x = species->particles->getPtrMomentum(0);
    double* momentum_y = species->particles->getPtrMomentum(1);
    double* momentum_z = species->particles->getPtrMomentum(2);
    double* weight     = species->particles->getPtrWeight();
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target is_device_ptr( position, momentum_x, momentum_y, momentum_z, weight ) map( tofrom : change_in_energy )
    #pragma omp teams distribute parallel for reduction( + : change_in_energy )
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel deviceptr(position,momentum_x,momentum_y,momentum_z,weight)
    #pragma acc loop gang worker vector reduction(+ : change_in_energy)
#else
    #pragma omp simd reduction(+ : change_in_energy)
#endif
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        if ( position[ ipart ] >= limit_sup) {
            double LorentzFactor = sqrt( 1. + momentum_x[ipart] * momentum_x[ipart] + momentum_y[ipart] * momentum_y[ipart] + momentum_z[ipart] * momentum_z[ipart] );
            change_in_energy += weight[ ipart ]*( LorentzFactor-1.0 ); // energy lost REDUCTION
            position[ ipart ] = 2.*limit_sup - position[ ipart ];
            momentum_x[ ipart ] = 0.;
            momentum_y[ ipart ] = 0.;
            momentum_z[ ipart ] = 0.;
        }
    }
    energy_change = change_in_energy;
}

void stop_particle_wall( Species *species, int imin, int imax, int direction, double wall_position, double dt, std::vector<double> &invgf, Random * /*rand*/, double &energy_change )
{
    energy_change = 0;
    double* position = species->particles->getPtrPosition(direction);
    double* momentum = species->particles->getPtrMomentum(direction);
    double* momentum_x = species->particles->getPtrMomentum(0);
    double* momentum_y = species->particles->getPtrMomentum(1);
    double* momentum_z = species->particles->getPtrMomentum(2);
    double* weight     = species->particles->getPtrWeight();
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        double particle_position     = position[ipart];
        double particle_position_old = particle_position - dt*invgf[ipart]*momentum[ipart];
        if ( ( wall_position-particle_position_old )*( wall_position-particle_position )<0 ) {
            double LorentzFactor = sqrt( 1. + momentum_x[ipart] * momentum_x[ipart] + momentum_y[ipart] * momentum_y[ipart] + momentum_z[ipart] * momentum_z[ipart] );
            energy_change += weight[ ipart ]*( LorentzFactor-1.0 ); // energy lost REDUCTION
            position[ ipart ] = 2.*wall_position - position[ ipart ];
            momentum_x[ ipart ] = 0.;
            momentum_y[ ipart ] = 0.;
            momentum_z[ ipart ] = 0.;
        }
    }
}

void stop_particle_AM( Species *species, int imin, int imax, int /*direction*/, double limit_sup, double /*dt*/, std::vector<double> &/*invgf*/, Random * /*rand*/, double &energy_change )
{
    double* position_y = species->particles->getPtrPosition(1);
    double* position_z = species->particles->getPtrPosition(2);
    double* momentum_x = species->particles->getPtrMomentum(0);
    double* momentum_y = species->particles->getPtrMomentum(1);
    double* momentum_z = species->particles->getPtrMomentum(2);
    double* weight     = species->particles->getPtrWeight();

    energy_change = 0;
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        double distance2ToAxis = position_y[ipart]*position_y[ipart]+position_z[ipart]*position_z[ipart];
        if ( distance2ToAxis >= limit_sup*limit_sup ) {
            double LorentzFactor = sqrt( 1. + momentum_x[ipart] * momentum_x[ipart] + momentum_y[ipart] * momentum_y[ipart] + momentum_z[ipart] * momentum_z[ipart] );
            energy_change += weight[ ipart ]*( LorentzFactor-1.0 ); // energy lost REDUCTION
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

void thermalize_particle_inf( Species *species, int imin, int imax, int direction, double limit_inf, double /*dt*/, std::vector<double> &/*invgf*/, Random * rand, double &energy_change )
{
    int nDim = species->nDim_particle;
    double* position = species->particles->getPtrPosition(direction);
    double* momentum = species->particles->getPtrMomentum(direction);
    double* momentumRefl_2D = species->particles->getPtrMomentum((direction+1)%nDim);
    double* momentumRefl_3D = species->particles->getPtrMomentum((direction+2)%nDim);
    double* momentum_x = species->particles->getPtrMomentum(0);
    double* momentum_y = species->particles->getPtrMomentum(1);
    double* momentum_z = species->particles->getPtrMomentum(2);
    double* weight     = species->particles->getPtrWeight();
#if defined( SMILEI_ACCELERATOR_GPU ) 
    uint32_t xorshift32_state = rand->xorshift32_state;
#endif
    double change_in_energy = 0.0;
    double thermal_momentum = species->thermal_momentum_[direction];
    double thermal_momentum1;
    double thermal_momentum2;
    double v0 = species->thermal_velocity_[0];
    if (nDim>1) {
        thermal_momentum1 = species->thermal_momentum_[(direction+1)%nDim];
        if (nDim>2) {
            thermal_momentum2 = species->thermal_momentum_[(direction+2)%nDim];
        }
    }
    double vx, vy, vz, v2, g, gm1, Lxx, Lyy, Lzz, Lxy, Lxz, Lyz;
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
    }
#if defined( SMILEI_ACCELERATOR_GPU) // GPU
        const int nchunks = (imax-imin)/32 + 1 ;
    #if defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp target is_device_ptr( position, momentum, momentumRefl_2D, momentumRefl_3D, momentum_x, momentum_y, momentum_z, weight ) map( tofrom : change_in_energy )
        #pragma omp teams distribute thread_limit(32) reduction( + : change_in_energy )
    #elif defined( SMILEI_ACCELERATOR_GPU_OACC )
        #pragma acc parallel loop gang vector_length(32)  reduction(+ : change_in_energy) independent deviceptr(position, momentum, momentumRefl_2D, momentumRefl_3D,momentum_x,momentum_y,momentum_z,weight)
    #endif
        for (int ichunk = 0 ; ichunk < nchunks ; ++ichunk ) {
                int chunk_size = (ichunk==nchunks-1) ? (imax-imin)%32 : 32;
                uint32_t xorshift32_state_local = xorshift32_state + ichunk;
                uint32_t xorshift32_state_array[32];
            #if defined( SMILEI_ACCELERATOR_GPU_OACC )
                #pragma acc loop seq
            #elif defined( SMILEI_ACCELERATOR_GPU_OMP )
                //#pragma omp single // does not work with rocm
            #endif        
                // Fill  xorshift32_state_array[...] with local state
                for( int i = 0; i < chunk_size; ++i ){
                    xorshift32_state_array[i] = Random_namespace::xorshift32(xorshift32_state_local);
                }
            #if defined( SMILEI_ACCELERATOR_GPU_OACC )
                #pragma acc loop vector
            #elif defined( SMILEI_ACCELERATOR_GPU_OMP )
                #pragma omp parallel for 
            #endif
                for( int i = 0; i < chunk_size ; ++i ){
                    int ipart = imin + ichunk * 32 + i;
#else //CPU
        #pragma omp simd reduction(+ : change_in_energy)
        for (int ipart = imin ; ipart < imax ; ++ipart ) {
#endif            
            if ( position[ ipart ] < limit_inf) {
                // checking the particle's velocity compared to the thermal one
                double p2 = momentum_x[ipart] * momentum_x[ipart] + momentum_y[ipart] * momentum_y[ipart] + momentum_z[ipart] * momentum_z[ipart];
                double LorentzFactor = sqrt( 1.+p2 );
                double v = sqrt( p2 )/LorentzFactor;

                // energy before thermalization
                double initial_energy = LorentzFactor - 1.0;
                // Apply bcs depending on the particle velocity
                // --------------------------------------------
                if( v > 3.0 * v0) {     //IF VELOCITY > 3*THERMAL VELOCITY THEN THERMALIZE IT

                        // velocity of the particle after thermalization/reflection
                        // change of velocity in the direction normal to the reflection plane
                        double sign_vel = -momentum[ ipart ]/std::abs( momentum[ ipart ] );
                    #if defined( SMILEI_ACCELERATOR_GPU ) 
                        momentum[ ipart ] = sign_vel * thermal_momentum * std::sqrt( -std::log( 1.0 - Random_namespace::uniform1(xorshift32_state_array[i]) ) );
                    #else
                        momentum[ ipart ] = sign_vel * thermal_momentum * std::sqrt( -std::log( 1.0 - rand->uniform1() ) );
                    #endif

                        // change of momentum in the direction(s) along the reflection plane
                        if (nDim>1) {
                            #if defined( SMILEI_ACCELERATOR_GPU ) 
                                momentumRefl_2D[ ipart ] = thermal_momentum1 * Random_namespace::perp_rand_dp(xorshift32_state_array[i]);
                            #else
                                momentumRefl_2D[ ipart ] = thermal_momentum1 * perp_rand( rand );
                            #endif
                            if (nDim>2) {
                                #if defined( SMILEI_ACCELERATOR_GPU ) 
                                    momentumRefl_3D[ ipart ] = thermal_momentum2 * Random_namespace::perp_rand_dp(xorshift32_state_array[i]);
                                #else
                                    momentumRefl_3D[ ipart ] = thermal_momentum2 * perp_rand( rand );
                                #endif
                            }
                        }  
                        // Adding the mean velocity (using relativistic composition)
                        double gp, px, py, pz;
                        if( v2>0. ) {
                            // Lorentz transformation of the momentum
                            gp = sqrt( 1.0 + momentum_x[ipart] * momentum_x[ipart] + momentum_y[ipart] * momentum_y[ipart] + momentum_z[ipart] * momentum_z[ipart] );
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
                LorentzFactor = sqrt( 1. + momentum_x[ipart] * momentum_x[ipart] + momentum_y[ipart] * momentum_y[ipart] + momentum_z[ipart] * momentum_z[ipart] );
                change_in_energy += weight[ ipart ] * ( initial_energy - LorentzFactor + 1.0 );


                // HERE IS AN ATTEMPT TO INTRODUCE A SPACE DEPENDENCE ON THE BCs
                // double val_min(params.dens_profile.vacuum_length[1]), val_max(params.dens_profile.vacuum_length[1]+params.dens_profile.length_params_y[0]);

                //if ( ( species->particles->position(1,ipart) >= val_min ) && ( species->particles->position(1,ipart) <= val_max ) ) {
                // nrj computed during diagnostics
                //species->particles->position(direction, ipart) = limit_pos - species->particles->position(direction, ipart);
                //species->particles->momentum(direction, ipart) = sqrt(params.thermal_velocity_[direction]) * tabFcts.erfinv( rand->uniform() );
                //}
                //else {
                //stop_particle( species->particles, ipart, direction, limit_pos, params, energy_change );
                //}
                
            }
        }
#if defined( SMILEI_ACCELERATOR_GPU ) 
        } //End for loop on chunks.
        xorshift32_state += 32;
        rand->xorshift32_state = xorshift32_state;
#endif
        energy_change = change_in_energy;

}

void thermalize_particle_sup( Species *species, int imin, int imax, int direction, double limit_sup, double /*dt*/, std::vector<double> &/*invgf*/, Random * rand, double &energy_change )
{
    int nDim = species->nDim_particle;
    double* position = species->particles->getPtrPosition(direction);
    double* momentum = species->particles->getPtrMomentum(direction);
    double* momentumRefl_2D = species->particles->getPtrMomentum((direction+1)%nDim);
    double* momentumRefl_3D = species->particles->getPtrMomentum((direction+2)%nDim);
    double* momentum_x = species->particles->getPtrMomentum(0);
    double* momentum_y = species->particles->getPtrMomentum(1);
    double* momentum_z = species->particles->getPtrMomentum(2);
    double* weight     = species->particles->getPtrWeight();
#if defined( SMILEI_ACCELERATOR_GPU ) 
    uint32_t xorshift32_state = rand->xorshift32_state;
#endif
    double change_in_energy = 0.0;
    double thermal_momentum = species->thermal_momentum_[direction];
    double thermal_momentum1;
    double thermal_momentum2;
    double v0 = species->thermal_velocity_[0];
    if (nDim>1) {
        thermal_momentum1 = species->thermal_momentum_[(direction+1)%nDim];
        if (nDim>2) {
            thermal_momentum2 = species->thermal_momentum_[(direction+2)%nDim];
        }
    }
    double vx, vy, vz, v2, g, gm1, Lxx, Lyy, Lzz, Lxy, Lxz, Lyz;
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
    }

#if defined( SMILEI_ACCELERATOR_GPU) // GPU
    const int nchunks = (imax-imin)/32 + 1 ;
    #if defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp target is_device_ptr( position, momentum, momentumRefl_2D, momentumRefl_3D, momentum_x, momentum_y, momentum_z, weight ) map( tofrom : change_in_energy )
        #pragma omp teams distribute thread_limit(32) reduction( + : change_in_energy )
    #elif defined( SMILEI_ACCELERATOR_GPU_OACC )
        #pragma acc parallel loop gang vector_length(32)  reduction(+ : change_in_energy) independent deviceptr(position, momentum, momentumRefl_2D, momentumRefl_3D,momentum_x,momentum_y,momentum_z,weight)
    #endif
        for (int ichunk = 0 ; ichunk < nchunks ; ++ichunk ) {
                int chunk_size = (ichunk==nchunks-1) ? (imax-imin)%32 : 32;
                uint32_t xorshift32_state_local = xorshift32_state + ichunk;
                uint32_t xorshift32_state_array[32];
            #if defined( SMILEI_ACCELERATOR_GPU_OACC )
                #pragma acc loop seq
            #elif defined( SMILEI_ACCELERATOR_GPU_OMP )
            //#pragma omp single // does not work with rocm
            #endif        
            // Fill  xorshift32_state_array[...] with local state
                for( int i = 0; i < chunk_size; ++i ){
                    xorshift32_state_array[i] = Random_namespace::xorshift32(xorshift32_state_local);
                }
            #if defined( SMILEI_ACCELERATOR_GPU_OACC )
                #pragma acc loop vector
            #elif defined( SMILEI_ACCELERATOR_GPU_OMP )
                #pragma omp parallel for 
            #endif
                for( int i = 0; i < chunk_size ; ++i ){
                    int ipart = imin + ichunk * 32 + i;
#else //CPU
        #pragma omp simd reduction(+ : change_in_energy)
        for (int ipart = imin ; ipart < imax ; ++ipart ) {
#endif            
            if ( position[ ipart ] >= limit_sup) {
                // checking the particle's velocity compared to the thermal one
                double p2 = momentum_x[ipart] * momentum_x[ipart] + momentum_y[ipart] * momentum_y[ipart] + momentum_z[ipart] * momentum_z[ipart];
                double LorentzFactor = sqrt( 1.+p2 );
                double v = sqrt( p2 )/LorentzFactor;

                // energy before thermalization
                double initial_energy = LorentzFactor - 1.0;

                // Apply bcs depending on the particle velocity
                // --------------------------------------------
                if( v > 3.0 * v0 ) {     //IF VELOCITY > 3*THERMAL VELOCITY THEN THERMALIZE IT

                    // velocity of the particle after thermalization/reflection
                    // change of velocity in the direction normal to the reflection plane
                    double sign_vel = -momentum[ ipart ]/std::abs( momentum[ ipart ] );
                    #if defined( SMILEI_ACCELERATOR_GPU ) 
                        momentum[ ipart ] = sign_vel * thermal_momentum * std::sqrt( -std::log( 1.0 - Random_namespace::uniform1(xorshift32_state_array[i]) ) );
                    #else
                        momentum[ ipart ] = sign_vel * thermal_momentum * std::sqrt( -std::log( 1.0 - rand->uniform1() ) );
                    #endif

                    // change of momentum in the direction(s) along the reflection plane
                        if (nDim>1) {
                            #if defined( SMILEI_ACCELERATOR_GPU ) 
                                momentumRefl_2D[ ipart ] = thermal_momentum1 * Random_namespace::perp_rand_dp(xorshift32_state_array[i]);
                            #else
                                momentumRefl_2D[ ipart ] = thermal_momentum1 * perp_rand( rand );
                            #endif
                            if (nDim>2) {
                                #if defined( SMILEI_ACCELERATOR_GPU ) 
                                    momentumRefl_3D[ ipart ] = thermal_momentum2 * Random_namespace::perp_rand_dp(xorshift32_state_array[i]);
                                #else
                                    momentumRefl_3D[ ipart ] = thermal_momentum2 * perp_rand( rand );
                                #endif
                            }
                        }

                    // Adding the mean velocity (using relativistic composition)
                    double gp, px, py, pz;
                    // mean-velocity
                    if( v2>0. ) {
                        // Lorentz transformation of the momentum
                        gp = sqrt( 1.0 + momentum_x[ipart] * momentum_x[ipart] + momentum_y[ipart] * momentum_y[ipart] + momentum_z[ipart] * momentum_z[ipart] );
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

                // The upper boundary does not belong to the domain so the reflection is done just before it.
                position[ ipart ] = 2*std::nextafter(limit_sup, 0) - position[ ipart ];

                // energy lost during thermalization
                LorentzFactor = sqrt( 1. + momentum_x[ipart] * momentum_x[ipart] + momentum_y[ipart] * momentum_y[ipart] + momentum_z[ipart] * momentum_z[ipart] );
                change_in_energy += weight[ ipart ] * ( initial_energy - LorentzFactor + 1.0 );

                /* HERE IS AN ATTEMPT TO INTRODUCE A SPACE DEPENDENCE ON THE BCs
                // double val_min(params.dens_profile.vacuum_length[1]), val_max(params.dens_profile.vacuum_length[1]+params.dens_profile.length_params_y[0]);

                if ( ( species->particles->position(1,ipart) >= val_min ) && ( species->particles->position(1,ipart) <= val_max ) ) {
                // nrj computed during diagnostics
                species->particles->position(direction, ipart) = limit_pos - species->particles->position(direction, ipart);
                species->particles->momentum(direction, ipart) = sqrt(params.thermal_velocity_[direction]) * tabFcts.erfinv( rand->uniform() );
                }
                else {
                stop_particle( species->particles, ipart, direction, limit_pos, params, energy_change );
                }
                */
            }
        }
#if defined( SMILEI_ACCELERATOR_GPU ) 
        } //End for loop on chunks.
        xorshift32_state += 32;
        rand->xorshift32_state = xorshift32_state;
#endif
    energy_change = change_in_energy;
}


void thermalize_particle_wall( Species *species, int imin, int imax, int direction, double wall_position, double dt, std::vector<double> &invgf, Random * rand, double &energy_change )
{
    int nDim = species->nDim_particle;
    double* position = species->particles->getPtrPosition(direction);
    double* momentum = species->particles->getPtrMomentum(direction);
    double* momentumRefl_2D = species->particles->getPtrMomentum((direction+1)%nDim);
    double* momentumRefl_3D = species->particles->getPtrMomentum((direction+2)%nDim);
    double* momentum_x = species->particles->getPtrMomentum(0);
    double* momentum_y = species->particles->getPtrMomentum(1);
    double* momentum_z = species->particles->getPtrMomentum(2);
    double* weight     = species->particles->getPtrWeight();

    energy_change = 0;
    for (int ipart=imin ; ipart<imax ; ipart++ ) {
        double particle_position     = position[ipart];
        double particle_position_old = particle_position - dt*invgf[ipart]*species->particles->Momentum[direction][ipart];
        if ( ( wall_position-particle_position_old )*( wall_position-particle_position )<0 ) {
            // checking the particle's velocity compared to the thermal one
            double p2 = momentum_x[ipart] * momentum_x[ipart] + momentum_y[ipart] * momentum_y[ipart] + momentum_z[ipart] * momentum_z[ipart];
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
                    *                             std::sqrt( -std::log( 1.0-rand->uniform1() ) );

                // change of momentum in the direction(s) along the reflection plane
                if (nDim>1) {
                    momentumRefl_2D[ ipart ] = species->thermal_momentum_[(direction+1)%nDim] * perp_rand( rand );
                    if (nDim>2) {
                        momentumRefl_3D[ ipart ] = species->thermal_momentum_[(direction+2)%nDim] * perp_rand( rand );
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
                    gp = sqrt( 1.0 + momentum_x[ipart] * momentum_x[ipart] + momentum_y[ipart] * momentum_y[ipart] + momentum_z[ipart] * momentum_z[ipart] );
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
            LorentzFactor = sqrt( 1. + momentum_x[ipart] * momentum_x[ipart] + momentum_y[ipart] * momentum_y[ipart] + momentum_z[ipart] * momentum_z[ipart] );
            energy_change += weight[ ipart ]*( initial_energy - LorentzFactor+1.0 );


            /* HERE IS AN ATTEMPT TO INTRODUCE A SPACE DEPENDENCE ON THE BCs
            // double val_min(params.dens_profile.vacuum_length[1]), val_max(params.dens_profile.vacuum_length[1]+params.dens_profile.length_params_y[0]);

            if ( ( species->particles->position(1,ipart) >= val_min ) && ( species->particles->position(1,ipart) <= val_max ) ) {
            // nrj computed during diagnostics
            species->particles->position(direction, ipart) = limit_pos - species->particles->position(direction, ipart);
            species->particles->momentum(direction, ipart) = sqrt(params.thermal_velocity_[direction]) * tabFcts.erfinv( rand->uniform() );
            }
            else {
            stop_particle( species->particles, ipart, direction, limit_pos, params, energy_change );
            }
            */
        }
    }
}


// ==== angle_threshold particle boundary condition ====

namespace {

// Trim whitespace
inline std::string strip_spaces( std::string s )
{
    s.erase(
        std::remove_if(
            s.begin(), s.end(),
            []( unsigned char c ) { return std::isspace( c ); }
        ),
        s.end()
    );
    return s;
}

// Parse boundary condition string:
//   "angle_threshold:3deg"  -> 3 degrees in radians
//   "angle_threshold:0.05"  -> 0.05 radians (default unit)
//   "angle_threshold:0.05rad" -> 0.05 radians
// If no parameter is provided ("angle_threshold"), defaults to pi/2 (≈ always thermalize).
inline double parse_angle0_rad_from_bc( const std::string &bc )
{
    const double pi = std::acos( -1.0 );
    const double half_pi = 0.5 * pi;

    auto pos = bc.find( ':' );
    if( pos == std::string::npos ) {
        return half_pi;
    }

    std::string s = strip_spaces( bc.substr( pos + 1 ) );
    bool is_deg = false;
    bool is_rad = false;

    if( s.size() >= 3 && s.compare( s.size() - 3, 3, "deg" ) == 0 ) {
        is_deg = true;
        s = s.substr( 0, s.size() - 3 );
    } else if( s.size() >= 3 && s.compare( s.size() - 3, 3, "rad" ) == 0 ) {
        is_rad = true;
        s = s.substr( 0, s.size() - 3 );
    }

    double val = std::atof( s.c_str() );
    if( val < 0.0 ) {
        val = -val;
    }

    if( is_deg ) {
        val *= pi / 180.0;
    } else if( is_rad ) {
        // already radians
    } else {
        // default: radians
    }

    // Clamp to [0, pi/2] for sanity
    if( val < 0.0 ) val = 0.0;
    if( val > half_pi ) val = half_pi;

    return val;
}

// Compute alpha = atan2(|p_perp|, |p_par|), where p_par is the momentum component along "direction"
// and p_perp is the magnitude of the momentum component(s) perpendicular to it.
#ifdef SMILEI_ACCELERATOR_GPU_OACC
#pragma acc routine seq
#endif
#ifdef SMILEI_ACCELERATOR_GPU_OMP
#pragma omp declare target
#endif

inline double compute_alpha( int direction,
                             const double *px, const double *py, const double *pz,
                             int ipart )
{
    // Use full 3-velocity (3V) information even in 1D/2D geometries (e.g., 1D3V, 2D3V).
    // direction is the boundary normal axis (0:x, 1:y, 2:z).
    const double p_par =
        ( direction == 0 ) ? px[ipart] :
        ( direction == 1 ) ? py[ipart] :
                             pz[ipart];

    double p_perp2 = 0.0;
    if( direction != 0 ) p_perp2 += px[ipart]*px[ipart];
    if( direction != 1 ) p_perp2 += py[ipart]*py[ipart];
    if( direction != 2 ) p_perp2 += pz[ipart]*pz[ipart];

    return std::atan2( std::sqrt( p_perp2 ), std::fabs( p_par ) );
}
#ifdef SMILEI_ACCELERATOR_GPU_OMP
#pragma omp end declare target
#endif


} // namespace

// angle_threshold BC (domain lower boundary):
// - if alpha < alpha0 : apply "thermalize" BC
// - else             : apply "reflective" BC
void angle_threshold_particle_inf( Species *species, int imin, int imax, int direction,
                                   double limit_inf, double dt, std::vector<double> &invgf,
                                   Random * rand, double &energy_change )
{
    // Decide alpha0 from the BC string on this side
    const std::string &bc = species->boundary_conditions_[direction][0];
    const double alpha0 = parse_angle0_rad_from_bc( bc );

    const int nDim = species->nDim_particle;

    double *pos_dir = species->particles->getPtrPosition( direction );
    double *p_dir   = species->particles->getPtrMomentum( direction );
    double *px      = species->particles->getPtrMomentum( 0 );
    double *py      = species->particles->getPtrMomentum( 1 );
    double *pz      = species->particles->getPtrMomentum( 2 );

    // First pass: reflect particles that are outside AND alpha >= alpha0
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc parallel deviceptr(pos_dir,p_dir,px,py,pz)
    #pragma acc loop gang worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target is_device_ptr(pos_dir,p_dir,px,py,pz)
    #pragma omp teams distribute parallel for
#endif
    for( int ipart=imin; ipart<imax; ipart++ ) {
        if( pos_dir[ipart] < limit_inf ) {
            const double alpha = compute_alpha( direction, px, py, pz, ipart );
            if( alpha >= alpha0 ) {
                // Reflect
                pos_dir[ipart] = 2.0 * limit_inf - pos_dir[ipart];
                p_dir[ipart]   = -p_dir[ipart];
            }
        }
    }

    // Second pass: remaining outside particles get the standard thermalize BC
    double thermal_energy_change = 0.0;
    thermalize_particle_inf( species, imin, imax, direction, limit_inf, dt, invgf, rand, thermal_energy_change );
    energy_change = thermal_energy_change;
}

// angle_threshold BC (domain upper boundary):
// - if alpha < alpha0 : apply "thermalize" BC
// - else             : apply "reflective" BC
void angle_threshold_particle_sup( Species *species, int imin, int imax, int direction,
                                   double limit_sup, double dt, std::vector<double> &invgf,
                                   Random * rand, double &energy_change )
{
    // Decide alpha0 from the BC string on this side
    const std::string &bc = species->boundary_conditions_[direction][1];
    const double alpha0 = parse_angle0_rad_from_bc( bc );

    const int nDim = species->nDim_particle;

    double *pos_dir = species->particles->getPtrPosition( direction );
    double *p_dir   = species->particles->getPtrMomentum( direction );
    double *px      = species->particles->getPtrMomentum( 0 );
    double *py      = species->particles->getPtrMomentum( 1 );
    double *pz      = species->particles->getPtrMomentum( 2 );

    // First pass: reflect particles that are outside AND alpha >= alpha0
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc parallel deviceptr(pos_dir,p_dir,px,py,pz)
    #pragma acc loop gang worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target is_device_ptr(pos_dir,p_dir,px,py,pz)
    #pragma omp teams distribute parallel for
#endif
    for( int ipart=imin; ipart<imax; ipart++ ) {
        if( pos_dir[ipart] >= limit_sup ) {
            const double alpha = compute_alpha( direction, px, py, pz, ipart );
            if( alpha >= alpha0 ) {
                // Reflect (upper boundary is outside the domain -> reflect just before it)
                pos_dir[ipart] = 2.0 * std::nextafter( limit_sup, 0.0 ) - pos_dir[ipart];
                p_dir[ipart]   = -p_dir[ipart];
            }
        }
    }

    // Second pass: remaining outside particles get the standard thermalize BC
    double thermal_energy_change = 0.0;
    thermalize_particle_sup( species, imin, imax, direction, limit_sup, dt, invgf, rand, thermal_energy_change );
    energy_change = thermal_energy_change;
}

// angle_threshold BC for internal walls.
// Implemented for completeness: it behaves like thermalize_particle_wall for alpha < alpha0,
// and like reflect_particle_wall for alpha >= alpha0.
void angle_threshold_particle_wall( Species *species, int imin, int imax, int direction,
                                    double wall_position, double dt, std::vector<double> &invgf,
                                    Random * rand, double &energy_change )
{
    const int nDim = species->nDim_particle;

    // Try to find an "angle_threshold:..." string on either side; if absent, default pi/2.
    const std::string &bc0 = species->boundary_conditions_[direction][0];
    const std::string &bc1 = species->boundary_conditions_[direction][1];
    double alpha0 = parse_angle0_rad_from_bc( bc0 );
    if( bc0.rfind( "angle_threshold", 0 ) != 0 ) {
        alpha0 = parse_angle0_rad_from_bc( bc1 );
    }

    double *pos_dir = species->particles->getPtrPosition( direction );
    double *p_dir   = species->particles->getPtrMomentum( direction );

    double *px      = species->particles->getPtrMomentum( 0 );
    double *py      = species->particles->getPtrMomentum( 1 );
    double *pz      = species->particles->getPtrMomentum( 2 );
    double *weight  = species->particles->getPtrWeight();

    energy_change = 0.0;

    for( int ipart=imin; ipart<imax; ipart++ ) {

        const double particle_position     = pos_dir[ipart];
        const double particle_position_old = particle_position - dt * invgf[ipart] * species->particles->Momentum[direction][ipart];

        // crossing test (same as reflect_particle_wall / thermalize_particle_wall)
        if( ( wall_position - particle_position_old ) * ( wall_position - particle_position ) < 0.0 ) {

            const double alpha = compute_alpha( direction, px, py, pz, ipart );

            if( alpha >= alpha0 ) {
                // Pure reflection (no energy loss)
                pos_dir[ipart] = 2.0 * wall_position - pos_dir[ipart];
                p_dir[ipart]   = -p_dir[ipart];
                continue;
            }

            // Otherwise: use the existing thermalize wall behavior (copied from thermalize_particle_wall),
            // which may thermalize or reflect depending on v compared to 3*v_th.
            // --- begin adapted block ---
            double p2 = px[ipart] * px[ipart] + py[ipart] * py[ipart] + pz[ipart] * pz[ipart];
            double LorentzFactor = std::sqrt( 1.0 + p2 );
            double v = std::sqrt( p2 ) / LorentzFactor;

            double initial_energy = LorentzFactor - 1.0;

            if( v > 3.0 * species->thermal_velocity_[0] ) {

                double sign_vel = -1.0;
                if( std::abs( p_dir[ipart] ) > 0.0 ) {
                    sign_vel = -p_dir[ipart] / std::abs( p_dir[ipart] );
                }
                p_dir[ipart] = sign_vel * species->thermal_momentum_[direction]
                    * std::sqrt( -std::log( 1.0 - rand->uniform1() ) );

                // tangential component(s)
                if( nDim > 1 ) {
                    // choose transverse index as in original code
                    int i2 = (direction + 1) % nDim;
                    double *p_t1 = species->particles->getPtrMomentum( i2 );
                    p_t1[ipart] = species->thermal_momentum_[i2] * perp_rand( rand );

                    if( nDim > 2 ) {
                        int i3 = (direction + 2) % nDim;
                        double *p_t2 = species->particles->getPtrMomentum( i3 );
                        p_t2[ipart] = species->thermal_momentum_[i3] * perp_rand( rand );
                    }
                }

            } else {
                // low velocity: reflect
                p_dir[ipart] = -p_dir[ipart];
            }

            pos_dir[ipart] = 2.0 * wall_position - pos_dir[ipart];

            LorentzFactor = std::sqrt( 1.0 + px[ipart] * px[ipart] + py[ipart] * py[ipart] + pz[ipart] * pz[ipart] );
            energy_change += weight[ipart] * ( initial_energy - LorentzFactor + 1.0 );
            // --- end adapted block ---
        }
    }
}
