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
#include "userFunctions.h"
#include "Random.h"

inline double perp_rand( Random * rand ) {
    double a = userFunctions::erfinv_dp( rand->uniform1() ); 
    if( rand->cointoss() ) { 
        a *= -1.;
    }
    return a;
}

void internal_inf( Species *species, int imin, int imax, int direction, double limit_inf, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );

void internal_sup( Species *species, int imin, int imax, int direction, double limit_sup, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );

void internal_inf_AM( Species *species, int imin, int imax, int direction, double limit_inf, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );

void internal_sup_AM( Species *species, int imin, int imax, int direction, double limit_sup, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );

void reflect_particle_inf( Species *species, int imin, int imax, int direction, double limit_inf, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );

void reflect_particle_sup( Species *species, int imin, int imax, int direction, double limit_sup, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );

void reflect_particle_wall( Species *species, int imin, int imax, int direction, double limit_sup, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );

// direction not used below, direction is "r"
void refl_particle_AM( Species *species, int imin, int imax, int direction, double limit_sup, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );

void remove_particle_inf( Species *species, int imin, int imax, int direction, double limit_inf, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );

void remove_particle_sup( Species *species, int imin, int imax, int direction, double limit_sup, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );

void remove_particle_wall( Species *species, int imin, int imax, int direction, double limit_sup, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );

void remove_particle_AM( Species *species, int imin, int imax, int direction, double limit_sup, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );

//! Delete photon (mass_==0) at the boundary and keep the energy for diagnostics
void remove_photon_inf( Species *species, int imin, int imax, int direction, double limit_inf, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );

void remove_photon_sup( Species *species, int imin, int imax, int direction, double limit_sup, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );

void stop_particle_inf( Species *species, int imin, int imax, int direction, double limit_inf, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );

void stop_particle_sup( Species *species, int imin, int imax, int direction, double limit_sup, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );

void stop_particle_wall( Species *species, int imin, int imax, int direction, double limit_sup, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );

void stop_particle_AM( Species *species, int imin, int imax, int direction, double limit_pos, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );

//!\todo (MG) at the moment the particle is thermalize whether or not there is a plasma initially at the boundary.
// ATTENTION: here the thermalization assumes a Maxwellian distribution, maybe we should add some checks on thermal_boundary_temperature (MG)!
void thermalize_particle_inf( Species *species, int imin, int imax, int direction, double limit_pos, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );

void thermalize_particle_sup( Species *species, int imin, int imax, int direction, double limit_pos, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );

void thermalize_particle_wall( Species *species, int imin, int imax, int direction, double limit_pos, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );


#endif
