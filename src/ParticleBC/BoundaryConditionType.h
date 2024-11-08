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
#include "Random.h"

inline double perp_rand( Random * rand ) {
    double a = userFunctions::erfinv( rand->uniform1() ); // to be switched to erfinv 3
    if( rand->cointoss() ) { 
        a *= -1.;
    }
    return a;
}

/**
 * copied from erfinv_DP_1.cu by Prof. Mike Giles.
 * https://people.maths.ox.ac.uk/gilesm/
 * https://people.maths.ox.ac.uk/gilesm/codes/erfinv/
 *
 * Original code is written for CUDA.
 * Mutsuo Saito modified original code for C++.
 */
inline double erfinv_v3(double x)
{
    double w, p;
    double sign;
    if (x > 0) {
        sign = 1.0;
    } else {
        sign = -1.0;
        x = abs(x);
    }
    w = - log((1.0-x)*(1.0+x));

    if ( w < 6.250000 ) {
        w = w - 3.125000;
        p =  -3.6444120640178196996e-21;
        p =   -1.685059138182016589e-19 + p*w;
        p =   1.2858480715256400167e-18 + p*w;
        p =    1.115787767802518096e-17 + p*w;
        p =   -1.333171662854620906e-16 + p*w;
        p =   2.0972767875968561637e-17 + p*w;
        p =   6.6376381343583238325e-15 + p*w;
        p =  -4.0545662729752068639e-14 + p*w;
        p =  -8.1519341976054721522e-14 + p*w;
        p =   2.6335093153082322977e-12 + p*w;
        p =  -1.2975133253453532498e-11 + p*w;
        p =  -5.4154120542946279317e-11 + p*w;
        p =    1.051212273321532285e-09 + p*w;
        p =  -4.1126339803469836976e-09 + p*w;
        p =  -2.9070369957882005086e-08 + p*w;
        p =   4.2347877827932403518e-07 + p*w;
        p =  -1.3654692000834678645e-06 + p*w;
        p =  -1.3882523362786468719e-05 + p*w;
        p =    0.0001867342080340571352 + p*w;
        p =  -0.00074070253416626697512 + p*w;
        p =   -0.0060336708714301490533 + p*w;
        p =      0.24015818242558961693 + p*w;
        p =       1.6536545626831027356 + p*w;
    }
    else if ( w < 16.000000 ) {
        w = sqrt(w) - 3.250000;
        p =   2.2137376921775787049e-09;
        p =   9.0756561938885390979e-08 + p*w;
        p =  -2.7517406297064545428e-07 + p*w;
        p =   1.8239629214389227755e-08 + p*w;
        p =   1.5027403968909827627e-06 + p*w;
        p =   -4.013867526981545969e-06 + p*w;
        p =   2.9234449089955446044e-06 + p*w;
        p =   1.2475304481671778723e-05 + p*w;
        p =  -4.7318229009055733981e-05 + p*w;
        p =   6.8284851459573175448e-05 + p*w;
        p =   2.4031110387097893999e-05 + p*w;
        p =   -0.0003550375203628474796 + p*w;
        p =   0.00095328937973738049703 + p*w;
        p =   -0.0016882755560235047313 + p*w;
        p =    0.0024914420961078508066 + p*w;
        p =   -0.0037512085075692412107 + p*w;
        p =     0.005370914553590063617 + p*w;
        p =       1.0052589676941592334 + p*w;
        p =       3.0838856104922207635 + p*w;
    }
    else {
        w = sqrt(w) - 5.000000;
        p =  -2.7109920616438573243e-11;
        p =  -2.5556418169965252055e-10 + p*w;
        p =   1.5076572693500548083e-09 + p*w;
        p =  -3.7894654401267369937e-09 + p*w;
        p =   7.6157012080783393804e-09 + p*w;
        p =  -1.4960026627149240478e-08 + p*w;
        p =   2.9147953450901080826e-08 + p*w;
        p =  -6.7711997758452339498e-08 + p*w;
        p =   2.2900482228026654717e-07 + p*w;
        p =  -9.9298272942317002539e-07 + p*w;
        p =   4.5260625972231537039e-06 + p*w;
        p =  -1.9681778105531670567e-05 + p*w;
        p =   7.5995277030017761139e-05 + p*w;
        p =  -0.00021503011930044477347 + p*w;
        p =  -0.00013871931833623122026 + p*w;
        p =       1.0103004648645343977 + p*w;
        p =       4.8499064014085844221 + p*w;
    }
    return sign * p * x;
}

inline uint32_t xorshift32(uint32_t xorshift32_state)
{
    // Algorithm "xor" from p. 4 of Marsaglia, "Xorshift RNGs" 
    xorshift32_state ^= xorshift32_state << 13;
    xorshift32_state ^= xorshift32_state >> 17;
    xorshift32_state ^= xorshift32_state << 5;
    return xorshift32_state;
}

static constexpr double xorshift32_invmax1 = (1.-1e-11)/4294967296.;

inline double uniform1(uint32_t xorshift32_state) {
    return xorshift32(xorshift32_state) * xorshift32_invmax1;
}

inline double perp_rand_gpu_v3(uint32_t xorshift32_state) {
    double a = erfinv_v3( uniform1(xorshift32_state) ); //userFunctions::
    // technically we could also use the erfinv() function fron cuda, it would require compiling with -cuda though ...
    // the study showed the gap in perf for BC thermal was not worth the added depend
    if( xorshift32(xorshift32_state) & 1  ) { // rand->cointoss()
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
