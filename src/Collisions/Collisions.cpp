
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <fstream>

#include "Collisions.h"

using namespace std;


// Constructor
Collisions::Collisions(
    Params &params,
    double coulomb_log,
    double coulomb_log_factor
) :
    coulomb_log_( coulomb_log ),
    coulomb_log_factor_( coulomb_log_factor ),
    coeff1_( 4.046650232e-21*params.reference_angular_frequency_SI ), // h*omega/(2*me*c^2)
    coeff2_( 2.817940327e-15*params.reference_angular_frequency_SI/299792458. ), // re omega / c
    coeff3_( coeff2_ * coulomb_log_factor_ ),
    coeff4_( 1. / cbrt( 3.*coeff2_ ) )
{
}

Collisions::Collisions() :
    coulomb_log_( -1. ),
    coulomb_log_factor_( 0 ),
    coeff1_( 0 ),
    coeff2_( 0 ),
    coeff3_( 0 ),
    coeff4_( 0 )
{
}

void Collisions::prepare()
{
    npairs_tot_  = 0.;
    smean_       = 0.;
    logLmean_    = 0.;
    #if defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc data copyin(npairs_tot_, smean_, logLmean_)
    #elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target update to(npairs_tot_, smean_, logLmean_)
    #endif
}

void Collisions::finish( Params &, Patch *, std::vector<Diagnostic *> &, bool, std::vector<unsigned int>, std::vector<unsigned int>, int )
{
    #if defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc data copyout(npairs_tot_, smean_, logLmean_)
    #elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target update from(npairs_tot_, smean_, logLmean_)
    #endif
    if( npairs_tot_>0. ) {
        smean_    /= npairs_tot_;
        logLmean_ /= npairs_tot_;
    }
}



