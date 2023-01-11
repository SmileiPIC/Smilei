// ----------------------------------------------------------------------------
//! \file RadiationCorrLandauLifshitz.cpp
//
//! \brief This class is for the continuous radiation loss.
//!        This model includes a quantum correction.
//
//! The implementation is adapted from the thesis results of M. Lobet
//! See http://www.theses.fr/2015BORD0361
//
// ----------------------------------------------------------------------------

#include "RadiationCorrLandauLifshitz.h"

#include <cstring>
#include <fstream>

#include <cmath>

// ---------------------------------------------------------------------------------------------------------------------
//! Constructor for RadiationCorrLandauLifshitz
//! Inherited from Radiation
// ---------------------------------------------------------------------------------------------------------------------
RadiationCorrLandauLifshitz::RadiationCorrLandauLifshitz( Params &params, Species *species, Random * rand )
    : Radiation( params, species, rand )
{
}

// ---------------------------------------------------------------------------------------------------------------------
//! Destructor for RadiationCorrLandauLifshitz
// ---------------------------------------------------------------------------------------------------------------------
RadiationCorrLandauLifshitz::~RadiationCorrLandauLifshitz()
{
}

// ---------------------------------------------------------------------------------------------------------------------
//! Overloading of the operator (): perform the Discontinuous radiation reaction
//! induced by the nonlinear inverse Compton scattering
//
//! \param particles   particle object containing the particle properties
//! \param photons     Particles object that will receive emitted photons
//! \param smpi        MPI properties
//! \param RadiationTables Cross-section data tables and useful functions
//                     for nonlinear inverse Compton scattering
//! \param istart      Index of the first particle
//! \param iend        Index of the last particle
//! \param ithread     Thread index
//! \param radiated_energy     overall energy radiated during the call to this method
// ---------------------------------------------------------------------------------------------------------------------
void RadiationCorrLandauLifshitz::operator()(
    Particles       &particles,
    Particles       */*photons*/,
    SmileiMPI       *smpi,
    RadiationTables &radiation_tables,
    double          &radiated_energy,
    int             istart,
    int             iend,
    int             ithread,
    int             /*ibin*/,
    int             ipart_ref)
{

    // _______________________________________________________________
    // Parameters
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    //std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);

    const int nparts = smpi->getBufferSize(ithread);
    const double *const __restrict__ Ex = &( ( *Epart )[0*nparts] );
    const double *const __restrict__ Ey = &( ( *Epart )[1*nparts] );
    const double *const __restrict__ Ez = &( ( *Epart )[2*nparts] );
    const double *const __restrict__ Bx = &( ( *Bpart )[0*nparts] );
    const double *const __restrict__ By = &( ( *Bpart )[1*nparts] );
    const double *const __restrict__ Bz = &( ( *Bpart )[2*nparts] );

    // 1/mass^2
    const double one_over_mass_square = one_over_mass_*one_over_mass_;

    // Minimum value of chi for the radiation
    const double minimum_chi_continuous = radiation_tables.getMinimumChiContinuous();

    // Momentum shortcut
    double *const __restrict__ momentum_x = particles.getPtrMomentum(0);
    double *const __restrict__ momentum_y = particles.getPtrMomentum(1);
    double *const __restrict__ momentum_z = particles.getPtrMomentum(2);

    // Charge shortcut
    const short *const __restrict__ charge = particles.getPtrCharge();

    // Weight shortcut
    const double *const __restrict__ weight = particles.getPtrWeight();

    // Optical depth for the Monte-Carlo process
    double *const __restrict__ chi = particles.getPtrChi();

    // Local vector to store the radiated energy
    double * rad_norm_energy = new double [iend-istart];
    // double * rad_norm_energy = (double*) aligned_alloc(64, (iend-istart)*sizeof(double));
    #pragma omp simd
    for( int ipart=0 ; ipart<iend-istart; ipart++ ) {
        rad_norm_energy[ipart] = 0;
    }

    double radiated_energy_loc = 0;

    // _______________________________________________________________
    // Computation

    #pragma omp simd
    for( int ipart=istart ; ipart<iend; ipart++ ) {
        
        // Charge divided by the square of the mass
        const double charge_over_mass_square = ( double )( charge[ipart] )*one_over_mass_square;

        // Temporary Lorentz factor
        const double gamma = sqrt( 1.0 + momentum_x[ipart]*momentum_x[ipart]
                      + momentum_y[ipart]*momentum_y[ipart]
                      + momentum_z[ipart]*momentum_z[ipart] );

        // Computation of the Lorentz invariant quantum parameter
        const double particle_chi = Radiation::computeParticleChi( charge_over_mass_square,
                       momentum_x[ipart], momentum_y[ipart], momentum_z[ipart],
                       gamma,
                       Ex[ipart-ipart_ref], Ey[ipart-ipart_ref], Ez[ipart-ipart_ref] ,
                       Bx[ipart-ipart_ref], By[ipart-ipart_ref], Bz[ipart-ipart_ref] );

        // Effect on the momentum
        // (Should be vectorized with masked instructions)
        if( gamma > 1.1 && particle_chi >= minimum_chi_continuous ) {

            // Radiated energy during the time step
            const double temp =
                radiation_tables.getRidgersCorrectedRadiatedEnergy( particle_chi, dt_ ) * gamma/( gamma*gamma - 1 );

            // Update of the momentum
            momentum_x[ipart] -= temp*momentum_x[ipart];
            momentum_y[ipart] -= temp*momentum_y[ipart];
            momentum_z[ipart] -= temp*momentum_z[ipart];

            // Exact energy loss due to the radiation
            rad_norm_energy[ipart - istart] = gamma - sqrt( 1.0
                                              + momentum_x[ipart]*momentum_x[ipart]
                                              + momentum_y[ipart]*momentum_y[ipart]
                                              + momentum_z[ipart]*momentum_z[ipart] );

        }
    }

    // _______________________________________________________________
    // Computation of the thread radiated energy

    #pragma omp simd reduction(+:radiated_energy_loc)
    for( int ipart=0 ; ipart<iend-istart; ipart++ ) {
        radiated_energy_loc += weight[ipart]*rad_norm_energy[ipart] ;
    }

    // _______________________________________________________________
    // Update of the quantum parameter

    #pragma omp simd
    for( int ipart=istart ; ipart<iend; ipart++ ) {
        const double charge_over_mass_square = ( double )( charge[ipart] )*one_over_mass_square;

        // Gamma
        const double gamma = sqrt( 1.0 + momentum_x[ipart]*momentum_x[ipart]
                      + momentum_y[ipart]*momentum_y[ipart]
                      + momentum_z[ipart]*momentum_z[ipart] );

        // Computation of the Lorentz invariant quantum parameter
        chi[ipart] = Radiation::computeParticleChi( charge_over_mass_square,
                       momentum_x[ipart], momentum_y[ipart], momentum_z[ipart],
                       gamma,
                       Ex[ipart-ipart_ref], Ey[ipart-ipart_ref], Ez[ipart-ipart_ref],
                       Bx[ipart-ipart_ref], By[ipart-ipart_ref], Bz[ipart-ipart_ref] );

    }

    // Add the local energy to the patch one
    radiated_energy += radiated_energy_loc;

    // _______________________________________________________________
    // Cleaning

    // free(rad_norm_energy);
    delete [] rad_norm_energy;

}
