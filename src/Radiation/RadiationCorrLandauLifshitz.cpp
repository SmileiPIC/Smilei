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
//! \param photon_species species that will receive emitted photons
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
    Species         *photon_species,
    SmileiMPI       *smpi,
    RadiationTables &RadiationTables,
    double          &radiated_energy,
    int             istart,
    int             iend,
    int             ithread,
    int             ipart_ref)
{

    // _______________________________________________________________
    // Parameters
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    //std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);

    int nparts = Epart->size()/3;
    double *Ex = &( ( *Epart )[0*nparts] );
    double *Ey = &( ( *Epart )[1*nparts] );
    double *Ez = &( ( *Epart )[2*nparts] );
    double *Bx = &( ( *Bpart )[0*nparts] );
    double *By = &( ( *Bpart )[1*nparts] );
    double *Bz = &( ( *Bpart )[2*nparts] );

    // Charge divided by the square of the mass
    double charge_over_mass_square;

    // 1/mass^2
    const double one_over_mass_square = one_over_mass_*one_over_mass_;

    // Temporary quantum parameter
    double particle_chi;

    // Temporary Lorentz factor
    double gamma;

    // Temporary double parameter
    double temp;

    // Momentum shortcut
    double* momentum_x = particles.getPtrMomentum(0);
    double* momentum_y = particles.getPtrMomentum(1);
    double* momentum_z = particles.getPtrMomentum(2);

    // Charge shortcut
    short *charge = particles.getPtrCharge();

    // Weight shortcut
    double *weight = particles.getPtrWeight();

    // Optical depth for the Monte-Carlo process
    double* chi = particles.getPtrChi();

    // Local vector to store the radiated energy
    double * rad_norm_energy = new double [iend-istart];
    #pragma omp simd
    for( int ipart=0 ; ipart<iend-istart; ipart++ ) {
        rad_norm_energy[ipart] = 0;
    }

    double radiated_energy_loc = 0;

    // _______________________________________________________________
    // Computation

    #ifndef _GPU
        #pragma omp simd
    #else
        int np = iend-istart;
        #pragma acc parallel \
            create(rad_norm_energy[0:np]) \
            present(Ex[istart:np],Ey[istart:np],Ez[istart:np],\
            Bx[istart:np],By[istart:np],Bz[istart:np],radiated_energy) \
            deviceptr(momentum_x,momentum_y,momentum_z,charge,weight,chi) \
            reduction(+:radiated_energy_loc)
    {
        #pragma acc loop gang worker vector
    #endif
    for( int ipart=istart ; ipart<iend; ipart++ ) {
        charge_over_mass_square = ( double )( charge[ipart] )*one_over_mass_square;

        // Gamma
        gamma = sqrt( 1.0 + momentum_x[ipart]*momentum_x[ipart]
                      + momentum_y[ipart]*momentum_y[ipart]
                      + momentum_z[ipart]*momentum_z[ipart] );

        // Computation of the Lorentz invariant quantum parameter
        particle_chi = Radiation::computeParticleChi( charge_over_mass_square,
                       momentum_x[ipart], momentum_y[ipart], momentum_z[ipart],
                       gamma,
                       ( *( Ex+ipart-ipart_ref ) ), ( *( Ey+ipart-ipart_ref ) ), ( *( Ez+ipart-ipart_ref ) ),
                       ( *( Bx+ipart-ipart_ref ) ), ( *( By+ipart-ipart_ref ) ), ( *( Bz+ipart-ipart_ref ) ) );

        // Effect on the momentum
        // (Should be vectorized with masked instructions)
        if( (gamma>1.) && (particle_chi >= RadiationTables.getMinimumChiContinuous()) ) {

            // Radiated energy during the time step
            temp =
                RadiationTables.getRidgersCorrectedRadiatedEnergy( particle_chi, dt_ );

            // Temporary factor
            temp *= gamma/( gamma*gamma - 1 );

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

    #ifndef _GPU
        #pragma omp simd reduction(+:radiated_energy_loc)
    #else
        #pragma acc loop reduction(+:radiated_energy_loc) gang worker vector
    #endif
    for( int ipart=0 ; ipart<iend-istart; ipart++ ) {
        radiated_energy_loc += weight[ipart]*rad_norm_energy[ipart] ;
    }

    // _______________________________________________________________
    // Update of the quantum parameter

    #ifndef _GPU
        #pragma omp simd
    #else
        #pragma acc loop gang worker vector
    #endif
    for( int ipart=istart ; ipart<iend; ipart++ ) {
        charge_over_mass_square = ( double )( charge[ipart] )*one_over_mass_square;

        // Gamma
        gamma = sqrt( 1.0 + momentum_x[ipart]*momentum_x[ipart]
                      + momentum_y[ipart]*momentum_y[ipart]
                      + momentum_z[ipart]*momentum_z[ipart] );

        // Computation of the Lorentz invariant quantum parameter
        chi[ipart] = Radiation::computeParticleChi( charge_over_mass_square,
                       momentum_x[ipart], momentum_y[ipart], momentum_z[ipart],
                       gamma,
                       ( *( Ex+ipart-ipart_ref ) ), ( *( Ey+ipart-ipart_ref ) ), ( *( Ez+ipart-ipart_ref ) ),
                       ( *( Bx+ipart-ipart_ref ) ), ( *( By+ipart-ipart_ref ) ), ( *( Bz+ipart-ipart_ref ) ) );

    }

    #ifdef _GPU
    } // end acc parallel
    //#pragma acc update self(radiated_energy)
    #endif

    // Add the local energy to the patch one
    radiated_energy += radiated_energy_loc;

    // _______________________________________________________________
    // Cleaning

    delete [] rad_norm_energy;

}
