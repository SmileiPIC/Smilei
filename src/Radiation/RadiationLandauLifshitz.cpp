// ----------------------------------------------------------------------------
//! \file RadiationLandauLifshitz.cpp
//
//! \brief This class is for the classical continuous radiation loss with
//!        the Landau-Lifshitz model.
//!        This model does not include a quantum correction.
//
//! References:
//! L. D. Landau and E. M. Lifshitz, The classical theory of fields, 1947
//! F. Niel et al., 2017
//
// ----------------------------------------------------------------------------

#include "RadiationLandauLifshitz.h"

// -----------------------------------------------------------------------------
//! Constructor for RadiationNLandauLifshitz
//! Inherited from Radiation
// -----------------------------------------------------------------------------
RadiationLandauLifshitz::RadiationLandauLifshitz( Params &params,
        Species *species, Random * rand  )
    : Radiation( params, species, rand )
{
}

// -----------------------------------------------------------------------------
//! Destructor for RadiationLandauLifshitz
// -----------------------------------------------------------------------------
RadiationLandauLifshitz::~RadiationLandauLifshitz()
{
}

// -----------------------------------------------------------------------------
//! Overloading of the operator (): perform the Landau-Lifshitz classical
//! radiation reaction
//
//! \param particles   particle object containing the particle properties
//! \param smpi        MPI properties
//! \param RadiationTables Cross-section data tables and useful functions
//                     for nonlinear inverse Compton scattering
//! \param istart      Index of the first particle
//! \param iend        Index of the last particle
//! \param ithread     Thread index
//! \param radiated_energy     overall energy radiated during the call to this method
// -----------------------------------------------------------------------------
void RadiationLandauLifshitz::operator()(
    Particles       &particles,
    Species         *photon_species,
    SmileiMPI       *smpi,
    RadiationTables &RadiationTables,
    double          &radiated_energy,
    int istart,
    int iend,
    int ithread,
    int ipart_ref )
{

    // _______________________________________________________________
    // Parameters
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    //std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);

    int nparts = Epart->size()/3;
    double * __restrict__ Ex = &( ( *Epart )[0*nparts] );
    double * __restrict__ Ey = &( ( *Epart )[1*nparts] );
    double * __restrict__ Ez = &( ( *Epart )[2*nparts] );
    double * __restrict__ Bx = &( ( *Bpart )[0*nparts] );
    double * __restrict__ By = &( ( *Bpart )[1*nparts] );
    double * __restrict__ Bz = &( ( *Bpart )[2*nparts] );

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
    double * __restrict__ momentum_x = particles.getPtrMomentum(0);
    double * __restrict__ momentum_y = particles.getPtrMomentum(1);
    double * __restrict__ momentum_z = particles.getPtrMomentum(2);

    // Charge shortcut
    short * __restrict__ charge = particles.getPtrCharge();

    // Weight shortcut
    double * __restrict__ weight = particles.getPtrWeight();

    // Optical depth for the Monte-Carlo process
    double * __restrict__ chi = particles.getPtrChi();

    // cumulative Radiated energy from istart to iend
    double radiated_energy_loc = 0;

#ifndef _GPU
    // Local vector to store the radiated energy

    // double * rad_norm_energy = new double [iend-istart];
    double  * rad_norm_energy = (double*) aligned_alloc(64, (iend-istart)*sizeof(double));
    #pragma omp simd
    for( int ipart=0 ; ipart<iend-istart; ipart++ ) {
        rad_norm_energy[ipart] = 0;
    }
#else
    double rad_norm_energy = 0;
#endif

    // _______________________________________________________________
    // Computation

    #ifndef _GPU
        #pragma omp simd aligned(rad_norm_energy:64)
    #else
        int np = iend-istart;
        #pragma acc parallel \
            present(Ex[istart:np],Ey[istart:np],Ez[istart:np],\
            Bx[istart:np],By[istart:np],Bz[istart:np]) \
            deviceptr(momentum_x,momentum_y,momentum_z,charge,weight,chi) \
            reduction(+:radiated_energy_loc)
    {
        #pragma acc loop reduction(+:radiated_energy_loc) gang worker vector private(rad_norm_energy)
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
                       Ex[ipart-ipart_ref], Ey[ipart-ipart_ref], Ez[ipart-ipart_ref] ,
                       Bx[ipart-ipart_ref], By[ipart-ipart_ref], Bz[ipart-ipart_ref] );

        // Effect on the momentum
        if( (gamma>1.) && (particle_chi >= RadiationTables.getMinimumChiContinuous()) ) {

            // Radiated energy during the time step
            temp =
                RadiationTables.getClassicalRadiatedEnergy( particle_chi, dt_ );

            // Effect on the momentum
            // Temporary factor
            temp *= gamma/( gamma*gamma-1. );

            // Update of the momentum
            momentum_x[ipart] -= temp*momentum_x[ipart];
            momentum_y[ipart] -= temp*momentum_y[ipart];
            momentum_z[ipart] -= temp*momentum_z[ipart];
                                              
    // _______________________________________________________________
    // Computation of the thread radiated energy
                                              
#ifndef _GPU

            // Exact energy loss due to the radiation
            rad_norm_energy[ipart] = gamma - std::sqrt( 1.0
                                              + momentum_x[ipart]*momentum_x[ipart]
                                              + momentum_y[ipart]*momentum_y[ipart]
                                              + momentum_z[ipart]*momentum_z[ipart] );
        }
    }

    #pragma omp simd reduction(+:radiated_energy_loc)
    for( int ipart=0 ; ipart<iend-istart; ipart++ ) {
        radiated_energy_loc += weight[ipart]*rad_norm_energy[ipart] ;
    }
#else
            radiated_energy_loc += weight[ipart]*(gamma - std::sqrt( 1.0
                                              + momentum_x[ipart]*momentum_x[ipart]
                                              + momentum_y[ipart]*momentum_y[ipart]
                                              + momentum_z[ipart]*momentum_z[ipart] );
#endif


    // _______________________________________________________________
    // Update of the quantum parameter

#ifndef _GPU
    #pragma omp simd private(gamma)
    for( int ipart=istart ; ipart<iend; ipart++ ) {
#endif

        charge_over_mass_square = ( double )( charge[ipart] )*one_over_mass_square;

        // Gamma
        gamma = sqrt( 1.0 + momentum_x[ipart]*momentum_x[ipart]
                      + momentum_y[ipart]*momentum_y[ipart]
                      + momentum_z[ipart]*momentum_z[ipart] );

        // Computation of the Lorentz invariant quantum parameter
        chi[ipart] = computeParticleChi( charge_over_mass_square,
                     momentum_x[ipart], momentum_y[ipart], momentum_z[ipart],
                     gamma,
                     Ex[ipart-ipart_ref], Ey[ipart-ipart_ref], Ez[ipart-ipart_ref],
                     Bx[ipart-ipart_ref], By[ipart-ipart_ref], Bz[ipart-ipart_ref] );

    #ifndef _GPU
    } // end loop ipart
    #else
            } // end if
        } // end loop ipart
    } // end acc parallel
    #endif

    radiated_energy += radiated_energy_loc;

    #ifndef _GPU
        // _______________________________________________________________
        // Cleaning

        delete [] rad_norm_energy;
    #endif

}
            
