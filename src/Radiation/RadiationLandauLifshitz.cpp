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
RadiationLandauLifshitz::RadiationLandauLifshitz(Params& params,
                                                 Species * species)
      : Radiation(params, species)
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
// -----------------------------------------------------------------------------
void RadiationLandauLifshitz::operator() (
        Particles &particles,
        Species * photon_species,
        SmileiMPI* smpi,
        RadiationTables &RadiationTables,
        int istart,
        int iend,
        int ithread, int ipart_ref)
{

    // _______________________________________________________________
    // Parameters
    std::vector<double> *Epart = &(smpi->dynamics_Epart[ithread]);
    std::vector<double> *Bpart = &(smpi->dynamics_Bpart[ithread]);
    //std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);

    int nparts = Epart->size()/3;
    double* Ex = &( (*Epart)[0*nparts] );
    double* Ey = &( (*Epart)[1*nparts] );
    double* Ez = &( (*Epart)[2*nparts] );
    double* Bx = &( (*Bpart)[0*nparts] );
    double* By = &( (*Bpart)[1*nparts] );
    double* Bz = &( (*Bpart)[2*nparts] );

    // Charge divided by the square of the mass
    double charge_over_mass2;

    // 1/mass^2
    const double one_over_mass_2 = pow(one_over_mass_,2.);

    // Temporary quantum parameter
    double chipa;

    // Temporary Lorentz factor
    double gamma;

    // Temporary double parameter
    double temp;

    // Momentum shortcut
    double* momentum[3];
    for ( int i = 0 ; i<3 ; i++ )
        momentum[i] =  &( particles.momentum(i,0) );

    // Charge shortcut
    short* charge = &( particles.charge(0) );

    // Weight shortcut
    double* weight = &( particles.weight(0) );

    // Optical depth for the Monte-Carlo process
    // double* chi = &( particles.chi(0));

    // Local vector to store the radiated energy
    std::vector <double> rad_norm_energy (iend-istart,0);

    // Reinitialize the cumulative radiated energy for the current thread
    this->radiated_energy = 0.;

    // _______________________________________________________________
    // Computation

    #pragma omp simd
    for (int ipart=istart ; ipart<iend; ipart++ ) {
        charge_over_mass2 = (double)(charge[ipart])*one_over_mass_2;

        // Gamma
        gamma = sqrt(1.0 + momentum[0][ipart]*momentum[0][ipart]
                             + momentum[1][ipart]*momentum[1][ipart]
                             + momentum[2][ipart]*momentum[2][ipart]);

        // Computation of the Lorentz invariant quantum parameter
        chipa = Radiation::compute_chipa(charge_over_mass2,
                     momentum[0][ipart],momentum[1][ipart],momentum[2][ipart],
                     gamma,
                     (*(Ex+ipart-ipart_ref)),(*(Ey+ipart-ipart_ref)),(*(Ez+ipart-ipart_ref)),
                     (*(Bx+ipart-ipart_ref)),(*(By+ipart-ipart_ref)),(*(Bz+ipart-ipart_ref)) );

        // Effect on the momentum
        if (chipa >= RadiationTables.get_chipa_radiation_threshold())
        {

            // Radiated energy during the time step
            temp =
            RadiationTables.get_classical_cont_rad_energy(chipa,dt);

            // Effect on the momentum
            // Temporary factor
            temp *= gamma/(gamma*gamma-1.);

            // Update of the momentum
            momentum[0][ipart] -= temp*momentum[0][ipart];
            momentum[1][ipart] -= temp*momentum[1][ipart];
            momentum[2][ipart] -= temp*momentum[2][ipart];

            // Exact energy loss due to the radiation
            rad_norm_energy[ipart - istart] = gamma - sqrt(1.0
                                         + momentum[0][ipart]*momentum[0][ipart]
                                         + momentum[1][ipart]*momentum[1][ipart]
                                         + momentum[2][ipart]*momentum[2][ipart]);
        }
    }

    // _______________________________________________________________
    // Computation of the thread radiated energy

    double radiated_energy_loc = 0;

    #pragma omp simd reduction(+:radiated_energy_loc)
    for (int ipart=istart ; ipart<iend; ipart++ )
    {
        radiated_energy_loc += weight[ipart]*rad_norm_energy[ipart - istart] ;
    }
    radiated_energy += radiated_energy_loc;
}
