// ----------------------------------------------------------------------------
//! \file RadiationNiel.cpp
//
//! \brief This file implements the class RadiationNiel.
//!        This class is for the semi-classical Fokker-Planck model of
//!        synchrotron-like radiation loss developped by Niel et al. as an
//!        extension of the classical Landau-Lifshitz model in the weak quantum
//!        regime.
//!        This model includew a quantum correction + stochastic diffusive
//!        operator.
//
//! \details See these references for more information.
//! F. Niel et al., 2017
//! L. D. Landau and E. M. Lifshitz, The classical theory of fields, 1947
// ----------------------------------------------------------------------------

#include "RadiationNiel.h"

// -----------------------------------------------------------------------------
//! Constructor for RadiationNLL
//! Inherited from Radiation
// -----------------------------------------------------------------------------
RadiationNiel::RadiationNiel(Params& params, Species * species)
      : Radiation(params, species)
{
}

// -----------------------------------------------------------------------------
//! Destructor for RadiationNiel
// -----------------------------------------------------------------------------
RadiationNiel::~RadiationNiel()
{
}

// -----------------------------------------------------------------------------
//! Overloading of the operator (): perform the corrected Landau-Lifshitz
//! classical radiation reaction + stochastic diffusive operator.
//
//! \param particles   particle object containing the particle properties
//! \param smpi        MPI properties
//! \param RadiationTables Cross-section data tables and useful functions
//                     for nonlinear inverse Compton scattering
//! \param istart      Index of the first particle
//! \param iend        Index of the last particle
//! \param ithread     Thread index
// -----------------------------------------------------------------------------
void RadiationNiel::operator() (Particles &particles,
        SmileiMPI* smpi,
        RadiationTables &RadiationTables,
        int istart,
        int iend,
        int ithread)
{

    // _______________________________________________________________
    // Parameters
    std::vector<LocalFields> *Epart = &(smpi->dynamics_Epart[ithread]);
    std::vector<LocalFields> *Bpart = &(smpi->dynamics_Bpart[ithread]);
    //std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);

    // Charge divided by the square of the mass
    double charge_over_mass2;

    // 1/mass^2
    double one_over_mass_2 = pow(one_over_mass_,2.);

    // Sqrt(dt), used intensively in these loops
    double sqrtdt = sqrt(dt);

    // Temporary quantum parameter
    double chipa;

    // Temporary Lorentz factor
    double gamma;

    // Temporary double parameter
    double temp;

    // Radiated energy
    double rad_energy;

    // Stochastic diffusive term fo Niel et al.
    double diffusion;

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

    // Reinitialize the cumulative radiated energy for the current thread
    this->radiated_energy = 0.;

    // _______________________________________________________________
    // Computation

    //#pragma omp simd
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
                     (*Epart)[ipart].x,(*Epart)[ipart].y,(*Epart)[ipart].z,
                     (*Bpart)[ipart].x,(*Bpart)[ipart].y,(*Bpart)[ipart].z);

         // Radiated energy during the time step
         rad_energy =
         RadiationTables.get_corrected_cont_rad_energy_Ridgers(chipa,dt);

        // Below chipa = 1E-3, radiation losses are negligible
        if (chipa > 1E-3)
        {

            // Diffusive stochastic component during dt
            diffusion = RadiationTables.get_Niel_stochastic_term(gamma,
                                                                 chipa,sqrtdt);
        }
        else
        {
            diffusion = 0;
        }

        // Effect on the momentum
        // Temporary factor
        temp = (rad_energy - diffusion)*gamma/(gamma*gamma-1.);

        // Update of the momentum
        momentum[0][ipart] -= temp*momentum[0][ipart];
        momentum[1][ipart] -= temp*momentum[1][ipart];
        momentum[2][ipart] -= temp*momentum[2][ipart];

        // Incrementation of the radiated energy cumulative parameter
        radiated_energy += weight[ipart]*(gamma - sqrt(1.0
                            + momentum[0][ipart]*momentum[0][ipart]
                            + momentum[1][ipart]*momentum[1][ipart]
                            + momentum[2][ipart]*momentum[2][ipart]));

    }

}

// -----------------------------------------------------------------------------
//! Overloading of the operator (): perform the corrected Landau-Lifshitz
//! classical radiation reaction + stochastic diffusive operator.
//
//! \param particles   particle object containing the particle properties
//! \param smpi        MPI properties
//! \param RadiationTables Cross-section data tables and useful functions
//                     for nonlinear inverse Compton scattering
//! \param istart      Index of the first particle
//! \param iend        Index of the last particle
//! \param ithread     Thread index
// -----------------------------------------------------------------------------
/*void RadiationNiel::operator() (Particles &particles,
        SmileiMPI* smpi,
        RadiationTables &RadiationTables,
        int istart,
        int iend,
        int ithread)
{

    // _______________________________________________________________
    // Parameters
    std::vector<LocalFields> *Epart = &(smpi->dynamics_Epart[ithread]);
    std::vector<LocalFields> *Bpart = &(smpi->dynamics_Bpart[ithread]);
    //std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);

    // Charge divided by the square of the mass
    double charge_over_mass2;

    // 1/mass^2
    double one_over_mass_2 = pow(one_over_mass_,2.);

    // Sqrt(dt), used intensively in these loops
    double sqrtdt = sqrt(dt);

    // Number of particles
    int nbparticles = iend-istart;

    // Temporary quantum parameter
    double * chipa = new double[nbparticles];

    // Temporary Lorentz factor
    double * gamma = new double[nbparticles];

    // Temporary double parameter
    double temp;

    // Particle id
    int ipart;

    // Radiated energy
    double rad_energy;

    // Stochastic diffusive term fo Niel et al.
    double * diffusion = new double[nbparticles];

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

    // Reinitialize the cumulative radiated energy for the current thread
    this->radiated_energy = 0.;

    // _______________________________________________________________
    // Computation

    #pragma omp simd
    for (int i=0 ; i < nbparticles; i++ ) {

        // Particle number
        ipart = istart + i;

        charge_over_mass2 = (double)(charge[ipart])*one_over_mass_2;

        // Gamma
        gamma[i] = sqrt(1.0 + momentum[0][ipart]*momentum[0][ipart]
                             + momentum[1][ipart]*momentum[1][ipart]
                             + momentum[2][ipart]*momentum[2][ipart]);

        // Computation of the Lorentz invariant quantum parameter
        chipa[i] = Radiation::compute_chipa(charge_over_mass2,
                     momentum[0][ipart],momentum[1][ipart],momentum[2][ipart],
                     gamma[i],
                     (*Epart)[ipart].x,(*Epart)[ipart].y,(*Epart)[ipart].z,
                     (*Bpart)[ipart].x,(*Bpart)[ipart].y,(*Bpart)[ipart].z);
    }

    // _______________________________________________________________
    // Computation of the diffusion coefficients

    for (int i=0 ; i < nbparticles; i++ )
    {
        // Below chipa = 1E-3, radiation losses are negligible
        if (chipa[i] > 1E-3)
        {

            // Diffusive stochastic component during dt
            diffusion[i] = RadiationTables.get_Niel_stochastic_term(gamma[i],
                                                            chipa[i],sqrtdt);
        }
        else
        {
            diffusion[i] = 0.;
        }
    }

    // _______________________________________________________________
    // Update of the momentum

    #pragma omp simd
    for (int i=0 ; i < nbparticles; i++ ) {

        // Particle number
        ipart = istart + i;

        // Radiated energy during the time step
        rad_energy =
        RadiationTables.get_corrected_cont_rad_energy_Ridgers(chipa[i],dt);

        // Effect on the momentum
        // Temporary factor
        temp = (rad_energy - diffusion[i])*gamma[i]/(gamma[i]*gamma[i]-1.);

        // Update of the momentum
        momentum[0][ipart] -= temp*momentum[0][ipart];
        momentum[1][ipart] -= temp*momentum[1][ipart];
        momentum[2][ipart] -= temp*momentum[2][ipart];

    }

    // _______________________________________________________________
    // Computation of the thread radiated energy

    double radiated_energy_loc = 0;

    #pragma omp simd reduction(+:radiated_energy_loc)
    for (int ipart=istart ; ipart<iend; ipart++ )
    {
        radiated_energy_loc += weight[ipart]*(gamma[ipart-istart] - sqrt(1.0
                            + momentum[0][ipart]*momentum[0][ipart]
                            + momentum[1][ipart]*momentum[1][ipart]
                            + momentum[2][ipart]*momentum[2][ipart]));
    }
    radiated_energy += radiated_energy_loc;

}*/
