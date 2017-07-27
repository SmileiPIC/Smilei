// ----------------------------------------------------------------------------
//! \file MultiphotonBreitWheeler.cpp
//
//! \brief This file contains the class methods for the generic class
//!  MultiphotonBreitWheeler for the photon decay into pairs via the
//!  mutliphoton Breit-Wheeler process.
//
// ----------------------------------------------------------------------------

#include "MultiphotonBreitWheeler.h"

// -----------------------------------------------------------------------------
//! Constructor for Radiation
// input: simulation parameters & Species index
//! \param params simulation parameters
//! \param species Species index
// -----------------------------------------------------------------------------
MultiphotonBreitWheeler::MultiphotonBreitWheeler(Params& params, Species * species)
{
    // Dimension position
    nDim_ = params.nDim_particle;

    // Time step
    dt    = params.timestep;

    // Normalized Schwinger Electric Field
    norm_E_Schwinger = params.electron_mass*params.c_vacuum*params.c_vacuum
                     / (params.red_planck_cst*params.referenceAngularFrequency_SI);

    // Inverse of norm_E_Schwinger
    inv_norm_E_Schwinger = 1./norm_E_Schwinger;

}

// -----------------------------------------------------------------------------
//! Destructor for MultiphotonBreitWheeler
// -----------------------------------------------------------------------------
MultiphotonBreitWheeler::~MultiphotonBreitWheeler()
{
}

// -----------------------------------------------------------------------------
//! Computation of the quantum parameter for the given
//! thread of photons
//! \param Particles class containg the particle property arrays
//! \param smpi class for mpi parameters
//! \param istart      Index of the first particle
//! \param iend        Index of the last particle
//! \param ithread     Thread index
// -----------------------------------------------------------------------------
void MultiphotonBreitWheeler::compute_thread_chiph(Particles &particles,
        SmileiMPI* smpi,
        int istart,
        int iend,
        int ithread)
{
    // _______________________________________________________________
    // Parameters
    std::vector<LocalFields> *Epart = &(smpi->dynamics_Epart[ithread]);
    std::vector<LocalFields> *Bpart = &(smpi->dynamics_Bpart[ithread]);

    // Temporary Lorentz factor
    double gamma;

    // Momentum shortcut
    double* momentum[3];
    for ( int i = 0 ; i<3 ; i++ )
        momentum[i] =  &( particles.momentum(i,0) );

    // Optical depth for the Monte-Carlo process
    double* chi = &( particles.chi(0));

    // _______________________________________________________________
    // Computation

    #pragma omp simd
    for (int ipart=istart ; ipart<iend; ipart++ )
    {

        // Gamma
        gamma = sqrt(momentum[0][ipart]*momentum[0][ipart]
                    + momentum[1][ipart]*momentum[1][ipart]
                    + momentum[2][ipart]*momentum[2][ipart]);

        // Computation of the Lorentz invariant quantum parameter
        chi[ipart] = compute_chiph(
                 momentum[0][ipart],momentum[1][ipart],momentum[2][ipart],
                 gamma,
                 (*Epart)[ipart].x,(*Epart)[ipart].y,(*Epart)[ipart].z,
                 (*Bpart)[ipart].x,(*Bpart)[ipart].y,(*Bpart)[ipart].z);

    }
}

// ---------------------------------------------------------------------------------------------------------------------
//! Overloading of the operator (): perform the pair generation
//! Monte-Carlo process for the multiphoton Breit-Wheeler
//
//! \param particles   particle object containing the particle properties
//! \param smpi        MPI properties
//! \param MultiphotonBreitWheelerTables Cross-section data tables and useful
//                     functions for multiphoton Breit-Wheeler
//! \param istart      Index of the first particle
//! \param iend        Index of the last particle
//! \param ithread     Thread index
// ---------------------------------------------------------------------------------------------------------------------
void MultiphotonBreitWheeler::operator() (
        Particles &particles,
        SmileiMPI* smpi,
        MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables,
        int istart,
        int iend,
        int ithread)
{
    // _______________________________________________________________
    // Parameters
    std::vector<LocalFields> *Epart = &(smpi->dynamics_Epart[ithread]);
    std::vector<LocalFields> *Bpart = &(smpi->dynamics_Bpart[ithread]);

    // Temporary quantum parameter
    double chiph;

    // Temporary Lorentz factor
    double gamma;

    // Temporary value
    double temp;

    // Time to event
    double event_time;

    // time spent in the iteration
    double local_it_time;

    // Number of Monte-Carlo iteration
    int mc_it_nb;

    // Momentum shortcut
    double* momentum[3];
    for ( int i = 0 ; i<3 ; i++ )
        momentum[i] =  &( particles.momentum(i,0) );

    // Position shortcut
    double* position[3];
    for ( int i = 0 ; i<nDim_ ; i++ )
        position[i] =  &( particles.position(i,0) );

    // Optical depth for the Monte-Carlo process
    double* tau = &( particles.tau(0));

    // _______________________________________________________________
    // Computation

    for (int ipart=istart ; ipart<iend; ipart++ )
    {

        // Gamma
        gamma = sqrt(momentum[0][ipart]*momentum[0][ipart]
                    + momentum[1][ipart]*momentum[1][ipart]
                    + momentum[2][ipart]*momentum[2][ipart]);

        // Computation of the Lorentz invariant quantum parameter
        chiph = MultiphotonBreitWheeler::compute_chiph(
                 momentum[0][ipart],momentum[1][ipart],momentum[2][ipart],
                 gamma,
                 (*Epart)[ipart].x,(*Epart)[ipart].y,(*Epart)[ipart].z,
                 (*Bpart)[ipart].x,(*Bpart)[ipart].y,(*Bpart)[ipart].z);

        // If the photon has enough energy
        // We also check that chiph > chiph_threshold,
        // else chiph is too low to induce a decay
        if ((gamma > 2.) && (chiph > 1E-2))
        {
            // Init local variables
            event_time = 0;
            local_it_time = 0;
            mc_it_nb = 0;

            // New even
            // If tau[ipart] <= 0, this is a new process
            if (tau[ipart] <= epsilon_tau)
            {
             // New final optical depth to reach for emision
             while (tau[ipart] <= epsilon_tau)
                tau[ipart] = -log(1.-Rand::uniform());

            }

            // Photon decay: emission under progress
            // If epsilon_tau > 0
            if (tau[ipart] > epsilon_tau)
            {
                // from the cross section
                temp = MultiphotonBreitWheelerTables.compute_dNBWdt(chiph,gamma);
            }

        }
    }
}
