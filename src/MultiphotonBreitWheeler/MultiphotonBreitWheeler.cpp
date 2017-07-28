// ----------------------------------------------------------------------------
//! \file MultiphotonBreitWheeler.cpp
//
//! \brief This file contains the class methods for the generic class
//!  MultiphotonBreitWheeler for the photon decay into pairs via the
//!  mutliphoton Breit-Wheeler process.
//
// ----------------------------------------------------------------------------

#include "MultiphotonBreitWheeler.h"
#include "Species.h"

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

    // Number of positrons and electrons generated per event
    this->mBW_pair_creation_sampling[0] = species->mBW_pair_creation_sampling[0];
    this->mBW_pair_creation_inv_sampling[0] = 1. / mBW_pair_creation_sampling[0];

    mBW_pair_creation_sampling[1] = species->mBW_pair_creation_sampling[1];
    mBW_pair_creation_inv_sampling[1] = 1. / mBW_pair_creation_sampling[1];

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

    // Momentum shortcut
    double* momentum[3];
    for ( int i = 0 ; i<3 ; i++ )
        momentum[i] =  &( particles.momentum(i,0) );

    // Position shortcut
    double* position[3];
    for ( int i = 0 ; i<nDim_ ; i++ )
        position[i] =  &( particles.position(i,0) );

    // Weight shortcut
    double* weight = &( particles.weight(0) );

    // Optical depth for the Monte-Carlo process
    double* tau = &( particles.tau(0));

    // Photon id
    // uint64_t * id = &( particles.id(0));

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
            else if (tau[ipart] > epsilon_tau)
            {
                // from the cross section
                temp = MultiphotonBreitWheelerTables.compute_dNBWdt(chiph,gamma);

                // Time to decay
                // If this time is above the remaining iteration time,
                // There is a synchronization at the end of the pic iteration
                // and the process continues at the next
                event_time = std::min(tau[ipart]/temp, dt);

                // Update of the optical depth
                tau[ipart] -= temp*event_time;

                // If the final optical depth is reached
                // The photon decays into pairs
                if (tau[ipart] <= epsilon_tau)
                {

                    // Update of the position
                    // Move the photons
#ifdef  __DEBUG
                    for ( int i = 0 ; i<nDim_ ; i++ )
                        position_old[i][ipart] = position[i][ipart];
#endif
                    for ( int i = 0 ; i<nDim_ ; i++ )
                        position[i][ipart]     += event_time*momentum[i][ipart]/gamma;

                    // Generation of the pairs
                    MultiphotonBreitWheeler::pair_emission(ipart,
                                           chiph,gamma,
                                           position,
                                           momentum,
                                           weight,
                                           MultiphotonBreitWheelerTables);

                    // Optical depth becomes negative meaning
                    // that a new drawing is possible
                    // at the next Monte-Carlo iteration
                    tau[ipart] = -1.;
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
//! Perform the creation of pairs from a photon
//! \param ipart              photon index
//! \param chipa              photon quantum parameter
//! \param gammapa            photon normalized energy
//! \param position           photon position
//! \param momentum           photon momentum
//! \param MultiphotonBreitWheelerTables    Cross-section data tables
//!                       and useful functions
//!                       for the multiphoton Breit-Wheeler process
// -----------------------------------------------------------------------------
void MultiphotonBreitWheeler::pair_emission(int ipart,
                double & chiph,
                double & gammaph,
                double * position[3],
                double * momentum[3],
                double * weight,
                MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables)
{

    // _______________________________________________
    // Parameters
    double u[3]; // propagation direction
    int nparticles; // Total number of particles in the temporary arrays

    // pair propagation direction // direction of the photon
    for ( int i = 0 ; i<nDim_ ; i++ ) {
        u[i] = momentum[i][ipart]/gammaph;
    }

    // _______________________________________________
    // Creation of an electron

    // Creation of new electrons in the temporary array new_pair[0]
    new_pair[0].create_particles(mBW_pair_creation_sampling[0]);

    // Final size
    nparticles = new_pair[0].size();

    // For all new electrons...
    for (int idNew=nparticles-mBW_pair_creation_sampling[0]; idNew<nparticles; idNew++)
    {
        // Positions
        for (int i=0; i<nDim_; i++) {
            new_pair[0].position(i,idNew)=position[i][ipart];
        }

        // Momentum
        for (int i=0; i<3; i++) {
            new_pair[0].momentum(i,idNew) =
            sqrt(pow(0.5*gammaph,2)-1)*u[i];
        }

        new_pair[0].weight(idNew)=weight[ipart]*mBW_pair_creation_inv_sampling[0];
        new_pair[0].charge(idNew)=-1.;

        if (new_pair[0].isQuantumParameter)
        {
            new_pair[0].chi(idNew) = 0.5*chiph;
        }

        if (new_pair[0].isMonteCarlo)
        {
            new_pair[0].tau(idNew) = -1.;
        }

    }

    // _______________________________________________
    // Creation of a positron

    // Creation of new positrons in the temporary array new_pair[1]
    new_pair[1].create_particles(mBW_pair_creation_sampling[1]);

    // Final size
    nparticles = new_pair[1].size();

    // For all new positrons...
    for (int idNew=nparticles-mBW_pair_creation_sampling[1]; idNew<nparticles; idNew++)
    {
        // Positions
        for (int i=0; i<nDim_; i++) {
            new_pair[1].position(i,idNew)=position[i][ipart];
        }

        // Momentum
        for (int i=0; i<3; i++) {
            new_pair[1].momentum(i,idNew) =
            sqrt(pow(0.5*gammaph,2)-1)*u[i];
        }

        new_pair[1].weight(idNew)=weight[ipart]*mBW_pair_creation_inv_sampling[1];
        new_pair[1].charge(idNew)=1.;

        if (new_pair[1].isQuantumParameter)
        {
            new_pair[1].chi(idNew) = 0.5*chiph;
        }

        if (new_pair[1].isMonteCarlo)
        {
            new_pair[1].tau(idNew) = -1.;
        }
    }

    // The photon with negtive weight will be deleted latter
    //weight[ipart] = -1;

}
