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

    //! Threshold under which pair creation is not considered
    chiph_threashold = 1E-2;

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
void MultiphotonBreitWheeler::operator() (Particles &particles,
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
    // We use dynamics_invgf to store gamma
    std::vector<double> * gamma = &(smpi->dynamics_invgf[ithread]);

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
    // double* weight = &( particles.weight(0) );

    // Optical depth for the Monte-Carlo process
    double* tau = &( particles.tau(0));

    // Quantum parameter
    double * chiph = &( particles.chi(0));

    // Photon id
    // uint64_t * id = &( particles.id(0));

    // Total energy converted into pairs for this species during this timestep
    this->pair_converted_energy = 0;

    // _______________________________________________________________
    // Computation

    // 1. Computation of gamma and chi
    //    Can be vectorized
    #pragma omp simd
    for (int ipart=istart ; ipart<iend; ipart++ )
    {
        // Gamma
        (*gamma)[ipart] = sqrt(momentum[0][ipart]*momentum[0][ipart]
                    + momentum[1][ipart]*momentum[1][ipart]
                    + momentum[2][ipart]*momentum[2][ipart]);

        // Computation of the Lorentz invariant quantum parameter
        chiph[ipart] = MultiphotonBreitWheeler::compute_chiph(
                 momentum[0][ipart],momentum[1][ipart],momentum[2][ipart],
                 (*gamma)[ipart],
                 (*Epart)[ipart].x,(*Epart)[ipart].y,(*Epart)[ipart].z,
                 (*Bpart)[ipart].x,(*Bpart)[ipart].y,(*Bpart)[ipart].z);
    }

    // 2. Monte-Carlo process
    //    No vectorized
    for (int ipart=istart ; ipart<iend; ipart++ )
    {

        // If the photon has enough energy
        // We also check that chiph > chiph_threshold,
        // else chiph is too low to induce a decay
        if (((*gamma)[ipart] > 2.) && (chiph[ipart] > chiph_threashold))
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
                temp = MultiphotonBreitWheelerTables.compute_dNBWdt(chiph[ipart],(*gamma)[ipart]);

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
                        particles.position_old(i,ipart) = position[i][ipart];
#endif
                    for ( int i = 0 ; i<nDim_ ; i++ )
                        position[i][ipart]     += event_time*momentum[i][ipart]/(*gamma)[ipart];


                    // Generation of the pairs
                    /*MultiphotonBreitWheeler::pair_emission(ipart,
                                           chiph[ipart],(*gamma)[ipart],
                                           position,
                                           momentum,
                                           weight,
                                           MultiphotonBreitWheelerTables);*/

                   MultiphotonBreitWheeler::pair_emission_2(ipart,
                                          particles,
                                          (*gamma)[ipart],
                                          dt - event_time,
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
//! \param remaining_dt       remaining time before the end of the iteration
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
                double remaining_dt,
                MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables)
{

    // _______________________________________________
    // Parameters

    double   u[3]; // propagation direction
    int      nparticles; // Total number of particles in the temporary arrays
    double * chi = new double[2];
    double   inv_chiph_gammaph;
    double   p;
    double   inv_gamma;

    inv_chiph_gammaph = (gammaph-2.)/chiph;

    // Get the pair quantum parameters to compute the energy
    chi = MultiphotonBreitWheelerTables.compute_pair_chi(chiph);

    // pair propagation direction // direction of the photon
    for ( int i = 0 ; i<3 ; i++ ) {
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

        // Momentum
        p = sqrt(pow(1.+chi[0]*inv_chiph_gammaph,2)-1);
        for (int i=0; i<3; i++) {
            new_pair[0].momentum(i,idNew) =
            p*u[i];
        }

        // gamma
        inv_gamma = 1./sqrt(1.+p*p);

        // Positions
        for (int i=0; i<nDim_; i++) {
            new_pair[0].position(i,idNew)=position[i][ipart]
               + new_pair[0].momentum(i,idNew)*remaining_dt*inv_gamma;
        }

        // Old positions
#ifdef  __DEBUG
        for (int i=0; i<nDim_; i++) {
            new_pair[0].position_old(i,idNew)=position[i][ipart];
        }
#endif

        new_pair[0].weight(idNew)=weight[ipart]*mBW_pair_creation_inv_sampling[0];
        new_pair[0].charge(idNew)=-1.;

        if (new_pair[0].isQuantumParameter)
        {
            new_pair[0].chi(idNew) = chi[0];
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

        // Momentum
        p = sqrt(pow(1.+chi[1]*inv_chiph_gammaph,2)-1);
        for (int i=0; i<3; i++) {
            new_pair[1].momentum(i,idNew) =
            p*u[i];
        }

        // gamma
        inv_gamma = 1./sqrt(1.+p*p);

        // Positions
        for (int i=0; i<nDim_; i++) {
            new_pair[1].position(i,idNew)=position[i][ipart]
               + new_pair[1].momentum(i,idNew)*remaining_dt*inv_gamma;
        }

        // Old positions
#ifdef  __DEBUG
        for (int i=0; i<nDim_; i++) {
            new_pair[1].position_old(i,idNew)=position[i][ipart];
        }
#endif

        new_pair[1].weight(idNew)=weight[ipart]*mBW_pair_creation_inv_sampling[1];
        new_pair[1].charge(idNew)=1.;

        if (new_pair[1].isQuantumParameter)
        {
            new_pair[1].chi(idNew) = chi[1];
        }

        if (new_pair[1].isMonteCarlo)
        {
            new_pair[1].tau(idNew) = -1.;
        }
    }

    // Total energy converted into pairs during the current timestep
    this->pair_converted_energy += weight[ipart]*gammaph;

    // The photon with negtive weight will be deleted latter
    weight[ipart] = -1;

}


// -----------------------------------------------------------------------------
//! Second version of pair_emission:
//! Perform the creation of pairs from a photon with particles as an argument
//! \param ipart              photon index
//! \param particles          object particles containing the photons and their properties
//! \param gammaph            photon normalized energy
//! \param remaining_dt       remaining time before the end of the iteration
//! \param MultiphotonBreitWheelerTables    Cross-section data tables
//!                       and useful functions
//!                       for the multiphoton Breit-Wheeler process
// -----------------------------------------------------------------------------
void MultiphotonBreitWheeler::pair_emission_2(int ipart,
                Particles & particles,
                double & gammaph,
                double remaining_dt,
                MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables)
{

    // _______________________________________________
    // Parameters

    int      nparticles;           // Total number of particles in the temporary arrays
    int      k,i;
    double   u[3];                 // propagation direction
    double * chi = new double[2];  // temporary quantum parameters
    double   inv_chiph_gammaph;    // (gamma_ph - 2) / chi
    double   p;
    double   inv_gamma;

    inv_chiph_gammaph = (gammaph-2.)/particles.chi(ipart);

    // Get the pair quantum parameters to compute the energy
    chi = MultiphotonBreitWheelerTables.compute_pair_chi( particles.chi(ipart));

    // pair propagation direction // direction of the photon
    for (k = 0 ; k<3 ; k++ ) {
        u[k] = particles.momentum(k,ipart)/gammaph;
    }

    // _______________________________________________
    // Electron (k=0) and positron (k=1) generation

    for (k=0 ; k < 2 ; k++)
    {

        // Creation of new electrons in the temporary array new_pair[0]
        new_pair[k].create_particles(mBW_pair_creation_sampling[k]);

        // Final size
        nparticles = new_pair[k].size();

        // For all new electrons...
        for (int idNew=nparticles-mBW_pair_creation_sampling[k]; idNew<nparticles; idNew++)
        {

            // Momentum
            p = sqrt(pow(1.+chi[k]*inv_chiph_gammaph,2)-1);
            for (i=0; i<3; i++) {
                new_pair[k].momentum(i,idNew) =
                p*u[i];
            }

            // gamma
            inv_gamma = 1./sqrt(1.+p*p);

            // Positions
            for (i=0; i<nDim_; i++) {
                new_pair[k].position(i,idNew)=particles.position(i,ipart)
               + new_pair[k].momentum(i,idNew)*remaining_dt*inv_gamma;
            }

            // Old positions
#ifdef  __DEBUG
            for (i=0; i<nDim_; i++) {
                new_pair[k].position_old(i,idNew)=particles.position(i,ipart) ;
            }
#endif

            new_pair[k].weight(idNew)=particles.weight(ipart)*mBW_pair_creation_inv_sampling[k];
            new_pair[k].charge(idNew)= k*2-1;

            if (new_pair[k].isQuantumParameter)
            {
                new_pair[k].chi(idNew) = chi[k];
            }

            if (new_pair[k].isMonteCarlo)
            {
                new_pair[k].tau(idNew) = -1.;
            }
        }
    }

    // Total energy converted into pairs during the current timestep
    this->pair_converted_energy += particles.weight(ipart)*gammaph;

    // The photon with negtive weight will be deleted latter
    particles.weight(ipart) = -1;

}

// -----------------------------------------------------------------------------
//! Clean photons that decayed into pairs (weight <= 0)
//! \param particles   particle object containing the particle
//!                    properties of the current species
//! \param istart      Index of the first particle
//! \param iend        Index of the last particle
//! \param ithread     Thread index
// -----------------------------------------------------------------------------
void MultiphotonBreitWheeler::decayed_photon_cleaning(
                Particles &particles,
                int ibin, int nbin,
                int * bmin,int * bmax)
{
    if (bmax[ibin] > bmin[ibin])
    {
        // Weight shortcut
        double* weight = &( particles.weight(0) );

        // Index of the last existing photon (weight > 0)
        int last_photon_index;
        int first_photon_index;
        int ii;
        int nb_deleted_photon;

        // Backward loop over the photons to fing the first existing photon
        last_photon_index = bmax[ibin]-1;
        first_photon_index = bmin[ibin];
        while ((last_photon_index >= bmin[ibin])
           && (weight[last_photon_index] <= 0))
        {
                last_photon_index--;
        }
        while ((first_photon_index < bmax[ibin])
           && (weight[first_photon_index] > 0))
        {
                first_photon_index++;
        }
        // At this level, last_photon_index is the position of the last photon
        // that will not be erased

        // Backward loop over the photons to fill holes in the photon particle array
        for (int ipart=last_photon_index-1 ; ipart>=bmin[ibin]; ipart-- )
        {
            if (weight[ipart] <= 0)
            {
                if (ipart < last_photon_index)
                {
                    // The last existing photon comes to the position of
                    // the deleted photon
                    particles.overwrite_part(last_photon_index,ipart);
                    last_photon_index --;
                }
            }
        }

        // Removal of the photons
        nb_deleted_photon = bmax[ibin]-last_photon_index-1;
        /*if (last_photon_index+1 <= bmax[ibin]-1)
        {
            std::cerr << "Photon cleaning: " << last_photon_index+1
                      << " bmin[ibin]:" << bmin[ibin]
                      << " bmax[ibin]-1: " << bmax[ibin]-1
                      << " nb deleted ph: " << nb_deleted_photon
                      << std::endl;
        }*/
        /*std::cerr << "Photon cleaning: " << last_photon_index+1
                  << " " << bmax[ibin]
                  << " " << nb_deleted_photon
                  << std::endl;*/
        if (nb_deleted_photon > 0)
        {
            particles.erase_particle(last_photon_index+1,nb_deleted_photon);
            bmax[ibin] = last_photon_index+1;
            for (ii=ibin+1; ii<nbin; ii++) {
                bmin[ii] -= nb_deleted_photon;
                bmax[ii] -= nb_deleted_photon;
            }
        }
    }
}
