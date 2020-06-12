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
MultiphotonBreitWheeler::MultiphotonBreitWheeler( Params &params, Species *species, Random * rand )
{
    // Dimension position
    n_dimensions_ = params.nDim_particle;

    // Time step
    dt_    = params.timestep;

    // Normalized Schwinger Electric Field
    norm_E_Schwinger_ = params.electron_mass*params.c_vacuum*params.c_vacuum
                        / ( params.red_planck_cst*params.reference_angular_frequency_SI );

    // Inverse of norm_E_Schwinger_
    inv_norm_E_Schwinger_ = 1./norm_E_Schwinger_;

    // Number of positrons and electrons generated per event
    mBW_pair_creation_sampling_[0] = species->mBW_pair_creation_sampling_[0];
    mBW_pair_creation_inv_sampling_[0] = 1. / mBW_pair_creation_sampling_[0];

    mBW_pair_creation_sampling_[1] = species->mBW_pair_creation_sampling_[1];
    mBW_pair_creation_inv_sampling_[1] = 1. / mBW_pair_creation_sampling_[1];

    // Threshold under which pair creation is not considered
    chiph_threshold_ = 1E-2;
    
    // Local random generator
    rand_ = rand;

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
void MultiphotonBreitWheeler::compute_thread_chiph( Particles &particles,
        SmileiMPI *smpi,
        int istart,
        int iend,
        int ithread, int ipart_ref )
{
    // _______________________________________________________________
    // Parameters
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );

    // Temporary Lorentz factor
    double gamma;

    // Momentum shortcut
    double *momentum[3];
    for( int i = 0 ; i<3 ; i++ ) {
        momentum[i] =  &( particles.momentum( i, 0 ) );
    }

    // Optical depth for the Monte-Carlo process
    double *chi = &( particles.chi( 0 ) );

    // _______________________________________________________________
    // Computation

    int nparts = Epart->size()/3;
    double *Ex = &( ( *Epart )[0*nparts] );
    double *Ey = &( ( *Epart )[1*nparts] );
    double *Ez = &( ( *Epart )[2*nparts] );
    double *Bx = &( ( *Bpart )[0*nparts] );
    double *By = &( ( *Bpart )[1*nparts] );
    double *Bz = &( ( *Bpart )[2*nparts] );
    #pragma omp simd
    for( int ipart=istart ; ipart<iend; ipart++ ) {

        // Gamma
        gamma = sqrt( momentum[0][ipart]*momentum[0][ipart]
                      + momentum[1][ipart]*momentum[1][ipart]
                      + momentum[2][ipart]*momentum[2][ipart] );

        // Computation of the Lorentz invariant quantum parameter
        chi[ipart] = compute_chiph(
                         momentum[0][ipart], momentum[1][ipart], momentum[2][ipart],
                         gamma,
                         ( *( Ex+ipart-ipart_ref ) ), ( *( Ey+ipart-ipart_ref ) ), ( *( Ez+ipart-ipart_ref ) ),
                         ( *( Bx+ipart-ipart_ref ) ), ( *( By+ipart-ipart_ref ) ), ( *( Bz+ipart-ipart_ref ) ) );

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
void MultiphotonBreitWheeler::operator()( Particles &particles,
        SmileiMPI *smpi,
        MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables,
        int istart,
        int iend,
        int ithread, int ipart_ref )
{
    // _______________________________________________________________
    // Parameters
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    // We use dynamics_invgf to store gamma
    std::vector<double> *gamma = &( smpi->dynamics_invgf[ithread] );

    int nparts = Epart->size()/3;
    double *Ex = &( ( *Epart )[0*nparts] );
    double *Ey = &( ( *Epart )[1*nparts] );
    double *Ez = &( ( *Epart )[2*nparts] );
    double *Bx = &( ( *Bpart )[0*nparts] );
    double *By = &( ( *Bpart )[1*nparts] );
    double *Bz = &( ( *Bpart )[2*nparts] );

    // Temporary value
    double temp;

    // Time to event
    double event_time;

    // Momentum shortcut
    double *momentum[3];
    for( int i = 0 ; i<3 ; i++ ) {
        momentum[i] =  &( particles.momentum( i, 0 ) );
    }

    // Position shortcut
    // Commented particles displasment while particles injection not managed  in a better way
    //    for now particles could be created outside of the local domain
    //    without been subject do boundary conditions (including domain exchange)
    //double* position[3];
    //for ( int i = 0 ; i<n_dimensions_ ; i++ )
    //    position[i] =  &( particles.position(i,0) );

    // Weight shortcut
    // double* weight = &( particles.weight(0) );

    // Optical depth for the Monte-Carlo process
    double *tau = &( particles.tau( 0 ) );

    // Quantum parameter
    double *photon_chi = &( particles.chi( 0 ) );

    // Photon id
    // uint64_t * id = &( particles.id(0));

    // Total energy converted into pairs for this species during this timestep
    this->pair_converted_energy_ = 0;

    // _______________________________________________________________
    // Computation

    // 1. Computation of gamma and chi
    //    Can be vectorized
    #pragma omp simd
    for( int ipart=istart ; ipart<iend; ipart++ ) {
        // Gamma
        ( *gamma )[ipart] = sqrt( momentum[0][ipart]*momentum[0][ipart]
                                  + momentum[1][ipart]*momentum[1][ipart]
                                  + momentum[2][ipart]*momentum[2][ipart] );

        // Computation of the Lorentz invariant quantum parameter
        photon_chi[ipart] = MultiphotonBreitWheeler::compute_chiph(
                                momentum[0][ipart], momentum[1][ipart], momentum[2][ipart],
                                ( *gamma )[ipart],
                                ( *( Ex+ipart-ipart_ref ) ), ( *( Ey+ipart-ipart_ref ) ), ( *( Ez+ipart-ipart_ref ) ),
                                ( *( Bx+ipart-ipart_ref ) ), ( *( By+ipart-ipart_ref ) ), ( *( Bz+ipart-ipart_ref ) ) );
    }

    // 2. Monte-Carlo process
    //    No vectorized
    for( int ipart=istart ; ipart<iend; ipart++ ) {

        // If the photon has enough energy
        // We also check that photon_chi > chiph_threshold,
        // else photon_chi is too low to induce a decay
        if( ( ( *gamma )[ipart] > 2. ) && ( photon_chi[ipart] > chiph_threshold_ ) ) {
            // Init local variables
            event_time = 0;

            // New even
            // If tau[ipart] <= 0, this is a new process
            if( tau[ipart] <= epsilon_tau_ ) {
                // New final optical depth to reach for emision
                while( tau[ipart] <= epsilon_tau_ ) {
                    //tau[ipart] = -log( 1.-Rand::uniform() );
                    tau[ipart] = -log( 1.-rand_->uniform() );
                }

            }

            // Photon decay: emission under progress
            // If epsilon_tau_ > 0
            else if( tau[ipart] > epsilon_tau_ ) {
                // from the cross section
                temp = MultiphotonBreitWheelerTables.computeBreitWheelerPairProductionRate( photon_chi[ipart], ( *gamma )[ipart] );
                
                // Time to decay
                // If this time is above the remaining iteration time,
                // There is a synchronization at the end of the pic iteration
                // and the process continues at the next
                event_time = std::min( tau[ipart]/temp, dt_ );

                // Update of the optical depth
                tau[ipart] -= temp*event_time;

                // If the final optical depth is reached
                // The photon decays into pairs
                if( tau[ipart] <= epsilon_tau_ ) {

                    // Update of the position
                    // Move the photons

//#ifdef  __DEBUG
//                    for ( int i = 0 ; i<n_dimensions_ ; i++ )
//                        particles.position_old(i,ipart) = position[i][ipart];
//#endif
//                    for ( int i = 0 ; i<n_dimensions_ ; i++ )
//                        position[i][ipart]     += event_time*momentum[i][ipart]/(*gamma)[ipart];

                    // Generation of the pairs
                    MultiphotonBreitWheeler::pair_emission( ipart,
                                                            particles,
                                                            ( *gamma )[ipart],
                                                            dt_ - event_time,
                                                            MultiphotonBreitWheelerTables );


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
void MultiphotonBreitWheeler::pair_emission( int ipart,
        Particles &particles,
        double &gammaph,
        double remaining_dt,
        MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables )
{

    // _______________________________________________
    // Parameters

    int      nparticles;           // Total number of particles in the temporary arrays
    int      k, i;
    double   u[3];                 // propagation direction
    double *chi = new double[2];   // temporary quantum parameters
    double   inv_chiph_gammaph;    // (gamma_ph - 2) / chi
    double   p;
    // Commented particles displasment while particles injection not managed  in a better way
    //    for now particles could be created outside of the local domain
    //    without been subject do boundary conditions (including domain exchange)
    //double   inv_gamma;

    inv_chiph_gammaph = ( gammaph-2. )/particles.chi( ipart );

    // Get the pair quantum parameters to compute the energy
    chi = MultiphotonBreitWheelerTables.computePairQuantumParameter( particles.chi( ipart ), rand_ );
    
    // pair propagation direction // direction of the photon
    for( k = 0 ; k<3 ; k++ ) {
        u[k] = particles.momentum( k, ipart )/gammaph;
    }

    // _______________________________________________
    // Electron (k=0) and positron (k=1) generation

    for( k=0 ; k < 2 ; k++ ) {

        // Creation of new electrons in the temporary array new_pair[0]
        new_pair[k].createParticles( mBW_pair_creation_sampling_[k] );
        
        // Final size
        nparticles = new_pair[k].size();

        // For all new electrons...
        for( int idNew=nparticles-mBW_pair_creation_sampling_[k]; idNew<nparticles; idNew++ ) {

            // Momentum
            p = sqrt( pow( 1.+chi[k]*inv_chiph_gammaph, 2 )-1 );
            for( i=0; i<3; i++ ) {
                new_pair[k].momentum( i, idNew ) =
                    p*u[i];
            }

            // gamma
            //inv_gamma = 1./sqrt(1.+p*p);

            // Positions
            for( i=0; i<n_dimensions_; i++ ) {
                new_pair[k].position( i, idNew )=particles.position( i, ipart );
//               + new_pair[k].momentum(i,idNew)*remaining_dt*inv_gamma;
            }

            // Old positions
#ifdef  __DEBUG
            for( i=0; i<n_dimensions_; i++ ) {
                new_pair[k].position_old( i, idNew )=particles.position( i, ipart ) ;
            }
#endif

            new_pair[k].weight( idNew )=particles.weight( ipart )*mBW_pair_creation_inv_sampling_[k];
            new_pair[k].charge( idNew )= k*2-1;

            if( new_pair[k].isQuantumParameter ) {
                new_pair[k].chi( idNew ) = chi[k];
            }

            if( new_pair[k].isMonteCarlo ) {
                new_pair[k].tau( idNew ) = -1.;
            }
        }
    }

    // Total energy converted into pairs during the current timestep
    pair_converted_energy_ += particles.weight( ipart )*gammaph;

    // The photon with negtive weight will be deleted latter
    particles.weight( ipart ) = -1;

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
    SmileiMPI *smpi,
    int ibin, int nbin,
    int *bmin, int *bmax, int ithread )
{
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    std::vector<double> *gamma = &( smpi->dynamics_invgf[ithread] );
    std::vector<int> *iold = &( smpi->dynamics_iold[ithread] );
    std::vector<double> *deltaold = &( smpi->dynamics_deltaold[ithread] );

    std::vector<double> *thetaold = NULL;
    if ( smpi->dynamics_thetaold.size() )
        thetaold = &( smpi->dynamics_thetaold[ithread] );

    int nparts = Epart->size()/3;


    if( bmax[ibin] > bmin[ibin] ) {
        // Weight shortcut
        double *weight = &( particles.weight( 0 ) );

        // Index of the last existing photon (weight > 0)
        int last_photon_index;
        int first_photon_index;
        int ii;
        int nb_deleted_photon;

        // Backward loop over the photons to fing the first existing photon
        last_photon_index = bmax[ibin]-1;
        first_photon_index = bmin[ibin];
        while( ( last_photon_index >= bmin[ibin] )
                && ( weight[last_photon_index] <= 0 ) ) {
            last_photon_index--;
        }
        while( ( first_photon_index < bmax[ibin] )
                && ( weight[first_photon_index] > 0 ) ) {
            first_photon_index++;
        }
        // At this level, last_photon_index is the position of the last photon
        // that will not be erased

        // Backward loop over the photons to fill holes in the photon particle array
        for( int ipart=last_photon_index-1 ; ipart>=bmin[ibin]; ipart-- ) {
            if( weight[ipart] <= 0 ) {
                if( ipart < last_photon_index ) {
                    // The last existing photon comes to the position of
                    // the deleted photon
                    particles.overwriteParticle( last_photon_index, ipart );
                    // Overwrite bufferised data
                    for ( int iDim=2 ; iDim>=0 ; iDim-- ) {
                        (*Epart)[iDim*nparts+ipart] = (*Epart)[iDim*nparts+last_photon_index];
                        (*Bpart)[iDim*nparts+ipart] = (*Bpart)[iDim*nparts+last_photon_index];
                    }
                    for ( int iDim=n_dimensions_-1 ; iDim>=0 ; iDim-- ) {
                        (*iold)[iDim*nparts+ipart] = (*iold)[iDim*nparts+last_photon_index];
                        (*deltaold)[iDim*nparts+ipart] = (*deltaold)[iDim*nparts+last_photon_index];
                    }
                    (*gamma)[0*nparts+ipart] = (*gamma)[0*nparts+last_photon_index];

                    if (thetaold) {
                        (*thetaold)[0*nparts+ipart] = (*thetaold)[0*nparts+last_photon_index];
                    }

                    last_photon_index --;


                }
            }
        }

        // Removal of the photons
        nb_deleted_photon = bmax[ibin]-last_photon_index-1;

        if( nb_deleted_photon > 0 ) {
            particles.eraseParticle( last_photon_index+1, nb_deleted_photon );
            // Erase bufferised data
            for ( int iDim=2 ; iDim>=0 ; iDim-- ) {
                Epart->erase(Epart->begin()+iDim*nparts+last_photon_index+1,Epart->begin()+iDim*nparts+last_photon_index+1+nb_deleted_photon);
                Bpart->erase(Bpart->begin()+iDim*nparts+last_photon_index+1,Bpart->begin()+iDim*nparts+last_photon_index+1+nb_deleted_photon);
            }
            for ( int iDim=n_dimensions_-1 ; iDim>=0 ; iDim-- ) {
                iold->erase(iold->begin()+iDim*nparts+last_photon_index+1,iold->begin()+iDim*nparts+last_photon_index+1+nb_deleted_photon);
                deltaold->erase(deltaold->begin()+iDim*nparts+last_photon_index+1,deltaold->begin()+iDim*nparts+last_photon_index+1+nb_deleted_photon);
            }
            gamma->erase(gamma->begin()+0*nparts+last_photon_index+1,gamma->begin()+0*nparts+last_photon_index+1+nb_deleted_photon);

            if (thetaold) {
                thetaold->erase(thetaold->begin()+0*nparts+last_photon_index+1,thetaold->begin()+0*nparts+last_photon_index+1+nb_deleted_photon);
            }



            bmax[ibin] = last_photon_index+1;
            for( ii=ibin+1; ii<nbin; ii++ ) {
                bmin[ii] -= nb_deleted_photon;
                bmax[ii] -= nb_deleted_photon;
            }
        }
    }
}
