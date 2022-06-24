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
    norm_E_Schwinger_ = params.electron_mass*params.c_vacuum_*params.c_vacuum_
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
void MultiphotonBreitWheeler::computeThreadPhotonChi( Particles &particles,
        SmileiMPI *smpi,
        int istart,
        int iend,
        int ithread, int ipart_ref )
{
    // _______________________________________________________________
    // Parameters
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );

    int nparts = Epart->size()/3;
    const double *const __restrict__ Ex = &( ( *Epart )[0*nparts] );
    const double *const __restrict__ Ey = &( ( *Epart )[1*nparts] );
    const double *const __restrict__ Ez = &( ( *Epart )[2*nparts] );
    const double *const __restrict__ Bx = &( ( *Bpart )[0*nparts] );
    const double *const __restrict__ By = &( ( *Bpart )[1*nparts] );
    const double *const __restrict__ Bz = &( ( *Bpart )[2*nparts] );
    
    // Particles Momentum shortcut
    const double *const __restrict__ momentum_x = particles.getPtrMomentum(0);
    const double *const __restrict__ momentum_y = particles.getPtrMomentum(1);
    const double *const __restrict__ momentum_z = particles.getPtrMomentum(2);

    // Quantum parameter
    double *const __restrict__ chi = particles.getPtrChi();

    // _______________________________________________________________
    // Computation

    
    #pragma omp simd
    for( int ipart=istart ; ipart<iend; ipart++ ) {

        // Gamma (Lorentz factor)
        const double gamma = std::sqrt( momentum_x[ipart]*momentum_x[ipart]
                     + momentum_y[ipart]*momentum_y[ipart]
                    + momentum_z[ipart]*momentum_z[ipart] );

        // Computation of the Lorentz invariant quantum parameter
        chi[ipart] = computePhotonChi(
                         momentum_x[ipart], momentum_y[ipart], momentum_z[ipart],
                         gamma,
                         Ex[ipart-ipart_ref], Ey[ipart-ipart_ref], Ez[ipart-ipart_ref],
                         Bx[ipart-ipart_ref], By[ipart-ipart_ref], Bz[ipart-ipart_ref] );

    }
}

// ---------------------------------------------------------------------------------------------------------------------
//! Overloading of the operator (): perform the pair generation
//! Monte-Carlo process for the multiphoton Breit-Wheeler
//
//! \param particles   particle object containing the particle properties
//! \param smpi        MPI properties
//! \param MultiphotonBreitWheelerTables Cross-section data tables and useful
//!                     functions for multiphoton Breit-Wheeler
//! \param pair_energy energy converted into pairs
//! \param istart      Index of the first particle
//! \param iend        Index of the last particle
//! \param ithread     Thread index
// ---------------------------------------------------------------------------------------------------------------------
void MultiphotonBreitWheeler::operator()( Particles &particles,
        SmileiMPI* smpi,
        Particles** new_pair,
        MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables,
        double & pair_energy,
        int istart,
        int iend,
        int ithread, int ipart_ref )
{
    // _______________________________________________________________
    // Parameters
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    
    // We use dynamics_invgf to store gamma
    double * const __restrict__ photon_gamma = &( smpi->dynamics_invgf[ithread][0] );

    int nparts = Epart->size()/3;
    
    const double *const __restrict__ Ex = &( ( *Epart )[0*nparts] );
    const double *const __restrict__ Ey = &( ( *Epart )[1*nparts] );
    const double *const __restrict__ Ez = &( ( *Epart )[2*nparts] );
    const double *const __restrict__ Bx = &( ( *Bpart )[0*nparts] );
    const double *const __restrict__ By = &( ( *Bpart )[1*nparts] );
    const double *const __restrict__ Bz = &( ( *Bpart )[2*nparts] );

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
    double *const __restrict__ position_x = particles.getPtrPosition( 0 );
    double *const __restrict__ position_y = n_dimensions_ > 1 ? particles.getPtrPosition( 1 ) : nullptr;
    double *const __restrict__ position_z = n_dimensions_ > 2 ? particles.getPtrPosition( 2 ) : nullptr;

    // Particles Momentum shortcut
    double *const __restrict__ momentum_x = particles.getPtrMomentum(0);
    double *const __restrict__ momentum_y = particles.getPtrMomentum(1);
    double *const __restrict__ momentum_z = particles.getPtrMomentum(2);

    // Weight shortcut
    double* weight = particles.getPtrWeight();

    // Optical depth for the Monte-Carlo process
    double *tau = &( particles.tau( 0 ) );

    // Quantum parameter
    double *const __restrict__ photon_chi = particles.getPtrChi();

    // Photon id
    // uint64_t * id = &( particles.id(0));

    // Reserve pair particles (else, pointer could become obsolete)
    double np = new_pair[0]->size();
    new_pair[0]->reserve( np + mBW_pair_creation_sampling_[0] * (iend - istart) );
    new_pair[1]->reserve( np + mBW_pair_creation_sampling_[1] * (iend - istart) );

    // _______________________________________________________________
    // Computation

    // 1. Computation of gamma and chi
    //    Can be vectorized
    #pragma omp simd
    for( int ipart=istart ; ipart<iend; ipart++ ) {
        // Gamma
        photon_gamma[ipart] = std::sqrt( momentum_x[ipart]*momentum_x[ipart]
                                  + momentum_y[ipart]*momentum_y[ipart]
                                  + momentum_z[ipart]*momentum_z[ipart] );

        // Computation of the Lorentz invariant quantum parameter
        photon_chi[ipart] = MultiphotonBreitWheeler::computePhotonChi(
                                momentum_x[ipart], momentum_y[ipart], momentum_z[ipart],
                                photon_gamma[ipart],
                                Ex[ipart-ipart_ref], Ey[ipart-ipart_ref], Ez[ipart-ipart_ref],
                                Bx[ipart-ipart_ref], By[ipart-ipart_ref], Bz[ipart-ipart_ref] );
    }

    // 2. Monte-Carlo process
    //    No vectorized
    for( int ipart=istart ; ipart<iend; ipart++ ) {

        // If the photon has enough energy
        // We also check that photon_chi > chiph_threshold,
        // else photon_chi is too low to induce a decay
        if( ( photon_gamma[ipart] > 2. ) && ( photon_chi[ipart] > chiph_threshold_ ) ) {
            // Init local variables
            event_time = 0;

            // New even
            // If tau[ipart] <= 0, this is a new process
            if( tau[ipart] <= epsilon_tau_ ) {
                // New final optical depth to reach for emision
                while( tau[ipart] <= epsilon_tau_ ) {
                    //tau[ipart] = -log( 1.-Rand::uniform() );
                    tau[ipart] = -std::log( 1.-rand_->uniform() );
                }

            }

            // Photon decay: emission under progress
            // If epsilon_tau_ > 0
            else if( tau[ipart] > epsilon_tau_ ) {
                // from the cross section
                temp = MultiphotonBreitWheelerTables.computeBreitWheelerPairProductionRate( photon_chi[ipart], photon_gamma [ipart] );

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
                    // pair_energy += MultiphotonBreitWheeler::pair_emission( ipart,
                    //                                         particles,
                    //                                         ( *gamma )[ipart],
                    //                                         dt_ - event_time,
                    //                                         MultiphotonBreitWheelerTables );

                    double inv_chiph_gammaph = ( photon_gamma[ipart]-2. ) / photon_chi[ipart];

                    double pair_chi[2];

                    // Get the pair quantum parameters to compute the energy
                    MultiphotonBreitWheelerTables.computePairQuantumParameter( photon_chi[ipart], pair_chi, rand_ );

                    // pair propagation direction // direction of the photon
                    double ux = momentum_x[ipart]/photon_gamma[ipart];
                    double uy = momentum_y[ipart]/photon_gamma[ipart];
                    double uz = momentum_z[ipart]/photon_gamma[ipart];

                    // Loop on the pair (2 particles)
                    for( int k=0 ; k < 2 ; k++ ) {

                        // Creation of new electrons in the temporary array new_pair[0]
                        new_pair[k]->createParticles( mBW_pair_creation_sampling_[k] );

                        // Final size
                        const int nparticles = new_pair[k]->size();

                        // For all new electrons...
                        for( int idNew=nparticles-mBW_pair_creation_sampling_[k]; idNew<nparticles; idNew++ ) {

                            // Momentum
                            const double p = std::sqrt( std::pow( 1.+pair_chi[k]*inv_chiph_gammaph, 2 )-1 );
                            new_pair[k]->momentum( 0, idNew ) = p*ux;
                            new_pair[k]->momentum( 1, idNew ) = p*uy;
                            new_pair[k]->momentum( 2, idNew ) = p*uz;

                            // gamma
                            //inv_gamma = 1./sqrt(1.+p*p);

                            // Positions
                            for( int i=0; i<n_dimensions_; i++ ) {
                                new_pair[k]->position( i, idNew )=particles.position( i, ipart );
                //               + new_pair[k].momentum(i,idNew)*remaining_dt*inv_gamma;
                            }

                            // Old positions
                            if( particles.Position_old.size() > 0 ) {
                                for( int i=0; i<n_dimensions_; i++ ) {
                                    new_pair[k]->position_old( i, idNew )=particles.position( i, ipart ) ;
                                }
                            }

                            new_pair[k]->weight( idNew )=weight[ipart]*mBW_pair_creation_inv_sampling_[k];
                            new_pair[k]->charge( idNew )= k*2-1;

                            if( new_pair[k]->isQuantumParameter ) {
                                new_pair[k]->chi( idNew ) = pair_chi[k];
                            }

                            if( new_pair[k]->isMonteCarlo ) {
                                new_pair[k]->tau( idNew ) = -1.;
                            }
                        }
                    }

                    // Total energy converted into pairs during the current timestep
                    pair_energy += weight[ipart]*photon_gamma[ipart];

                    // The photon with negtive weight will be deleted latter
                    weight[ipart] = -1;

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

    std::vector<std::complex<double>> *thetaold = NULL;
    if ( smpi->dynamics_eithetaold.size() )
        thetaold = &( smpi->dynamics_eithetaold[ithread] );

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
