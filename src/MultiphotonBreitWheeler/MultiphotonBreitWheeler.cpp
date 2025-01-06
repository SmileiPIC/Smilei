// ----------------------------------------------------------------------------
//! \file MultiphotonBreitWheeler.cpp
//
//! \brief This file contains the class methods for the generic class
//!  MultiphotonBreitWheeler for the photon decay into pairs via the
//!  multiphoton Breit-Wheeler process.
//
// ----------------------------------------------------------------------------

#include "MultiphotonBreitWheeler.h"
#include "Species.h"

#if defined(SMILEI_ACCELERATOR_GPU_OACC)
    #define __HIP_PLATFORM_NVCC__
    #define __HIP_PLATFORM_NVIDIA__
    #include "gpuRandom.h"
#elif defined(SMILEI_ACCELERATOR_GPU_OMP)
    #define __HIP_PLATFORM_HCC__
    #define __HIP_PLATFORM_AMD__
    #include "gpuRandom.h"
#endif

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

    // Get the respectiv charge
    // mBW_pair_species_charge_[0] = species->mBW_pair_species_[0]->max_charge_;
    // mBW_pair_species_charge_[1] = species->mBW_pair_species_[1]->max_charge_;

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

    int nparts = smpi->getBufferSize(ithread);
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
//! \param particles        particle object containing the particle properties
//! \param smpi             MPI properties
//! \param new_pair         Particles object containing the new generated pairs
//! \param new_pair_species Species object representing the pair species
//! \param mBW_tables       Cross-section data tables and useful
//!                         functions for multiphoton Breit-Wheeler
//! \param pair_energy      Energy converted into pairs
//! \param istart           Index of the first particle
//! \param iend             Index of the last particle
//! \param ithread          Thread index
// ---------------------------------------------------------------------------------------------------------------------
void MultiphotonBreitWheeler::operator()( Particles &particles,
        SmileiMPI *                   smpi,
        Particles **                  new_pair,
        Species **                     new_pair_species,
        MultiphotonBreitWheelerTables &mBW_tables,
        double &                      pair_energy,
        int                           istart,
        int                           iend,
        int ithread, int ibin, int ipart_ref )
{
    // _______________________________________________________________
    // Parameters
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );

    // We use dynamics_invgf to store gamma
    double * const __restrict__ photon_gamma = &( smpi->dynamics_invgf[ithread][0] );

    const int nparts = smpi->getBufferSize(ithread);
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

    // Position shortcut
    double *const __restrict__ position_x = particles.getPtrPosition( 0 );
    double *const __restrict__ position_y = n_dimensions_ > 1 ? particles.getPtrPosition( 1 ) : nullptr;
    double *const __restrict__ position_z = n_dimensions_ > 2 ? particles.getPtrPosition( 2 ) : nullptr;

    // Particles Momentum shortcut
    double *const __restrict__ momentum_x = particles.getPtrMomentum(0);
    double *const __restrict__ momentum_y = particles.getPtrMomentum(1);
    double *const __restrict__ momentum_z = particles.getPtrMomentum(2);

    // Weight shortcut
    double *const __restrict__ weight = particles.getPtrWeight();

    // Optical depth for the Monte-Carlo process
    double *const __restrict__ tau =  particles.getPtrTau();

    // Quantum parameter
    double *const __restrict__ photon_chi = particles.getPtrChi();

    // Photon id
    // uint64_t * id = &( particles.id(0));

    // Reserve pair particles (else, pointer could become obsolete)
    // double np = new_pair[0]->size();
    // new_pair[0]->reserve( np + mBW_pair_creation_sampling_[0] * (iend - istart) );
    // new_pair[1]->reserve( np + mBW_pair_creation_sampling_[1] * (iend - istart) );

    // Pair shortcut
    double *const __restrict__ pair0_position_x = new_pair[0]->getPtrPosition( 0 );
    double *const __restrict__ pair0_position_y = (n_dimensions_ > 1 ? new_pair[0]->getPtrPosition( 1 ) : nullptr) ;
    double *const __restrict__ pair0_position_z = (n_dimensions_ > 2 ? new_pair[0]->getPtrPosition( 2 ) : nullptr) ;

    double *const __restrict__ pair0_position_old_x = particles.keepOldPositions() ? new_pair[0]->getPtrPositionOld( 0 ) : nullptr;
    double *const __restrict__ pair0_position_old_y = (particles.Position_old.size() > 1 ? new_pair[0]->getPtrPositionOld( 1 ) : nullptr) ;
    double *const __restrict__ pair0_position_old_z = (particles.Position_old.size() > 2 ? new_pair[0]->getPtrPositionOld( 2 ) : nullptr) ;

    double *const __restrict__ pair0_momentum_x = new_pair[0]->getPtrMomentum( 0 );
    double *const __restrict__ pair0_momentum_y = new_pair[0]->getPtrMomentum( 1 );
    double *const __restrict__ pair0_momentum_z = new_pair[0]->getPtrMomentum( 2 );

    double *const __restrict__ pair0_weight = new_pair[0]->getPtrWeight();
    short *const __restrict__ pair0_charge = new_pair[0]->getPtrCharge();

    double *const __restrict__ pair0_chi = new_pair[0]->has_quantum_parameter ? new_pair[0]->getPtrChi() : nullptr;
    double *const __restrict__ pair0_tau = new_pair[0]->has_Monte_Carlo_process ? new_pair[0]->getPtrTau() : nullptr;

    double *const __restrict__ pair1_position_x = new_pair[1]->getPtrPosition( 0 );
    double *const __restrict__ pair1_position_y = (n_dimensions_ > 1 ? new_pair[1]->getPtrPosition( 1 ) : nullptr);
    double *const __restrict__ pair1_position_z = (n_dimensions_ > 2 ? new_pair[1]->getPtrPosition( 2 ) : nullptr);

    double *const __restrict__ pair1_position_old_x = particles.keepOldPositions() ? new_pair[1]->getPtrPositionOld( 0 ) : nullptr;
    double *const __restrict__ pair1_position_old_y = (particles.Position_old.size() > 1 ? new_pair[1]->getPtrPositionOld( 1 ) : nullptr) ;
    double *const __restrict__ pair1_position_old_z = (particles.Position_old.size() > 2 ? new_pair[1]->getPtrPositionOld( 2 ) : nullptr) ;

    double *const __restrict__ pair1_momentum_x = new_pair[1]->getPtrMomentum( 0 );
    double *const __restrict__ pair1_momentum_y = new_pair[1]->getPtrMomentum( 1 );
    double *const __restrict__ pair1_momentum_z = new_pair[1]->getPtrMomentum( 2 );

    double *const __restrict__ pair1_weight = new_pair[1]->getPtrWeight();
    short *const __restrict__ pair1_charge = new_pair[1]->getPtrCharge();

    double *const __restrict__ pair1_chi = new_pair[1]->has_quantum_parameter ? new_pair[1]->getPtrChi() : nullptr;
    double *const __restrict__ pair1_tau = new_pair[1]->has_Monte_Carlo_process ? new_pair[1]->getPtrTau() : nullptr;

#ifdef SMILEI_ACCELERATOR_GPU_OACC
    // Parameters for random generator
    unsigned long long seed;
    unsigned long long seq;
    unsigned long long offset;
    // curandState_t state_1;
    // curandState_t state_2;
    // hiprandState_t state_1;
    // hiprandState_t state_2;
    
    seed = 12345ULL;
    seq = 0ULL;
    offset = 0ULL;

    // Parameters for linear alleatory number generator
    const int a = 1664525;
    const int c = 1013904223;
    const int m = std::pow(2,32);

    // Initialize initial seed for linear generator
    double initial_seed_1 = rand_->uniform();
    double initial_seed_2 = rand_->uniform();

    int seed_curand_1;
    int seed_curand_2;

    // _______________________________________________________________
    // Computation

    #pragma acc parallel \
    present(photon_gamma[0:nparts], \
    Ex[0:nparts],Ey[0:nparts],Ez[0:nparts],\
    Bx[0:nparts],By[0:nparts],Bz[0:nparts], \
    mBW_tables.T_.data_[0:mBW_tables.T_.size_], \
    mBW_tables.xi_.data_[0:mBW_tables.xi_.size_], \
    mBW_tables.xi_.axis1_min_[0:mBW_tables.xi_.dim_size_[0]]) \
    deviceptr(position_x, position_y, position_z, \
            momentum_x,momentum_y,momentum_z,charge,weight,tau,photon_chi, \
            pair0_position_x, pair0_position_y, pair0_position_z, \
            pair0_momentum_x, pair0_momentum_y, pair0_momentum_z, \
            pair0_weight, pair0_charge, pair0_chi, pair0_tau, \
            pair1_position_x, pair1_position_y, pair1_position_z, \
            pair1_momentum_x, pair1_momentum_y, pair1_momentum_z, \
            pair1_weight, pair1_charge, pair1_chi, pair1_tau \
    )
    {

        smilei::tools::gpu::Random prng_state_1;
        smilei::tools::gpu::Random prng_state_2;
        //curandState_t state_1;
        //curandState_t state_2;
        
        
        #pragma acc loop gang worker vector \
        private(seed_curand_1, seed_curand_2)

#else

    // 1. Computation of gamma and chi
    //    Can be vectorized
    #pragma omp simd
    
#endif
    
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
                                
#ifndef SMILEI_ACCELERATOR_GPU_OACC
                                
    }

    // 2. Monte-Carlo process
    //    No vectorized
    for( int ipart=istart ; ipart<iend; ipart++ ) {

#endif

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
                    
#ifndef SMILEI_ACCELERATOR_GPU_OACC
                    tau[ipart] = -std::log( 1.-rand_->uniform() );
#else
                    
                    seed_curand_1 = (int) (ipart+1)*(initial_seed_1+1); //Seed for linear generator
                    //seed_curand_1 = std::fmod(a * seed_curand_1 + c, m); //Linear generator
                    seed_curand_1 = (a * seed_curand_1 + c) % m; //Linear generator
           		
                    prng_state_1.init( seed_curand_1, seq, offset ); //Cuda generator initialization
                    //hiprand_init(seed_curand_1, seq, offset, &state_1); //Cuda generator initialization
                    //curand_init(seed_curand_1, seq, offset, &state_1); //Cuda generator initialization
                    
                    const double random_number = prng_state_1.uniform(); //Generating number
                    
                    tau[ipart] = -std::log( 1.-random_number );
                    initial_seed_1 = random_number;
#endif
                }

            }

            // Photon decay: emission under progress
            // If epsilon_tau_ > 0
            else if( tau[ipart] > epsilon_tau_ ) {
                // from the cross section
                temp = mBW_tables.computeBreitWheelerPairProductionRate( photon_chi[ipart], photon_gamma [ipart] );

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
                    // pair_energy += MultiphotonBreitWheeler::pair_emission( ipart,
                    //                                         particles,
                    //                                         ( *gamma )[ipart],
                    //                                         dt_ - event_time,
                    //                                         mBW_tables );

                    double inv_chiph_gammaph = ( photon_gamma[ipart]-2. ) / photon_chi[ipart];

                    // pair quantum parameters
                    double pair_chi[2];

                    // Draw random number in [0,1[
#ifndef SMILEI_ACCELERATOR_GPU_OACC
                    const double random_number = rand_->uniform();
#else
                    seed_curand_2 = (int) (ipart + 1)*(initial_seed_2 + 1); //Seed for linear generator
                    //seed_curand_2 = std::fmod(a * seed_curand_2 + c, m); //Linear generator
                    seed_curand_2 = (a * seed_curand_2 + c) % m; //Linear generator
                        
                    prng_state_2.init( seed_curand_2, seq, offset ); //Random generator initialization
	
                    const double random_number = prng_state_2.uniform(); //Generating number
                    //random_number = curand_uniform(&state_2); //Generating number
#endif

                    // Get the pair quantum parameters to compute the energy
                    mBW_tables.computePairQuantumParameter( photon_chi[ipart], &pair_chi[0], random_number );

                    // pair propagation direction // direction of the photon
                    double ux = momentum_x[ipart]/photon_gamma[ipart];
                    double uy = momentum_y[ipart]/photon_gamma[ipart];
                    double uz = momentum_z[ipart]/photon_gamma[ipart];

                    SMILEI_UNUSED( ibin );
                    // Creation of new electrons in the temporary array new_pair[0]
                    new_pair[0]->createParticles( mBW_pair_creation_sampling_[0] );
#ifndef SMILEI_ACCELERATOR_GPU_OACC
                    // Final size
                    int nparticles = new_pair[0]->size();

                    // Start index in the pair buffer
                    int i_pair_start = nparticles - mBW_pair_creation_sampling_[0];
#else
                    int i_pair_start = (istart + ipart)*mBW_pair_creation_sampling_[0];
#endif

                    // For all new paticles
#ifndef SMILEI_ACCELERATOR_GPU_OACC
                    #pragma omp simd
#endif
                    for( int ipair=i_pair_start; ipair < i_pair_start+mBW_pair_creation_sampling_[0]; ipair++ ) {

                        // Momentum
                        const double p = std::sqrt( ( 1.+pair_chi[0]*inv_chiph_gammaph)*( 1.+pair_chi[0]*inv_chiph_gammaph) - 1 );
                        pair0_momentum_x[ipair] = p*ux;
                        pair0_momentum_y[ipair] = p*uy;
                        pair0_momentum_z[ipair] = p*uz;

                        // gamma
                        //inv_gamma = 1./sqrt(1.+p*p);

                        // Positions
                        pair0_position_x[ipair]=position_x[ipart];
                        if (n_dimensions_>1) {
                            pair0_position_y[ipair]=position_y[ipart];
                            if (n_dimensions_>2) {
                                pair0_position_z[ipair]=position_z[ipart];
                            }
                        }
            //               + new_pair[k].momentum(i,ipair)*remaining_dt*inv_gamma;

#ifndef SMILEI_ACCELERATOR_GPU_OACC
                        // Old positions
                        if( particles.keepOldPositions() ) {
                            pair0_position_old_x[ipair]=position_x[ipart] ;
                            if (n_dimensions_>1) {
                                pair0_position_old_y[ipair]=position_y[ipart] ;
                            }
                            if (n_dimensions_>2) {
                                pair0_position_old_z[ipair]=position_z[ipart] ;
                            }
                        }
#endif

                        pair0_weight[ipair]=weight[ipart]*mBW_pair_creation_inv_sampling_[0];
                        pair0_charge[ipair]=new_pair_species[0]->max_charge_;

                        if( new_pair[0]->has_quantum_parameter ) {
                            pair0_chi[ipair] = pair_chi[0];
                        }

                        if( new_pair[0]->has_Monte_Carlo_process ) {
                            pair0_tau[ipair] = -1.;
                        }
                    }

                    // Create particle for the second pair species
                    new_pair[1]->createParticles( mBW_pair_creation_sampling_[1] );

#ifndef SMILEI_ACCELERATOR_GPU_OACC
                    // Final size
                    nparticles = new_pair[1]->size();

                    // Start index in the pair buffer
                    i_pair_start = nparticles - mBW_pair_creation_sampling_[1];
#else
                    i_pair_start = (istart + ipart)*mBW_pair_creation_sampling_[1];
#endif

                    // For all new paticles
#ifndef SMILEI_ACCELERATOR_GPU_OACC
                    #pragma omp simd
#endif
                    for( auto ipair=i_pair_start; ipair < i_pair_start + mBW_pair_creation_sampling_[1]; ipair++ ) {

                        // Momentum
                        const double p = std::sqrt( ( 1.+pair_chi[1]*inv_chiph_gammaph) * ( 1.+pair_chi[1]*inv_chiph_gammaph) - 1 );
                        pair1_momentum_x[ipair] = p*ux;
                        pair1_momentum_y[ipair] = p*uy;
                        pair1_momentum_z[ipair] = p*uz;

                        // gamma
                        //inv_gamma = 1./sqrt(1.+p*p);

                        // Positions

                        pair1_position_x[ipair]=position_x[ipart];
                        if (n_dimensions_>1) {
                            pair1_position_y[ipair]=position_y[ipart];
                            if (n_dimensions_>2) {
                                pair1_position_z[ipair]=position_z[ipart];
                            }
                        }
            //               + new_pair[k].momentum(i,ipair)*remaining_dt*inv_gamma;

#ifndef SMILEI_ACCELERATOR_GPU_OACC
                        // Old positions
                        if( particles.keepOldPositions() ) {
                            pair1_position_old_x[ipair]=position_x[ipart] ;
                            if (n_dimensions_>1) {
                                pair1_position_old_y[ipair]=position_y[ipart] ;
                            }
                            if (n_dimensions_>2) {
                                pair1_position_old_z[ipair]=position_z[ipart] ;
                            }
                        }
#endif

                        pair1_weight[ipair]=weight[ipart]*mBW_pair_creation_inv_sampling_[1];
                        pair1_charge[ipair]=new_pair_species[1]->max_charge_;

                        if( new_pair[1]->has_quantum_parameter ) {
                            pair1_chi[ipair] = pair_chi[1];
                        }

                        if( new_pair[1]->has_Monte_Carlo_process ) {
                            pair1_tau[ipair] = -1.;
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
    } // end ipart loop
    
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    }
#endif
}

// -----------------------------------------------------------------------------
//! Clean photons that decayed into pairs (weight <= 0)
//! \param particles   particle object containing the particle
//!                    properties of the current species
//! \param istart      Index of the first particle
//! \param iend        Index of the last particle
//! \param ithread     Thread index
// -----------------------------------------------------------------------------
void MultiphotonBreitWheeler::removeDecayedPhotons(
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

    const int nparts = smpi->getBufferSize(ithread);

    if( bmax[ibin] > bmin[ibin] ) {
        // Weight shortcut
        double *weight = &( particles.weight( 0 ) );

        // Backward loop over the photons to fing the first existing photon
        int last_photon_index = bmax[ibin]-1; // Index of the last existing photon (weight > 0)
        int first_photon_index = bmin[ibin];
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

                    int iDim_max_iold = n_dimensions_-1;
                    if (thetaold) iDim_max_iold = 1; // in AMcylindrical iold and deltaold are defined in 2D
                    for ( int iDim=iDim_max_iold ; iDim>=0 ; iDim-- ) {
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
        const unsigned int nb_deleted_photon = bmax[ibin]-last_photon_index-1;

        if( nb_deleted_photon > 0 ) {
            particles.eraseParticle( last_photon_index+1, nb_deleted_photon );
            // Erase bufferised data
            for ( int iDim=2 ; iDim>=0 ; iDim-- ) {
                Epart->erase(Epart->begin()+iDim*nparts+last_photon_index+1,Epart->begin()+iDim*nparts+last_photon_index+1+nb_deleted_photon);
                Bpart->erase(Bpart->begin()+iDim*nparts+last_photon_index+1,Bpart->begin()+iDim*nparts+last_photon_index+1+nb_deleted_photon);
            }

            int iDim_max_iold = n_dimensions_-1;
            if (thetaold) iDim_max_iold = 1; // in AMcylindrical iold and deltaold are defined in 2D
            for ( int iDim=iDim_max_iold ; iDim>=0 ; iDim-- ) {
                iold->erase(iold->begin()+iDim*nparts+last_photon_index+1,iold->begin()+iDim*nparts+last_photon_index+1+nb_deleted_photon);
                deltaold->erase(deltaold->begin()+iDim*nparts+last_photon_index+1,deltaold->begin()+iDim*nparts+last_photon_index+1+nb_deleted_photon);
            }
            gamma->erase(gamma->begin()+0*nparts+last_photon_index+1,gamma->begin()+0*nparts+last_photon_index+1+nb_deleted_photon);

            if (thetaold) {
                thetaold->erase(thetaold->begin()+0*nparts+last_photon_index+1,thetaold->begin()+0*nparts+last_photon_index+1+nb_deleted_photon);
            }

            bmax[ibin] = last_photon_index+1;
            for( int ii=ibin+1; ii<nbin; ii++ ) {
                bmin[ii] -= nb_deleted_photon;
                bmax[ii] -= nb_deleted_photon;
            }
        }
    }
}

// -----------------------------------------------------------------------------
    //! Clean photons that decayed into pairs (weight <= 0) and resize each bin
    //! But keeping the space between bins (so called no compression)
    //! \param particles   particle object containing the particle
    //!                    properties of the current species
    //! \param smpi        MPI properties
    //! \param ibin        Index of the current bin
    //! \param bmin        Pointer toward the first particle index of the bin in the Particles object
    //! \param bmax        Pointer toward the last particle index of the bin in the Particles object
    //! \param ithread     Thread index
// -----------------------------------------------------------------------------

void MultiphotonBreitWheeler::removeDecayedPhotonsWithoutBinCompression(
    Particles &particles,
    SmileiMPI *smpi,
    int ibin,
    int *bmin, int *bmax, int ithread )
{

    double *const Epart     = smpi->dynamics_Epart[ithread].data();
    double *const Bpart     = smpi->dynamics_Bpart[ithread].data();
    double *const gamma     = smpi->dynamics_invgf[ithread].data();
    int *const iold         = smpi->dynamics_iold[ithread].data();
    double *const deltaold  = smpi->dynamics_deltaold[ithread].data();

    const int nparts = smpi->getBufferSize(ithread);

    std::complex<double> * thetaold = NULL;
    if ( smpi->dynamics_eithetaold.size() )
        thetaold = smpi->dynamics_eithetaold[ithread].data();

    if( bmax[ibin] > bmin[ibin] ) {
        // Weight shortcut
        double *const weight = particles.getPtrWeight();
        //int nb_deleted_photon;

        // Backward loop over the photons to find the first existing photon
        int last_photon_index = bmax[ibin]-1; // Index of the last existing photon (weight > 0)
        int first_photon_index = bmin[ibin]; // Index of the first photon
        while( ( last_photon_index >= bmin[ibin] )
                && ( weight[last_photon_index] <= 0 ) ) {
            last_photon_index--;
        }
        while( ( first_photon_index < bmax[ibin] )
                && ( weight[first_photon_index] > 0 ) ) {
            first_photon_index++;
        }
        // At this level, last_photon_index is the position of the last still-existing photon (weight > 0)
        // that will not be erased

        // Backward loop over the photons to fill holes in the photon particle array (at the bin level only)
        for( int ipart=last_photon_index-1 ; ipart>=bmin[ibin]; ipart-- ) {
            if( weight[ipart] <= 0 ) {
                if( ipart < last_photon_index ) {
                    // The last existing photon comes to the position of
                    // the deleted photon
#ifndef SMILEI_ACCELERATOR_GPU_OACC
                    particles.overwriteParticle( last_photon_index, ipart );
#else
#endif
                    // Overwrite bufferised data
                    for ( int iDim=2 ; iDim>=0 ; iDim-- ) {
                        Epart[iDim*nparts+ipart] = Epart[iDim*nparts+last_photon_index];
                        Bpart[iDim*nparts+ipart] = Bpart[iDim*nparts+last_photon_index];
                    }
                    for ( int iDim=n_dimensions_-1 ; iDim>=0 ; iDim-- ) {
                        iold[iDim*nparts+ipart] = iold[iDim*nparts+last_photon_index];
                        deltaold[iDim*nparts+ipart] = deltaold[iDim*nparts+last_photon_index];
                    }
                    gamma[ipart] = gamma[0*nparts+last_photon_index];

                    if (thetaold) {
                        thetaold[0*nparts+ipart] = thetaold[0*nparts+last_photon_index];
                    }
                    last_photon_index --;
                }
            }
        } // end for ipart

        // Update of the bin boundaries
        // const unsigned int nb_deleted_photon = bmax[ibin]-last_photon_index-1;

        // We suppress the deleted photons at the end of the bin
        if( last_photon_index + 1 < bmax[ibin] ) {
            bmax[ibin] = last_photon_index+1;

            // std::cerr
            //         << " ibin: " << ibin
            //         << " - bmin: " << bmin[ibin]
            //         << " - bmax: " << bmax[ibin]
            //         << " - nb_deleted_photon: " << nb_deleted_photon
            //         << std::endl;
        }

    } // if bmax[ibin] > bmin[ibin]
}

// void MultiphotonBreitWheeler::joinNewElectronPositronPairs(Particles **new_pair, unsigned int Nbins)
// {

//     for( int k=0 ; k < 2 ; k++ ) {
//        for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
//            // number of particles to add from the bin
//            unsigned int nparticles_to_add = new_pair_per_bin[ibin][k].size();
//            new_pair[k]->createParticles(nparticles_to_add);

//            for (unsigned int ipart = 0; ipart < nparticles_to_add ; ipart++){

//               int idNew = (new_pair[k]->size() - nparticles_to_add) + ipart;

//               // momenta
//               for( int i=0; i<3; i++ ) {
//                   new_pair[k]->momentum( i, idNew ) = new_pair_per_bin[ibin][k].momentum( i, ipart );
//               }

//               // positions
//               for( int i=0; i<n_dimensions_; i++ ) {
//                   new_pair[k]->position( i, idNew ) = new_pair_per_bin[ibin][k].position( i, ipart );
//               }

//               // old positions
//               if( new_pair[k]->Position_old.size() > 0 ) {
//                   for( int i=0; i<n_dimensions_; i++ ) {
//                       new_pair[k]->position_old( i, idNew ) = new_pair_per_bin[ibin][k].position_old( i, ipart );
//                   }
//               }

//               // weight
//               new_pair[k]->weight( idNew ) = new_pair_per_bin[ibin][k].weight( ipart );

//               // charge
//               new_pair[k]->charge( idNew ) = new_pair_per_bin[ibin][k].charge( ipart );

//               // chi
//               if( new_pair[k]->has_quantum_parameter ) {
//                   new_pair[k]->chi( idNew ) = new_pair_per_bin[ibin][k].chi( ipart );
//               }

//               //tau
//               if( new_pair[k]->has_Monte_Carlo_process ) {
//                   new_pair[k]->tau( idNew ) = new_pair_per_bin[ibin][k].tau( ipart );
//               }
//            } // end ipart
//            new_pair_per_bin[ibin][k].clear();
//        } // end ibin loop
//     } // end k loop (different species loop)

// } // end joinNewElectronPositronPairs
