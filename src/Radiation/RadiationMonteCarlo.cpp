// ----------------------------------------------------------------------------
//! \file RadiationMonteCarlo.cpp
//
//! \brief This class performs the Nonlinear Inverse Compton Scattering
//! on particles using a Monte-Carlo approach.
//
//! The implementation is adapted from the thesis results of M. Lobet
//! See http://www.theses.fr/2015BORD0361
//
// ----------------------------------------------------------------------------

#include "RadiationMonteCarlo.h"

#include <cstring>
#include <fstream>

#if defined(SMILEI_ACCELERATOR_GPU_OACC)
    #define __HIP_PLATFORM_NVCC__
    #define __HIP_PLATFORM_NVIDIA__
    #include "gpuRandom.h"
#elif defined(SMILEI_ACCELERATOR_GPU_OMP)
    #define __HIP_PLATFORM_HCC__
    #define __HIP_PLATFORM_AMD__
    #include "gpuRandom.h"
#endif


// ---------------------------------------------------------------------------------------------------------------------
//! Constructor for RadiationMonteCarlo
//! Inherit from Radiation
// ---------------------------------------------------------------------------------------------------------------------
RadiationMonteCarlo::RadiationMonteCarlo( Params &params, Species *species, Random * rand  )
    : Radiation( params, species, rand )
{
    radiation_photon_sampling_        = species->radiation_photon_sampling_;
    max_photon_emissions_             = species->radiation_max_emissions_;
    radiation_photon_gamma_threshold_ = species->radiation_photon_gamma_threshold_;
    inv_radiation_photon_sampling_    = 1. / radiation_photon_sampling_;
}

// ---------------------------------------------------------------------------------------------------------------------
//! Destructor for RadiationMonteCarlo
// ---------------------------------------------------------------------------------------------------------------------
RadiationMonteCarlo::~RadiationMonteCarlo()
{
}

// ---------------------------------------------------------------------------------------------------------------------
//! Overloading of the operator (): perform the Discontinuous radiation reaction
//! induced by the nonlinear inverse Compton scattering
//
//! \param particles   particle object containing the particle properties
//! \param photons     Particles object that will receive emitted photons
//! \param smpi        MPI properties
//! \param radiation_tables Cross-section data tables and useful functions
//                     for nonlinear inverse Compton scattering
//! \param istart      Index of the first particle
//! \param iend        Index of the last particle
//! \param ithread     Thread index
//! \param radiated_energy     overall energy radiated during the call to this method
// ---------------------------------------------------------------------------------------------------------------------
void RadiationMonteCarlo::operator()(
    Particles       &particles,
    Particles       *photons,
    SmileiMPI       *smpi,
    RadiationTables &radiation_tables,
    double          &radiated_energy,
    int             istart,
    int             iend,
    int             ithread,
    int             ibin,
    int             ipart_ref)
{
    SMILEI_UNUSED( ibin );
    // _______________________________________________________________
    // Parameters

    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    //std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);

    // Total number of particles
    const int nparts = smpi->getBufferSize(ithread);

    const double *const __restrict__ Ex = &( ( *Epart )[0*nparts] );
    const double *const __restrict__ Ey = &( ( *Epart )[1*nparts] );
    const double *const __restrict__ Ez = &( ( *Epart )[2*nparts] );
    const double *const __restrict__ Bx = &( ( *Bpart )[0*nparts] );
    const double *const __restrict__ By = &( ( *Bpart )[1*nparts] );
    const double *const __restrict__ Bz = &( ( *Bpart )[2*nparts] );


    // 1 / mass^2
    const double one_over_mass_square = one_over_mass_*one_over_mass_;

    // Radiated energy
    double cont_rad_energy;

    // Temporary double parameter
    double temp;

#ifdef SMILEI_ACCELERATOR_GPU_OACC
    unsigned long long seed; // Parameters for CUDA generator
    unsigned long long seq;
    unsigned long long offset;
    // curandState_t state_1;
    // curandState_t state_2;
    // hiprandState_t state_1;
    // hiprandState_t state_2;
    
    seed = 12345ULL;
    seq = 0ULL;
    offset = 0ULL;
#endif

    // Parameter to store the local radiated energy
    double radiated_energy_loc = 0;

    //random temporary number
    double random_number;

    // Particle properties ----------------------------------------------------------------

    // Particles position shortcut
    double *const __restrict__ position_x = particles.getPtrPosition( 0 );
    double *const __restrict__ position_y = nDim_ > 1 ? particles.getPtrPosition( 1 ) : nullptr;
    double *const __restrict__ position_z = nDim_ > 2 ? particles.getPtrPosition( 2 ) : nullptr;

    // Particles Momentum shortcut
    double *const __restrict__ momentum_x = particles.getPtrMomentum(0);
    double *const __restrict__ momentum_y = particles.getPtrMomentum(1);
    double *const __restrict__ momentum_z = particles.getPtrMomentum(2);

    // Charge shortcut
    const short *const __restrict__ charge = particles.getPtrCharge();

    // Weight shortcut
    const double *const __restrict__ weight = particles.getPtrWeight();

    // Optical depth for the Monte-Carlo process
    double *const __restrict__ tau = particles.getPtrTau();

    // Quantum parameter
    double *const __restrict__ chi = particles.getPtrChi();

    
    // Photon properties ----------------------------------------------------------------

    // Number of photons
    int nphotons;
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    int nphotons_start;
#endif
    
    // Buffer size for each particle
    const double photon_buffer_size_per_particle = radiation_photon_sampling_ * max_photon_emissions_;
    
    if (photons) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            // We reserve a large number of potential photons on device since we can't reallocate
            nphotons_start = photons->deviceSize();
            //static_cast<nvidiaParticles*>(photons)->deviceReserve( nphotons + (iend - istart) * photon_buffer_size_per_particle );
            static_cast<nvidiaParticles*>(photons)->createParticles( (iend - istart) * photon_buffer_size_per_particle );
            //std::cerr << "photons size: " << static_cast<nvidiaParticles*>(photons)->deviceSize()
            //          << " new: " << (iend - istart)*photon_buffer_size_per_particle
            //          << std::endl;

#else
            nphotons = photons->size();
            // We reserve a large number of photons
            photons->reserve( nphotons + (iend - istart) * photon_buffer_size_per_particle );
#endif
    } else {
        nphotons = 0;
    }

    // Photon position shortcut
    double *const __restrict__ photon_position_x = photons ? photons->getPtrPosition( 0 ) : nullptr;
    double *const __restrict__ photon_position_y = photons ? (nDim_ > 1 ? photons->getPtrPosition( 1 ) : nullptr) : nullptr;
    double *const __restrict__ photon_position_z = photons ? (nDim_ > 2 ? photons->getPtrPosition( 2 ) : nullptr) : nullptr;

    // Particles Momentum shortcut
    double *const __restrict__ photon_momentum_x = photons ? photons->getPtrMomentum(0) : nullptr;
    double *const __restrict__ photon_momentum_y = photons ? photons->getPtrMomentum(1) : nullptr;
    double *const __restrict__ photon_momentum_z = photons ? photons->getPtrMomentum(2) : nullptr;

    // Charge shortcut
    short *const __restrict__ photon_charge = photons ? photons->getPtrCharge() : nullptr;

    // Weight shortcut
    double *const __restrict__ photon_weight = photons ? photons->getPtrWeight() : nullptr;

    // Quantum Parameter
    double *const __restrict__ photon_chi_array = photons ? (photons->has_quantum_parameter ? photons->getPtrChi() : nullptr) : nullptr;

    double *const __restrict__ photon_tau = photons ? (photons->has_Monte_Carlo_process ? photons->getPtrTau() : nullptr) : nullptr;

#ifdef SMILEI_ACCELERATOR_GPU_OACC
    // Cell keys as a mask
    int *const __restrict__ photon_cell_keys = photons ? photons->getPtrCellKeys() : nullptr;
#endif

    // Table properties ----------------------------------------------------------------
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    // Size of tables
    // int size_of_Table_integfochi = RadiationTables.integfochi_.size_particle_chi_;
    // int size_of_Table_min_photon_chi = RadiationTables.xi_.size_particle_chi_;
    // int size_of_Table_xi = RadiationTables.xi_.size_particle_chi_*
    //                        RadiationTables.xi_.size_photon_chi_;
#endif


    // Tables for MC
    // const double *const table_integfochi = &(radiation_tables.integfochi_.data_[0]);
    // const double *const table_min_photon_chi = &(radiation_tables.xi_.axis1_min_[0]);
    // double * table_xi = &(radiation_tables.xi_.table_[0]);

    // _______________________________________________________________
    // Computation
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    // Management of the data on GPU though this data region
    int np = iend-istart;
    
    // Initialize initial seed for linear generator
    double initial_seed_1 = rand_->uniform();
    double initial_seed_2 = rand_->uniform();

    // Parameters for linear alleatory number generator
    const int a = 1664525;
    const int c = 1013904223;
    const int m = std::pow(2,32);

    // Variable to save seed for CUDA generators
    int seed_curand_1;
    int seed_curand_2;
    
    #pragma acc data present(Ex[istart:np],Ey[istart:np],Ez[istart:np],\
            Bx[istart:np],By[istart:np],Bz[istart:np], \
            radiation_tables.integfochi_.data_[0:radiation_tables.integfochi_.size_], \
            radiation_tables.xi_.data_[0:radiation_tables.xi_.size_], \
            radiation_tables.xi_.axis1_min_[0:radiation_tables.xi_.dim_size_[0]]) \
            deviceptr(momentum_x,momentum_y,momentum_z,position_x, \
            position_y,position_z,charge,weight,tau,chi, \
            photon_position_x, \
            photon_position_y, \
            photon_position_z, \
            photon_momentum_x, \
            photon_momentum_y, \
            photon_momentum_z, \
            photon_weight, \
            photon_charge, \
            photon_chi_array, \
            photon_tau, \
            photon_cell_keys \
            )
    {
    #pragma acc parallel \
    present(Ex[istart:np],Ey[istart:np],Ez[istart:np],\
    Bx[istart:np],By[istart:np],Bz[istart:np], \
    radiation_tables.integfochi_.data_[0:radiation_tables.integfochi_.size_], \
    radiation_tables.xi_.data_[0:radiation_tables.xi_.size_], \
    radiation_tables.xi_.axis1_min_[0:radiation_tables.xi_.dim_size_[0]]) \
    deviceptr(position_x, \
            position_y, \
            position_z, \
            momentum_x,momentum_y,momentum_z,charge,weight,tau,chi, \
            photon_position_x, \
            photon_position_y, \
            photon_position_z, \
            photon_momentum_x, \
            photon_momentum_y, \
            photon_momentum_z, \
            photon_weight, \
            photon_charge, \
            photon_chi_array, \
            photon_tau, \
            photon_cell_keys, \
    )
    {
        
        smilei::tools::gpu::Random prng_state_1;
        smilei::tools::gpu::Random prng_state_2;
        //curandState_t state_1;
        //curandState_t state_2;
        
        #pragma acc loop gang worker vector \
        private(random_number, seed_curand_1, seed_curand_2) \
        reduction(+:radiated_energy_loc)

#endif

    for( int ipart=istart ; ipart<iend; ipart++ ) {

        // charge / mass^2
        const double charge_over_mass_square = ( double )( charge[ipart] )*one_over_mass_square;

        // Time to emission
        double emission_time = 0;

        // time spent in the iteration
        double local_it_time = 0;

        // Number of Monte-Carlo iteration
        int mc_it_nb = 0;

        // Number of emitted photons per particles
        int i_photon_emission = 0;

        // Monte-Carlo Manager inside the time step
        while( ( local_it_time < dt_ )
                &&( mc_it_nb < max_monte_carlo_iterations_ ) ) {

            // Gamma
            const double particle_gamma = std::sqrt( 1.0 + momentum_x[ipart]*momentum_x[ipart]
                          + momentum_y[ipart]*momentum_y[ipart]
                          + momentum_z[ipart]*momentum_z[ipart] );
            // does not apply the MC routine for particles with 0 kinetic energy
            if( particle_gamma < 1.1 ){
                break;
            }

            // Computation of the Lorentz invariant quantum parameter
            const double particle_chi = Radiation::computeParticleChi( charge_over_mass_square,
                           momentum_x[ipart], momentum_y[ipart], momentum_z[ipart],
                           particle_gamma,
                           Ex[ipart-ipart_ref], Ey[ipart-ipart_ref], Ez[ipart-ipart_ref],
                           Bx[ipart-ipart_ref], By[ipart-ipart_ref], Bz[ipart-ipart_ref] );

            // Update the quantum parameter in species
            // chi[ipart] = particle_chi;

            // Discontinuous emission: New emission
            // If tau[ipart] <= 0, this is a new emission
            // We also check that particle_chi > chipa_threshold,
            // else particle_chi is too low to induce a discontinuous emission
            if( ( particle_chi > radiation_tables.getMinimumChiDiscontinuous() )
                    && ( tau[ipart] <= epsilon_tau_ ) ) {
                // New final optical depth to reach for emision
                while( tau[ipart] <= epsilon_tau_ ) {
                    //tau[ipart] = -log( 1.-Rand::uniform() );
                    #ifndef SMILEI_ACCELERATOR_GPU_OACC
                        tau[ipart] = -std::log( 1.-rand_->uniform() );
                    #else
                        seed_curand_1 = (int) (ipart+1)*(initial_seed_1+1); //Seed for linear generator
                        seed_curand_1 = (a * seed_curand_1 + c) % m; //Linear generator
               		
                        prng_state_1.init( seed_curand_1, seq, offset ); //Cuda generator initialization
                        //hiprand_init(seed_curand_1, seq, offset, &state_1); //Cuda generator initialization
                        //curand_init(seed_curand_1, seq, offset, &state_1); //Cuda generator initialization
                        
                        random_number = prng_state_1.uniform(); //Generating number
                        //random_number = hiprand_uniform(&state_1); //Generating number
                        //random_number = curand_uniform(&state_1); //Generating number
                        
                        tau[ipart] = -std::log( 1.- random_number );
                        initial_seed_1 = random_number;
                    #endif
                }

            }

            // Discontinuous emission: emission under progress
            // If epsilon_tau_ > 0
            if( tau[ipart] > epsilon_tau_ ) {

                // from the cross section
                temp = radiation_tables.computePhotonProductionYield(
                                              particle_chi,
                                              particle_gamma);

                // Time to discontinuous emission
                // If this time is > the remaining iteration time,
                // we have a synchronization
                emission_time = std::min( tau[ipart]/temp, dt_ - local_it_time );

                // Update of the optical depth
                tau[ipart] -= temp*emission_time;

                // If the final optical depth is reached, photons are emitted
                if( tau[ipart] <= epsilon_tau_ ) {


                    // Draw random number in [0,1[
                    #ifndef SMILEI_ACCELERATOR_GPU_OACC
                        random_number = rand_->uniform();
                    #else
                        seed_curand_2 = (int) (ipart + 1)*(initial_seed_2 + 1); //Seed for linear generator
                        seed_curand_2 = (a * seed_curand_2 + c) % m; //Linear generator
                        
                        prng_state_2.init( seed_curand_2, seq, offset ); //Random generator initialization
	
                        random_number = prng_state_2.uniform(); //Generating number
                        //random_number = curand_uniform(&state_2); //Generating number
                    #endif

                    // Emission of a photon without tasks
                    // Radiated energy is incremented only if the macro-photon is not created

                    // Get the photon quantum parameter from the table xip
                    // photon_chi = RadiationTables.computeRandomPhotonChi( particle_chi );
                    double photon_chi = radiation_tables.computeRandomPhotonChiWithInterpolation( particle_chi, random_number);

                    // compute the photon gamma factor
                    double photon_gamma = photon_chi/particle_chi*( particle_gamma-1.0 );

                    // *****************************************************************
                    // Creation of the new photon

                    // Update of the particle properties
                    // direction d'emission // direction de l'electron (1/gamma << 1)
                    // With momentum conservation
                    double inv_old_norm_p = photon_gamma/std::sqrt( particle_gamma*particle_gamma - 1.0 );
                    momentum_x[ipart] -= momentum_x[ipart]*inv_old_norm_p;
                    momentum_y[ipart] -= momentum_y[ipart]*inv_old_norm_p;
                    momentum_z[ipart] -= momentum_z[ipart]*inv_old_norm_p;

                    // With energy conservation
                    /*inv_old_norm_p = 1./sqrt(particle_gamma*particle_gamma - 1.0);
                    particle_gamma -= photon_gamma;
                    new_norm_p = sqrt(particle_gamma*particle_gamma - 1.0);
                    px *= new_norm_p * inv_old_norm_p;
                    py *= new_norm_p * inv_old_norm_p;
                    pz *= new_norm_p * inv_old_norm_p;*/

                    // Creation of macro-photons if requested
                    // Check that the photons is defined and the threshold on the energy
                    if(          photons
                            && ( photon_gamma >= radiation_photon_gamma_threshold_ )
                            && ( i_photon_emission < max_photon_emissions_)) {
                                
// CPU implementation (non-threaded implementation)
#ifndef SMILEI_ACCELERATOR_GPU_OACC

                        // Creation of new photons in the temporary array photons
                        photons->createParticles( radiation_photon_sampling_ );

                        // New number of photons
                        nphotons += radiation_photon_sampling_;

                        // Inverse of the momentum norm
                        inv_old_norm_p = 1./std::sqrt( momentum_x[ipart]*momentum_x[ipart]
                                                  + momentum_y[ipart]*momentum_y[ipart]
                                                  + momentum_z[ipart]*momentum_z[ipart] );

                        // For all new photons
                        #pragma omp simd
                        for( auto iphoton=nphotons-radiation_photon_sampling_; iphoton<nphotons; iphoton++ ) {

                            // std::cerr  << photons << " "
                            //            << iphoton << " "
                            //            << photons->size() << " "
                            //            << radiation_photon_sampling_ << " "
                            //            << ipart << " "
                            //            << photon_position_x << " "
                            //            << photons->getPtrPosition( 0 ) << " "
                            //            << std::endl;

                            photon_position_x[iphoton]=position_x[ipart];
                            if (nDim_>1) {
                                photon_position_y[iphoton]=position_y[ipart];
                                if (nDim_>2) {
                                    photon_position_z[iphoton]=position_z[ipart];
                                }
                            }

                            photon_momentum_x[iphoton] =
                                photon_gamma*momentum_x[ipart]*inv_old_norm_p;
                            photon_momentum_y[iphoton] =
                                photon_gamma*momentum_y[ipart]*inv_old_norm_p;
                            photon_momentum_z[iphoton] =
                                photon_gamma*momentum_z[ipart]*inv_old_norm_p;

                            photon_weight[iphoton] = weight[ipart]*inv_radiation_photon_sampling_;
                            photon_charge[iphoton] = 0;

                            if( photons->has_quantum_parameter ) {
                                photon_chi_array[iphoton] = photon_chi;
                            }

                            if( photons->has_Monte_Carlo_process ) {
                                photon_tau[iphoton] = -1.;
                            }

                        } // end for iphoton

                        // Number of emitted photons

                        i_photon_emission += 1;
                        
// GPU optimized implementation (SIMT compatible implementation)
// Each particle has a buffer of `max_photon_emissions_` possible emission
// meaning `max_photon_emissions_*radiation_photon_sampling_` photons
#else

                        // Inverse of the momentum norm
                        inv_old_norm_p = 1./std::sqrt( momentum_x[ipart]*momentum_x[ipart]
                                                  + momentum_y[ipart]*momentum_y[ipart]
                                                  + momentum_z[ipart]*momentum_z[ipart] );

                        const int iphoton_start = nphotons_start // initial number of photons
                                  + (ipart - istart) * photon_buffer_size_per_particle // beginning of the buffer for ipart
                                  + i_photon_emission * radiation_photon_sampling_; // already emitted photons (i.e. buffer usage)


                        // For all new photons
                        for( auto iphoton=0; iphoton<radiation_photon_sampling_; iphoton++ ) {
                            
                            // Photon positions at particle positions
                            photon_position_x[iphoton_start + iphoton]=position_x[ipart];
                            if (nDim_>1) {
                                photon_position_y[iphoton_start + iphoton]=position_y[ipart];
                                if (nDim_>2) {
                                    photon_position_z[iphoton_start + iphoton]=position_z[ipart];
                                }
                            }
                            
                            // Photon momentum
                            photon_momentum_x[iphoton_start + iphoton] =
                                photon_gamma*momentum_x[ipart]*inv_old_norm_p;
                            photon_momentum_y[iphoton_start + iphoton] =
                                photon_gamma*momentum_y[ipart]*inv_old_norm_p;
                            photon_momentum_z[iphoton_start + iphoton] =
                                photon_gamma*momentum_z[ipart]*inv_old_norm_p;
                                
                            photon_weight[iphoton_start + iphoton] = weight[ipart]*inv_radiation_photon_sampling_;
                            photon_charge[iphoton_start + iphoton] = 0;

                            if( photons->has_quantum_parameter ) {
                                photon_chi_array[iphoton_start + iphoton] = photon_chi;
                            }

                            if( photons->has_Monte_Carlo_process ) {
                                photon_tau[iphoton_start + iphoton] = -1.;
                            }
                                
                            // cell_keys plays the role of a mask here
                            // If the photon is created, then cell_keys is True
                            photon_cell_keys[iphoton_start + iphoton] = 1;
                            
                        }
                        
                        i_photon_emission += 1;
                        
#endif

                    }
                    // If no emission of a macro-photon:
                    // Addition of the emitted energy in the cumulating parameter
                    // for the scalar diagnostics
                    else {
                        photon_gamma = particle_gamma - std::sqrt( 1.0 + momentum_x[ipart]*momentum_x[ipart]
                                                         + momentum_y[ipart]*momentum_y[ipart]
                                                         + momentum_z[ipart]*momentum_z[ipart] );
                        radiated_energy_loc += weight[ipart]*photon_gamma;
                    }

                    // Optical depth becomes negative meaning
                    // that a new drawing is possible
                    // at the next Monte-Carlo iteration
                    tau[ipart] = -1.;
                } // end if tau

                // Incrementation of the Monte-Carlo iteration counter
                mc_it_nb ++;
                // Update of the local time
                local_it_time += emission_time;

            // Continuous emission
            // particle_chi needs to be below the discontinuous threshold
            // particle_chi needs to be above the continuous thresholdF
            // No discontiuous emission is in progress:
            // tau[ipart] <= epsilon_tau_
            }
            else if( particle_chi <=  radiation_tables.getMinimumChiDiscontinuous()
                     && tau[ipart] <= epsilon_tau_
                     && particle_chi >  radiation_tables.getMinimumChiContinuous() ) {


                // Remaining time of the iteration
                emission_time = dt_ - local_it_time;

                // Radiated energy during emission_time
                cont_rad_energy =
                    radiation_tables.getRidgersCorrectedRadiatedEnergy( particle_chi,
                            emission_time );

                // Effect on the momentum
                temp = cont_rad_energy*particle_gamma/( particle_gamma*particle_gamma-1. );
                momentum_x[ipart] -= temp*momentum_x[ipart];
                momentum_y[ipart] -= temp*momentum_y[ipart];
                momentum_z[ipart] -= temp*momentum_z[ipart];

                // Incrementation of the radiated energy cumulative parameter
                radiated_energy_loc += weight[ipart]*( particle_gamma - std::sqrt( 1.0
                                                    + momentum_x[ipart]*momentum_x[ipart]
                                                    + momentum_y[ipart]*momentum_y[ipart]
                                                    + momentum_z[ipart]*momentum_z[ipart] ) );

                // End for this particle
                local_it_time = dt_;
            }
            // No emission since particle_chi is too low
            else { // if (particle_chi < radiation_tables.getMinimumChiContinuous())
                local_it_time = dt_;

            } // end if
        } // end while
    } // end for

#ifdef SMILEI_ACCELERATOR_GPU_OACC
    } // end acc parallel
#endif

    //if (photons) std::cerr << photons->deviceSize()  << std::endl;

    // Remove extra space to save memory
#ifndef SMILEI_ACCELERATOR_GPU_OACC
    if (photons) {
        photons->shrinkToFit( true );
    }
#endif

    // Update the patch radiated energy
    radiated_energy += radiated_energy_loc;
    //std::cerr << " " << radiated_energy << std::endl;
    
    // ____________________________________________________
    // Update of the quantum parameter chi

#ifndef SMILEI_ACCELERATOR_GPU_OACC
        #pragma omp simd
#else
    int np = iend-istart;
    #pragma acc parallel present(Ex[istart:np],Ey[istart:np],Ez[istart:np],\
    Bx[istart:np],By[istart:np],Bz[istart:np]) \
    deviceptr(momentum_x,momentum_y,momentum_z, charge,weight,chi) \
    private(gamma)
    {

        #pragma acc loop gang worker vector
#endif
        for( int ipart=istart ; ipart<iend; ipart++ ) {
            const double charge_over_mass_square = ( double )( charge[ipart] )*one_over_mass_square;

            // Gamma
            const double particle_gamma = std::sqrt( 1.0 + momentum_x[ipart]*momentum_x[ipart]
                          + momentum_y[ipart]*momentum_y[ipart]
                          + momentum_z[ipart]*momentum_z[ipart] );

            // Computation of the Lorentz invariant quantum parameter
            chi[ipart] = Radiation::computeParticleChi( charge_over_mass_square,
                         momentum_x[ipart], momentum_y[ipart], momentum_z[ipart],
                         particle_gamma,
                         Ex[ipart-ipart_ref], Ey[ipart-ipart_ref], Ez[ipart-ipart_ref],
                         Bx[ipart-ipart_ref], By[ipart-ipart_ref], Bz[ipart-ipart_ref] );

        }

    #ifdef SMILEI_ACCELERATOR_GPU_OACC
    } // end acc parallel
    #endif

#ifdef SMILEI_ACCELERATOR_GPU_OACC
    }   // end acc data
#endif

}
