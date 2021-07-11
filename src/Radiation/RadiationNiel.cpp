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
RadiationNiel::RadiationNiel( Params &params, Species *species, Random * rand  )
    : Radiation( params, species, rand )
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
//! **Vectorized version** But needs to be improved
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
void RadiationNiel::operator()(
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

    int nparts = Epart->size()/3;
    double *Ex = &( ( *Epart )[0*nparts] );
    double *Ey = &( ( *Epart )[1*nparts] );
    double *Ez = &( ( *Epart )[2*nparts] );
    double *Bx = &( ( *Bpart )[0*nparts] );
    double *By = &( ( *Bpart )[1*nparts] );
    double *Bz = &( ( *Bpart )[2*nparts] );

    // Used to store gamma directly
    double *gamma = &( smpi->dynamics_invgf[ithread][0] );

    // Charge divided by the square of the mass
    double charge_over_mass_square = 0.;

    // 1/mass^2
    const double one_over_mass_square = one_over_mass_*one_over_mass_;

    // Sqrt(dt_), used intensively in these loops
    const double sqrtdt = std::sqrt( dt_ );

    // Number of particles
    const int nbparticles = iend-istart;

    // Temporary double parameter
    double temp;

    // Particle id
    int ipart;
    int ipartz = 0; //Starting from zero

    // Radiated energy
    double rad_energy;

    // Stochastic diffusive term fo Niel et al.
    double diffusion[nbparticles];

    // Random Number
    double random_numbers[nbparticles];

    // Momentum shortcut
    double* momentum_x = particles.getPtrMomentum(0);
    double* momentum_y = particles.getPtrMomentum(1);
    double* momentum_z = particles.getPtrMomentum(2);

    // Charge shortcut
    short *charge = particles.getPtrCharge();

    // Weight shortcut
    double *weight = particles.getPtrWeight();

    // Quantum parameter
    double* particle_chi = particles.getPtrChi();

    // Niel table
    double* table = &(RadiationTables.niel_.table_[0]);

    const double minimum_chi_continuous_ = RadiationTables.getMinimumChiContinuous();
    const double factor_classical_radiated_power_      = RadiationTables.getFactorClassicalRadiatedPower();
    const int niel_computation_method = RadiationTables.getNielHComputationMethodIndex();
    const int size_of_table_Niel = RadiationTables.niel_.size_particle_chi_;

    // Parameter to store the local radiated energy
    double radiated_energy_loc = 0;
    double new_gamma = 0;
    double ipa = 0;

    #ifdef _GPU
    unsigned long long seed;
    unsigned long long seq;
    unsigned long long offset;
    #endif

    #ifdef _GPU
    // Management of the data on GPU though this data region
    int np = iend-istart;
    #pragma acc data create(random_numbers[0:nbparticles], diffusion[0:nbparticles]) \
            present(Ex[istart:np],Ey[istart:np],Ez[istart:np],\
            Bx[istart:np],By[istart:np],Bz[istart:np],gamma[istart:np], table[0:size_of_table_Niel]) \
            deviceptr(momentum_x,momentum_y,momentum_z,charge,weight,particle_chi)
    {
    #endif

    // _______________________________________________________________
    // Computation
    // 1) Vectorized computation of gamma and the particle quantum parameter
        #ifndef _GPU
            #pragma omp simd
        #else
            #pragma acc parallel \
            present(Ex[istart:np],Ey[istart:np],Ez[istart:np],\
            Bx[istart:np],By[istart:np],Bz[istart:np],gamma[istart:np], \
            random_numbers[0:nbparticles], diffusion[0:nbparticles], table[0:size_of_table_Niel]) \
            deviceptr(momentum_x,momentum_y,momentum_z,charge,weight,particle_chi) \
            private(temp,ipartp, rad_energy,new_gamma) reduction(+:radiated_energy_loc)
        {
            curandState_t state;
            seed = 12345ULL;
            seq = 0ULL;
            offset = 0ULL;
            curand_init(seed, seq, offset, &state);
            #pragma acc loop gang worker vector
    #endif
        for( ipart=istart ; ipart< iend; ipart++ ) {

            charge_over_mass_square = ( double )( charge[ipart] )*one_over_mass_square;

            // Gamma
            gamma[ipart-ipart_ref] = sqrt( 1.0 + momentum_x[ipart]*momentum_x[ipart]
                                 + momentum_y[ipart]*momentum_y[ipart]
                                 + momentum_z[ipart]*momentum_z[ipart] );

            // Computation of the Lorentz invariant quantum parameter
            particle_chi[ipart] = Radiation::computeParticleChi( charge_over_mass_square,
                                  momentum_x[ipart], momentum_y[ipart], momentum_z[ipart],
                                  gamma[ipart],
                                  ( *( Ex+ipart-ipart_ref ) ), ( *( Ey+ipart-ipart_ref ) ), ( *( Ez+ipart-ipart_ref ) ),
                                  ( *( Bx+ipart-ipart_ref ) ), ( *( By+ipart-ipart_ref ) ), ( *( Bz+ipart-ipart_ref ) ) );

    //double t1 = MPI_Wtime();

    // 2) Computation of the random number


    #ifdef _GPU

    /*#pragma acc parallel num_gangs(1) present(random_numbers[0:nbparticles]) private(state, p,temp)
    {
        #pragma acc loop seq*/
        /*for( ipart=0 ; ipart < nbparticles; ipart++ ) {*/
            if( particle_chi[ipart] > minimum_chi_continuous_ ) {

                random_numbers[ipart] = std::sqrt( 2. )*sqrtdt*curand_normal(&state);
            
             /*}
         }*/
    
    //}

    #else

    // Vectorized computation of the random number in a uniform distribution
    // #pragma omp simd
    }
    
    for( ipart=0 ; ipart < nbparticles; ipart++ ) {

        // Below particle_chi = minimum_chi_continuous_, radiation losses are negligible
        if( particle_chi[ipart+istart] > minimum_chi_continuous_ ) {

            random_numbers[ipart] = 2.*rand_->uniform() -1.;
        
        }
    }

    // Vectorized computation of the random number in a normal distribution
    double p;
    #pragma omp simd private(p,temp)
    for( ipart=0 ; ipart < nbparticles; ipart++ ) {
        // Below particle_chi = minimum_chi_continuous_, radiation losses are negligible
        if( particle_chi[ipart+istart] > minimum_chi_continuous_ ) {
            temp = -std::log( ( 1.0-random_numbers[ipart] )*( 1.0+random_numbers[ipart] ) );

            if( temp < 5.000000 ) {
                temp = temp - 2.500000;
                p = +2.81022636000e-08      ;
                p = +3.43273939000e-07 + p*temp;
                p = -3.52338770000e-06 + p*temp;
                p = -4.39150654000e-06 + p*temp;
                p = +0.00021858087e+00 + p*temp;
                p = -0.00125372503e+00 + p*temp;
                p = -0.00417768164e+00 + p*temp;
                p = +0.24664072700e+00 + p*temp;
                p = +1.50140941000e+00 + p*temp;
            } else {
                temp = std::sqrt( temp ) - 3.000000;
                p = -0.000200214257      ;
                p = +0.000100950558 + p*temp;
                p = +0.001349343220 + p*temp;
                p = -0.003673428440 + p*temp;
                p = +0.005739507730 + p*temp;
                p = -0.007622461300 + p*temp;
                p = +0.009438870470 + p*temp;
                p = +1.001674060000 + p*temp;
                p = +2.832976820000 + p*temp;
            }

            random_numbers[ipart] *= p*sqrtdt*std::sqrt( 2. );

            //random_numbers[ipart] = userFunctions::erfinv2(random_numbers[ipart])*sqrtdt*sqrt(2.);
        }
    }
    #pragma omp simd private(temp,rad_energy,new_gamma) reduction(+:radiated_energy_loc)
    for( ipartz=0 ; ipartz < nbparticles; ipartz++ ) {    

        if( particle_chi[ipartz+istart] > minimum_chi_continuous_ ) {
    #endif

    //double t2 = MPI_Wtime();

    // 3) Computation of the diffusion coefficients

    new_gamma = 0;

            if( niel_computation_method == 0 ) { // Using the table (non-vectorized)
         
                temp = RadiationTables.getHNielFromTable( particle_chi[ipartz+istart], table );
                
                diffusion[ipartz] = sqrt( factor_classical_radiated_power_*gamma[ipartz+istart-ipart_ref]*temp )*random_numbers[ipartz];

            }else if( niel_computation_method == 1 ) { // Using the fit at order 5 (vectorized)

                temp = RadiationTools::getHNielFitOrder5( particle_chi[ipartz + istart] );

                diffusion[ipartz] = sqrt( factor_classical_radiated_power_*gamma[ipartz + istart-ipart_ref]*temp )*random_numbers[ipartz];

            }else if( niel_computation_method == 2 ) { // Using the fit at order 10 (vectorized)
 
                temp = RadiationTools::getHNielFitOrder10( particle_chi[ipartz + istart] );

                diffusion[ipartz] = sqrt( factor_classical_radiated_power_*gamma[ipartz + istart-ipart_ref]*temp )*random_numbers[ipartz];
 
           }else if( niel_computation_method == 3) { // Using Ridgers

                temp = RadiationTools::getHNielFitRidgers( particle_chi[ipartz + istart] );

                diffusion[ipartz] = sqrt( factor_classical_radiated_power_*gamma[ipartz + istart-ipart_ref]*temp )*random_numbers[ipartz];
            }
                //double t3 = MPI_Wtime();

                 // 4) Vectorized update of the momentum
                // Radiated energy during the time step
                rad_energy = RadiationTables.getRidgersCorrectedRadiatedEnergy( particle_chi[ipart], dt_ );

                // Effect on the momentum
                // Temporary factor
                temp = ( rad_energy - diffusion[ipart])
                       * gamma[ipart-ipart_ref]/( gamma[ipart-ipart_ref]*gamma[ipart-ipart_ref]-1. );

                // Update of the momentum
                momentum_x[ipart] -= temp*momentum_x[ipart];
                momentum_y[ipart] -= temp*momentum_y[ipart];
                momentum_z[ipart] -= temp*momentum_z[ipart];

        }

        charge_over_mass_square = ( double )( charge[ipart] )*one_over_mass_square;

        new_gamma = sqrt( 1.0
                    + momentum_x[ipart]*momentum_x[ipart]
                    + momentum_y[ipart]*momentum_y[ipart]
                    + momentum_z[ipart]*momentum_z[ipart] );

        radiated_energy_loc += weight[ipart]*( gamma[ipart] - new_gamma );

        particle_chi[ipart] = Radiation::computeParticleChi( charge_over_mass_square,
                    momentum_x[ipart], momentum_y[ipart], momentum_z[ipart],
                    new_gamma,
                    ( *( Ex+ipart-ipart_ref ) ), ( *( Ey+ipart-ipart_ref ) ), ( *( Ez+ipart-ipart_ref ) ),
                    ( *( Bx+ipart-ipart_ref ) ), ( *( By+ipart-ipart_ref ) ), ( *( Bz+ipart-ipart_ref ) ) );

    #ifdef _GPU
    ipartz++;
    #endif

    }

    #ifdef _GPU
    } // end acc parallel loop
    #endif

    //double t5 = MPI_Wtime();

    #ifdef _GPU
    }   // end acc data
    #endif

    // Update the patch radiated energy
    radiated_energy += radiated_energy_loc;

    //std::cerr << "" << std::endl;
    //std::cerr << "" << istart << " " << nbparticles << " " << ithread << std::endl;
    //std::cerr << "Computation of chi: " << t1 - t0 << std::endl;
    //std::cerr << "Computation of random numbers: " << t2 - t1 << std::endl;
    //std::cerr << "Computation of the diffusion: " << t3 - t2 << std::endl;
    //std::cerr << "Computation of the momentum: " << t4 - t3 << std::endl;
    //std::cerr << "Computation of the radiated energy: " << t5 - t4 << std::endl;

}