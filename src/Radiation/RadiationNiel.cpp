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
    double charge_over_mass_square;

    // 1/mass^2
    const double one_over_mass_square = pow( one_over_mass_, 2. );

    // Sqrt(dt_), used intensively in these loops
    const double sqrtdt = sqrt( dt_ );

    // Number of particles
    const int nbparticles = iend-istart;

    // Temporary double parameter
    double temp;

    // Particle id
    int ipart;

    // Radiated energy
    double rad_energy;

    // Stochastic diffusive term fo Niel et al.
    double diffusion[nbparticles];

    // Random Number
    double random_numbers[nbparticles];

    // Momentum shortcut
    double *momentum[3];
    for( int i = 0 ; i<3 ; i++ ) {
        momentum[i] =  &( particles.momentum( i, istart ) );
    }

    // Charge shortcut
    short *charge = &( particles.charge( istart ) );

    // Weight shortcut
    double *weight = &( particles.weight( istart ) );

    // Quantum parameter
    double *particle_chi = &( particles.chi( istart ) );

    const double minimum_chi_continuous_ = RadiationTables.getMinimumChiContinuous();
    const double factor_classical_radiated_power_      = RadiationTables.getFactorClassicalRadiatedPower();
    const std::string niel_computation_method = RadiationTables.getNielHComputationMethod();
    
    // _______________________________________________________________
    // Computation

    //double t0 = MPI_Wtime();

    // Vectorized computation of gamma and the particle quantum parameter
    #pragma omp simd
    for( ipart=0 ; ipart< nbparticles; ipart++ ) {

        charge_over_mass_square = ( double )( charge[ipart] )*one_over_mass_square;

        // Gamma
        gamma[ipart] = sqrt( 1.0 + momentum[0][ipart]*momentum[0][ipart]
                             + momentum[1][ipart]*momentum[1][ipart]
                             + momentum[2][ipart]*momentum[2][ipart] );

        // Computation of the Lorentz invariant quantum parameter
        particle_chi[ipart] = Radiation::computeParticleChi( charge_over_mass_square,
                              momentum[0][ipart], momentum[1][ipart], momentum[2][ipart],
                              gamma[ipart],
                              ( *( Ex+ipart-ipart_ref ) ), ( *( Ey+ipart-ipart_ref ) ), ( *( Ez+ipart-ipart_ref ) ),
                              ( *( Bx+ipart-ipart_ref ) ), ( *( By+ipart-ipart_ref ) ), ( *( Bz+ipart-ipart_ref ) ) );
    }

    //double t1 = MPI_Wtime();

    // Non-vectorized computation of the random number
    /*for (ipart=0 ; ipart < nbparticles; ipart++ )
    {

        // Below particle_chi = minimum_chi_continuous_, radiation losses are negligible
        if (particle_chi[ipart] > minimum_chi_continuous_)
        {

          // Pick a random number in the normal distribution of standard
          // deviation sqrt(dt_) (variance dt_)
          random_numbers[ipart] = Rand::normal(sqrtdt);
        }
    }*/

    // Vectorized computation of the random number in a uniform distribution
    // #pragma omp simd
    for( ipart=0 ; ipart < nbparticles; ipart++ ) {

        // Below particle_chi = minimum_chi_continuous_, radiation losses are negligible
        if( particle_chi[ipart] > minimum_chi_continuous_ ) {

            // Pick a random number in the normal distribution of standard
            // deviation sqrt(dt_) (variance dt_)
            //random_numbers[ipart] = 2.*Rand::uniform() -1.;
            //random_numbers[ipart] = 2.*drand48() -1.;
            random_numbers[ipart] = 2.*rand_->uniform() -1.;
        }
        // else {
        //     random_numbers[ipart] = 0;
        // }
    }

    // Vectorized computation of the random number in a normal distribution
    double p;
    #pragma omp simd private(p,temp)
    for( ipart=0 ; ipart < nbparticles; ipart++ ) {
        // Below particle_chi = minimum_chi_continuous_, radiation losses are negligible
        if( particle_chi[ipart] > minimum_chi_continuous_ ) {
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

            random_numbers[ipart] *= p*sqrtdt*sqrt( 2. );

            //random_numbers[ipart] = userFunctions::erfinv2(random_numbers[ipart])*sqrtdt*sqrt(2.);
        }
    }

    //double t2 = MPI_Wtime();

    // Computation of the diffusion coefficients
    // Using the table (non-vectorized)
    if( niel_computation_method == "table" ) {
        // #pragma omp simd
        for( ipart=0 ; ipart < nbparticles; ipart++ ) {

            // Below particle_chi = minimum_chi_continuous_, radiation losses are negligible
            if( particle_chi[ipart] > minimum_chi_continuous_ ) {

                //h = RadiationTables.getHNielFitOrder10(particle_chi[ipart]);
                //h = RadiationTables.getHNielFitOrder5(particle_chi[ipart]);
                temp = RadiationTables.getHNielFromTable( particle_chi[ipart] );

                diffusion[ipart] = sqrt( factor_classical_radiated_power_*gamma[ipart]*temp )*random_numbers[ipart];
            }
        }
    }
    // Using the fit at order 5 (vectorized)
    else if( niel_computation_method == "fit5" ) {
        #pragma omp simd private(temp)
        for( ipart=0 ; ipart < nbparticles; ipart++ ) {

            // Below particle_chi = minimum_chi_continuous_, radiation losses are negligible
            if( particle_chi[ipart] > minimum_chi_continuous_ ) {

                temp = RadiationTools::getHNielFitOrder5( particle_chi[ipart] );

                diffusion[ipart] = sqrt( factor_classical_radiated_power_*gamma[ipart]*temp )*random_numbers[ipart];
            }
        }
    }
    // Using the fit at order 10 (vectorized)
    else if( niel_computation_method == "fit10" ) {
        #pragma omp simd private(temp)
        for( ipart=0 ; ipart < nbparticles; ipart++ ) {

            // Below particle_chi = minimum_chi_continuous_, radiation losses are negligible
            if( particle_chi[ipart] > minimum_chi_continuous_ ) {

                temp = RadiationTools::getHNielFitOrder10( particle_chi[ipart] );

                diffusion[ipart] = sqrt( factor_classical_radiated_power_*gamma[ipart]*temp )*random_numbers[ipart];
            }
        }
    }
    // Using Ridgers
    else if( niel_computation_method == "ridgers" ) {
    
        #pragma omp simd private(temp)
        for( ipart=0 ; ipart < nbparticles; ipart++ ) {

            // Below particle_chi = minimum_chi_continuous_, radiation losses are negligible
            if( particle_chi[ipart] > minimum_chi_continuous_ ) {

                temp = RadiationTools::getHNielFitRidgers( particle_chi[ipart] );

                diffusion[ipart] = sqrt( factor_classical_radiated_power_*gamma[ipart]*temp )*random_numbers[ipart];
            }
        }
    }
    //double t3 = MPI_Wtime();

    // Vectorized update of the momentum
    #pragma omp simd private(temp,rad_energy)
    for( ipart=0 ; ipart<nbparticles; ipart++ ) {
        // Below particle_chi = minimum_chi_continuous_, radiation losses are negligible
        if( particle_chi[ipart] > minimum_chi_continuous_ ) {

            // Radiated energy during the time step
            rad_energy =
                RadiationTables.getRidgersCorrectedRadiatedEnergy( particle_chi[ipart], dt_ );

            // Effect on the momentum
            // Temporary factor
            temp = ( rad_energy - diffusion[ipart] )
                   * gamma[ipart]/( gamma[ipart]*gamma[ipart]-1. );

            // Update of the momentum
            momentum[0][ipart] -= temp*momentum[0][ipart];
            momentum[1][ipart] -= temp*momentum[1][ipart];
            momentum[2][ipart] -= temp*momentum[2][ipart];

        }
    }

    //double t4 = MPI_Wtime();

    // ____________________________________________________
    // Vectorized computation of the thread radiated energy
    // and update of the quantum parameter
    
    double radiated_energy_loc = 0;
    double new_gamma = 0;

    #pragma omp simd private(new_gamma) reduction(+:radiated_energy_loc)
    for( int ipart=0 ; ipart<nbparticles; ipart++ ) {
        
        new_gamma = sqrt( 1.0
                       + momentum[0][ipart]*momentum[0][ipart]
                       + momentum[1][ipart]*momentum[1][ipart]
                       + momentum[2][ipart]*momentum[2][ipart] );
        
        radiated_energy_loc += weight[ipart]*( gamma[ipart] - new_gamma );
        
        particle_chi[ipart] = Radiation::computeParticleChi( charge_over_mass_square,
                     momentum[0][ipart], momentum[1][ipart], momentum[2][ipart],
                     new_gamma,
                     ( *( Ex+ipart-ipart_ref ) ), ( *( Ey+ipart-ipart_ref ) ), ( *( Ez+ipart-ipart_ref ) ),
                     ( *( Bx+ipart-ipart_ref ) ), ( *( By+ipart-ipart_ref ) ), ( *( Bz+ipart-ipart_ref ) ) );
                     
    }
    radiated_energy += radiated_energy_loc;

    //double t5 = MPI_Wtime();

    //std::cerr << "" << std::endl;
    //std::cerr << "" << istart << " " << nbparticles << " " << ithread << std::endl;
    //std::cerr << "Computation of chi: " << t1 - t0 << std::endl;
    //std::cerr << "Computation of random numbers: " << t2 - t1 << std::endl;
    //std::cerr << "Computation of the diffusion: " << t3 - t2 << std::endl;
    //std::cerr << "Computation of the momentum: " << t4 - t3 << std::endl;
    //std::cerr << "Computation of the radiated energy: " << t5 - t4 << std::endl;
}
