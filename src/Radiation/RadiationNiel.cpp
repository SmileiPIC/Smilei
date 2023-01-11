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
//! \param radiation_tables Cross-section data tables and useful functions
//                     for nonlinear inverse Compton scattering
//! \param istart      Index of the first particle
//! \param iend        Index of the last particle
//! \param ithread     Thread index
//! \param radiated_energy     overall energy radiated during the call to this method
// -----------------------------------------------------------------------------
void RadiationNiel::operator()(
    Particles       &particles,
    Particles       */*photons*/,
    SmileiMPI       *smpi,
    RadiationTables &radiation_tables,
    double          &radiated_energy,
    int istart,
    int iend,
    int ithread,
    int /*ibin*/,
    int ipart_ref )
{

    // _______________________________________________________________
    // Parameters
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );

    const int nparts = smpi->getBufferSize(ithread);
    const double *const __restrict__ Ex = &( ( *Epart )[0*nparts] );
    const double *const __restrict__ Ey = &( ( *Epart )[1*nparts] );
    const double *const __restrict__ Ez = &( ( *Epart )[2*nparts] );
    const double *const __restrict__ Bx = &( ( *Bpart )[0*nparts] );
    const double *const __restrict__ By = &( ( *Bpart )[1*nparts] );
    const double *const __restrict__ Bz = &( ( *Bpart )[2*nparts] );

    // Used to store gamma directly
    double *const __restrict__ gamma = &( smpi->dynamics_invgf[ithread][0] );

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

    // Radiated energy
    double rad_energy;

    // Stochastic diffusive term for Niel et al.
    double * diffusion = new double [nbparticles];

    // Random Number
    double * random_numbers = new double [nbparticles];

    // Momentum shortcut
    double*const __restrict__ momentum_x = particles.getPtrMomentum(0);
    double*const __restrict__ momentum_y = particles.getPtrMomentum(1);
    double*const __restrict__ momentum_z = particles.getPtrMomentum(2);

    // Charge shortcut
    const short*const __restrict__ charge = particles.getPtrCharge();

    // Weight shortcut
    const double*const __restrict__ weight = particles.getPtrWeight();

    // Quantum parameter
    double*const __restrict__ particle_chi = particles.getPtrChi();

    const double minimum_chi_continuous           = radiation_tables.getMinimumChiContinuous();
    const double factor_classical_radiated_power  = radiation_tables.getFactorClassicalRadiatedPower();
    const int niel_computation_method             = radiation_tables.getNielHComputationMethodIndex();

    // _______________________________________________________________
    // Computation

    //double t0 = MPI_Wtime();

    // Vectorized computation of gamma and the particle quantum parameter
    #pragma omp simd
    for( ipart=istart ; ipart< iend; ipart++ ) {

        const double charge_over_mass_square = ( double )( charge[ipart] )*one_over_mass_square;

        // Gamma
        gamma[ipart-ipart_ref] = std::sqrt( 1.0 + momentum_x[ipart]*momentum_x[ipart]
                             + momentum_y[ipart]*momentum_y[ipart]
                             + momentum_z[ipart]*momentum_z[ipart] );

        // Computation of the Lorentz invariant quantum parameter
        particle_chi[ipart] = Radiation::computeParticleChi( charge_over_mass_square,
                              momentum_x[ipart], momentum_y[ipart], momentum_z[ipart],
                              gamma[ipart],
                              Ex[ipart-ipart_ref], Ey[ipart-ipart_ref], Ez[ipart-ipart_ref],
                              Bx[ipart-ipart_ref], By[ipart-ipart_ref], Bz[ipart-ipart_ref] );
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
        if( particle_chi[ipart+istart] > minimum_chi_continuous ) {

            // Pick a random number in the normal distribution of standard
            // deviation sqrt(dt_) (variance dt_)
            //random_numbers[ipart] = 2.*Rand::uniform() -1.;
            //random_numbers[ipart] = 2.*drand48() -1.;
            random_numbers[ipart] = rand_->uniform2();
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
        if( particle_chi[ipart+istart] > minimum_chi_continuous ) {
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

    //double t2 = MPI_Wtime();

    // Computation of the diffusion coefficients
    // Using the table (non-vectorized)
    if( niel_computation_method == 0 ) {
        // #pragma omp simd
        for( ipart=0 ; ipart < nbparticles; ipart++ ) {

            // Below particle_chi = minimum_chi_continuous_, radiation losses are negligible
            if( particle_chi[ipart+istart] > minimum_chi_continuous ) {

                //h = radiation_tables.getHNielFitOrder10(particle_chi[ipart]);
                //h = radiation_tables.getHNielFitOrder5(particle_chi[ipart]);
                //temp = radiation_tables.getHNielFromTable( particle_chi[ipart+istart] );
                temp = radiation_tables.niel_.get( particle_chi[ipart+istart] );

                diffusion[ipart] = std::sqrt( factor_classical_radiated_power*gamma[ipart+istart-ipart_ref]*temp )*random_numbers[ipart];
            }
        }
    }
    // Using the fit at order 5 (vectorized)
    else if( niel_computation_method == 1 ) {
        #pragma omp simd private(temp)
        for( ipart=0 ; ipart < nbparticles; ipart++ ) {

            int ipartp = ipart + istart;

            // Below particle_chi = minimum_chi_continuous_, radiation losses are negligible
            if( particle_chi[ipartp] > minimum_chi_continuous ) {

                temp = RadiationTools::getHNielFitOrder5( particle_chi[ipartp] );

                diffusion[ipart] = std::sqrt( factor_classical_radiated_power*gamma[ipartp-ipart_ref]*temp )*random_numbers[ipart];
            }
        }
    }
    // Using the fit at order 10 (vectorized)
    else if( niel_computation_method == 2 ) {
        #pragma omp simd private(temp)
        for( ipart=0 ; ipart < nbparticles; ipart++ ) {

            int ipartp = ipart + istart;

            // Below particle_chi = minimum_chi_continuous_, radiation losses are negligible
            if( particle_chi[ipartp] > minimum_chi_continuous ) {

                temp = RadiationTools::getHNielFitOrder10( particle_chi[ipartp] );

                diffusion[ipart] = std::sqrt( factor_classical_radiated_power*gamma[ipartp-ipart_ref]*temp )*random_numbers[ipart];
            }
        }
    }
    // Using Ridgers
    else if( niel_computation_method == 3) {

        #pragma omp simd private(temp)
        for( ipart=0 ; ipart < nbparticles; ipart++ ) {

            int ipartp = ipart + istart;

            // Below particle_chi = minimum_chi_continuous_, radiation losses are negligible
            if( particle_chi[ipartp] > minimum_chi_continuous ) {

                temp = RadiationTools::getHNielFitRidgers( particle_chi[ipartp] );

                diffusion[ipart] = std::sqrt( factor_classical_radiated_power*gamma[ipartp-ipart_ref]*temp )*random_numbers[ipart];
            }
        }
    }
    //double t3 = MPI_Wtime();

    // Vectorized update of the momentum
    #pragma omp simd private(temp,rad_energy)
    for( ipart=istart ; ipart<iend; ipart++ ) {
        // Below particle_chi = minimum_chi_continuous_, radiation losses are negligible
        if( gamma[ipart-ipart_ref] > 1.1 && particle_chi[ipart] > minimum_chi_continuous ) {

            // Radiated energy during the time step
            rad_energy =
                radiation_tables.getRidgersCorrectedRadiatedEnergy( particle_chi[ipart], dt_ );

            // Effect on the momentum
            // Temporary factor
            temp = ( rad_energy - diffusion[ipart-istart] )
                   * gamma[ipart-ipart_ref]/( gamma[ipart-ipart_ref]*gamma[ipart-ipart_ref]-1. );

            // Update of the momentum
            momentum_x[ipart] -= temp*momentum_x[ipart];
            momentum_y[ipart] -= temp*momentum_y[ipart];
            momentum_z[ipart] -= temp*momentum_z[ipart];

        }
    }

    //double t4 = MPI_Wtime();

    // ____________________________________________________
    // Vectorized computation of the thread radiated energy
    // and update of the quantum parameter

    double radiated_energy_loc = 0;
    double new_gamma = 0;

    #pragma omp simd private(new_gamma) reduction(+:radiated_energy_loc)
    for( int ipart=istart ; ipart<iend; ipart++ ) {

    const double charge_over_mass_square = ( double )( charge[ipart] )*one_over_mass_square;

        new_gamma = std::sqrt( 1.0
                       + momentum_x[ipart]*momentum_x[ipart]
                       + momentum_y[ipart]*momentum_y[ipart]
                       + momentum_z[ipart]*momentum_z[ipart] );

        radiated_energy_loc += weight[ipart]*( gamma[ipart] - new_gamma );

        particle_chi[ipart] = Radiation::computeParticleChi( charge_over_mass_square,
                     momentum_x[ipart], momentum_y[ipart], momentum_z[ipart],
                     new_gamma,
                     Ex[ipart-ipart_ref], Ey[ipart-ipart_ref], Ez[ipart-ipart_ref],
                     Bx[ipart-ipart_ref], By[ipart-ipart_ref], Bz[ipart-ipart_ref] );

    }
    radiated_energy += radiated_energy_loc;

    // Destruction
    delete [] diffusion;
    delete [] random_numbers;

    //double t5 = MPI_Wtime();

    //std::cerr << "" << std::endl;
    //std::cerr << "" << istart << " " << nbparticles << " " << ithread << std::endl;
    //std::cerr << "Computation of chi: " << t1 - t0 << std::endl;
    //std::cerr << "Computation of random numbers: " << t2 - t1 << std::endl;
    //std::cerr << "Computation of the diffusion: " << t3 - t2 << std::endl;
    //std::cerr << "Computation of the momentum: " << t4 - t3 << std::endl;
    //std::cerr << "Computation of the radiated energy: " << t5 - t4 << std::endl;
}
