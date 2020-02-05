// ----------------------------------------------------------------------------
//! \file RadiationTables.h
//
//! \brief This class contains the tables and the functions to generate them
//! for the Nonlinear Inverse Compton Scattering
//
//! \details This header contains the definition of the class RadiationTables.
//! The implementation is adapted from the thesis results of M. Lobet
//! See http://www.theses.fr/2015BORD0361
// ----------------------------------------------------------------------------

#ifndef RADIATIONTABLES_H
#define RADIATIONTABLES_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <iomanip>
#include <cmath>
#include "userFunctions.h"
#include "Params.h"
#include "H5.h"

//------------------------------------------------------------------------------
//! RadiationTables class: holds parameters, tables and functions to compute
//! cross-sections,
//! optical depths and other useful parameters for the Compton Monte-Carlo
//! pusher.
//------------------------------------------------------------------------------
class RadiationTables
{


public:

    //! Constructor for RadiationTables
    RadiationTables();

    //! Destructor for RadiationTables
    ~RadiationTables();

    //! Initialization of the parmeters for the nonlinear
    //! inverse Compton scattering
    void initializeParameters( Params &params , SmileiMPI *smpi );

    // ---------------------------------------------------------------------
    // PHYSICAL COMPUTATION
    // ---------------------------------------------------------------------
    
    //! Computation of the photon production yield dNph/dt which is
    //! also the cross-section for the Monte-Carlo
    double computePhotonProductionYield( double particle_chi, double particle_gamma );

    //! Determine randomly a photon quantum parameter photon_chi
    //! for an emission process
    //! from a particle chi value (particle_chi) and
    //! using the tables xip and chiphmin
    //! \param particle_chi particle quantum parameter
    double computeRandomPhotonChi( double particle_chi );

    //! Return the value of the function h(particle_chi) of Niel et al.
    //! Use an integration of Gauss-Legendre
    //
    //! \param particle_chi particle quantum parameter
    //! \param nb_iterations number of iterations for the Gauss-Legendre integration
    //! \param eps epsilon for the modified bessel function
    double computeHNiel( double particle_chi, int nb_iterations, double eps );

    //! Return the value of the function h(particle_chi) of Niel et al.
    //! from the computed table niel_.table
    //! \param particle_chi particle quantum parameter
    double getHNielFromTable( double particle_chi );

    //! Return the stochastic diffusive component of the pusher
    //! of Niel et al.
    //! \param gamma particle Lorentz factor
    //! \param particle_chi particle quantum parameter
    //! \param dt time step
    double getNielStochasticTerm( double gamma,
                                  double particle_chi,
                                  double dt );

    //! Computation of the corrected continuous quantum radiated energy
    //! during dt from the quantum parameter particle_chi using the Ridgers
    //! formulae.
    //! \param particle_chi particle quantum parameter
    //! \param dt time step
    //#pragma omp declare simd
    double inline getRidgersCorrectedRadiatedEnergy( double particle_chi,
            double dt )
    {
        return computeRidgersFit( particle_chi )*dt*particle_chi*particle_chi*factor_classical_radiated_power_;
    };

    //! Get of the classical continuous radiated energy during dt
    //! \param particle_chi particle quantum parameter
    //! \param dt time step
    double inline getClassicalRadiatedEnergy( double particle_chi, double dt )
    {
        return dt*particle_chi*particle_chi*factor_classical_radiated_power_;
    };

    //! Return the minimum_chi_discontinuous_ value
    //! Under this value, no discontinuous radiation reaction
    double inline getMinimumChiDiscontinuous()
    {
        return minimum_chi_discontinuous_;
    }

    //! Return the minimum_chi_continuous_ value
    //! Under this value, no continuous radiation reaction
    double inline getMinimumChiContinuous()
    {
        return minimum_chi_continuous_;
    }

    //! Computation of the function g of Erber using the Ridgers
    //! approximation formulae
    //! \param particle_chi particle quantum parameter
    //#pragma omp declare simd
    double inline computeRidgersFit( double particle_chi )
    {
        return std::pow( 1.0 + 4.8*( 1.0+particle_chi )*std::log( 1.0 + 1.7*particle_chi )
                    + 2.44*particle_chi*particle_chi, -2.0/3.0 );
    };

    std::string inline getNielHComputationMethod()
    {
        return this->niel_.computation_method;
    }

    // -----------------------------------------------------------------------------
    //! Return the value of the function h(particle_chi) of Niel et al.
    //! from a polynomial numerical fit at order 10
    //! Valid between particle_chi in 1E-3 and 1E1
    //! \param particle_chi particle quantum parameter
    // -----------------------------------------------------------------------------
    double inline getHNielFitOrder10( double particle_chi )
    {
        // Max relative error ~2E-4
        return exp( -3.231764974833856e-08 * pow( std::log( particle_chi ), 10 )
                    -7.574417415366786e-07 * pow( std::log( particle_chi ), 9 )
                    -5.437005218419013e-06 * pow( std::log( particle_chi ), 8 )
                    -4.359062260446135e-06 * pow( std::log( particle_chi ), 7 )
                    + 5.417842511821415e-05 * pow( std::log( particle_chi ), 6 )
                    -1.263905701127627e-04 * pow( std::log( particle_chi ), 5 )
                    + 9.899812622393002e-04 * pow( std::log( particle_chi ), 4 )
                    + 1.076648497464146e-02 * pow( std::log( particle_chi ), 3 )
                    -1.624860613422593e-01 * pow( std::log( particle_chi ), 2 )
                    + 1.496340836237785e+00 * std::log( particle_chi )
                    -2.756744141581370e+00 );
    }

    // -----------------------------------------------------------------------------
    //! Return the value of the function h(particle_chi) of Niel et al.
    //! from a polynomial numerical fit at order 5
    //! Valid between particle_chi in 1E-3 and 1E1
    //! \param particle_chi particle quantum parameter
    // -----------------------------------------------------------------------------
    double inline getHNielFitOrder5( double particle_chi )
    {
        // Max relative error ~0.02
        return exp( 1.399937206900322e-04 * pow( std::log( particle_chi ), 5 )
                    + 3.123718241260330e-03 * pow( std::log( particle_chi ), 4 )
                    + 1.096559086628964e-02 * pow( std::log( particle_chi ), 3 )
                    -1.733977278199592e-01 * pow( std::log( particle_chi ), 2 )
                    + 1.492675770100125e+00 * std::log( particle_chi )
                    -2.748991631516466e+00 );
    }

    // -----------------------------------------------------------------------------
    //! Return the value of the function h(particle_chi) of Niel et al.
    //! using the numerical fit of Ridgers in
    //! Ridgers et al., ArXiv 1708.04511 (2017)
    //! \param particle_chi particle quantum parameter
    // -----------------------------------------------------------------------------
    double inline getHNielFitRidgers( double particle_chi )
    {
        return pow( particle_chi, 3 )*1.9846415503393384*pow( 1. +
                ( 1. + 4.528*particle_chi )*std::log( 1.+12.29*particle_chi ) + 4.632*pow( particle_chi, 2 ), -7./6. );
    }

    // -----------------------------------------------------------------------------
    //! Return the classical power factor factor_classical_radiated_power_.
    // -----------------------------------------------------------------------------
    double inline getFactorClassicalRadiatedPower()
    {
        return factor_classical_radiated_power_;
    }

    // ---------------------------------------------------------------------
    // TABLE READING
    // ---------------------------------------------------------------------

    //! Read the external table h
    //! \param smpi Object of class SmileiMPI containing MPI properties
    void readHTable( SmileiMPI *smpi );

    //! Read the external table integfochi
    //! \param smpi Object of class SmileiMPI containing MPI properties
    void readIntegfochiTable( SmileiMPI *smpi );

    //! Read the external table xip_chiphmin and xip
    //! \param smpi Object of class SmileiMPI containing MPI properties
    void readXiTable( SmileiMPI *smpi );

    //! Read the external all external tables for the radiation
    //! \param smpi Object of class SmileiMPI containing MPI properties
    void readTables( Params &params, SmileiMPI *smpi );

    // ---------------------------------------------------------------------
    // TABLE COMMUNICATIONS
    // ---------------------------------------------------------------------

    //! Bcast of the external table h
    //! \param smpi Object of class SmileiMPI containing MPI properties
    void bcastHTable( SmileiMPI *smpi );

    //! Bcast of the external table integfochi
    //! \param smpi Object of class SmileiMPI containing MPI properties
    void bcastIntegfochiTable( SmileiMPI *smpi );

    //! Bcast of the external table xip_chiphmin and xip
    //! \param smpi Object of class SmileiMPI containing MPI properties
    void bcastTableXi( SmileiMPI *smpi );

private:

    // ---------------------------------------------
    // General parameters
    // ---------------------------------------------

    //! Output format of the tables
    std::string output_format_;

    //! Path to the tables
    std::string table_path_;

    //! Flag that activate the table computation
    bool compute_table_;

    //! Minimum threshold above which the Monte-Carlo algorithm is working
    //! This avoids using the Monte-Carlo algorithm when particle_chi is too low
    double minimum_chi_discontinuous_;

    //! Under this value, no radiation loss
    double minimum_chi_continuous_;

    // ---------------------------------------------
    // Table h for the
    // stochastic diffusive operator of Niel et al.
    // ---------------------------------------------
    
    struct Niel {
        
        //! Array containing tabulated values of the function h for the
        //! stochastic diffusive operator of Niel et al.
        std::vector<double > table;
        
        //! Minimum boundary of the table h
        double chipa_min;
        
        //! Maximum boundary of the table h
        double chipa_max;
        
        //! Inverse delta chi for the table h
        double chipa_inv_delta;
        
        //! Delta chi for the table h
        double chipa_delta;
        
        //! Log10 of the minimum boundary of the table h
        double log10_chipa_min;
        
        //! Method to be used to get the h values (table, fit5, fit10)
        std::string computation_method;
        
        //! Dimension of the array h
        int size_particle_chi;
        
    };
    
    struct Niel niel_;

    // ---------------------------------------------
    // Table integfochi
    // ---------------------------------------------

    struct IntegrationFoverChi {
        
        //! Array containing tabulated values for the computation
        //! of the photon production rate dN_{\gamma}/dt
        //! (which is also the optical depth for the Monte-Carlo process).
        //! This table is the integration of the Synchrotron emissivity
        //! refers to as F over the quantum parameter Chi.
        std::vector<double > table;
        
        //! Minimum boundary of the table integfochi_table
        double chipa_min;
        
        //! Maximum boundary of the table integfochi_table
        double chipa_max;
        
        //! Minimum boundary of the table integfochi_table
        double chipa_dim;
        
        //! Log10 of the minimum boundary of the table integfochi_table
        double log10_chipa_min;

        //! Delta chi for the table integfochi_table
        double chipa_delta;

        //! Inverse delta chi for the table integfochi_table
        double chipa_inv_delta;
        
    };
    
    struct IntegrationFoverChi integfochi;

    // ---------------------------------------------
    // Structure for min_photon_chi_for_xi and xi
    // ---------------------------------------------

    struct Xi {
        
        //! Table containing the cumulative distribution function \f$P(0 \rightarrow \chi_{\gamma})\f$
        //! that gives gives the probability for a photon emission in the range \f$[0, \chi_{\gamma}]\f$
        std::vector<double> table ;
        
        //! Table containing the photon_chi min values
        //! Under this value, photon energy is
        //! considered negligible
        std::vector<double > min_photon_chi_table;
        
        //! Logarithm of the minimum boundary for particle_chi in the table xip
        //! and xip_chiphmin
        double log10_chipa_min;
        
        //! Maximum boundary for particle_chi in the table xip and xip_chiphmin
        double chipa_max;
        
        //! Minimum boundary for particle_chi in the table xip and xip_chiphmin
        double chipa_min;
        
        //! Delta for the particle_chi discretization  in the table xip and xip_chiphmin
        double chipa_delta;
        
        //! Inverse of the delta for the particle_chi discretization
        //! in the table xip and xip_chiphmin
        double chipa_inv_delta;
        
        //! Dimension of the discretized parameter particle_chi
        int chipa_dim;
        
        //! Dimension of the discretized parameter photon_chi
        int chiph_dim;
        
        //! 1/(xi_.chiph_dim - 1)
        double inv_chiph_dim_minus_one;

        //! xip power
        double power;

        //! xip threshold
        double threshold;
        
    };
    
    struct Xi xi_;

    // ---------------------------------------------
    // Factors
    // ---------------------------------------------

    //! Factor for the computation of dNphdt
    double factor_dNph_dt_;

    //! Factor for the Classical radiated power
    //! 2.*params.fine_struct_cst/(3.*normalized_Compton_wavelength_);
    double factor_classical_radiated_power_;

    //! Normalized reduced Compton wavelength
    double normalized_Compton_wavelength_;

};

#endif
