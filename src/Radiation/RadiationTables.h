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
#include "Random.h"
#include "Table.h"
#include "Table2D.h"
#include "RadiationTools.h"

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
    void initialization( Params &params , SmileiMPI *smpi );

    // ---------------------------------------------------------------------
    // PHYSICAL COMPUTATION
    // ---------------------------------------------------------------------

    //! Computation of the photon production yield dNph/dt which is
    //! also the cross-section for the Monte-Carlo
    //! param[in] particle_chi particle quantum parameter
    //! param[in] particle_gamma particle Lorentz factor
    //! param[in] integfochi_table table of the discretized integrated f/chi function for Photon production yield computation
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc routine seq
#endif
    double computePhotonProductionYield( const double particle_chi,
                                         const double particle_gamma);

    //! Determine randomly a photon quantum parameter photon_chi
    //! for an emission process
    //! from a particle chi value (particle_chi) and
    //! using the tables xip and chiphmin
    //! \param particle_chi particle quantum parameter
    // double computeRandomPhotonChi( double particle_chi );

    //! Computation of the photon quantum parameter photon_chi for emission
    //! ramdomly and using the tables xi and chiphmin
    //! \param[in] particle_chi particle quantum parameter
    //! \param[in] xi
    //! \param[in] table_min_photon_chi
    //! \param[in] table_xi
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc routine seq
#endif
    double computeRandomPhotonChiWithInterpolation( double particle_chi,
                                                    double xi);

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
    
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc routine seq
#endif
    double getHNielFromTable( double particle_chi, double * tableNiel);

    //! Return the stochastic diffusive component of the pusher
    //! of Niel et al.
    //! \param gamma particle Lorentz factor
    //! \param particle_chi particle quantum parameter
    //! \param dt time step
    // double getNielStochasticTerm( double gamma,
    //                               double particle_chi,
    //                               double dt,
    //                               Random * rand);

    //! Computation of the corrected continuous quantum radiated energy
    //! during dt from the quantum parameter particle_chi using the Ridgers
    //! formulae.
    //! \param particle_chi particle quantum parameter
    //! \param dt time step
    //#pragma omp declare simd
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc routine seq
#endif
    inline double __attribute__((always_inline)) getRidgersCorrectedRadiatedEnergy( const double particle_chi,
            const double dt )
    {
        return computeRidgersFit( particle_chi )*dt*particle_chi*particle_chi*factor_classical_radiated_power_;
    };

    //! Computation of the function g of Erber using the Ridgers
    //! approximation formulae
    //! \param particle_chi particle quantum parameter
    //#pragma omp declare simd
    static inline double __attribute__((always_inline)) computeRidgersFit( double particle_chi )
    {
        double a = 1.0 + 4.8 * ( 1.0 + particle_chi )*std::log( 1.0 + 1.7 * particle_chi ) 
                    + 2.44 * particle_chi * particle_chi;
        return 1.0 / std::cbrt( a * a );
    };

    //! Get of the classical continuous radiated energy during dt
    //! \param particle_chi particle quantum parameter
    //! \param dt time step
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc routine seq
#endif
    inline double __attribute__((always_inline)) getClassicalRadiatedEnergy( double particle_chi, double dt )
    {
        return dt*particle_chi*particle_chi*factor_classical_radiated_power_;
    };

    //! Return the minimum_chi_discontinuous_ value
    //! Under this value, no discontinuous radiation reaction
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc routine seq
#endif
    inline double __attribute__((always_inline)) getMinimumChiDiscontinuous()
    {
        return minimum_chi_discontinuous_;
    }

    //! Return the minimum_chi_continuous_ value
    //! Under this value, no continuous radiation reaction
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc routine seq
#endif
    inline double __attribute__((always_inline)) getMinimumChiContinuous()
    {
        return minimum_chi_continuous_;
    }

    inline std::string __attribute__((always_inline)) getNielHComputationMethod()
    {
        return this->niel_computation_method_;
    }

    inline int __attribute__((always_inline)) getNielHComputationMethodIndex()
    {
        return this->niel_computation_method_index_;
    }

    // -----------------------------------------------------------------------------
    //! Return the classical power factor factor_classical_radiated_power_.
    // -----------------------------------------------------------------------------
    inline double __attribute__((always_inline)) getFactorClassicalRadiatedPower()
    {
        return factor_classical_radiated_power_;
    }

    // ---------------------------------------------------------------------
    // TABLE READING
    // ---------------------------------------------------------------------

    //! Read the external table h
    //! \param smpi Object of class SmileiMPI containing MPI properties
    void readHTable( SmileiMPI *smpi );

    //! Read the external table integfochi_
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

    //! Bcast of the external table xip_chiphmin and xip
    //! \param smpi Object of class SmileiMPI containing MPI properties
    void bcastTableXi( SmileiMPI *smpi );


    void needNielTables();

    // ---------------------------------------------
    // Table h for the
    // stochastic diffusive operator of Niel et al.
    // ---------------------------------------------
    
    // 1d array
    // axe 0 = particle_chi
    Table niel_;

    // ---------------------------------------------
    // Table integfochi
    // ---------------------------------------------

    // 1d array
    // axe 0 = particle_chi
    Table integfochi_;

    // ---------------------------------------------
    // Structure for min_photon_chi_for_xi and xi
    // ---------------------------------------------

    // 2d array
    // axe0: particle_chi
    // axe1: photon_chi
    Table2D xi_;
    
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

    //! Method to be used to get the h values (table, fit5, fit10)
    std::string niel_computation_method_;

    //! Index for the computational method
    int niel_computation_method_index_;

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
