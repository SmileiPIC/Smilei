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
#include "RadiationTools.h"
#include "H5.h"
#include "Random.h"

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
    double computePhotonProductionYield( double particle_chi, double particle_gamma );

    //! Determine randomly a photon quantum parameter photon_chi
    //! for an emission process
    //! from a particle chi value (particle_chi) and
    //! using the tables xip and chiphmin
    //! \param particle_chi particle quantum parameter
    // double computeRandomPhotonChi( double particle_chi );

    //! Computation of the photon quantum parameter photon_chi for emission
    //! ramdomly and using the tables xi and chiphmin
    //! \param particle_chi particle quantum parameter
    double computeRandomPhotonChiWithInterpolation( double particle_chi, Random * rand);

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
    
    #pragma acc routine seq
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
    inline double getRidgersCorrectedRadiatedEnergy( double particle_chi,
            double dt )
    {
        return RadiationTools::computeRidgersFit( particle_chi )*dt*particle_chi*particle_chi*factor_classical_radiated_power_;
    };

    //! Get of the classical continuous radiated energy during dt
    //! \param particle_chi particle quantum parameter
    //! \param dt time step
    inline double getClassicalRadiatedEnergy( double particle_chi, double dt )
    {
        return dt*particle_chi*particle_chi*factor_classical_radiated_power_;
    };

    //! Return the minimum_chi_discontinuous_ value
    //! Under this value, no discontinuous radiation reaction
    inline double getMinimumChiDiscontinuous()
    {
        return minimum_chi_discontinuous_;
    }

    //! Return the minimum_chi_continuous_ value
    //! Under this value, no continuous radiation reaction
    inline double getMinimumChiContinuous()
    {
        return minimum_chi_continuous_;
    }

    inline std::string getNielHComputationMethod()
    {
        return this->niel_.computation_method_;
    }

    inline int getNielHComputationMethodIndex()
    {
        return this->niel_.computation_method_index_;
    }

    // -----------------------------------------------------------------------------
    //! Return the classical power factor factor_classical_radiated_power_.
    // -----------------------------------------------------------------------------
    inline double getFactorClassicalRadiatedPower()
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

    //! Bcast of the external table h
    //! \param smpi Object of class SmileiMPI containing MPI properties
    void bcastHTable( SmileiMPI *smpi );

    //! Bcast of the external table integfochi_
    //! \param smpi Object of class SmileiMPI containing MPI properties
    void bcastIntegfochiTable( SmileiMPI *smpi );

    //! Bcast of the external table xip_chiphmin and xip
    //! \param smpi Object of class SmileiMPI containing MPI properties
    void bcastTableXi( SmileiMPI *smpi );


    void needNielTables();

    // ---------------------------------------------
    // Table h for the
    // stochastic diffusive operator of Niel et al.
    // ---------------------------------------------

    struct Niel {

        //! Array containing tabulated values of the function h for the
        //! stochastic diffusive operator of Niel et al.
        std::vector<double > table_;

        //! Minimum boundary of the table h
        double min_particle_chi_;

        //! Maximum boundary of the table h
        double max_particle_chi_;

        //! Inverse delta chi for the table h
        double inv_particle_chi_delta_;

        //! Delta chi for the table h
        double particle_chi_delta_;

        //! Log10 of the minimum boundary of the table h
        double log10_min_particle_chi_;

        //! Method to be used to get the h values (table, fit5, fit10)
        std::string computation_method_;

        //! Index for the computational method
        int computation_method_index_;

        //! Dimension of the array h
        int size_particle_chi_;

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
        std::vector<double > table_;

        //! Minimum boundary of the table integfochi_table
        double min_particle_chi_;

        //! Maximum boundary of the table integfochi_table
        double max_particle_chi_;

        //! Minimum boundary of the table integfochi_table
        int size_particle_chi_;

        //! Log10 of the minimum boundary of the table integfochi_table
        double log10_min_particle_chi_;

        //! Delta chi for the table integfochi_table
        double particle_chi_delta_;

        //! Inverse delta chi for the table integfochi_table
        double inv_particle_chi_delta_;

    };

    struct IntegrationFoverChi integfochi_;

    // ---------------------------------------------
    // Structure for min_photon_chi_for_xi and xi
    // ---------------------------------------------

    struct Xi {

        //! Table containing the cumulative distribution function \f$P(0 \rightarrow \chi_{\gamma})\f$
        //! that gives gives the probability for a photon emission in the range \f$[0, \chi_{\gamma}]\f$
        std::vector<double> table_ ;

        //! Table containing the photon_chi min values
        //! Under this value, photon energy is
        //! considered negligible
        std::vector<double > min_photon_chi_table_;

        //! Logarithm of the minimum boundary for particle_chi in the table xip
        //! and xip_chiphmin
        double log10_min_particle_chi_;

        //! Maximum boundary for particle_chi in the table xip and xip_chiphmin
        double max_particle_chi_;

        //! Minimum boundary for particle_chi in the table xip and xip_chiphmin
        double min_particle_chi_;

        //! Delta for the particle_chi discretization  in the table xip and xip_chiphmin
        double particle_chi_delta_;

        //! Inverse of the delta for the particle_chi discretization
        //! in the table xip and xip_chiphmin
        double inv_particle_chi_delta_;

        //! Dimension of the discretized parameter particle_chi
        int size_particle_chi_;

        //! Dimension of the discretized parameter photon_chi
        int size_photon_chi_;

        //! 1/(xi_.size_photon_chi_ - 1)
        double inv_size_photon_chi_minus_one_;

        //! xip power
        // double power_;

        //! xip threshold
        // double threshold_;

    };

    struct Xi xi_;
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
