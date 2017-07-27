// ----------------------------------------------------------------------------
//! \file MultiphotonBreitWheelerTables.h
//
//! \brief This class contains the mathods and tools to generate and manage
//! the physical tables for the multiphoton Breit-wheeler process.
//
//! \details This header contains the definition of the class
//! MultiphotonBreitWheelerTables.
//! The implementation is adapted from the thesis results of M. Lobet
//! See http://www.theses.fr/2015BORD0361
// ----------------------------------------------------------------------------

#ifndef MBWTABLES_H
#define MBWTABLES_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "Params.h"
#include "H5.h"
#include "userFunctions.h"

//------------------------------------------------------------------------------
//! MutliphotonBreitWheelerTables class: holds parameters, tables and
//! functions to compute cross-sections,
//! optical depths and other useful parameters for the pair creation Monte-Carlo
//! process.
//------------------------------------------------------------------------------
class MultiphotonBreitWheelerTables
{

    public:

        //! Constructor for MutliphotonBreitWheeler
        MultiphotonBreitWheelerTables();

        //! Destructor for MutliphotonBreitWheeler
        ~MultiphotonBreitWheelerTables();

        //! Initialization of the parmeters
        //! \param params Object Params for the parameters from the input script
        void initialization(Params& params);

        // ---------------------------------------------------------------------
        // PHYSICAL COMPUTATION
        // ---------------------------------------------------------------------

        //! Computation of the production rate of pairs per photon
        //! \param chiph photon quantum parameter
        //! \param gamma photon normalized energy
        double compute_dNBWdt(double chiph, double gamma);

        //! Computation of the value T(chiph) using the approximated
        //! formula of Erber
        //! \param chiph photon quantum parameter
        //! \param nbit number of iteration for the Bessel evaluation
        //! \param eps epsilon for the Bessel evaluation
        double compute_Erber_T(double chiph,int nbit,
                           double eps);

        //! Computation of the value T(chiph) using the formula of Ritus
        //! \param chiph photon quantum parameter
        //! \param nbit number of iteration for the Gauss-Legendre integration
        //! \param eps epsilon for the Bessel evaluation
        double compute_Ritus_T(double chiph,
                               int nbit,double eps);

       //! Computation of the value T(chiph) using the formula of Ritus
       //! \param chiph photon quantum parameter
       //! \param nbit number of iteration for the Gauss-Legendre integration
       //! \param eps epsilon for the Bessel evaluation
       double compute_Ritus_dTdchi(double chiph,
                             double chipa,int nbit,double eps);

        // ---------------------------------------------------------------------
        // TABLE COMPUTATION
        // ---------------------------------------------------------------------

        //! Computation of the table T_table that is a discetization of the
        //! T function for the multiphoton Breit-Wheeler process
        //! \param smpi Object of class SmileiMPI containing MPI properties
        void compute_T_table(SmileiMPI *smpi);

        //! Output the computed tables so that thay can be read at the next run.
        //! \param params list of simulation parameters
        //! \param smpi MPI parameters
        void compute_tables(Params& params,
                            SmileiMPI *smpi);

        // ---------------------------------------------------------------------
        // TABLE OUTPUTS
        // ---------------------------------------------------------------------

        //! Ouput in a file of the table values of T for the
        //! mutliphoton Breit-Wheeler process
        void output_T_table();

        //! Output the computed tables so that thay can be read at the next run.
        //! Table output by the master MPI rank
        //! \param smpi Object of class SmileiMPI containing MPI properties
        void output_tables(SmileiMPI *smpi);

        // ---------------------------------------------------------------------
        // TABLE READING
        // ---------------------------------------------------------------------

        //! Read the external table T
        //! \param smpi Object of class SmileiMPI containing MPI properties
        bool read_T_table(SmileiMPI *smpi);

        // ---------------------------------------------------------------------
        // TABLE COMMUNICATIONS
        // ---------------------------------------------------------------------

        //! Bcast of the external table T
        //! \param smpi Object of class SmileiMPI containing MPI properties
        void bcast_T_table(SmileiMPI *smpi);

    private:

        // ---------------------------------------------
        // General parameters
        // ---------------------------------------------

        //! Output format of the tables
        std::string output_format;

        //! Path to the tables
        std::string table_path;

        // ---------------------------------------------
        // Table T for the
        // pair creation Monte-Carlo process
        // ---------------------------------------------

        //! Array containing tabulated values of the function T
        std::vector<double > T_table;

        //! Minimum boundary of the table T
        double T_chiph_min;

        //! Log10 of the minimum boundary of the table T
        double T_log10_chiph_min;

        //! Maximum boundary of the table T
        double T_chiph_max;

        //! Delta chi for the table T
        double T_chiph_delta;

        //! Inverse delta chi for the table h
        double T_chiph_inv_delta;

        //! Dimension of the array T
        int T_dim;

        //! This variable is true if the table is computed, false if read
        bool T_computed;

        // ---------------------------------------------
        // Factors
        // ---------------------------------------------

        //! Factor for the computation of dN_{BW} / dt from T
        double factor_dNBWdt;

        //! Normalized reduced Compton wavelength
        double norm_lambda_compton;

};

#endif
