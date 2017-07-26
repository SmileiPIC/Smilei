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

        //! Computation of the value T(chiph) using the approximated
        //! formula of Erber
        //! \param chiph photon quantum parameter
        //! \param nbit number of iteration for the Bessel evaluation
        //! \param eps epsilon for the Bessel evaluation
        double computeErberT(double chiph,int nbit,
                           double eps);

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
