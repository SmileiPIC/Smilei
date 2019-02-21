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
        //! \param photon_chi photon quantum parameter
        //! \param gamma photon normalized energy
        double compute_dNBWdt(double photon_chi, double gamma);

        //! Computation of the value T(photon_chi) using the approximated
        //! formula of Erber
        //! \param photon_chi photon quantum parameter
        //! \param nb_iterations number of iteration for the Bessel evaluation
        //! \param eps epsilon for the Bessel evaluation
        double compute_Erber_T(double photon_chi,int nb_iterations,
                           double eps);

        //! Computation of the value T(photon_chi) using the formula of Ritus
        //! \param photon_chi photon quantum parameter
        //! \param particle_chi particle quantum parameter for integration (=0.5*photon_chi for full integration)
        //! \param nb_iterations number of iteration for the Gauss-Legendre integration
        //! \param eps epsilon for the Bessel evaluation
        double compute_integration_Ritus_dTdchi(double photon_chi,
                               double particle_chi,
                               int nb_iterations,
                               double eps);

       //! Computation of the value T(photon_chi) using the formula of Ritus
       //! \param photon_chi photon quantum parameter
       //! \param nb_iterations number of iteration for the Gauss-Legendre integration
       //! \param eps epsilon for the Bessel evaluation
       double compute_Ritus_dTdchi(double photon_chi,
                             double particle_chi,int nb_iterations,double eps);

        //! Computation of the electron and positron quantum parameters for
        //! the multiphoton Breit-Wheeler pair creation
        //! \param photon_chi photon quantum parameter
        double * compute_pair_chi(double photon_chi);

        // ---------------------------------------------------------------------
        // TABLE COMPUTATION
        // ---------------------------------------------------------------------

        //! Computation of the table T_table that is a discetization of the
        //! T function for the multiphoton Breit-Wheeler process
        //! \param smpi Object of class SmileiMPI containing MPI properties
        void compute_T_table(SmileiMPI *smpi);


        //! Computation of the minimum particle quantum parameter chipamin
        //! for the photon xip array and computation of the photon xip array.
        //! \details Under the minimum particle_chi value, the particle kinetic energy is
        //! considered negligible. All energy goes to the other.
        //! \param smpi Object of class SmileiMPI containing MPI properties
        void compute_xip_table(SmileiMPI *smpi);

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

        //! File output of xip_chipamin_table and xip_table
        void output_xip_table();

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

        //! Read the external table xip_chipamin and xip
        //! \param smpi Object of class SmileiMPI containing MPI properties
        bool read_xip_table(SmileiMPI *smpi);

        // ---------------------------------------------------------------------
        // TABLE COMMUNICATIONS
        // ---------------------------------------------------------------------

        //! Bcast of the external table T
        //! \param smpi Object of class SmileiMPI containing MPI properties
        void bcast_T_table(SmileiMPI *smpi);

        //! Bcast of the external table xip_chipamin and xip
        //! \param smpi Object of class SmileiMPI containing MPI properties
        void bcast_xip_table(SmileiMPI *smpi);

    private:

        // ---------------------------------------------
        // General parameters
        // ---------------------------------------------

        //! Output format of the tables
        std::string output_format_;

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
        // Table particle_chi min for xip table
        // ---------------------------------------------

        //! Table containing the particle_chi min values
        //! Under this value, electron kinetic energy of the pair is
        //! considered negligible
        std::vector<double > xip_chipamin_table;

        // ---------------------------------------------
        // Table xip
        // ---------------------------------------------

        //! Table containing the cumulative distribution function \f$P(0 \rightarrow \chi_{e^-})\f$
        //! that gives gives the probability for a photon to decay into pair
        //! with an electron of energy in the range \f$[0, \chi_{e^-}]\f$
        //! This enables to compute the energy repartition between the electron and the positron
        std::vector<double> xip_table;

        //! Minimum boundary for photon_chi in the table xip and xip_chipamin
        double xip_chiph_min;

        //! Logarithm of the minimum boundary for photon_chi in the table xip
        //! and xip_chipamin
        double xip_log10_chiph_min;

        //! Maximum boundary for photon_chi in the table xip and xip_chipamin
        double xip_chiph_max;

        //! Delta for the photon_chi discretization in the table xip and xip_chipamin
        double xip_chiph_delta;

        //! Inverse of the delta for the photon_chi discretization
        //! in the table xip and xip_chipamin
        double xip_chiph_inv_delta;

        //! Dimension of the discretized parameter photon_chi
        int xip_chiph_dim;

        //! Dimension of the discretized parameter particle_chi
        int xip_chipa_dim;

        //! 1/(xip_chipa_dim - 1)
        double xip_inv_chipa_dim_minus_one;

        //! xip power
        double xip_power;

        //! xip threshold
        double xip_threshold;

        //! This variable is true if the table is computed, false if read
        bool xip_computed;

        // ---------------------------------------------
        // Factors
        // ---------------------------------------------

        //! Factor for the computation of dN_{BW} / dt from T
        double factor_dNBWdt;

        //! Normalized reduced Compton wavelength
        double normalized_Compton_wavelength_;

};

#endif
