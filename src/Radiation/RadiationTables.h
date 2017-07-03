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
        void initParams(Params& params);

        // ---------------------------------------------------------------------
        // PHYSICAL COMPUTATION
        // ---------------------------------------------------------------------

        //! Synchrotron emissivity from Ritus
        static double compute_sync_emissivity_ritus(double chie,
                double chiph,
                int nbit,
                double eps);

        //! Computation of the cross-section dNph/dt
        double compute_dNphdt(double chipa,double gfpa);

        //! Compute integration of F/chi between
        //! using Gauss-Legendre for a given chie value
        static double compute_integfochi(double chie,
                double chipmin,
                double chipmax,
                int nbit,
                double eps);

        //! Computation of the photon quantum parameter chiph for emission
        //! ramdomly and using the tables xip and chiphmin
        //! \param chipa particle quantum parameter
        double compute_chiph_emission(double chipa);

        // -----------------------------------------------------------------------------
        //! Return the value of the function h(chipa) of Niel et al.
        //! Use an integration of Gauss-Legendre
        //
        //! \param chipa particle quantum parameter
        //! \param nbit number of iterations for the Gauss-Legendre integration
        //! \param eps epsilon for the modified bessel function
        // -----------------------------------------------------------------------------
        double get_h_Niel(double chipa,int nbit, double eps);

        //! Computation of the corrected continuous quantum radiated energy
        //! during dt from the quantum parameter chipa using the Ridgers
        //! formulae.
        //! \param chipa particle quantum parameter
        //! \param dt time step
        //#pragma omp declare simd
        double inline get_corrected_cont_rad_energy_Ridgers(double chipa,
                                                            double dt)
        {
            return compute_g_Ridgers(chipa)*dt*chipa*chipa*factor_cla_rad_power;
        };

        //! Get of the classical continuous radiated energy during dt
        //! \param chipa particle quantum parameter
        //! \param dt time step
        double inline get_classical_cont_rad_energy(double chipa, double dt)
        {
            return dt*chipa*chipa*factor_cla_rad_power;
        };

        //! Return the get_chipa_disc_min_threshold value
        double inline get_chipa_disc_min_threshold()
        {
            return chipa_disc_min_threshold;
        }

        //! Computation of the function g of Erber using the Ridgers
        //! approximation formulae
        //! \param chipa particle quantum parameter
        //#pragma omp declare simd
        double inline compute_g_Ridgers(double chipa)
        {
            return pow(1. + 4.8*(1.+chipa)*log(1. + 1.7*chipa)
                          + 2.44*chipa*chipa,-2./3.);
        };

        // ---------------------------------------------------------------------
        // TABLE COMPUTATION
        // ---------------------------------------------------------------------

        //! Computation of the table h that is a discetization of the h function
        //! in the stochastic model of Niel.
        //! \param smpi Object of class SmileiMPI containing MPI properties
        void compute_h_table(SmileiMPI *smpi);

        //! Generate table values for Integration of F/chi: Integfochi
        //! \param smpi Object of class SmileiMPI containing MPI properties
        void compute_integfochi_table(SmileiMPI *smpi);

        //! Computation of the minimum photon quantum parameter for the array
        //! xip (xip_chiphmin) and computation of the xip array.
        //! \param smpi Object of class SmileiMPI containing MPI properties
        void compute_xip_table(SmileiMPI *smpi);

        //! Compute all the tables
        void compute_tables(Params& params, SmileiMPI *smpi);

        // ---------------------------------------------------------------------
        // TABLE OUTPUTS
        // ---------------------------------------------------------------------

        //! Write in a file table values of the h table
        void output_h_table();

        //! Write in a file table values for Integration of F/chi: Integfochi
        void output_integfochi_table();

        //! Write in a file the table xip_chiphmin and xip
        void output_xip_table();

        //! Output all computed tables so that they can be
        //! read at the next run
        //! Table output by the master MPI rank
        //! \param smpi Object of class SmileiMPI containing MPI properties
        void output_tables(SmileiMPI *smpi);

        // ---------------------------------------------------------------------
        // TABLE READING
        // ---------------------------------------------------------------------

        //! Read the external table h
        //! \param smpi Object of class SmileiMPI containing MPI properties
        bool read_h_table(SmileiMPI *smpi);

        //! Read the external table integfochi
        //! \param smpi Object of class SmileiMPI containing MPI properties
        bool read_integfochi_table(SmileiMPI *smpi);

        //! Read the external table xip_chiphmin and xip
        //! \param smpi Object of class SmileiMPI containing MPI properties
        bool read_xip_table(SmileiMPI *smpi);

        // ---------------------------------------------------------------------
        // TABLE COMMUNICATIONS
        // ---------------------------------------------------------------------

        //! Bcast of the external table h
        //! \param smpi Object of class SmileiMPI containing MPI properties
        void bcast_h_table(SmileiMPI *smpi);

        //! Bcast of the external table integfochi
        //! \param smpi Object of class SmileiMPI containing MPI properties
        void bcast_integfochi_table(SmileiMPI *smpi);

        //! Bcast of the external table xip_chiphmin and xip
        //! \param smpi Object of class SmileiMPI containing MPI properties
        void bcast_xip_table(SmileiMPI *smpi);

    private:

        // ---------------------------------------------
        // General parameters
        // ---------------------------------------------

        //! Output format of the tables
        std::string output_format;

        //! Path to the tables
        std::string table_path;

        //! Minimum threashold above which the Monte-Carlo algorithm is working
        //! This avoids using the Monte-Carlo algorithm when chipa is too low
        double chipa_disc_min_threshold;

        // ---------------------------------------------
        // Table h for the
        // stochastic diffusive operator of Niel et al.
        // ---------------------------------------------

        //! Array containing tabulated values of the function h for the
        //! stochastic diffusive operator of Niel et al.
        std::vector<double > h_table;

        //! Minimum boundary of the table h
        double chipa_h_min;

        //! Log10 of the minimum boundary of the table h
        double log10_chipa_h_min;

        //! Maximum boundary of the table h
        double chipa_h_max;

        //! Delta chi for the table h
        double chipa_h_delta;

        //! Dimension of the array h
        int chipa_h_dim;

        //! Inverse delta chi for the table h
        double chipa_h_inv_delta;

        //! This variable is true if the table is computed, false if read
        bool h_computed;

        // ---------------------------------------------
        // Table Integfochi
        // ---------------------------------------------

        //! Array containing tabulated values for the computation
        //! of the photon production rate dN_{\gamma}/dt
        //! (which is also the optical depth for the Monte-Carlo process).
        //! This table is the integration of the Synchrotron emissivity
        //! refers to as F over the quantum parameter Chi.
        std::vector<double > Integfochi;

        //! Minimum boundary of the table Integfochi
        double chipa_integfochi_min;

        //! Log10 of the minimum boundary of the table Integfochi
        double log10_chipa_integfochi_min;

        //! Maximum boundary of the table Integfochi
        double chipa_integfochi_max;

        //! Delta chi for the table Integfochi
        double delta_chipa_integfochi;

        //! Inverse delta chi for the table Integfochi
        double inv_delta_chipa_integfochi;

        //! Dimension of the array Integfochi
        int dim_integfochi;

        //! This variable is true if the table is computed, false if read
        bool integfochi_computed;

        // ---------------------------------------------
        // Table chiph min for xip table
        // ---------------------------------------------

        //! Table containing the chiph min values
        //! Under this value, photon energy is
        //! considered negligible
        std::vector<double > xip_chiphmin_table;

        // ---------------------------------------------
        // Table xip
        // ---------------------------------------------

        //! Table containing the cumulative distribution function \f$P(0 \rightarrow \chi_{\gamma})\f$
        //! that gives gives the probability for a photon emission in the range \f$[0, \chi_{\gamma}]\f$
        std::vector<double> xip_table;

        //! Minimum boundary for chipa in the table xip and xip_chiphmin
        double chipa_xip_min;

        //! Logarithm of the minimum boundary for chipa in the table xip
        //! and xip_chiphmin
        double log10_chipa_xip_min;

        //! Maximum boundary for chipa in the table xip and xip_chiphmin
        double chipa_xip_max;

        //! Delta for the chipa discretization  in the table xip and xip_chiphmin
        double chipa_xip_delta;

        //! Inverse of the delta for the chipa discretization
        //! in the table xip and xip_chiphmin
        double inv_chipa_xip_delta;

        //! Dimension of the discretized parameter chipa
        int chipa_xip_dim;

        //! Dimension of the discretized parameter chiph
        int chiph_xip_dim;

        //! 1/(chiph_xip_dim - 1)
        double inv_chiph_xip_dim_minus_one;

        //! xip power
        double xip_power;

        //! xip threshold
        double xip_threshold;

        //! This variable is true if the table is computed, false if read
        bool xip_computed;

        // ---------------------------------------------
        // Factors
        // ---------------------------------------------

        //! Factor for the computation of dNphdt
        double factor_dNphdt;

        //! Factor Classical radiated power
        double factor_cla_rad_power;

        //! Fine structure constant
        const double fine_struct_cst = 7.2973525698e-3;

        //! Reduced Planck Constant (J.s)
        const double red_planck_cst = 1.054571628E-34;

        //! Electron mass
        const double electron_mass = 9.109382616e-31;

        //! Speed of light in vacuum (m/s)
        const double c_vacuum = 299792458;

        //! Normalized reduced Compton wavelength
        double norm_lambda_compton;

};

#endif
