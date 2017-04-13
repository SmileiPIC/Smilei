// ----------------------------------------------------------------------------
//! \file NLICompton.h
//
//! \brief Nonlinear Inverse Compton Scattering
//
//! \details This header contains the definition of the class NLICompton.
//
// ----------------------------------------------------------------------------

#ifndef NLICOMPTON_H
#define NLICOMPTON_H


#include <iostream>
#include <fstream>
#include <vector>

//----------------------------------------------------------------------------------------------------------------------
//! NLICompton class: holds parameters, tables and functions to compute cross-sections, optical depths and other usefull parameters for the Compton Monte-Carlo pusher.
//----------------------------------------------------------------------------------------------------------------------
class NLICompton
{

    public:

        //! Constructor for NLICompton
        NLICompton();

        //! Destructor for NLICompton
        ~NLICompton(){}

        //! Compute integration of F/chi between 
        //! using Gauss-Legendre for a given chie value
        static double compute_integfochi(double chie,
                double chipmin,
                double chipmax,
                int nbit,
                double eps);

        //! Synchrotron emissivity from Ritus
        static double compute_sync_emissivity_ritus(double chie,
                double chiph, 
                int nbit, 
                double eps);

        //! Generate table values for Integration of F/chi: Integfochi
        void compute_integfochi_table();

        //! Write in a file table values for Integration of F/chi: Integfochi
        void output_integfochi_table();
 
    private:


        //! Array containing tabulated values for the computation of the photon production rate dN_{\gamma}/dt (which is also the optical depth for the Monte-Carlo process). This table is the integration of the Synchrotron emissivity refers to as F over the quantum parameter Chi.
        std::vector<double > Integfochi;

        //! Minimum boundary of the table Integfochi
        double chie_integfochi_min;

        //! Minimum boundary of the table Integfochi
        double chie_integfochi_max; 

        //! Delta chi for the table integfochi
        double delta_chie_integfochi;

        //! Dimension of the array Integfochi
        unsigned int dim_integfochi;

        
};

#endif
