// ----------------------------------------------------------------------------
//! \file NonlinearInverseComptonScattering.h
//
//! \brief This class contains methods to generate the nonlinear inverse Compton scattering tables
//
// ----------------------------------------------------------------------------

#ifndef COMPTON_SCATTERING_H
#define COMPTON_SCATTERING_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include "Tools.h"
#include <mpi.h>
#include "H5.h"

class NonlinearComptonScattering
{

    public:
        
    //! Creation of the tables
    static void createTables(int argc, std::string * arguments);
        
    //! Compute the integration of the synchrotron emissivity S/chi
    //! refered to as K in the documentation
    //! between min_photon_chi and max_photon_chi
    //! using Gauss-Legendre for a given particle_chi value
    //! \param nb_iterations number of iteration for the Gauss-Legendre
    //! \param eps relative error on the integration
    static double integrateSynchrotronEmissivity( double particle_chi,
            double min_photon_chi,
            double max_photon_chi,
            int nb_iterations,
            double eps );

    //! Synchrotron emissivity from Ritus
    //! \param particle_chi particle quantum parameter
    //! \param photon_chi photon quantum parameter
    //! \param nb_iterations number of iterations for the Gauss-Legendre integration
    //! \param eps epsilon for the modified bessel function
    static double computeRitusSynchrotronEmissivity( double particle_chi,
            double photon_chi,
            int nb_iterations,
            double eps );

    //! Return the value of the function h(particle_chi) of Niel et al.
    //! Use an integration of Gauss-Legendre
    //
    //! \param particle_chi particle quantum parameter
    //! \param nb_iterations number of iterations for the Gauss-Legendre integration
    //! \param eps epsilon for the modified bessel function
    static double computeHNiel( double particle_chi,
            int nb_iterations, double eps );

};

#endif
