// ----------------------------------------------------------------------------
//! \file MultiphotonBreitWheeler.h
//
//! \brief This class contains the methods and tools to generate and manage
//! the physical tables for the multiphoton Breit-wheeler process.
//
//! \details
//! The implementation is adapted from the following results:
//! - Niel et al.
//! - M. Lobet (http://www.theses.fr/2015BORD0361)
// ----------------------------------------------------------------------------

#ifndef MULTIPHOTON_BREIT_WHEELER_H
#define MULTIPHOTON_BREIT_WHEELER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include "Tools.h"
#include <mpi.h>
#include "H5.h"

class MultiphotonBreitWheeler
{

    public:
        
    //! Creation of the tables
    static void createTables(int argc, std::string * arguments);
        
    // -----------------------------------------------------------------------------
    //! Computation of the value of the integration of  dT(photon_chi)/dhicpa
    //! using the formula of Ritus.
    //! \param photon_chi photon quantum parameter
    //! \param particle_chi particle quantum parameter for integration (=0.5*photon_chi for full integration)
    //! \param nb_iterations number of iteration for the Gauss-Legendre integration
    //! \param eps epsilon for the Bessel evaluation
    // -----------------------------------------------------------------------------
    static double computeIntegrationRitusDerivative( double photon_chi,
            double particle_chi,
            int nb_iterations,
            double eps );
        
    // -----------------------------------------------------------------------------
    //! Computation of the value dT/dchipa(photon_chi) using the formula of Ritus
    //! \param photon_chi photon quantum parameter
    //! \param particle_chi particle quantum parameter
    //! \param nb_iterations number of iteration for the Gauss-Legendre integration
    //! \param eps epsilon for the Bessel evaluation
    // -----------------------------------------------------------------------------
    static double computeRitusDerivative( double photon_chi,
            double particle_chi, int nb_iterations, double eps );
    
    // -----------------------------------------------------------------------------
    //! Computation of the value T(photon_chi) using the approximated
    //! formula of Erber
    //! \param photon_chi photon quantum parameter
    //! \param nb_iterations number of iteration for the Bessel evaluation
    //! \param eps epsilon for the Bessel evaluation
    // -----------------------------------------------------------------------------
    static double computeErberT( double photon_chi);
            
};

#endif
