// ----------------------------------------------------------------------------
//! \file MultiphotonBreitWheelerTools.h
//
//! \brief This class contains tools
//! for the Multiphoton Breit-Wheeler pair generation
//
// ----------------------------------------------------------------------------

#ifndef MULTIPHOTONBREITWHEELERTOOLS_H
#define MULTIPHOTONBREITWHEELERTOOLS_H

#include <cmath>

#include "MultiphotonBreitWheelerTables.h"

class MultiphotonBreitWheelerTools {

    public:

        // -----------------------------------------------------------------------------
        //! Computation of the production rate of pairs per photon
        //! \param photon_chi photon quantum parameter
        //! \param gamma photon normalized energy
        // -----------------------------------------------------------------------------
        static double computeBreitWheelerPairProductionRate( 
            double photon_chi, 
            double gamma,
            MultiphotonBreitWheelerTables * mBW_tables);

    private:


};
#endif
