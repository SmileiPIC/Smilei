// ----------------------------------------------------------------------------
//! \file MultiphotonBreitWheelerTools.cpp
//
//! \brief This class contains tools
//! for the Multiphoton Breit-Wheeler pair generation
//
// ----------------------------------------------------------------------------

#include "MultiphotonBreitWheelerTools.h"


// -----------------------------------------------------------------------------
// PHYSICAL COMPUTATION
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//! Computation of the production rate of pairs per photon
//! \param photon_chi photon quantum parameter
//! \param gamma photon normalized energy
// -----------------------------------------------------------------------------
double MultiphotonBreitWheelerTools::computeBreitWheelerPairProductionRate( 
    double photon_chi,
    double gamma,
    MultiphotonBreitWheelerTables * mBW_tables )
{
    // ________________________________________
    // Parameters

    double logchiphm;
    double logchiphp;
    // Index
    int ichiph;
    // final value
    double dNBWdt;

    // ________________________________________
    // Computation

    // Log of the photon quantum parameter particle_chi
    double logchiph = std::log10( photon_chi );

    // Lower index for interpolation in the table integfochi
    ichiph = int( std::floor( ( logchiph-mBW_tables->T_.log10_min_photon_chi_ )
                         *mBW_tables->T_.photon_chi_inv_delta_ ) );

    // If photon_chi is below the lower bound of the table
    // An asymptotic approximation is used
    if( ichiph < 0 ) {
        ichiph = 0;
        // 0.2296 * sqrt(3) * pi [MG/correction by Antony]
        dNBWdt = 1.2493450020845291*std::exp( -8.0/( 3.0*photon_chi ) ) * photon_chi*photon_chi;
    }
    // If photon_chi is above the upper bound of the table
    // An asymptotic approximation is used
    else if( ichiph >= mBW_tables->T_.size_photon_chi_-1 ) {
        ichiph = mBW_tables->T_.size_photon_chi_-2;
        dNBWdt = 2.067731275227008*pow( photon_chi, 5.0/3.0 );
    } else {
        // Upper and lower values for linear interpolation
        logchiphm = ichiph*mBW_tables->T_.photon_chi_delta_ + mBW_tables->T_.log10_min_photon_chi_;
        logchiphp = logchiphm + mBW_tables->T_.photon_chi_delta_;

        // Interpolation
        dNBWdt = ( mBW_tables->T_.table_[ichiph+1]*fabs( logchiph-logchiphm ) +
                   mBW_tables->T_.table_[ichiph]*fabs( logchiphp - logchiph ) )*mBW_tables->T_.photon_chi_inv_delta_;
    }
    return mBW_tables->getFactorDNdWdt()*dNBWdt/(photon_chi*gamma);
}