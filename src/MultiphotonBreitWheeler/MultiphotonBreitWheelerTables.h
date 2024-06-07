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
#include <iomanip>

#include "Params.h"
#include "userFunctions.h"
#include "Random.h"
#include "Table.h"
#include "Table2D.h"

//------------------------------------------------------------------------------
//! MultiphotonBreitWheelerTables class: holds parameters, tables and
//! functions to compute cross-sections,
//! optical depths and other useful parameters for the pair creation Monte-Carlo
//! process.
//------------------------------------------------------------------------------
class MultiphotonBreitWheelerTables
{

public:

    //! Constructor for MultiphotonBreitWheeler
    MultiphotonBreitWheelerTables();

    //! Destructor for MultiphotonBreitWheeler
    ~MultiphotonBreitWheelerTables();

    //! Initialization of the parmeters
    //! \param params Object Params for the parameters from the input script
    void initialization( Params &params, SmileiMPI *smpi );

    // ---------------------------------------------------------------------
    // PHYSICAL COMPUTATION
    // ---------------------------------------------------------------------

    //! Computation of the electron and positron quantum parameters for
    //! the multiphoton Breit-Wheeler pair creation
    //! \param photon_chi photon quantum parameter
    //! \param[out] pair_chi quantum parameters of the pair
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc routine seq
#endif
    void computePairQuantumParameter( const double photon_chi, 
                                      double * pair_chi,
                                      const double xip );

    //! Return factor for dN / dWdt computation
    inline double  __attribute__((always_inline)) getFactorDNdWdt(void) {
        return factor_dNBW_dt_;
    }

    // -----------------------------------------------------------------------------
    //! Computation of the production rate of pairs per photon
    //! \param photon_chi photon quantum parameter
    //! \param gamma photon normalized energy
    // -----------------------------------------------------------------------------
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc routine seq
#endif
    double computeBreitWheelerPairProductionRate( 
        const double photon_chi, 
        const double photon_gamma);

    // ---------------------------------------------------------------------
    // TABLE READING
    // ---------------------------------------------------------------------

    //! Read the external table T
    //! \param smpi Object of class SmileiMPI containing MPI properties
    void readTableT( SmileiMPI *smpi );

    //! Read the external table xip_chipamin and xip
    //! \param smpi Object of class SmileiMPI containing MPI properties
    void readTableXi( SmileiMPI *smpi );

    //! Read all external tables
    //! \param smpi Object of class SmileiMPI containing MPI properties
    void readTables( Params &params, SmileiMPI *smpi );

    // ---------------------------------------------
    // Structure for Table T used for the
    // pair creation Monte-Carlo process
    // ---------------------------------------------

    // 1d array
    // axe 0 : photon_chi
    Table T_;

    // ---------------------------------------------
    // Structure for xi and particle_chi min for xip table
    // ---------------------------------------------

    // 2d array:
    // - axe0: photon_chi
    // - axe1: particle_chi
    Table2D xi_;
    
private:

    // ---------------------------------------------
    // General parameters
    // ---------------------------------------------

    //! Path to the tables
    std::string table_path_;

    // ---------------------------------------------
    // Factors
    // ---------------------------------------------

    //! Factor for the computation of dN_{BW} / dt from T
    double factor_dNBW_dt_;

    //! Normalized reduced Compton wavelength
    double normalized_Compton_wavelength_;

};

#endif
