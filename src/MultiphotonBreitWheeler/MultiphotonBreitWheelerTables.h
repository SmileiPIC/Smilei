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
    void initialization( Params &params, SmileiMPI *smpi );

    // ---------------------------------------------------------------------
    // PHYSICAL COMPUTATION
    // ---------------------------------------------------------------------

    //! Computation of the production rate of pairs per photon
    //! \param photon_chi photon quantum parameter
    //! \param gamma photon normalized energy
    double computeBreitWheelerPairProductionRate( double photon_chi, double gamma );

    //! Computation of the electron and positron quantum parameters for
    //! the multiphoton Breit-Wheeler pair creation
    //! \param photon_chi photon quantum parameter
    double *computePairQuantumParameter( double photon_chi );


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

    // ---------------------------------------------------------------------
    // TABLE COMMUNICATIONS
    // ---------------------------------------------------------------------

    //! Bcast of the external table T
    //! \param smpi Object of class SmileiMPI containing MPI properties
    void bcastTableT( SmileiMPI *smpi );

    //! Bcast of the external table xip_chipamin and xip
    //! \param smpi Object of class SmileiMPI containing MPI properties
    void bcastTableXi( SmileiMPI *smpi );

private:

    // ---------------------------------------------
    // General parameters
    // ---------------------------------------------

    //! Path to the tables
    std::string table_path_;

    // ---------------------------------------------
    // Structure for Table T used for the
    // pair creation Monte-Carlo process
    // ---------------------------------------------

    struct T {
        
        //! Array containing tabulated values of the function T
        std::vector<double > table_;
        
        //! Minimum boundary of the table T
        double chiph_min;
        
        //! Log10 of the minimum boundary of the table T
        double log10_chiph_min;
        
        //! Maximum boundary of the table T
        double chiph_max;
        
    };

    struct T T_;

    //! Delta chi for the table T
    double T_chiph_delta;

    //! Inverse delta chi for the table h
    double T_chiph_inv_delta;

    //! Dimension of the array T
    int T_chiph_dim;

    // ---------------------------------------------
    // Structure for xi and particle_chi min for xip table
    // ---------------------------------------------

    struct Xi {
        
        //! Table containing the cumulative distribution function \f$P(0 \rightarrow \chi_{e^-})\f$
        //! that gives gives the probability for a photon to decay into pair
        //! with an electron of energy in the range \f$[0, \chi_{e^-}]\f$
        //! This enables to compute the energy repartition between the electron and the positron
        std::vector<double> table;
        
        //! Table containing the particle_chi min values
        //! Under this value, electron kinetic energy of the pair is
        //! considered negligible
        std::vector<double > min_particle_chi_;
        
        //! Minimum boundary for photon_chi in the table xi and xi_.chipamin
        double chiph_min;
        
        //! Logarithm of the minimum boundary for photon_chi in the table xi
        //! and xi_.chipamin
        double log10_chiph_min;
        
        //! Maximum boundary for photon_chi in the table xip and xip_chipamin
        double chiph_max;
        
        //! Delta for the photon_chi discretization in the table xip and xip_chipamin
        double chiph_delta;
        
        //! Inverse of the delta for the photon_chi discretization
        //! in the table xip and xip_chipamin
        double chiph_inv_delta;

        //! Dimension of the discretized parameter photon_chi
        int chiph_dim;

        //! Dimension of the discretized parameter particle_chi
        int chipa_dim;
        
        //! xip power
        double power;

        //! 1/(xi_.chipa_dim - 1)
        double inv_chipa_dim_minus_one;

        //! xip threshold
        double threshold;
        
    };
    
    struct Xi xi_;


    // ---------------------------------------------
    // Factors
    // ---------------------------------------------

    //! Factor for the computation of dN_{BW} / dt from T
    double factor_dNBW_dt_;

    //! Normalized reduced Compton wavelength
    double normalized_Compton_wavelength_;

};

#endif
