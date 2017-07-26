// ----------------------------------------------------------------------------
//! \file MultiphotonBreitWheelerTables.cpp
//
//! \brief This file is for the methods of the class
//! MultiphotonBreitWheelerTables.
//! This class contains the methods and tools to generate and manage
//! the physical tables for the multiphoton Breit-wheeler process.
//
//! \details
//! The implementation is adapted from the thesis results of M. Lobet
//! See http://www.theses.fr/2015BORD0361
// ----------------------------------------------------------------------------

#include "MultiphotonBreitWheelerTables.h"

// -----------------------------------------------------------------------------
// INITILIZATION AND DESTRUCTION
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Constructor for MutliphotonBreitWheelerTables
// -----------------------------------------------------------------------------
MultiphotonBreitWheelerTables::MultiphotonBreitWheelerTables()
{
    // T table
    T_table.resize(0);
    T_computed = false;
}

// -----------------------------------------------------------------------------
// Destructor for MutliphotonBreitWheelerTables
// -----------------------------------------------------------------------------
MultiphotonBreitWheelerTables::~MultiphotonBreitWheelerTables()
{
}

// -----------------------------------------------------------------------------
//
//! Initialization of the parameters for the multiphoton Breit-Wheeler pair
//! creation
//
//! \param params Object Params for the parameters from the input script
// -----------------------------------------------------------------------------
void MultiphotonBreitWheelerTables::initialization(Params& params)
{
    if (params.hasMultiphotonBreitWheeler)
    {
        TITLE("Initializing mutliphoton Breit-Wheeler")

        // Preliminary checks
        if (params.referenceAngularFrequency_SI <= 0.)
            ERROR("The parameter `referenceAngularFrequency_SI` needs "
                  << "to be defined and positive to compute radiation losses");

    }

    // If the namelist for Nonlinear Inverse Compton Scattering exists
    // We read the properties
    if(PyTools::nComponents("MultiphotonBreitWheeler"))
    {

        // Extraction of the parameter from the input file
        PyTools::extract("T_chiph_min", T_chiph_min, "MultiphotonBreitWheeler");
        PyTools::extract("T_chiph_max", T_chiph_max, "MultiphotonBreitWheeler");
        PyTools::extract("T_dim", T_dim, "MultiphotonBreitWheeler");

        T_log10_chiph_min = log10(T_chiph_min);

        // Format of the tables
        PyTools::extract("output_format", output_format, "MultiphotonBreitWheeler");

        // Path to the databases
        PyTools::extract("table_path", table_path, "MultiphotonBreitWheeler");
    }

    // Messages and checks
    if (params.hasMultiphotonBreitWheeler)
    {
        MESSAGE( "        table path: " << table_path);

        if (T_chiph_min >= T_chiph_max)
        {
            ERROR("T_chiph_min (" << T_chiph_min
                    << ") >= T_chiph_max (" << T_chiph_max << ")")
        }
    }

    MESSAGE("")
}
