// ----------------------------------------------------------------------------
//! \file RadiationNiel.h
//
//! \brief This header desribes the class RadiationNiel.
//!        This class is for the semi-classical Fokker-Planck model of
//!        synchrotron-like radiation loss developped by Niel et al. as an
//!        extension of the classical Landau-Lifshitz model in the weak quantum
//!        regime.
//!        This model includew a quantum correction + stochastic diffusive
//!        operator.
//
//! \details See these references for more information.
//! F. Niel et al., 2017
//! L. D. Landau and E. M. Lifshitz, The classical theory of fields, 1947
// ----------------------------------------------------------------------------

#ifndef RADIATIONNIEL_H
#define RADIATIONNIEL_H

#include "RadiationTables.h"
#include "Radiation.h"
#include "RadiationTools.h"
#include "userFunctions.h"

#include <cstring>
#include <fstream>
#include <cmath>

//------------------------------------------------------------------------------
//! RadiationLL class: holds parameters and functions to apply the
//! Landau-Lifshitz continuous radiation loss on Particles.
//------------------------------------------------------------------------------
class RadiationNiel : public Radiation
{

public:

    //! Constructor for RadiationLL
    RadiationNiel( Params &params, Species *species, Random * rand  );
    
    //! Destructor for RadiationLL
    ~RadiationNiel();
    
    // ---------------------------------------------------------------------
    //! Overloading of () operator: perform the Landau-Lifshitz
    //! radiation loops.
    //! \param particles   particle object containing the particle
    //!                    properties
    //! \param smpi        MPI properties
    //! \param nlicsTables Cross-section data tables and useful functions
    //                     for nonlinear inverse Compton scattering
    //! \param istart      Index of the first particle
    //! \param iend        Index of the last particle
    //! \param ithread     Thread index
    //! \param radiated_energy     overall energy radiated during the call to this method
    // ---------------------------------------------------------------------
    virtual void operator()(
        Particles &particles,
        Species *photon_species,
        SmileiMPI       *smpi,
        RadiationTables &RadiationTables,
        double          &radiated_energy,
        int             istart,
        int             iend,
        int             ithread,
        int             ipart_ref = 0
        );
        
protected:

    // ________________________________________
    // General parameters
    
    
private:

};

#endif
