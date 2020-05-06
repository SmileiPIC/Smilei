// ----------------------------------------------------------------------------
//! \file RadiationCorrLandauLifshitz.h
//
//! \brief This class is for the continuous radiation loss.
//!        This model includes a quantum correction.
//
//! \details This header contains the definition of the class RadiationCorrLandauLifshitz.
//! The implementation is adapted from the thesis results of M. Lobet
//! See http://www.theses.fr/2015BORD0361
// ----------------------------------------------------------------------------

#ifndef RADIATIONNLICSCONT_H
#define RADIATIONNLICSCONT_H

#include "RadiationTables.h"
#include "Radiation.h"
#include "userFunctions.h"

//------------------------------------------------------------------------------
//! RadiationCorrLandauLifshitz class: holds parameters and functions to apply the
//! continuous radiation loss on particles.
//------------------------------------------------------------------------------
class RadiationCorrLandauLifshitz : public Radiation
{

public:

    //! Constructor for RadiationCorrLandauLifshitz
    RadiationCorrLandauLifshitz( Params &params, Species *species, Random * rand  );
    
    //! Destructor for RadiationCorrLandauLifshitz
    ~RadiationCorrLandauLifshitz();
    
    // ---------------------------------------------------------------------
    //! Overloading of () operator: perform the discontinuous radiation
    //! reaction induced by the nonlinear inverse Compton scattering
    //! \param particles   particle object containing the particle
    //!                    properties
    //! \param photon_species species that will receive emitted photons
    //! \param smpi        MPI properties
    //! \param nlicsTables Cross-section data tables and useful functions
    //                     for nonlinear inverse Compton scattering
    //! \param istart      Index of the first particle
    //! \param iend        Index of the last particle
    //! \param ithread     Thread index
    //! \param radiated_energy     overall energy radiated during the call to this method
    // ---------------------------------------------------------------------
    virtual void operator()(
        Particles       &particles,
        Species         *photon_species,
        SmileiMPI       *smpi,
        RadiationTables &RadiationTables,
        double          &radiated_energy,
        int             istart,
        int             iend,
        int             ithread,
        int             ipart_ref = 0);
        
protected:

    // ________________________________________
    // General parameters
    
    //! Under this value, no radiation loss
    const double minimum_chi_continuous_ = 1e-5;
    
private:

};

#endif
