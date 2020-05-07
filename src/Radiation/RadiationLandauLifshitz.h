// ----------------------------------------------------------------------------
//! \file RadiationLandauLifshitz.h
//
//! \brief This class is for the classical continuous radiation loss with
//!        the Landau-Lifshitz model.
//!        This model does not include a quantum correction.
//
//! \details This header contains the definition of the class RadiationLandauLifshitz.
//! L. D. Landau and E. M. Lifshitz, The classical theory of fields, 1947
//! F. Niel et al., 2017
// ----------------------------------------------------------------------------

#ifndef RADIATIONLANDAULIFSHITZ_H
#define RADIATIONLANDAULIFSHITZ_H

#include "RadiationTables.h"
#include "Radiation.h"
#include "userFunctions.h"

#include <cstring>
#include <fstream>
#include <cmath>

//------------------------------------------------------------------------------
//! RadiationLandauLifshitz class: holds parameters and functions to apply the
//! Landau-Lifshitz continuous radiation loss on Particles.
//------------------------------------------------------------------------------
class RadiationLandauLifshitz : public Radiation
{

public:

    //! Constructor for RadiationLandauLifshitz
    RadiationLandauLifshitz( Params &params, Species *species, Random * rand  );
    
    //! Destructor for RadiationLandauLifshitz
    ~RadiationLandauLifshitz();
    
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
        Particles       &particles,
        Species         *photon_species,
        SmileiMPI       *smpi,
        RadiationTables &RadiationTables,
        double          &radiated_energy,
        int istart,
        int iend,
        int ithread,
        int ipart_ref = 0);
        
protected:

    // ________________________________________
    // General parameters
    
    
private:

};

#endif
