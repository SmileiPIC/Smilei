// ----------------------------------------------------------------------------
//! \file RadiationLL.h
//
//! \brief This class is for the classical continuous radiation loss with
//!        the Landau-Lifshitz model.
//!        This model does not include a quantum correction.
//
//! \details This header contains the definition of the class RadiationLL.
//! L. D. Landau and E. M. Lifshitz, The classical theory of fields, 1947
//! F. Niel et al., 2017
// ----------------------------------------------------------------------------

#ifndef RADIATIONLL_H
#define RADIATIONLL_H

#include "RadiationTables.h"
#include "Radiation.h"
#include "userFunctions.h"

#include <cstring>
#include <fstream>
#include <cmath>

//------------------------------------------------------------------------------
//! RadiationLL class: holds parameters and functions to apply the
//! Landau-Lifshitz continuous radiation loss on Particles.
//------------------------------------------------------------------------------
class RadiationLL : public Radiation {

    public:

        //! Constructor for RadiationLL
        RadiationLL(Params& params, Species * species);

        //! Destructor for RadiationLL
        ~RadiationLL();

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
        // ---------------------------------------------------------------------
        virtual void operator() (Particles &particles,
                SmileiMPI* smpi,
                RadiationTables &RadiationTables,
                int istart,
                int iend,
                int ithread);

    protected:

        // ________________________________________
        // General parameters


    private:

};

#endif
