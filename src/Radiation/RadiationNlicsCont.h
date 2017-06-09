// ----------------------------------------------------------------------------
//! \file RadiationNlicsCont.h
//
//! \brief This class is for the continuous radiation loss.
//!        This model includes a quantum correction.
//
//! \details This header contains the definition of the class RadiationNlicsCont.
//! The implementation is adapted from the thesis results of M. Lobet
//! See http://www.theses.fr/2015BORD0361
// ----------------------------------------------------------------------------

#ifndef RADIATIONNLICSCONT_H
#define RADIATIONNLICSCONT_H

#include "RadiationTables.h"
#include "Radiation.h"
#include "userFunctions.h"

//------------------------------------------------------------------------------
//! RadiationNlicsCont class: holds parameters and functions to apply the
//! continuous radiation loss on Particles,
//------------------------------------------------------------------------------
class RadiationNlicsCont : public Radiation {

    public:

        //! Constructor for RadiationNlicsCont
        RadiationNlicsCont(Params& params, Species * species);

        //! Destructor for RadiationNlicsCont
        ~RadiationNlicsCont();

        // ---------------------------------------------------------------------
        //! Overloading of () operator: perform the Discontinuous radiation
        //! reaction induced by the nonlinear inverse Compton scattering
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

        //! Under this value, no radiation loss
        const double chipa_cont_threshold = 1e-5;

    private:

};

#endif
