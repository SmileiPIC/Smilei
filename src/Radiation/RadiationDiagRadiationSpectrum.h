// ----------------------------------------------------------------------------
//! \file RadiationDiagRadiationSpectrum.h
// ----------------------------------------------------------------------------

#ifndef RADIATIONDIAGRADIATIONSPECTRUM_H
#define RADIATIONDIAGRADIATIONSPECTRUM_H

#include "RadiationTables.h"
#include "Radiation.h"
#include "userFunctions.h"

#include <cstring>
#include <fstream>
#include <cmath>

//------------------------------------------------------------------------------
//! RadiationDiagRadiationSpectrum class: holds parameters and functions to apply the
//! Landau-Lifshitz continuous radiation loss on Particles.
//------------------------------------------------------------------------------
class RadiationDiagRadiationSpectrum : public Radiation {

    public:

        //! Constructor for RadiationDiagRadiationSpectrum
        RadiationDiagRadiationSpectrum(Params& params, Species * species, Random * rand );

        //! Destructor for RadiationDiagRadiationSpectrum
        ~RadiationDiagRadiationSpectrum();

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
        virtual void operator()(
            Particles &particles,
            Species *photon_species,
            SmileiMPI *smpi,
            RadiationTables &RadiationTables,
            double          &radiated_energy,
            int istart,
            int iend,
            int ithread, int ipart_ref = 0 );

    protected:

        // ________________________________________
        // General parameters


    private:

};

#endif
