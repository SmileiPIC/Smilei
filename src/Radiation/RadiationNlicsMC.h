// ----------------------------------------------------------------------------
//! \file RadiationNlicsMC.h
//
//! \brief This class performs the Nonlinear Inverse Compton Scattering
//! on particles.
//
//! \details This header contains the definition of the class RadiationNlicsMC.
//! The implementation is adapted from the thesis results of M. Lobet
//! See http://www.theses.fr/2015BORD0361
// ----------------------------------------------------------------------------

#ifndef RADIATIONNLICSMC_H
#define RADIATIONNLICSMC_H

#include "RadiationTables.h"
#include "Radiation.h"
#include "userFunctions.h"

//----------------------------------------------------------------------------------------------------------------------
//! RadiationNlicsMC class: holds parameters and functions to apply the
//! nonlinear inverse Compton scattering on Particles,
//----------------------------------------------------------------------------------------------------------------------
class RadiationNlicsMC : public Radiation {

    public:

        //! Constructor for RadiationNlicsMC
        RadiationNlicsMC(Params& params, Species * species);

        //! Destructor for RadiationNlicsMC
        ~RadiationNlicsMC();

        // ---------------------------------------------------------------------
        //! Overloading of () operator: perform the Discontinuous radiation
        //! reaction induced by the nonlinear inverse Compton scattering
        //! \param particles   particle object containing the particle
        //!                    properties
        //! \param smpi        MPI properties
        //! \param RadiationTables Cross-section data tables and useful functions
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

        // ---------------------------------------------------------------------
        //! Perform the phoon emission (creation of a super-photon
        //! and slow down of the emitting particle)
        //! \param chipa          particle quantum parameter
        //! \param gammapa          particle gamma factor
        //! \param px             particle momentum in x
        //! \param py             particle momentum in y
        //! \param pz             particle momentum in z
        //! \param RadiationTables    Cross-section data tables and useful functions
        //                        for nonlinear inverse Compton scattering
        // ---------------------------------------------------------------------
        void photon_emission(double &chipa,
                             double & gammapa,
                             double & px,
                             double & py,
                             double & pz,
                             RadiationTables &RadiationTables);

    protected:

        // ________________________________________
        // General parameters

        //! Max number of Monte-Carlo iteration
        const int mc_it_nb_max = 100;

        //! Espilon to check when tau is near 0
        const double epsilon_tau = 1e-100;

        //! Under this value, the emission is considered continuous
        const double chipa_disc_threshold = 1e-2;

        //! Under this value, no radiation loss
        const double chipa_cont_threshold = 1e-5;

    private:

};

#endif
