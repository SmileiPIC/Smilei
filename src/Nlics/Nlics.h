// ----------------------------------------------------------------------------
//! \file Nlics.h
//
//! \brief This class performs the Nonlinear Inverse Compton Scattering
//! on particles.
//
//! \details This header contains the definition of the class Nlics.
//! The implementation is adapted from the thesis results of M. Lobet
//! See http://www.theses.fr/2015BORD0361
// ----------------------------------------------------------------------------

#ifndef NLICS_H
#define NLICS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "Params.h"
#include "Particles.h"
#include "Species.h"
#include "NlicsTables.h"

//----------------------------------------------------------------------------------------------------------------------
//! Nlics class: holds parameters and functions to apply the nonlinear inverse
//! Compton scattering on Particles,
//----------------------------------------------------------------------------------------------------------------------
class Nlics
{

    public:

        //! Constructor for Nlics
        Nlics(Params& params, Species * species);

        //! Destructor for Nlics
        ~Nlics();

        //! Overloading of () operator
        void operator() (Particles &particles,
                SmileiMPI* smpi,
                NlicsTables &nlicsTables,
                int istart,
                int iend,
                int ithread);


        //! Computation of the Lorentz invariant quantum parameter for particles
        //! \param charge_over_mass2 charge divided by the square of the mass
        //! \param px particle x momentum
        //! \param py particle y momentum
        //! \param pz particle z momentum
        //! \param gamma particle Lorentz factor
        //! \param Ex x component of the particle electric field
        //! \param Bx x component of the particle magnetic field
        double compute_chipa(double & charge_over_mass2,
                                     double & px, double & py, double & pz,
                                     double & gamma,
                                     double & Ex, double & Ey, double & Ez,
                                     double & Bx, double & By, double & Bz);

        // ---------------------------------------------------------------------------------------------------------------------
        //! Perform the phoon emission (creation of a super-photon
        //! and slow down of the emitting particle)
        //! \param chipa          particle quantum parameter
        //! \param gammapa          particle gamma factor
        //! \param px             particle momentum in x
        //! \param py             particle momentum in y
        //! \param pz             particle momentum in z
        //! \param nlicsTables    Cross-section data tables and useful functions
        //                        for nonlinear inverse Compton scattering
        // ---------------------------------------------------------------------------------------------------------------------
        void photon_emission(double &chipa,
                             double & gammapa,
                             double & px,
                             double & py,
                             double & pz,
                             NlicsTables &nlicsTables);

    private:

        // ________________________________________
        // General parameters

        //! Dimension of position
        int nDim_;

        //! Inverse species mass
        double one_over_mass_;

        //! Time step
        double dt;

        //! Max number of Monte-Carlo iteration
        const int mc_it_nb_max = 100;

        //! Espilon to check when tau is near 0
        const double epsilon_tau = 1e-100;

        //! Under this value, the emission is considered continuous
        const double chipa_disc_threshold = 1e-3;

        //! Under this value, no radiation loss
        const double chipa_cont_threshold = 1e-5;

        // _________________________________________
        // Factors

        //! Fine structure constant
        const double fine_struct_cst = 7.2973525698e-3;

        //! Reduced Planck Constant (J.s)
        const double red_planck_cst = 1.054571628E-34;

        //! Electron mass (kg)
        const double electron_mass = 9.109382616e-31;

        //! Speed of light in vacuum (m/s)
        const double c_vacuum = 299792458;

        //! Normalized Schwinger Electric field
        double norm_E_Schwinger;

        //! Inverse Normalized Schwinger Electric field
        double inv_norm_E_Schwinger;

};

#endif
