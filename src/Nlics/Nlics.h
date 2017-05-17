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
                int istart,
                int iend,
                int ithread);


        //! Computation of the Lorentz invariant quantum parameter for particles
        //! \param charge_over_mass2 charge divided by the square of the mass
        //! \param px particle x momentum
        //! \param gamma particle Lorentz factor
        //! \param Ex x component of the particle electric field
        //! \param Bx x component of the particle magnetic field
        double compute_chipa(double & charge_over_mass2,
                                     double & px, double & py, double & pz,
                                     double & gamma,
                                     double & Ex, double & Ey, double & Ez,
                                     double & Bx, double & By, double & Bz);

    private:

        // ________________________________________
        // General parameters

        //! Dimension of position
        int nDim_;

        //! Inverse species mass
        double one_over_mass_;

        //! Time step
        double dt;

        // Max number of Monte-Carlo iteration
        const int mc_it_nb_max = 100;

        // _________________________________________
        // Factors

        //! Fine structure constant
        const double fine_struct_cst = 7.2973525698e-3;

        //! Reduced Planck Constant (J.s)
        const double red_planck_cst = 1.054571628E-34;

        //! Electron mass
        const double electron_mass = 9.109382616e-31;

        //! Speed of light in vacuum (m/s)
        const double c_vacuum = 299792458;

};

#endif
