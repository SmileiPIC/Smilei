// ----------------------------------------------------------------------------
//! \file Radiation.h
//
//! \brief This file contains the header for the generic class Radiation
//   for the particle radiation losses.
//
// ----------------------------------------------------------------------------

#ifndef RADIATION_H
#define RADIATION_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "Params.h"
#include "Particles.h"
#include "Species.h"
#include "NlicsTables.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class Radiation
//  --------------------------------------------------------------------------------------------------------------------
class Radiation
{

public:
    //! Creator for Radiation
    Radiation(Params& params, Species *species);
    virtual ~Radiation();

    //! Overloading of () operator
    virtual void operator() (Particles &particles,
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
    //! \param Ey y component of the particle electric field
    //! \param Ez z component of the particle electric field
    //! \param Bx x component of the particle magnetic field
    //! \param By y component of the particle magnetic field
    //! \param Bz z component of the particle magnetic field
    double compute_chipa(double & charge_over_mass2,
                                 double & px, double & py, double & pz,
                                 double & gamma,
                                 double & Ex, double & Ey, double & Ez,
                                 double & Bx, double & By, double & Bz);

protected:

private:

    // ________________________________________
    // General parameters

    //! Dimension of position
    int nDim_;

    //! Inverse species mass
    double one_over_mass_;

    //! Time step
    double dt;

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

};//END class

#endif
