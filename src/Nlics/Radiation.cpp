// ----------------------------------------------------------------------------
//! \file Radiation.cpp
//
//! \brief This file contains the class functions for the generic class
//  Radiation for the particle radiation losses.
//
// ----------------------------------------------------------------------------

#include "Radiation.h"

// ---------------------------------------------------------------------------------------------------------------------
//! Constructor for Radiation
// input: simulation parameters & Species index
//! \param params simulation parameters
//! \param species Species index
// ---------------------------------------------------------------------------------------------------------------------
Radiation::Radiation(Params& params, Species * species)
{
    // Dimension position
    nDim_ = params.nDim_particle;

    // Time step
    dt    = params.timestep;

    // Inverse of the species mass
    one_over_mass_ = 1./species->mass;

    // Normalized Schwinger Electric Field
    norm_E_Schwinger = electron_mass*c_vacuum*c_vacuum
                     / (red_planck_cst*params.referenceAngularFrequency_SI);

    // Inverse of norm_E_Schwinger
    inv_norm_E_Schwinger = 1./norm_E_Schwinger;
}

// ---------------------------------------------------------------------------------------------------------------------
//! Computation of the Lorentz invariant quantum parameter for particles
//
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
// ---------------------------------------------------------------------------------------------------------------------
double Nlics::compute_chipa(double & charge_over_mass2,
                             double & px, double & py, double & pz,
                             double & gamma,
                             double & Ex, double & Ey, double & Ez,
                             double & Bx, double & By, double & Bz)
{
    double chipa;

    chipa = fabs(charge_over_mass2)*inv_norm_E_Schwinger
          * sqrt( fabs( pow(Ex*px + Ey*py + Ez*pz,2)
          - pow(gamma*Ex - By*pz + Bz*py,2)
          - pow(gamma*Ey - Bz*px + Bx*pz,2)
          - pow(gamma*Ez - Bx*py + By*px,2)));

    /*std::cerr << charge_over_mass2 << " " << inv_norm_E_Schwinger << " "
              << "P: " << px << " " << py << " " << pz << " " << " "
              << "Gamma: "<< gamma << " "
              << "E: " << Ex << " " << Ey << " " << Ez << " "
              << "B: " <<  Bx << " " << By << " " << Bz << " "
              << "chipa: " << chipa << std::endl;*/

    return chipa;
}
