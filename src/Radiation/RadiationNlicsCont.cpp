// ----------------------------------------------------------------------------
//! \file RadiationNlicsCont.cpp
//
//! \brief This class is for the continuous radiation loss.
//!        This model includes a quantum correction.
//
//! The implementation is adapted from the thesis results of M. Lobet
//! See http://www.theses.fr/2015BORD0361
//
// ----------------------------------------------------------------------------

#include "RadiationNlicsCont.h"

#include <cstring>
#include <fstream>

#include <cmath>

// ---------------------------------------------------------------------------------------------------------------------
//! Constructor for RadiationNlicsCont
//! Inherited from Radiation
// ---------------------------------------------------------------------------------------------------------------------
RadiationNlicsMC::RadiationNlicsMC(Params& params, Species * species)
      : Radiation(params, species)
{
}

// ---------------------------------------------------------------------------------------------------------------------
//! Destructor for RadiationNlicsCont
// ---------------------------------------------------------------------------------------------------------------------
RadiationNlicsMC::~RadiationNlicsMC()
{
}

// ---------------------------------------------------------------------------------------------------------------------
//! Overloading of the operator (): perform the Discontinuous radiation reaction
//! induced by the nonlinear inverse Compton scattering
//
//! \param particles   particle object containing the particle properties
//! \param smpi        MPI properties
//! \param nlicsTables Cross-section data tables and useful functions
//                     for nonlinear inverse Compton scattering
//! \param istart      Index of the first particle
//! \param iend        Index of the last particle
//! \param ithread     Thread index
// ---------------------------------------------------------------------------------------------------------------------
void RadiationNlicsMC::operator() (Particles &particles,
        SmileiMPI* smpi,
        NlicsTables &nlicsTables,
        int istart,
        int iend,
        int ithread)
{

    // _______________________________________________________________
    // Parameters
    std::vector<LocalFields> *Epart = &(smpi->dynamics_Epart[ithread]);
    std::vector<LocalFields> *Bpart = &(smpi->dynamics_Bpart[ithread]);
    std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);

    // Charge divided by the square of the mass
    double charge_over_mass2;

    // Temporary quantum parameter
    double chipa;

    // Temporary Lorentz factor
    double gamma;

    // Radiated energy
    double rad_energy;

    // Temporary double parameter
    double temp;

    // Momentum shortcut
    double* momentum[3];
    for ( int i = 0 ; i<3 ; i++ )
        momentum[i] =  &( particles.momentum(i,0) );

    // Charge shortcut
    short* charge = &( particles.charge(0) );

    // Optical depth for the Monte-Carlo process
    double* chi = &( particles.chi(0));

    // _______________________________________________________________
    // Computation

    for (int ipart=istart ; ipart<iend; ipart++ ) {
        charge_over_mass2 = (double)(charge[ipart])*pow(one_over_mass_,2.);

        // Gamma
        gamma = sqrt(1.0 + momentum[0][ipart]*momentum[0][ipart]
                             + momentum[1][ipart]*momentum[1][ipart]
                             + momentum[2][ipart]*momentum[2][ipart]);

        // Computation of the Lorentz invariant quantum parameter
        chipa = RadiationNlicsMC::compute_chipa(charge_over_mass2,
                     momentum[0][ipart],momentum[1][ipart],momentum[2][ipart],
                     gamma,
                     (*Epart)[ipart].x,(*Epart)[ipart].y,(*Epart)[ipart].z,
                     (*Bpart)[ipart].x,(*Bpart)[ipart].y,(*Bpart)[ipart].z);

        // Update the quantum parameter in species
        //chi[ipart] = chipa

        // Radiated energy during the time step
        rad_energy =
        nlicsTables.compute_cont_rad_energy_Ridgers(chipa,dt);

        // Effect on the momentum
        // Temporary factor
        temp = rad_energy*(*invgf)[ipart];
        // Update of the momentum
        momentum[0][ipart] -= temp*momentum[0][ipart];
        momentum[1][ipart] -= temp*momentum[1][ipart];
        momentum[2][ipart] -= temp*momentum[2][ipart];

    }
}
