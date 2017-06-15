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
RadiationNlicsCont::RadiationNlicsCont(Params& params, Species * species)
      : Radiation(params, species)
{
}

// ---------------------------------------------------------------------------------------------------------------------
//! Destructor for RadiationNlicsCont
// ---------------------------------------------------------------------------------------------------------------------
RadiationNlicsCont::~RadiationNlicsCont()
{
}

// ---------------------------------------------------------------------------------------------------------------------
//! Overloading of the operator (): perform the Discontinuous radiation reaction
//! induced by the nonlinear inverse Compton scattering
//
//! \param particles   particle object containing the particle properties
//! \param smpi        MPI properties
//! \param RadiationTables Cross-section data tables and useful functions
//                     for nonlinear inverse Compton scattering
//! \param istart      Index of the first particle
//! \param iend        Index of the last particle
//! \param ithread     Thread index
// ---------------------------------------------------------------------------------------------------------------------
void RadiationNlicsCont::operator() (Particles &particles,
        SmileiMPI* smpi,
        RadiationTables &RadiationTables,
        int istart,
        int iend,
        int ithread)
{

    // _______________________________________________________________
    // Parameters
    std::vector<LocalFields> *Epart = &(smpi->dynamics_Epart[ithread]);
    std::vector<LocalFields> *Bpart = &(smpi->dynamics_Bpart[ithread]);
    //std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);

    // Charge divided by the square of the mass
    double charge_over_mass2;

    // Temporary quantum parameter
    double chipa;

    // Temporary Lorentz factor
    double gamma;

    // Temporary double parameter
    double temp;

    // Momentum shortcut
    double* momentum[3];
    for ( int i = 0 ; i<3 ; i++ )
        momentum[i] =  &( particles.momentum(i,0) );

    // Charge shortcut
    short* charge = &( particles.charge(0) );

    // Weight shortcut
    double* weight = &( particles.weight(0) );

    // Optical depth for the Monte-Carlo process
    // double* chi = &( particles.chi(0));

    // Local vector to store the radiated energy
    std::vector <double> rad_norm_energy (iend-istart);

    // Reinitialize the cumulative radiated energy for the current thread
    this->radiated_energy = 0.;

    // _______________________________________________________________
    // Computation

    #pragma omp simd
    for (int ipart=istart ; ipart<iend; ipart++ ) {
        charge_over_mass2 = (double)(charge[ipart])*pow(one_over_mass_,2.);

        // Gamma
        gamma = sqrt(1.0 + momentum[0][ipart]*momentum[0][ipart]
                             + momentum[1][ipart]*momentum[1][ipart]
                             + momentum[2][ipart]*momentum[2][ipart]);

        // Computation of the Lorentz invariant quantum parameter
        chipa = Radiation::compute_chipa(charge_over_mass2,
                     momentum[0][ipart],momentum[1][ipart],momentum[2][ipart],
                     gamma,
                     (*Epart)[ipart].x,(*Epart)[ipart].y,(*Epart)[ipart].z,
                     (*Bpart)[ipart].x,(*Bpart)[ipart].y,(*Bpart)[ipart].z);

        // Update the quantum parameter in species
        // chi[ipart] = chipa;

        // Radiated energy during the time step
        rad_norm_energy[ipart - istart] =
        RadiationTables.compute_cont_rad_energy_Ridgers(chipa,dt);

        /*std::cerr << "rad_norm_energy: " << rad_norm_energy << " "
                  << "gamma: " << gamma << " "
                  << std::endl;*/

        // Effect on the momentum
        // Temporary factor
        temp = rad_norm_energy[ipart - istart]/gamma;
        // Update of the momentum
        momentum[0][ipart] -= temp*momentum[0][ipart];
        momentum[1][ipart] -= temp*momentum[1][ipart];
        momentum[2][ipart] -= temp*momentum[2][ipart];

        // Exact energy loss due to the radiation
        rad_norm_energy[ipart - istart] = gamma - sqrt(1.0
                                     + momentum[0][ipart]*momentum[0][ipart]
                                     + momentum[1][ipart]*momentum[1][ipart]
                                     + momentum[2][ipart]*momentum[2][ipart]);
    }

    // _______________________________________________________________
    // Computation of the thread radiated energy

    double radiated_energy_loc = 0;

    #pragma omp simd reduction(+:radiated_energy_loc)
    for (int ipart=istart ; ipart<iend; ipart++ )
    {
        radiated_energy_loc += weight[ipart]*rad_norm_energy[ipart - istart] ;
        /*std::cerr << weight[ipart]
                  << " " << rad_norm_energy[ipart - istart]
                  << std::endl;*/
    }
    radiated_energy += radiated_energy_loc;
}
