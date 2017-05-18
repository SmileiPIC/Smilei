// ----------------------------------------------------------------------------
//! \file Nlics.cpp
//
//! \brief This class performs the Nonlinear Inverse Compton Scattering
//! on particles.
//
//! The implementation is adapted from the thesis results of M. Lobet
//! See http://www.theses.fr/2015BORD0361
//
// ----------------------------------------------------------------------------

#include "Nlics.h"

#include <cstring>
#include <iostream>
#include <fstream>

#include <cmath>

#include "userFunctions.h"

// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Nlics
// ---------------------------------------------------------------------------------------------------------------------
Nlics::Nlics(Params& params, Species * species)
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
// Destructor for Nlics
// ---------------------------------------------------------------------------------------------------------------------
Nlics::~Nlics()
{
}


// ---------------------------------------------------------------------------------------------------------------------
//! Overloading of the operator (): perform the Discontinuous radiation reaction
//! induced by the nonlinear inverse Compton scattering
//
//! \param particles particle object containing the particle properties
//! \param smpi
//! \param nlicsTables
//! \param istart
//! \param iend
//! \param ithread
// ---------------------------------------------------------------------------------------------------------------------
void Nlics::operator() (Particles &particles,
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

    // Time to emission
    double emission_time;

    // time spent in the iteration
    double local_it_time;

    // Number of Monte-Carlo iteration
    int mc_it_nb;

    // Momentum shortcut
    double* momentum[3];
    for ( int i = 0 ; i<3 ; i++ )
        momentum[i] =  &( particles.momentum(i,0) );

    // Charge shortcut
    short* charge = &( particles.charge(0) );

    // Optical depth for the Monte-Carlo process
    double* tau = &( particles.tau(0));

    // Optical depth for the Monte-Carlo process
    double* chi = &( particles.chi(0));

    // _______________________________________________________________
    // Computation

    for (int ipart=istart ; ipart<iend; ipart++ ) {
        charge_over_mass2 = (double)(charge[ipart])*pow(one_over_mass_,2.);

        // Init local variables
        emission_time = 0;
        local_it_time = 0;
        mc_it_nb = 0;

        // Monte-Carlo Manager inside the time step
        while ((local_it_time < dt)
             &&(mc_it_nb < mc_it_nb_max))
        {

            // Gamma
            gamma = 1./(*invgf)[ipart];

            // Computation of the Lorentz invariant quantum parameter
            chipa = Nlics::compute_chipa(charge_over_mass2,
                     momentum[0][ipart],momentum[1][ipart],momentum[2][ipart],
                     gamma,
                     (*Epart)[ipart].x,(*Epart)[ipart].y,(*Epart)[ipart].z,
                     (*Bpart)[ipart].x,(*Bpart)[ipart].y,(*Bpart)[ipart].z);

            // Discontinuous emission: New emission
            // If tau[ipart] <= 0, this is a new emission
            // We also check that chipa > chipa_threshold,
            // else chipa is too low to induce a discontinuous emission
            if ((chipa > chipa_disc_threshold)
            && (tau[ipart] <= epsilon_tau) )
            {
                // New final optical depth to reach for emision
                tau[ipart] = -log( -Rand::uniform() + 1);
            }

            // Discontinuous emission: emission under progress
            // If epsilon_tau > 0
            if (tau[ipart] > epsilon_tau)
            {
                // from the cross section
                temp = nlicsTables.compute_dNphdt(chipa,gamma);

                // Time to discontinuous emission
                // If this time is > the remaining iteration time,
                // we have a synchronization
                emission_time = std::min(tau[ipart]/temp, dt - local_it_time);

                // Update of the optical depth
                tau[ipart] -= temp*emission_time;

                // If the final optical depth is reached
                if (tau[ipart] <= epsilon_tau)
                {

                    // Emission of a photon
                    Nlics::photon_emission();

                    // Optical depth becomes negative meaning
                    // that a new drawing is possible
                    // at the next Monte-Carlo iteration
                    tau[ipart] = -1;
                }

                // Incrementation of the Monte-Carlo iteration counter
                mc_it_nb ++;
                // Update of the local time
                local_it_time += emission_time;

            }

            // Continuous emission
            // chipa needs to be below the discontinuous threshold
            // chipa needs to be above the continuous threshold
            // No discontunous emission is in progress:
            // tau[ipart] <= epsilon_tau
            else if ((chipa <= chipa_disc_threshold)
            &&  (tau[ipart] <= epsilon_tau)
            &&  (chipa > chipa_cont_threshold))
            {

                // Remaining time of the iteration
                emission_time = dt - local_it_time;

                // Radiated energy during emission_time
                rad_energy =
                nlicsTables.compute_cont_rad_energy_Ridgers(chipa,
                                                             emission_time);

                // Effect on the momentum
                temp = rad_energy*(*invgf)[ipart];
                for ( int i = 0 ; i<3 ; i++ )
                    momentum[i][ipart] -= temp*momentum[i][ipart];

                local_it_time = dt;
            }

        }

    }
}

// ---------------------------------------------------------------------------------------------------------------------
//! Computation of the Lorentz invariant quantum parameter for particles
//
//! \param charge_over_mass2 charge divided by the square of the mass
//! \param px particle x momentum
//! \param gamma particle Lorentz factor
//! \param Ex x component of the particle electric field
//! \param Bx x component of the particle magnetic field
// ---------------------------------------------------------------------------------------------------------------------
double Nlics::compute_chipa(double & charge_over_mass2,
                             double & px, double & py, double & pz,
                             double & gamma,
                             double & Ex, double & Ey, double & Ez,
                             double & Bx, double & By, double & Bz)
{
    double chipa;

    chipa = charge_over_mass2*inv_norm_E_Schwinger
          * sqrt( pow(Ex*px + Ey*py + Ez*pz,2)
          - pow(gamma*Ex - By*pz + Bz*py,2)
          - pow(gamma*Ey - Bz*px + Bx*pz,2)
          - pow(gamma*Ez - Bx*py + By*px,2));

    return chipa;
}

void Nlics::photon_emission()
{
    
}
