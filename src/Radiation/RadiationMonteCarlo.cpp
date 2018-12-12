// ----------------------------------------------------------------------------
//! \file RadiationMonteCarlo.cpp
//
//! \brief This class performs the Nonlinear Inverse Compton Scattering
//! on particles using a Monte-Carlo approach.
//
//! The implementation is adapted from the thesis results of M. Lobet
//! See http://www.theses.fr/2015BORD0361
//
// ----------------------------------------------------------------------------

#include "RadiationMonteCarlo.h"

#include <cstring>
#include <fstream>

#include <cmath>

// ---------------------------------------------------------------------------------------------------------------------
//! Constructor for RadiationMonteCarlo
//! Inherit from Radiation
// ---------------------------------------------------------------------------------------------------------------------
RadiationMonteCarlo::RadiationMonteCarlo(Params& params, Species * species)
      : Radiation(params, species)
{
    this->radiation_photon_sampling_ = species->radiation_photon_sampling_;
    this->radiation_photon_gamma_threshold_ = species->radiation_photon_gamma_threshold_;
    this->inv_radiation_photon_sampling_ = 1. / this->radiation_photon_sampling_;
}

// ---------------------------------------------------------------------------------------------------------------------
//! Destructor for RadiationMonteCarlo
// ---------------------------------------------------------------------------------------------------------------------
RadiationMonteCarlo::~RadiationMonteCarlo()
{
}

// ---------------------------------------------------------------------------------------------------------------------
//! Overloading of the operator (): perform the Discontinuous radiation reaction
//! induced by the nonlinear inverse Compton scattering
//
//! \param particles   particle object containing the particle properties
//! \param photon_species species that will receive emitted photons
//! \param smpi        MPI properties
//! \param RadiationTables Cross-section data tables and useful functions
//                     for nonlinear inverse Compton scattering
//! \param istart      Index of the first particle
//! \param iend        Index of the last particle
//! \param ithread     Thread index
// ---------------------------------------------------------------------------------------------------------------------
void RadiationMonteCarlo::operator() (
        Particles &particles,
        Species * photon_species,
        SmileiMPI* smpi,
        RadiationTables &RadiationTables,
        int istart,
        int iend,
        int ithread, int ipart_ref)
{

    // _______________________________________________________________
    // Parameters
    std::vector<double> *Epart = &(smpi->dynamics_Epart[ithread]);
    std::vector<double> *Bpart = &(smpi->dynamics_Bpart[ithread]);
    //std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);

    int nparts = Epart->size()/3;
    double* Ex = &( (*Epart)[0*nparts] );
    double* Ey = &( (*Epart)[1*nparts] );
    double* Ez = &( (*Epart)[2*nparts] );
    double* Bx = &( (*Bpart)[0*nparts] );
    double* By = &( (*Bpart)[1*nparts] );
    double* Bz = &( (*Bpart)[2*nparts] );

    // Charge divided by the square of the mass
    double charge_over_mass2;

    const double one_over_mass_2 = pow(one_over_mass_,2.);

    // Temporary quantum parameter
    double chipa;

    // Temporary Lorentz factor
    double gamma;

    // Radiated energy
    double cont_rad_energy;

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

    // Position shortcut
    double* position[3];
    for ( int i = 0 ; i<nDim_ ; i++ )
        position[i] =  &( particles.position(i,0) );

    // Charge shortcut
    short* charge = &( particles.charge(0) );

    // Weight shortcut
    double* weight = &( particles.weight(0) );

    // Optical depth for the Monte-Carlo process
    double* tau = &( particles.tau(0));

    // Optical depth for the Monte-Carlo process
    // double* chi = &( particles.chi(0));

    // Reinitialize the cumulative radiated energy for the current thread
    this->radiated_energy = 0.;

    // _______________________________________________________________
    // Computation

    for (int ipart=istart ; ipart<iend; ipart++ ) {
        charge_over_mass2 = (double)(charge[ipart])*one_over_mass_2;

        // Init local variables
        emission_time = 0;
        local_it_time = 0;
        mc_it_nb = 0;

        // Monte-Carlo Manager inside the time step
        while ((local_it_time < dt)
             &&(mc_it_nb < max_monte_carlo_iterations_))
        {

            // Gamma
            gamma = sqrt(1.0 + momentum[0][ipart]*momentum[0][ipart]
                             + momentum[1][ipart]*momentum[1][ipart]
                             + momentum[2][ipart]*momentum[2][ipart]);

            // Computation of the Lorentz invariant quantum parameter
            chipa = Radiation::compute_chipa(charge_over_mass2,
                     momentum[0][ipart],momentum[1][ipart],momentum[2][ipart],
                     gamma,
                     (*(Ex+ipart-ipart_ref)),(*(Ey+ipart-ipart_ref)),(*(Ez+ipart-ipart_ref)),
                     (*(Bx+ipart-ipart_ref)),(*(By+ipart-ipart_ref)),(*(Bz+ipart-ipart_ref)) );

            // Update the quantum parameter in species
            // chi[ipart] = chipa;

            // Discontinuous emission: New emission
            // If tau[ipart] <= 0, this is a new emission
            // We also check that chipa > chipa_threshold,
            // else chipa is too low to induce a discontinuous emission
            if ((chipa > RadiationTables.get_chipa_disc_min_threshold())
            && (tau[ipart] <= epsilon_tau_) )
            {
                // New final optical depth to reach for emision
                while (tau[ipart] <= epsilon_tau_)
                   tau[ipart] = -log(1.-Rand::uniform());

            }

            // Discontinuous emission: emission under progress
            // If epsilon_tau_ > 0
            if (tau[ipart] > epsilon_tau_)
            {

                // from the cross section
                temp = RadiationTables.compute_dNphdt(chipa,gamma);

                // Time to discontinuous emission
                // If this time is > the remaining iteration time,
                // we have a synchronization
                emission_time = std::min(tau[ipart]/temp, dt - local_it_time);

                // Update of the optical depth
                tau[ipart] -= temp*emission_time;

                // If the final optical depth is reached
                if (tau[ipart] <= epsilon_tau_)
                {

                    // Emission of a photon
                    RadiationMonteCarlo::photon_emission(ipart,
                                           chipa,gamma,
                                           position,
                                           momentum,
                                           weight,
                                           photon_species,
                                           RadiationTables);

                    // Optical depth becomes negative meaning
                    // that a new drawing is possible
                    // at the next Monte-Carlo iteration
                    tau[ipart] = -1.;
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
            // tau[ipart] <= epsilon_tau_
            else if ((chipa <= RadiationTables.get_chipa_disc_min_threshold())
            &&  (tau[ipart] <= epsilon_tau_)
            &&  (chipa > RadiationTables.get_chipa_radiation_threshold())
            &&  (gamma > 1.))
            {

                // Remaining time of the iteration
                emission_time = dt - local_it_time;

                // Radiated energy during emission_time
                cont_rad_energy =
                RadiationTables.get_corrected_cont_rad_energy_Ridgers(chipa,
                                                             emission_time);

                // Effect on the momentum
                temp = cont_rad_energy*gamma/(gamma*gamma-1.);
                for ( int i = 0 ; i<3 ; i++ )
                {
                    momentum[i][ipart] -= temp*momentum[i][ipart];
                }

                // Incrementation of the radiated energy cumulative parameter
                radiated_energy += weight[ipart]*(gamma - sqrt(1.0
                                    + momentum[0][ipart]*momentum[0][ipart]
                                    + momentum[1][ipart]*momentum[1][ipart]
                                    + momentum[2][ipart]*momentum[2][ipart]));

                // End for this particle
                local_it_time = dt;
            }
            // No emission since chipa is too low
            else // if (chipa < RadiationTables.get_chipa_radiation_threshold())
            {
                local_it_time = dt;
            }

        }

    }

}

// ---------------------------------------------------------------------------------------------------------------------
//! Perform the photon emission (creation of a super-photon
//! and slow down of the emitting particle)
//! \param ipart              particle index
//! \param chipa              particle quantum parameter
//! \param gammapa            particle gamma factor
//! \param position           particle position
//! \param momentum           particle momentum
//! \param RadiationTables    Cross-section data tables and useful functions
//                        for nonlinear inverse Compton scattering
// ---------------------------------------------------------------------------------------------------------------------
void RadiationMonteCarlo::photon_emission(int ipart,
                            double &chipa,
                            double & gammapa,
                            double * position[3],
                            double * momentum[3],
                            double * weight,
                            Species * photon_species,
                            RadiationTables &RadiationTables)
{
    // ____________________________________________________
    // Parameters
    double chiph;      // Photon quantum parameter
    double gammaph;    // Photon gamma factor
    double inv_old_norm_p;
    //double new_norm_p;

    // Get the photon quantum parameter from the table xip
    chiph = RadiationTables.compute_chiph_emission(chipa);

    // compute the photon gamma factor
    gammaph = chiph/chipa*(gammapa-1.0);

    // ____________________________________________________
    // Creation of the new photon

    // ____________________________________________________
    // Update of the particle properties
    // direction d'emission // direction de l'electron (1/gamma << 1)
    // With momentum conservation
    inv_old_norm_p = gammaph/sqrt(gammapa*gammapa - 1.0);
    momentum[0][ipart] -= momentum[0][ipart]*inv_old_norm_p;
    momentum[1][ipart] -= momentum[1][ipart]*inv_old_norm_p;
    momentum[2][ipart] -= momentum[2][ipart]*inv_old_norm_p;

    // With energy conservation
    /*inv_old_norm_p = 1./sqrt(gammapa*gammapa - 1.0);
    gammapa -= gammaph;
    new_norm_p = sqrt(gammapa*gammapa - 1.0);
    px *= new_norm_p * inv_old_norm_p;
    py *= new_norm_p * inv_old_norm_p;
    pz *= new_norm_p * inv_old_norm_p;*/

    // Creation of macro-photons if requested
    // Check that the photon_species is defined and the threshold on the energy
    if (photon_species
        && (gammaph >= radiation_photon_gamma_threshold_))
    {
        /* ---------------------------------------------------------------------
        // First method: emission of a single photon

        // Creation of the new photon in the temporary array new_photons
        new_photons.create_particle();

        int idNew = new_photons.size() - 1;

        for (int i=0; i<nDim_; i++) {
            new_photons.position(i,idNew)=position[i][ipart];
        }

        inv_old_norm_p = 1./sqrt(momentum[0][ipart]*momentum[0][ipart]
                                + momentum[1][ipart]*momentum[1][ipart]
                                + momentum[2][ipart]*momentum[2][ipart]);

        for (unsigned int i=0; i<3; i++) {
            new_photons.momentum(i,idNew) =
            gammaph*momentum[i][ipart]*inv_old_norm_p;
        }

        new_photons.weight(idNew)=weight[ipart];
        new_photons.charge(idNew)=0;
        --------------------------------------------------------------------- */

        // Second method: emission of several photons for statistics following
        // the parameter radiation_photon_sampling_

        // Creation of new photons in the temporary array new_photons
        new_photons.create_particles(radiation_photon_sampling_);

        // Final size
        int npart = new_photons.size();

        // Inverse of the momentum norm
        inv_old_norm_p = 1./sqrt(momentum[0][ipart]*momentum[0][ipart]
                                + momentum[1][ipart]*momentum[1][ipart]
                                + momentum[2][ipart]*momentum[2][ipart]);

        // For all new photons...
        for (int idNew=npart-radiation_photon_sampling_; idNew<npart; idNew++)
        {
            for (int i=0; i<nDim_; i++) {
                new_photons.position(i,idNew)=position[i][ipart];
            }

            for (int i=0; i<3; i++) {
                new_photons.momentum(i,idNew) =
                gammaph*momentum[i][ipart]*inv_old_norm_p;
            }

            new_photons.weight(idNew)=weight[ipart]*inv_radiation_photon_sampling_;
            new_photons.charge(idNew)=0;

            if (new_photons.isQuantumParameter)
            {
                new_photons.chi(idNew) = chiph;
            }

            if (new_photons.isMonteCarlo)
            {
                new_photons.tau(idNew) = -1.;
            }

        }

    }
    // Addition of the emitted energy in the cumulating parameter
    // for the scalar diagnostics
    else
    {

        gammaph = gammapa - sqrt(1.0 + momentum[0][ipart]*momentum[0][ipart]
                                     + momentum[1][ipart]*momentum[1][ipart]
                                     + momentum[2][ipart]*momentum[2][ipart]);
        radiated_energy += weight[ipart]*gammaph;
    }
}
