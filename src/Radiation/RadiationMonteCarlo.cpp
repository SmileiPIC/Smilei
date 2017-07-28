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
    radiation_photon_sampling = species->radiation_photon_sampling;
    inv_radiation_photon_sampling = 1. / radiation_photon_sampling;
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
        int ithread)
{

    // _______________________________________________________________
    // Parameters
    std::vector<LocalFields> *Epart = &(smpi->dynamics_Epart[ithread]);
    std::vector<LocalFields> *Bpart = &(smpi->dynamics_Bpart[ithread]);
    //std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);

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
             &&(mc_it_nb < mc_it_nb_max))
        {

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

             /*if (mc_it_nb > 1) {
             std::cerr << "ipart: " << ipart
                       << " mc_it_nb: " << mc_it_nb
                       << " local_it_time: " << local_it_time << "/" << dt
                       << " gamma: " << gamma
                       << " px: " << momentum[0][ipart]
                       << " py: " << momentum[1][ipart]
                       << " chi: " << chipa << "(" << chipa_cont_threshold << ")"
                       << " Tau: " << tau[ipart] << "/" << epsilon_tau
                       << std::endl;
                       if (std::isnan(gamma))
                       {
                           ERROR("stop")
                       }
             }*/

            // Update the quantum parameter in species
            // chi[ipart] = chipa;

            // Discontinuous emission: New emission
            // If tau[ipart] <= 0, this is a new emission
            // We also check that chipa > chipa_threshold,
            // else chipa is too low to induce a discontinuous emission
            if ((chipa > RadiationTables.get_chipa_disc_min_threshold())
            && (tau[ipart] <= epsilon_tau) )
            {
                // New final optical depth to reach for emision
                while (tau[ipart] <= epsilon_tau)
                   tau[ipart] = -log(1.-Rand::uniform());

            }

            // Discontinuous emission: emission under progress
            // If epsilon_tau > 0
            if (tau[ipart] > epsilon_tau)
            {

                /*std::cerr << "Continue discontinuous emission - "
                        << "tau: " << tau[ipart]
                        << std::endl;*/

                // from the cross section
                temp = RadiationTables.compute_dNphdt(chipa,gamma);

                // Time to discontinuous emission
                // If this time is > the remaining iteration time,
                // we have a synchronization
                emission_time = std::min(tau[ipart]/temp, dt - local_it_time);

                // Update of the optical depth
                tau[ipart] -= temp*emission_time;

                /*if (ipart == 1)
                {
                    std::cerr << "tau: " << tau[ipart] << " "
                              << "temp: " << temp << " "
                              << "tau[ipart]/temp: " << tau[ipart]/temp << " "
                              << "emission_time: " << emission_time << " "
                              << "dt - local_it_time: " << dt - local_it_time << " "
                              << std::endl;
                }*/

                // If the final optical depth is reached
                if (tau[ipart] <= epsilon_tau)
                {

                    /*std::cerr << "Photon emission"
                            << std::endl;*/

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
            // tau[ipart] <= epsilon_tau
            else if ((chipa <= RadiationTables.get_chipa_disc_min_threshold())
            &&  (tau[ipart] <= epsilon_tau)
            &&  (chipa > chipa_cont_threshold))
            {

                /*std::cerr << "Continuous - "
                          << "chipa: " << chipa << std::endl;*/

                // Remaining time of the iteration
                emission_time = dt - local_it_time;

                // Radiated energy during emission_time
                cont_rad_energy =
                RadiationTables.get_corrected_cont_rad_energy_Ridgers(chipa,
                                                             emission_time);

                // Effect on the momentum
                temp = cont_rad_energy*gamma/(gamma*gamma-1.);
                for ( int i = 0 ; i<3 ; i++ )
                    momentum[i][ipart] -= temp*momentum[i][ipart];

                // Incrementation of the radiated energy cumulative parameter
                radiated_energy += weight[ipart]*(gamma - sqrt(1.0
                                    + momentum[0][ipart]*momentum[0][ipart]
                                    + momentum[1][ipart]*momentum[1][ipart]
                                    + momentum[2][ipart]*momentum[2][ipart]));

                // End for this particle
                local_it_time = dt;
            }
            // No emission since chipa is too low
            else if (chipa < chipa_cont_threshold)
            {
                local_it_time = dt;
            }

        }

        /*std::cerr << "end ipart: " << ipart << std::endl;*/

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
    if (photon_species)
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
        // the parameter radiation_photon_sampling

        // Creation of new photons in the temporary array new_photons
        new_photons.create_particles(radiation_photon_sampling);

        // Final size
        int npart = new_photons.size();

        // For all new photons...
        for (int idNew=npart-radiation_photon_sampling; idNew<npart; idNew++)
        {
            for (int i=0; i<nDim_; i++) {
                new_photons.position(i,idNew)=position[i][ipart];
            }

            inv_old_norm_p = 1./sqrt(momentum[0][ipart]*momentum[0][ipart]
                                    + momentum[1][ipart]*momentum[1][ipart]
                                    + momentum[2][ipart]*momentum[2][ipart]);

            for (int i=0; i<3; i++) {
                new_photons.momentum(i,idNew) =
                gammaph*momentum[i][ipart]*inv_old_norm_p;
            }

            new_photons.weight(idNew)=weight[ipart]*inv_radiation_photon_sampling;
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

    // Debugging
    /*std::cerr << "chipa: " << chipa << " "
              << "chiph: " << chiph << " "
              << "gammapa: " << gammapa << " "
              << "gammaph: " << gammaph << " "
              << "inv_old_norm_p: " << inv_old_norm_p << " "
              //<< "new_norm_p: " << new_norm_p << " "
              << "" << sqrt(1 + px*px + py*py + pz*pz) << " "
              << std::endl;*/

}
