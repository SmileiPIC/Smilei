#include "IonizationTunnel.h"

#include <cmath>
#include <functional>

#include "IonizationTables.h"
#include "Particles.h"
#include "Species.h"

using namespace std;

IonizationTunnel::IonizationTunnel(Params &params, Species *species) : Ionization(params, species)
{
    DEBUG("Creating the Tunnel Ionizaton class");

    // Ionization potential & quantum numbers (all in atomic units 1 au = 27.2116 eV)
    for (unsigned int Z = 0; Z < atomic_number_; Z++)
    {
        DEBUG("Z : " << Z);

        Potential[Z] = IonizationTables::ionization_energy(atomic_number_, Z) * eV_to_au;
        Azimuthal_quantum_number[Z] = IonizationTables::azimuthal_atomic_number(atomic_number_, Z);

        DEBUG("Potential: " << Potential[Z] << " Az.q.num: " << Azimuthal_quantum_number[Z]);

        double cst = ((double)Z + 1.0) * sqrt(2.0 / Potential[Z]);
        alpha_tunnel[Z] = cst - 1.0;
        beta_tunnel[Z] = pow(2, alpha_tunnel[Z]) * (8. * Azimuthal_quantum_number[Z] + 4.0) / (cst * tgamma(cst)) *
                         Potential[Z] * au_to_w0;
        gamma_tunnel[Z] = 2.0 * pow(2.0 * Potential[Z], 1.5);
    }

    DEBUG("Finished Creating the Tunnel Ionizaton class");
}

void IonizationTunnel::operator()(Particles *particles, unsigned int ipart_min, unsigned int ipart_max,
                                  vector<double> *Epart, Patch *patch, Projector *Proj, int ipart_ref)
{
    unsigned int Z;
    electricFields E;
    vector<double> IonizRate_tunnel(atomic_number_), Dnom_tunnel(atomic_number_, 0.);

    int nparts = Epart->size() / 3;
    E.x = &((*Epart)[0 * nparts]);
    E.y = &((*Epart)[1 * nparts]);
    E.z = &((*Epart)[2 * nparts]);

    for (unsigned int ipart = ipart_min; ipart < ipart_max; ipart++)
    {
        // Current charge state of the ion
        Z = (unsigned int)(particles->charge(ipart));

        // If ion already fully ionized then skip
        if (Z == atomic_number_)
        {
            continue;
        }

        // Absolute value of the electric field normalized in atomic units
        E.abs = EC_to_au * sqrt(pow(*(E.x + ipart - ipart_ref), 2) + pow(*(E.y + ipart - ipart_ref), 2) +
                                pow(*(E.z + ipart - ipart_ref), 2));
        if (E.abs < 1e-10)
        {
            continue;
        }
        E.inv = 1. / E.abs;

        tunnelMonteCarloRoutine(particles, ipart, Epart, patch, Proj, Z, E, IonizRate_tunnel, Dnom_tunnel, &IonizationTunnel::ionizationRateTunnel);
    }  // Loop on particles
}

double IonizationTunnel::ionizationRateTunnel(const int Z, const electricFields E)
{
    double delta = gamma_tunnel[Z] * E.inv;
    return beta_tunnel[Z] * exp(-delta * one_third + alpha_tunnel[Z] * log(delta));
}

void IonizationTunnel::ionizationTunnelWithTasks(Particles *particles, unsigned int ipart_min, unsigned int ipart_max,
                                                 vector<double> *Epart, Patch *patch, Projector *Proj, int ibin,
                                                 int bin_shift, double *b_Jx, double *b_Jy, double *b_Jz, int ipart_ref)
{
    unsigned int Z, Zp1, newZ, k_times;
    double TotalIonizPot, E, invE, factorJion, delta, ran_p, Mult, D_sum, P_sum, Pint_tunnel;
    vector<double> IonizRate_tunnel(atomic_number_), Dnom_tunnel(atomic_number_);
    LocalFields Jion;
    double factorJion_0 = au_to_mec2 * EC_to_au * EC_to_au * invdt;

    int nparts = Epart->size() / 3;
    double *Ex = &((*Epart)[0 * nparts]);
    double *Ey = &((*Epart)[1 * nparts]);
    double *Ez = &((*Epart)[2 * nparts]);

    for (unsigned int ipart = ipart_min; ipart < ipart_max; ipart++)
    {
        // Current charge state of the ion
        Z = (unsigned int)(particles->charge(ipart));

        // If ion already fully ionized then skip
        if (Z == atomic_number_)
        {
            continue;
        }

        // Absolute value of the electric field normalized in atomic units
        E = EC_to_au * sqrt(pow(*(Ex + ipart - ipart_ref), 2) + pow(*(Ey + ipart - ipart_ref), 2) +
                            pow(*(Ez + ipart - ipart_ref), 2));
        if (E < 1e-10)
        {
            continue;
        }

        // --------------------------------
        // Start of the Monte-Carlo routine
        // --------------------------------

        invE = 1. / E;
        factorJion = factorJion_0 * invE * invE;
        delta = gamma_tunnel[Z] * invE;
        ran_p = patch->rand_->uniform();
        IonizRate_tunnel[Z] = beta_tunnel[Z] * exp(-delta * one_third + alpha_tunnel[Z] * log(delta));

        // Total ionization potential (used to compute the ionization current)
        TotalIonizPot = 0.0;

        // k_times will give the nb of ionization events
        k_times = 0;
        Zp1 = Z + 1;

        if (Zp1 == atomic_number_)
        {
            // if ionization of the last electron: single ionization
            // -----------------------------------------------------
            if (ran_p < 1.0 - exp(-IonizRate_tunnel[Z] * dt))
            {
                TotalIonizPot += Potential[Z];
                k_times = 1;
            }
        }
        else
        {
            // else : multiple ionization can occur in one time-step
            //        partial & final ionization are decoupled (see Nuter Phys. Plasmas)
            // -------------------------------------------------------------------------

            // initialization
            Mult = 1.0;
            Dnom_tunnel[0] = 1.0;
            Pint_tunnel = exp(-IonizRate_tunnel[Z] * dt);  // cummulative prob.

            // multiple ionization loop while Pint_tunnel < ran_p and still partial ionization
            while ((Pint_tunnel < ran_p) and (k_times < atomic_number_ - Zp1))
            {
                newZ = Zp1 + k_times;
                delta = gamma_tunnel[newZ] * invE;
                IonizRate_tunnel[newZ] = beta_tunnel[newZ] * exp(-delta * one_third + alpha_tunnel[newZ] * log(delta));
                D_sum = 0.0;
                P_sum = 0.0;
                Mult *= IonizRate_tunnel[Z + k_times];
                for (unsigned int i = 0; i < k_times + 1; i++)
                {
                    Dnom_tunnel[i] = Dnom_tunnel[i] / (IonizRate_tunnel[newZ] - IonizRate_tunnel[Z + i]);
                    D_sum += Dnom_tunnel[i];
                    P_sum += exp(-IonizRate_tunnel[Z + i] * dt) * Dnom_tunnel[i];
                }
                Dnom_tunnel[k_times + 1] = -D_sum;  // bug fix
                P_sum = P_sum + Dnom_tunnel[k_times + 1] * exp(-IonizRate_tunnel[newZ] * dt);
                Pint_tunnel = Pint_tunnel + P_sum * Mult;

                TotalIonizPot += Potential[Z + k_times];
                k_times++;
            }  // END while

            // final ionization (of last electron)
            if (((1.0 - Pint_tunnel) > ran_p) && (k_times == atomic_number_ - Zp1))
            {
                TotalIonizPot += Potential[atomic_number_ - 1];
                k_times++;
            }
        }  // END Multiple ionization routine

        // Compute ionization current
        if (b_Jx != NULL)
        {  // For the moment ionization current is not accounted for in AM geometry
            factorJion *= TotalIonizPot;
            Jion.x = factorJion * *(Ex + ipart);
            Jion.y = factorJion * *(Ey + ipart);
            Jion.z = factorJion * *(Ez + ipart);

            Proj->ionizationCurrentsForTasks(b_Jx, b_Jy, b_Jz, *particles, ipart, Jion, bin_shift);
        }

        // Creation of the new electrons
        // (variable weights are used)
        // -----------------------------
        if (k_times != 0)
        {
            new_electrons_per_bin[ibin].createParticle();
            int idNew = new_electrons_per_bin[ibin].size() -
                        1;  // cout<<"ibin "<<ibin<<"size "<<new_electrons_per_bin[ibin].size()<<"capacity "<<new_electrons_per_bin[ibin].capacity()<<"\n"<<endl;
            for (unsigned int i = 0; i < new_electrons_per_bin[ibin].dimension(); i++)
            {
                new_electrons_per_bin[ibin].position(i, idNew) = particles->position(i, ipart);
            }
            for (unsigned int i = 0; i < 3; i++)
            {
                new_electrons_per_bin[ibin].momentum(i, idNew) = particles->momentum(i, ipart) * ionized_species_invmass;
            }
            new_electrons_per_bin[ibin].weight(idNew) = double(k_times) * particles->weight(ipart);
            new_electrons_per_bin[ibin].charge(idNew) = -1;

            if (save_ion_charge_)
            {
                ion_charge_per_bin_[ibin].push_back(particles->charge(ipart));
            }

            // // Increase the charge of the particle
            particles->charge(ipart) += k_times;
        }

    }  // Loop on particles
}
