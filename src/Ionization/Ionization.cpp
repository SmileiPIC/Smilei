#include "Ionization.h"

#include <vector>

#include "Species.h"

using namespace std;

Ionization::Ionization(Params &params, Species *species)
{
    reference_angular_frequency_SI = params.reference_angular_frequency_SI;

    dt = params.timestep;
    invdt = 1. / dt;
    nDim_field = params.nDim_field;
    nDim_particle = params.nDim_particle;
    ionized_species_invmass = 1. / species->mass_;

    // Normalization constant from Smilei normalization to/from atomic units
    eV_to_au = 1.0 / 27.2116;
    au_to_mec2 = 27.2116 / 510.998e3;
    EC_to_au = 3.314742578e-15 * reference_angular_frequency_SI;  // hbar omega / (me c^2 alpha^3)
    au_to_w0 = 4.134137172e+16 / reference_angular_frequency_SI;  // alpha^2 me c^2 / (hbar omega)

    atomic_number_ = species->atomic_number_;
    Potential.resize(atomic_number_);
    Azimuthal_quantum_number.resize(atomic_number_);

    one_third = 1.0 / 3.0;

    alpha_tunnel.resize(atomic_number_);
    beta_tunnel.resize(atomic_number_);
    gamma_tunnel.resize(atomic_number_);

#ifdef _OMPTASKS
    new_electrons_per_bin = new Particles[species->Nbins];
    ion_charge_per_bin_.resize(species->Nbins);
#endif
}

Ionization::~Ionization() {}

void Ionization::joinNewElectrons(unsigned int Nbins)
{
    // if tasks on bins are used for ionization, join the lists of new electrons
    // created in each bin, to have the list of new electrons for this species and patch

    size_t start = new_electrons.size();

    // Resize new_electrons
    size_t total_n_new = 0;
    for (size_t ibin = 0; ibin < Nbins; ibin++)
    {
        total_n_new += new_electrons_per_bin[ibin].size();
    }
    new_electrons.createParticles(total_n_new);
    // Also resize ion_charge_ if necessary
    if (save_ion_charge_)
    {
        ion_charge_.resize(start + total_n_new);
    }

    // Move each new_electrons_per_bin into new_electrons
    for (size_t ibin = 0; ibin < Nbins; ibin++)
    {
        size_t n_new = new_electrons_per_bin[ibin].size();
        for (size_t i = 0; i < new_electrons.dimension(); i++)
        {
            copy(&new_electrons_per_bin[ibin].position(i, 0), &new_electrons_per_bin[ibin].position(i, n_new),
                 &new_electrons.position(i, start));
        }
        for (size_t i = 0; i < new_electrons.dimension(); i++)
        {
            copy(&new_electrons_per_bin[ibin].momentum(i, 0), &new_electrons_per_bin[ibin].momentum(i, n_new),
                 &new_electrons.momentum(i, start));
        }
        copy(&new_electrons_per_bin[ibin].weight(0), &new_electrons_per_bin[ibin].weight(n_new), &new_electrons.weight(start));
        copy(&new_electrons_per_bin[ibin].charge(0), &new_electrons_per_bin[ibin].charge(n_new), &new_electrons.charge(start));
        new_electrons_per_bin[ibin].clear();
        if (save_ion_charge_)
        {
            copy(ion_charge_per_bin_[ibin].begin(), ion_charge_per_bin_[ibin].end(), ion_charge_.begin() + start);
            ion_charge_per_bin_[ibin].clear();
        }
        start += n_new;
    }
}

void Ionization::tunnelMonteCarloRoutine(Particles *particles, unsigned int ipart, vector<double> *Epart,
                                                      Patch *patch, Projector *Proj, const unsigned int Z,
                                                      const electricFields E, vector<double> &IonizRate_tunnel,
                                                      vector<double> &Dnom_tunnel,
                                                      function<double(const int, const electricFields)> ionizationRate)
{
    double TotalIonizPot, factorJion, ran_p, Mult, D_sum, P_sum, Pint_tunnel;
    LocalFields Jion;
    const double factorJion_0 = au_to_mec2 * EC_to_au * EC_to_au * invdt;
    unsigned int newZ, Zp1, k_times;
    factorJion = factorJion_0 * E.inv * E.inv;
    ran_p = patch->rand_->uniform();
    IonizRate_tunnel[Z] = ionizationRate(Z, E);

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
        //        partial & final ionization are decoupled (see Nuter Phys.
        //        Plasmas)
        // -------------------------------------------------------------------------

        // initialization
        Mult = 1.0;
        Dnom_tunnel[0] = 1.0;
        Pint_tunnel = exp(-IonizRate_tunnel[Z] * dt);  // cummulative prob.

        // multiple ionization loop while Pint_tunnel < ran_p and still partial
        // ionization
        while ((Pint_tunnel < ran_p) and (k_times < atomic_number_ - Zp1))
        {
            newZ = Zp1 + k_times;
            IonizRate_tunnel[newZ] = ionizationRate(newZ, E);
            D_sum = 0.0;
            P_sum = 0.0;
            Mult *= IonizRate_tunnel[Z + k_times];
            for (unsigned int i = 0; i < k_times + 1; i++)
            {
                Dnom_tunnel[i] = Dnom_tunnel[i] / (IonizRate_tunnel[newZ] - IonizRate_tunnel[Z + i]);
                D_sum += Dnom_tunnel[i];
                P_sum += exp(-IonizRate_tunnel[Z + i] * dt) * Dnom_tunnel[i];
            }
            Dnom_tunnel[k_times + 1] -= D_sum;
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
    if (patch->EMfields->Jx_ != NULL)
    {  // For the moment ionization current is
       // not accounted for in AM geometry
        factorJion *= TotalIonizPot;
        Jion.x = factorJion * *(E.x + ipart);
        Jion.y = factorJion * *(E.y + ipart);
        Jion.z = factorJion * *(E.z + ipart);

        Proj->ionizationCurrents(patch->EMfields->Jx_, patch->EMfields->Jy_, patch->EMfields->Jz_, *particles, ipart, Jion);
    }

    // Creation of the new electrons
    // (variable weights are used)
    // -----------------------------

    if (k_times != 0)
    {
        new_electrons.createParticle();
        int idNew = new_electrons.size() - 1;
        for (unsigned int i = 0; i < new_electrons.dimension(); i++)
        {
            new_electrons.position(i, idNew) = particles->position(i, ipart);
        }
        for (unsigned int i = 0; i < 3; i++)
        {
            new_electrons.momentum(i, idNew) = particles->momentum(i, ipart) * ionized_species_invmass;
        }
        new_electrons.weight(idNew) = double(k_times) * particles->weight(ipart);
        new_electrons.charge(idNew) = -1;

        if (save_ion_charge_)
        {
            ion_charge_.push_back(particles->charge(ipart));
        }

        // Increase the charge of the particle
        particles->charge(ipart) += k_times;
    }
}

