#ifndef IONIZATIONTUNNEL_H
#define IONIZATIONTUNNEL_H

#include <cmath>
#include <functional>
#include <vector>

#include "Ionization.h"
#include "Particles.h"
#include "Species.h"
#include "Tools.h"

class Particles;

template <int Model>
class IonizationTunnel : public Ionization
{
   public:
    inline IonizationTunnel(Params &params, Species *species);

    inline void operator()(Particles *, unsigned int, unsigned int, std::vector<double> *, Patch *, Projector *,
                           int ipart_ref = 0) override;

   private:
    template <int place>
    inline double ionizationRate(const int Z, const double E);

    double one_third;
    unsigned int atomic_number_;
    std::vector<double> Potential, Azimuthal_quantum_number;
    std::vector<double> alpha_tunnel, beta_tunnel, gamma_tunnel;

    // To be conditionally prepared
    // FullPPT
    std::vector<double> Magnetic_quantum_number;

    // Tong&Lin
    double ionization_tl_parameter_;
    std::vector<double> lambda_tunnel;

    // BSI
    const double IH = 13.598434005136;
    int rate_formula;
};

template <int Model>
inline void IonizationTunnel<Model>::operator()(Particles *particles, unsigned int ipart_min, unsigned int ipart_max,
                                                vector<double> *Epart, Patch *patch, Projector *Proj, int ipart_ref)
{
    unsigned int Z, Zp1, newZ, k_times;
    double TotalIonizPot, E, invE, factorJion, ran_p, Mult, D_sum, P_sum, Pint_tunnel;
    vector<double> IonizRate_tunnel(atomic_number_), Dnom_tunnel(atomic_number_);
    LocalFields Jion;
    double factorJion_0 = au_to_mec2 * EC_to_au * EC_to_au * invdt;

    int nparts = Epart->size() / 3;
    double *Ex = &((*Epart)[0 * nparts]);
    double *Ey = &((*Epart)[1 * nparts]);
    double *Ez = &((*Epart)[2 * nparts]);

    for (unsigned int ipart = ipart_min; ipart < ipart_max; ipart++) {
        // Current charge state of the ion
        Z = (unsigned int)(particles->charge(ipart));

        // If ion already fully ionized then skip
        if (Z == atomic_number_) {
            continue;
        }

        // Absolute value of the electric field normalized in atomic units
        E = EC_to_au * sqrt(*(Ex + ipart - ipart_ref) * *(Ex + ipart - ipart_ref) +
                            *(Ey + ipart - ipart_ref) * *(Ey + ipart - ipart_ref) +
                            *(Ez + ipart - ipart_ref) * *(Ez + ipart - ipart_ref));
        if (E < 1e-10) {
            continue;
        }

        // --------------------------------
        // Start of the Monte-Carlo routine
        // --------------------------------

        invE = 1. / E;
        factorJion = factorJion_0 * invE * invE;
        ran_p = patch->rand_->uniform();
        IonizRate_tunnel[Z] = ionizationRate<1>(Z, E);

        // Total ionization potential (used to compute the ionization current)
        TotalIonizPot = 0.0;

        // k_times will give the nb of ionization events
        k_times = 0;
        Zp1 = Z + 1;

        if (Zp1 == atomic_number_) {
            // if ionization of the last electron: single ionization
            // -----------------------------------------------------
            if (ran_p < 1.0 - exp(-IonizRate_tunnel[Z] * dt)) {
                TotalIonizPot += Potential[Z];
                k_times = 1;
            }

        } else {
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
            while ((Pint_tunnel < ran_p) and (k_times < atomic_number_ - Zp1)) {
                newZ = Zp1 + k_times;
                IonizRate_tunnel[newZ] = ionizationRate<2>(newZ, E);
                D_sum = 0.0;
                P_sum = 0.0;
                Mult *= IonizRate_tunnel[Z + k_times];
                for (unsigned int i = 0; i < k_times + 1; i++) {
                    Dnom_tunnel[i] = Dnom_tunnel[i] / (IonizRate_tunnel[newZ] - IonizRate_tunnel[Z + i]);
                    D_sum += Dnom_tunnel[i];
                    P_sum += exp(-IonizRate_tunnel[Z + i] * dt) * Dnom_tunnel[i];
                }
                Dnom_tunnel[k_times + 1] = -D_sum;
                P_sum = P_sum + Dnom_tunnel[k_times + 1] * exp(-IonizRate_tunnel[newZ] * dt);
                Pint_tunnel = Pint_tunnel + P_sum * Mult;

                TotalIonizPot += Potential[Z + k_times];
                k_times++;
            }  // END while

            // final ionization (of last electron)
            if (((1.0 - Pint_tunnel) > ran_p) && (k_times == atomic_number_ - Zp1)) {
                TotalIonizPot += Potential[atomic_number_ - 1];
                k_times++;
            }
        }  // END Multiple ionization routine

        // Compute ionization current
        if (patch->EMfields->Jx_ != NULL) {  // For the moment ionization current is
                                             // not accounted for in AM geometry
            factorJion *= TotalIonizPot;
            Jion.x = factorJion * *(Ex + ipart);
            Jion.y = factorJion * *(Ey + ipart);
            Jion.z = factorJion * *(Ez + ipart);

            Proj->ionizationCurrents(patch->EMfields->Jx_, patch->EMfields->Jy_, patch->EMfields->Jz_, *particles, ipart, Jion);
        }

        // Creation of the new electrons
        // (variable weights are used)
        // -----------------------------

        if (k_times != 0) {
            new_electrons.createParticle();
            int idNew = new_electrons.size() - 1;
            for (unsigned int i = 0; i < new_electrons.dimension(); i++) {
                new_electrons.position(i, idNew) = particles->position(i, ipart);
            }
            for (unsigned int i = 0; i < 3; i++) {
                new_electrons.momentum(i, idNew) = particles->momentum(i, ipart) * ionized_species_invmass;
            }
            new_electrons.weight(idNew) = double(k_times) * particles->weight(ipart);
            new_electrons.charge(idNew) = -1;

            if (save_ion_charge_) {
                ion_charge_.push_back(particles->charge(ipart));
            }

            // Increase the charge of the particle
            particles->charge(ipart) += k_times;
        }

    }  // Loop on particles
}

template <int Model>
template <int place>
inline double IonizationTunnel<Model>::ionizationRate(const int Z, const double E)
{
    double delta = gamma_tunnel[Z] / E;
    return beta_tunnel[Z] * exp(-delta * one_third + alpha_tunnel[Z] * log(delta));
}

// IonizationTunnel : 0
template <>
IonizationTunnel<0>::IonizationTunnel(Params &params, Species *species);

// IonizationTunnelFullPPT: 1
template <>
IonizationTunnel<1>::IonizationTunnel(Params &params, Species *species);

// Tong&Ling: 2
template <>
IonizationTunnel<2>::IonizationTunnel(Params &params, Species *species);

template <>
template <int place>
inline double IonizationTunnel<2>::ionizationRate(const int Z, const double E);

// BSI: 3
template <>
IonizationTunnel<3>::IonizationTunnel(Params &params, Species *species);

template <>
template <>
double IonizationTunnel<3>::ionizationRate<1>(const int Z, const double E);

template <>
template <>
double IonizationTunnel<3>::ionizationRate<2>(const int Z, const double E);

#endif
