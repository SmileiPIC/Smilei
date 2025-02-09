#ifndef IONIZATIONTUNNEL_H
#define IONIZATIONTUNNEL_H

#include <cmath>
#include <functional>
#include <vector>

#include "Ionization.h"
#include "IonizationTables.h"
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
    inline double ionizationRate(const int Z, const double E);

    static constexpr double one_third = 1. / 3.;
    unsigned int atomic_number_;
    std::vector<double> Potential, Azimuthal_quantum_number;
    std::vector<double> alpha_tunnel, beta_tunnel, gamma_tunnel;

    // Tong&Lin
    std::vector<double> lambda_tunnel;
};

template <int Model>
IonizationTunnel<Model>::IonizationTunnel(Params &params, Species *species) : Ionization(params, species)
{
    DEBUG("Creating the Tunnel Ionizaton class");
    double abs_m = 0;
    double ionization_tl_parameter_;

    // Ionization potential & quantum numbers (all in atomic units 1 au = 27.2116 eV)
    atomic_number_ = species->atomic_number_;
    Potential.resize(atomic_number_);
    Azimuthal_quantum_number.resize(atomic_number_);

    alpha_tunnel.resize(atomic_number_);
    beta_tunnel.resize(atomic_number_);
    gamma_tunnel.resize(atomic_number_);

    if (Model == 2) {
        ionization_tl_parameter_ = species->ionization_tl_parameter_;  // species->ionization_tl_parameter_ is
                                                                       // double Varies from 6 to 9. This is
                                                                       // the alpha parameter in Tong-Lin
                                                                       // exponential, see Eq. (6) in [M F
                                                                       // Ciappina and S V Popruzhenko 2020
                                                                       // Laser Phys. Lett. 17 025301 2020].
        lambda_tunnel.resize(atomic_number_);
    }

    for (unsigned int Z = 0; Z < atomic_number_; Z++) {
        DEBUG("Z : " << Z);

        if (Model == 1) {
            abs_m = abs(IonizationTables::magnetic_atomic_number(atomic_number_, Z));
        }

        Potential[Z] = IonizationTables::ionization_energy(atomic_number_, Z) * eV_to_au;
        Azimuthal_quantum_number[Z] = IonizationTables::azimuthal_atomic_number(atomic_number_, Z);

        DEBUG("Potential: " << Potential[Z] << " Az.q.num: " << Azimuthal_quantum_number[Z]);

        double cst = ((double)Z + 1.0) * sqrt(2.0 / Potential[Z]);
        alpha_tunnel[Z] = cst - 1.0 - abs_m;
        beta_tunnel[Z] = pow(2, alpha_tunnel[Z]) * (8. * Azimuthal_quantum_number[Z] + 4.0) / (cst * tgamma(cst)) *
                         Potential[Z] * au_to_w0 * tgamma(Azimuthal_quantum_number[Z] + abs_m + 1) /
                         (tgamma(abs_m + 1) * tgamma(Azimuthal_quantum_number[Z] - abs_m + 1));
        gamma_tunnel[Z] = 2.0 * sqrt(2.0 * Potential[Z] * 2.0 * Potential[Z] * 2.0 * Potential[Z]);
        if (Model == 2) {
            lambda_tunnel[Z] = ionization_tl_parameter_ * cst * cst / gamma_tunnel[Z];
        }
    }

    DEBUG("Finished Creating the Tunnel Ionizaton class");
}

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
        IonizRate_tunnel[Z] = ionizationRate(Z, E);

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
                IonizRate_tunnel[newZ] = ionizationRate(newZ, E);
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
inline double IonizationTunnel<Model>::ionizationRate(const int Z, const double E)
{
    double delta = gamma_tunnel[Z] / E;
    return beta_tunnel[Z] * exp(-delta * one_third + alpha_tunnel[Z] * log(delta));
}

// Tong&Ling: 2
template <>
inline double IonizationTunnel<2>::ionizationRate(const int Z, const double E)
{
    const double delta = gamma_tunnel[Z] / E;
    return beta_tunnel[Z] * exp(-delta * one_third + alpha_tunnel[Z] * log(delta) - E * lambda_tunnel[Z]);
}

// BSI: 3
template <>
inline double IonizationTunnel<3>::ionizationRate(const int Z, const double E)
{
    constexpr double IH = 13.598434005136;
    double ratio_of_IPs = IH / IonizationTables::ionization_energy(atomic_number_, Z);

    double BSI_rate_quadratic = 2.4 * (E * E) * ratio_of_IPs * ratio_of_IPs * au_to_w0;
    double BSI_rate_linear = 0.8 * E * sqrt(ratio_of_IPs) * au_to_w0;
    double delta = gamma_tunnel[Z] / E;
    double Tunnel_rate = beta_tunnel[Z] * exp(-delta / 3.0 + alpha_tunnel[Z] * log(delta));

    if (BSI_rate_quadratic >= BSI_rate_linear) {
        return BSI_rate_linear;
    } else if (std::min(Tunnel_rate, BSI_rate_quadratic) == BSI_rate_quadratic) {
        return BSI_rate_quadratic;
    } else {
        return Tunnel_rate;
    }
}

#endif
