#include "IonizationTunnelTL.h"

#include <cmath>

#include "IonizationTables.h"
#include "Particles.h"
#include "Species.h"

using namespace std;

IonizationTunnelTL::IonizationTunnelTL(Params &params, Species *species) : Ionization(params, species)
{
    DEBUG("Creating the Tong-Lin Tunnel Ionizaton class");

    ionization_tl_parameter_ =
        species->ionization_tl_parameter_;  // species->ionization_tl_parameter_ is double
                                            // Varies from 6 to 9. This is the alpha parameter in Tong-Lin exponential, see Eq. (6) in [M F Ciappina and S V Popruzhenko 2020 Laser Phys. Lett. 17 025301 2020].
    lambda_tunnel.resize(atomic_number_);

    // Ionization potential & quantum numbers (all in atomic units 1 au = 27.2116 eV)
    for (int Z = 0; Z < (int)atomic_number_; Z++)
    {
        DEBUG("Z : " << Z);
        Potential[Z] = IonizationTables::ionization_energy(atomic_number_, Z) * eV_to_au;
        Azimuthal_quantum_number[Z] = IonizationTables::azimuthal_atomic_number(atomic_number_, Z);

        DEBUG("potential: " << Potential[Z] << " Az.q.num: " << Azimuthal_quantum_number[Z]);
    
        double cst = ((double)Z + 1.0) * sqrt(2.0 / Potential[Z]);
        alpha_tunnel[Z] = cst - 1.0;
        beta_tunnel[Z] = pow(2, alpha_tunnel[Z]) * (8. * Azimuthal_quantum_number[Z] + 4.0) / (cst * tgamma(cst)) *
                         Potential[Z] * au_to_w0;
        gamma_tunnel[Z] = 2.0 * pow(2.0 * Potential[Z], 1.5);
        lambda_tunnel[Z] = ionization_tl_parameter_ * pow(cst, 2) / gamma_tunnel[Z];
    }

    DEBUG("Finished Creating the Tong-Lin Tunnel Ionizaton class");
}

void IonizationTunnelTL::operator()(Particles *particles, unsigned int ipart_min, unsigned int ipart_max,
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

        tunnelMonteCarloRoutine(particles, ipart, Epart, patch, Proj, Z, E, IonizRate_tunnel, Dnom_tunnel, &IonizationTunnelTL::ionizationRateTunnelTL);
    }  // Loop on particles
}

double IonizationTunnelTL::ionizationRateTunnelTL(int Z, double delta, double E)
{
    return beta_tunnel[Z] * exp(-delta * one_third + alpha_tunnel[Z] * log(delta) - E * lambda_tunnel[Z]);
}
