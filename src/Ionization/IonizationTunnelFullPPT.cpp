#include "IonizationTunnelFullPPT.h"

#include <cmath>

#include "IonizationTables.h"
#include "Particles.h"
#include "Species.h"

using namespace std;

IonizationTunnelFullPPT::IonizationTunnelFullPPT(Params &params, Species *species) : Ionization(params, species)
{
    DEBUG("Creating the Tunnel Ionizaton class");

    // Ionization potential & quantum numbers (all in atomic units 1 au = 27.2116 eV)
    Magnetic_quantum_number.resize(atomic_number_);

    for (unsigned int Z = 0; Z < atomic_number_; Z++)
    {
        DEBUG("Z : " << Z);

        Potential[Z] = IonizationTables::ionization_energy(atomic_number_, Z) * eV_to_au;
        Azimuthal_quantum_number[Z] = IonizationTables::azimuthal_atomic_number(atomic_number_, Z);
        Magnetic_quantum_number[Z] = IonizationTables::magnetic_atomic_number(atomic_number_, Z);

        DEBUG("Potential: " << Potential[Z] << " Az.q.num: " << Azimuthal_quantum_number[Z]
                            << " M.q.num: " << Magnetic_quantum_number[Z]);

        double cst = ((double)Z + 1.0) * sqrt(2.0 / Potential[Z]);
        double abs_m = abs(Magnetic_quantum_number[Z]);
        alpha_tunnel[Z] = cst - 1.0 - abs_m;
        beta_tunnel[Z] = pow(2, alpha_tunnel[Z]) * (8. * Azimuthal_quantum_number[Z] + 4.0) / (cst * tgamma(cst)) *
                         Potential[Z] * au_to_w0 * tgamma(Azimuthal_quantum_number[Z] + abs_m + 1) /
                         (tgamma(abs_m + 1) * tgamma(Azimuthal_quantum_number[Z] - abs_m + 1));
        gamma_tunnel[Z] = 2.0 * pow(2.0 * Potential[Z], 1.5);
    }

    DEBUG("Finished Creating the Tunnel Ionizaton class");
}

void IonizationTunnelFullPPT::operator()(Particles *particles, unsigned int ipart_min, unsigned int ipart_max,
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

        tunnelMonteCarloRoutine(particles, ipart, Epart, patch, Proj, Z, E, IonizRate_tunnel, Dnom_tunnel, &IonizationTunnelFullPPT::ionizationRateTunnel);
    }  // Loop on particles
}

double IonizationTunnelFullPPT::ionizationRateTunnel(const int Z, const electricFields E)
{
    double delta = gamma_tunnel[Z] * E.inv;
    return beta_tunnel[Z] * exp(-delta * one_third + alpha_tunnel[Z] * log(delta));
}

