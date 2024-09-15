#include "IonizationTunnel.h"
#include "IonizationTables.h"

template <>
IonizationTunnel<1>::IonizationTunnel(Params &params, Species *species) : Ionization(params, species)
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


