#include "IonizationTunnel.h"

#include "IonizationTables.h"

// Classic: 1
template <>
IonizationTunnel<0>::IonizationTunnel(Params &params, Species *species)
    : Ionization(params, species) {
  DEBUG("Creating the Tunnel Ionizaton class");

  atomic_number_ = species->atomic_number_;
  Potential.resize(atomic_number_);
  Azimuthal_quantum_number.resize(atomic_number_);

  one_third = 1.0 / 3.0;

  alpha_tunnel.resize(atomic_number_);
  beta_tunnel.resize(atomic_number_);
  gamma_tunnel.resize(atomic_number_);

  // Ionization potential & quantum numbers (all in atomic units 1 au = 27.2116
  // eV)
  for (unsigned int Z = 0; Z < atomic_number_; Z++) {
    DEBUG("Z : " << Z);

    Potential[Z] =
        IonizationTables::ionization_energy(atomic_number_, Z) * eV_to_au;
    Azimuthal_quantum_number[Z] =
        IonizationTables::azimuthal_atomic_number(atomic_number_, Z);

    DEBUG("Potential: " << Potential[Z]
                        << " Az.q.num: " << Azimuthal_quantum_number[Z]);

    double cst = ((double)Z + 1.0) * sqrt(2.0 / Potential[Z]);
    alpha_tunnel[Z] = cst - 1.0;
    beta_tunnel[Z] = pow(2, alpha_tunnel[Z]) *
                     (8. * Azimuthal_quantum_number[Z] + 4.0) /
                     (cst * tgamma(cst)) * Potential[Z] * au_to_w0;
    gamma_tunnel[Z] = 2.0 * sqrt(2.0 * Potential[Z] * 2.0 * Potential[Z] * 2.0 *
                                 Potential[Z]);
  }

  DEBUG("Finished Creating the Tunnel Ionizaton class");
}
