#include "IonizationTables.h"
#include "IonizationTunnel.h"

template <>
IonizationTunnel<3>::IonizationTunnel(Params &params, Species *species) : Ionization(params, species)
{
    DEBUG("Creating the Tunnel BSI Ionizaton class");

    atomic_number_ = species->atomic_number_;
    Potential.resize(atomic_number_);
    Azimuthal_quantum_number.resize(atomic_number_);

    one_third = 1.0 / 3.0;

    alpha_tunnel.resize(atomic_number_);
    beta_tunnel.resize(atomic_number_);
    gamma_tunnel.resize(atomic_number_);

    // Ionization potential & quantum numbers (all in atomic units 1 au = 27.2116 eV)
    for (int Z = 0; Z < (int)atomic_number_; Z++) {
        DEBUG("Z : " << Z);

        Potential[Z] = IonizationTables::ionization_energy(atomic_number_, Z) * eV_to_au;
        Azimuthal_quantum_number[Z] = IonizationTables::azimuthal_atomic_number(atomic_number_, Z);

        DEBUG("potential: " << Potential[Z] << " Az.q.num: " << Azimuthal_quantum_number[Z]);

        double cst = ((double)Z + 1.0) * sqrt(2.0 / Potential[Z]);
        alpha_tunnel[Z] = cst - 1.0;
        beta_tunnel[Z] = pow(2, alpha_tunnel[Z]) * (8. * Azimuthal_quantum_number[Z] + 4.0) / (cst * tgamma(cst)) *
                         Potential[Z] * au_to_w0;
        gamma_tunnel[Z] = 2.0 * pow(2.0 * Potential[Z], 1.5);
    }
    DEBUG("Finished Creating the BSI Tunnel Ionizaton class");
}

template <>
double IonizationTunnel<3>::ionizationRate1(const int Z, const double E)
{
    double ratio_of_IPs = IH / IonizationTables::ionization_energy(atomic_number_, Z);

    double BSI_rate_quadratic = 2.4 * (pow(E, 2))*pow(ratio_of_IPs, 2) * au_to_w0;
    double BSI_rate_linear = 0.8 * E * pow(ratio_of_IPs, 0.5) * au_to_w0;
    double delta = gamma_tunnel[Z] / E;
    double Tunnel_rate = beta_tunnel[Z] * exp(-delta / 3.0 + alpha_tunnel[Z] * log(delta));

    if (std::min(Tunnel_rate, BSI_rate_quadratic) == BSI_rate_quadratic) {
        rate_formula = 1;
        return BSI_rate_quadratic;
    } else if (BSI_rate_quadratic >= BSI_rate_linear) {
        rate_formula = 2;
        return BSI_rate_quadratic;
    } else {
        rate_formula = 0;
        return Tunnel_rate;
    }
}

template <>
double IonizationTunnel<3>::ionizationRate2(const int newZ, const double E)
{
    double ratio_of_IPs_newZ = IH / IonizationTables::ionization_energy(atomic_number_, newZ);
    double delta = gamma_tunnel[newZ] / E;
    if (rate_formula == 1) {
        return au_to_w0 * (2.4 * (pow(E, 2))*pow(ratio_of_IPs_newZ, 2));
    } else if (rate_formula == 2) {
        return au_to_w0 * (0.8 * E * pow(ratio_of_IPs_newZ, 0.5));
    } else {
        return beta_tunnel[newZ] * exp(-delta * one_third + alpha_tunnel[newZ] * log(delta));
    }
}
