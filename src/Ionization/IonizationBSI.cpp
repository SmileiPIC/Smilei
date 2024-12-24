#include "IonizationTunnel.h"
#include "IonizationTables.h"

template <>
IonizationTunnel<3>::IonizationTunnel(Params &params, Species *species) : Ionization(params, species)
{
    DEBUG("Creating the Tunnel BSI Ionizaton class");

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
    }
    DEBUG("Finished Creating the BSI Tunnel Ionizaton class");
}

template <>
double IonizationTunnel<3>::ionizationRate(const int Z, const double E, const int oldZ)
{
    auto normal = [this](const int Z, const double E) -> double {
        double delta = gamma_tunnel[Z] / E;
        return beta_tunnel[Z] * exp(-delta * one_third + alpha_tunnel[Z] * log(delta));
    };

    auto linear = [this](const int Z, const double E) -> double {
        const double ratio_of_IPs = IH / IonizationTables::ionization_energy(atomic_number_, Z);
        return au_to_w0 * (0.8 * E * pow(ratio_of_IPs, 0.5));
    };
    auto quadratic = [this](const int Z, const double E) -> double {
        const double ratio_of_IPs = IH / IonizationTables::ionization_energy(atomic_number_, Z);
        return au_to_w0 * (2.4 * (pow(E, 2)) * pow(ratio_of_IPs, 2));
    };

    double BSI_rate_quadratic = quadratic(oldZ+1, E);
    double BSI_rate_linear = linear(oldZ+1, E);
    double Tunnel_rate = normal(oldZ, E);

    if (BSI_rate_quadratic >= BSI_rate_linear)
    {
        return linear(Z, E);
    }
    else if (std::min(Tunnel_rate, BSI_rate_quadratic) == BSI_rate_quadratic)
    {
        return quadratic(Z, E);
    }
    else
    {
        return normal(Z, E);
    }

}
