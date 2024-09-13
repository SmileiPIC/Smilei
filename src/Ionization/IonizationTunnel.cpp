#include "IonizationTunnel.h"
#include "IonizationTables.h"


// Classic: 1
template <>
IonizationTunnel<0>::IonizationTunnel(Params &params, Species *species) : Ionization(params, species)
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

//  IonizationFullPPT: 1
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

// Tong&Ling: 2
template <>
IonizationTunnel<2>::IonizationTunnel(Params &params, Species *species) : Ionization(params, species)
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

template <>
double IonizationTunnel<2>::ionizationRate(const int Z, const double E, const double oldZ)
{
    const double delta = gamma_tunnel[Z] / E;
    return beta_tunnel[Z] * exp(-delta * one_third + alpha_tunnel[Z] * log(delta) - E * lambda_tunnel[Z]);
}

// BSI: 3
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
double IonizationTunnel<3>::ionizationRate(const int Z, const double E, const double oldZ)
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

    double BSI_rate_quadratic = quadratic(oldZ, E);
    double BSI_rate_linear = linear(oldZ, E);
    double Tunnel_rate = normal(oldZ, E);

    if (BSI_rate_quadratic >= BSI_rate_linear)
    {
        return quadratic(Z, E);
    }
    else if (std::min(Tunnel_rate, BSI_rate_quadratic) == BSI_rate_quadratic)
    {
        return linear(Z, E);
    }
    else
    {
        return normal(Z, E);
    }

}




