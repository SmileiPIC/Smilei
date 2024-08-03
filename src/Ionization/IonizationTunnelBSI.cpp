#include "IonizationTunnelBSI.h"

#include <cmath>
#include <vector>

#include "IonizationTables.h"
#include "Particles.h"
#include "Species.h"
using namespace std;

IonizationTunnelBSI::IonizationTunnelBSI(Params &params, Species *species) : Ionization(params, species)
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

void IonizationTunnelBSI::operator()(Particles *particles, unsigned int ipart_min, unsigned int ipart_max,
                                     vector<double> *Epart, Patch *patch, Projector *Proj, int ipart_ref)
{
    unsigned int Z;
    electricFields E;
    vector<double> IonizRate_tunnel(atomic_number_),
        Dnom_tunnel(atomic_number_);  // For MonteCarlo procedure using the tunnel
                                      // ioniz. rate.
    vector<double> IonizRate_BSI_linear(atomic_number_),
        Dnom_BSI_linear(atomic_number_);  // For MonteCarlo procedure using the BSI
                                          // ioniz. rate.
    vector<double> IonizRate_BSI_quadratic(atomic_number_), Dnom_BSI_quadratic(atomic_number_);

    double factorJion_0 = au_to_mec2 * EC_to_au * EC_to_au * invdt;  // Will be used to calc. ionization current.

    int nparts = Epart->size() / 3;  // Over 3 because there are 3 Dims. Epart is like [part1_Ex, part2_Ex,
                                     // ... partnparts_Ex, part1_Ey, part2_Ey,...partnparts_Ey, part1_Ez,
                                     // part2_Ez, partnparts_Ez]}
    E.x = &((*Epart)[0 * nparts]);
    E.y = &((*Epart)[1 * nparts]);
    E.z = &((*Epart)[2 * nparts]);

    for (unsigned int ipart = ipart_min; ipart < ipart_max; ipart++)
    {
        Z = (unsigned int)(particles->charge(ipart));  // current charge state of the particle number ipart of species
                                                       // called ''species'' with atomic no = atomic_number_
        double E_cr = pow(Potential[Z], 2) / (4 * (Z + 1));  //  Formula from DonnaStrickland, 1991, page5, eq(5) [Laser
                                                             //  ionization of noble gases by Coulomb-barrier suppression]
        // There's absolutely no usage of this E_cr in the code.

        // If ion already fully ionized then skip and go to next particle
        if (Z == atomic_number_)
        {
            continue;
        }  // atom fully ionized already, no electrons remained to be stripped away.

        // Absolute value of the electric field in AtomicUnits (i.e. normalized to
        // (in units of) E_atomic = 5.1422*10^11 V/m)
        E.abs = EC_to_au * sqrt(pow(*(E.x + ipart - ipart_ref),
                                    2)  // (Ex+ipart-ipart_ref) points to the location in the
                                        // container Epart at which E-field_x for particle ipart
                                        // sits.
                                + pow(*(E.y + ipart - ipart_ref),
                                      2)  // Similarly for y. Dereferencing it via *() means we
                                          // access and use the value sitting at that location.
                                + pow(*(E.z + ipart - ipart_ref),
                                      2));  // EC_to_au transforms from SMILEI-normalized units to
                                            // AtomicUnits.
        if (E.abs < 1e-10)
        {  // E is in Atomic Units
            continue;
        }

        // Used to be call to continuity_tool. After some shortening, the call was manually inlined
        double BSI_rate_quadratic = ionizationRateBSIQuadratic(Z, E);
        double BSI_rate_linear = ionizationRateBSILinear(Z, E);
        double Tunnel_rate = ionizationRateTunnel(Z, E);

        if (BSI_rate_quadratic >= BSI_rate_linear)
        {
            tunnelMonteCarloRoutine(particles, ipart, Epart, patch, Proj, Z, E, IonizRate_BSI_quadratic,
                                    Dnom_BSI_quadratic, &IonizationTunnelBSI::ionizationRateBSIQuadratic);
        }
        else if (std::min(Tunnel_rate, BSI_rate_quadratic) == BSI_rate_quadratic)
        {
            tunnelMonteCarloRoutine(particles, ipart, Epart, patch, Proj, Z, E, IonizRate_BSI_linear, Dnom_BSI_linear,
                                    &IonizationTunnelBSI::ionizationRateBSILinear);
        }
        else
        {
            tunnelMonteCarloRoutine(particles, ipart, Epart, patch, Proj, Z, E, IonizRate_tunnel, Dnom_tunnel,
                                    &IonizationTunnelBSI::ionizationRateTunnel);
        }
    }  // END loop on particles

}  // void IonizationTunnelBSI::operator()(arg1, arg2, ...) scope end.

double IonizationTunnelBSI::ionizationRateBSIQuadratic(const int Z, const electricFields E)  // Z is charge state
{
    double ratio_of_IPs = IH / IonizationTables::ionization_energy(atomic_number_, Z);
    return au_to_w0 * (2.4 * (pow(E.abs, 2)) * pow(ratio_of_IPs, 2));
}

double IonizationTunnelBSI::ionizationRateBSILinear(const int Z, const electricFields E)  // Z is charge state
{
    double ratio_of_IPs = IH / IonizationTables::ionization_energy(atomic_number_, Z);
    return au_to_w0 * (0.8 * E.abs * pow(ratio_of_IPs, 0.5));
}

double IonizationTunnelBSI::ionizationRateTunnel(const int Z, const electricFields E)
{
    double delta = gamma_tunnel[Z] * E.inv;
    return beta_tunnel[Z] * exp(-delta * one_third + alpha_tunnel[Z] * log(delta));
}

