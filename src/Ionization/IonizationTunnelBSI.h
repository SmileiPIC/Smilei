#ifndef IONIZATIONTUNNELBSI_H
#define IONIZATIONTUNNELBSI_H

#include <cmath>
#include <vector>

#include "Ionization.h"
#include "Tools.h"

class Particles;


class IonizationTunnelBSI : public Ionization {
    public:
        IonizationTunnelBSI(Params &params, Species *species);
        void operator()(Particles*, unsigned int, unsigned int, std::vector<double>*, Patch*, Projector*, int ipart_ref = 0) override;

    private:
        unsigned int atomic_number_;
        std::vector<double> Potential, Azimuthal_quantum_number;
        double one_third;
        std::vector<double> alpha_tunnel, beta_tunnel, gamma_tunnel; 

        double BSI_rate_quadratic_function(unsigned int Zp1, double E, unsigned int atomic_number_, double au_to_w0);
        double BSI_rate_linear_function(unsigned int Zp1, double E, unsigned int atomic_number_, double au_to_w0);
        double Tunnel_rate_function(unsigned int Z, double E, double alpha, double beta, double gamma);

        int continuity_tool(unsigned int Zp1, double E, double alpha, double beta, double gamma, double E_cr, double Potential, unsigned int atomic_number_, double au_to_w0);

        double au_to_eV=27.2116;
 
};

#endif
