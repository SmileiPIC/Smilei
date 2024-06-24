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

};

#endif