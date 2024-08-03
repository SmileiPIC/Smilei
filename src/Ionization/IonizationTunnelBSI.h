#ifndef IONIZATIONTUNNELBSI_H
#define IONIZATIONTUNNELBSI_H

#include <cmath>
#include <vector>

#include "Ionization.h"
#include "Tools.h"

class Particles;

class IonizationTunnelBSI : public Ionization
{
   public:
    IonizationTunnelBSI(Params &params, Species *species);
    void operator()(Particles *, unsigned int, unsigned int, std::vector<double> *, Patch *, Projector *, int ipart_ref = 0) override;

   private:
    double ionizationRateBSIQuadratic(int Z, electricFields E);
    double ionizationRateBSILinear(int Z, electricFields E);
    double ionizationRateTunnel(const int Z, const electricFields E);

    double au_to_eV = 27.2116;
    double IH = 13.598434005136;
};

#endif
