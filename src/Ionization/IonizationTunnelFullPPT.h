#ifndef IONIZATIONTUNNELFULLPPT_H
#define IONIZATIONTUNNELFULLPPT_H

#include <cmath>
#include <vector>

#include "Ionization.h"
#include "Tools.h"

class Particles;

//! calculate the particle tunnel ionization
class IonizationTunnelFullPPT : public Ionization
{
   public:
    //! Constructor for IonizationTunnelFullPPT: with no input argument
    IonizationTunnelFullPPT(Params &params, Species *species);
    void operator()(Particles *particles, unsigned int ipart_min, unsigned int ipart_max,
                                  std::vector<double> *Epart, Patch *patch, Projector *Proj, int ipart_ref = 0) override;

   private:
    std::vector<double> Magnetic_quantum_number;
    double ionizationRateTunnel(const int Z, const electricFields E);
};


#endif
