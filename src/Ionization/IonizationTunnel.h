#ifndef IONIZATIONTUNNEL_H
#define IONIZATIONTUNNEL_H

#include <cmath>
#include <functional>
#include <vector>

#include "Ionization.h"
#include "Tools.h"

class Particles;
struct electricFields;

//! calculate the particle tunnel ionization
class IonizationTunnel : public Ionization
{
   public:
    IonizationTunnel(Params &params, Species *species);

    void operator()(Particles *, unsigned int, unsigned int, std::vector<double> *, Patch *, Projector *,
                    int ipart_ref = 0) override;

    //! method for tunnel ionization with tasks
    void ionizationTunnelWithTasks(Particles *, unsigned int, unsigned int, std::vector<double> *, Patch *, Projector *,
                                   int, int, double *b_Jx, double *b_Jy, double *b_Jz, int ipart_ref = 0) override;

   private:
    double ionizationRateTunnel(const int Z, const electricFields E);
};

#endif
