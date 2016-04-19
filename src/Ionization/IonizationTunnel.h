#ifndef IONIZATIONTUNNEL_H
#define IONIZATIONTUNNEL_H

#include <cmath>

#include <vector>

#include "Ionization.h"
#include "Tools.h"

class Particles;

//! calculate the particle tunnel ionization
class IonizationTunnel : public Ionization
{

public:
    //! Constructor for IonizationTunnel: with no input argument
    IonizationTunnel(Params& params, Species * species);

    //! apply the Tunnel Ionization model to the species (without ionization current)
    virtual void operator() (Particles &particles, int ipart, LocalFields Epart);

    //! apply the Tunnel Ionization model to the species (with ionization current)
    virtual void operator() (Particles &particles, int ipart, LocalFields Epart, LocalFields Jion);

    double one_third;
    std::vector<double> alpha_tunnel;
    std::vector<double> beta_tunnel;
    std::vector<double> gamma_tunnel;

private:
};


#endif
