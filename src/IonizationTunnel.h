#ifndef IONIZATIONTUNNEL_H
#define IONIZATIONTUNNEL_H

#include "Ionization.h"
#include "Tools.h"

#include <vector>
#include <cmath>

class Particle;

//! calculate the particle tunnel ionization
class IonizationTunnel : public Ionization
{
    
public:
	//! Constructor for IonizationTunnel: with no input argument
	IonizationTunnel(PicParams *params, int ispec);
    
	//! apply the Tunnel Ionization model to the species (without ionization current)
	virtual void operator() (Particle* part, LocalFields Epart);
	
    //! apply the Tunnel Ionization model to the species (with ionization current)
	virtual void operator() (Particle* part, LocalFields Epart, LocalFields Jion);
    
    double one_third;
    std::vector<double> alpha_tunnel;
    std::vector<double> beta_tunnel;
    std::vector<double> gamma_tunnel;
    
private:
};


#endif
