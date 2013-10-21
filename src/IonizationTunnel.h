#ifndef IONIZATIONTUNNEL_H
#define IONIZATIONTUNNEL_H

#include "Ionization.h"
#include "Tools.h"

#include <vector>
#include <cmath>


//! class Ionization1D used to defined a 1d vector
class IonizationTunnel : public Ionization
{
    
public:
	//! Constructor for Ionization1D: with no input argument
	IonizationTunnel(PicParams *params, int ispec);
    
	//! apply the Tunnel Ionization model to the specie
	virtual void operator() (Particle* part, LocalFields Epart);
	
    
private:
};


#endif
