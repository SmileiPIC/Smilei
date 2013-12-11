#ifndef PROJECTOR_H
#define PROJECTOR_H

#include "PicParams.h"
#include "SmileiMPI.h"

class ElectroMagn;
class Field;
class Particle;


//----------------------------------------------------------------------------------------------------------------------
//! class Projector: contains the virtual operators used during the current projection
//----------------------------------------------------------------------------------------------------------------------
class Projector {
    
public: 
	//! Creator for the Projector
	Projector(PicParams*, SmileiMPI*) {};
	virtual ~Projector() {};
	//! \todo Comment more on this overloading of the () operator
	//! overloading of the () operator
	virtual void operator() (ElectroMagn* champs, Particle* part, double gf) = 0;
	//! \todo Comment more on this overloading of the () operator
	//! overloading of the () operator
	virtual void operator() (Field* rho, Particle* part) = 0;
    
private:
};

#endif

