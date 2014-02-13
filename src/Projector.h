#ifndef PROJECTOR_H
#define PROJECTOR_H

#include "PicParams.h"
#include "SmileiMPI.h"
#include "Field.h"

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
    //! overloading of the () operator
	virtual void operator() (Field* Jx, Field* Jy, Field* Jz, Field* rho, Particle* part, double gf) = 0;
    //! overloading of the () operator
	virtual void operator() (Field* Jx, Field* Jy, Field* Jz, Particle* part, LocalFields Jion) = 0;
	//! overloading of the () operator
	virtual void operator() (double* Jx, double* Jy, double* Jz, Particle* part, double gf, unsigned int bin, unsigned int b_dim0) = 0;
private:
};

#endif

