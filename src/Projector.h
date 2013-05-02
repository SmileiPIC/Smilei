#ifndef PROJECTOR_H
#define PROJECTOR_H

#include "PicParams.h"

class ElectroMagn;
class Field;
class Particle;


//----------------------------------------------------------------------------------------------------------------------
//! class Projector: contains the virtual operators used during the current projection
//----------------------------------------------------------------------------------------------------------------------
class Projector {
    
public: 
	//! Creator for the Projector
    Projector(PicParams*){};
    
    //! \todo Comment more on this overloading of the () operator
    //! overloading of the () operator
	virtual void operator() (ElectroMagn*, Particle*, double) = 0;
    
    //! \todo Comment more on this overloading of the () operator
    //! overloading of the () operator
	virtual void operator() (Field*, Particle*) = 0;
private:
};

#endif

