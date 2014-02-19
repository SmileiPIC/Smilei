#ifndef PROJECTOR_H
#define PROJECTOR_H

#include "PicParams.h"
#include "SmileiMPI.h"
#include "Field.h"

class ElectroMagn;
class Field;
class Particles;


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
    virtual void operator() (ElectroMagn* champs, Particles &particles, int ipart, double gf) = 0;
    //! \todo Comment more on this overloading of the () operator
    //! overloading of the () operator
    virtual void operator() (Field* rho, Particles &particles, int ipart) = 0;
    //! overloading of the () operator
    virtual void operator() (Field* Jx, Field* Jy, Field* Jz, Field* rho, Particles &particles, int ipart, double gf) = 0;
    //! overloading of the () operator
    virtual void operator() (Field* Jx, Field* Jy, Field* Jz, Particles &particles, int ipart, LocalFields Jion) = 0;
    //! overloading of the () operator
    virtual void operator() (double* Jx, double* Jy, double* Jz, Particles &particles, int ipart, double gf, unsigned int bin, unsigned int b_dim0) = 0;
private:
};

#endif

