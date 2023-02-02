
#ifndef ELECTROMAGNBC3D_H
#define ELECTROMAGNBC3D_H

#include "ElectroMagnBC.h"

class Params;
class Patch;
class ElectroMagn;

class ElectroMagnBC3D : public ElectroMagnBC
{
public:

    ElectroMagnBC3D( Params &params, Patch *patch, unsigned int i_boundary );
    virtual ~ElectroMagnBC3D();
    
    virtual void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) = 0;
    
protected:
    
};

#endif
