
#ifndef ELECTROMAGNBC2D_H
#define ELECTROMAGNBC2D_H

#include "ElectroMagnBC.h"

class Params;
class Patch;
class ElectroMagn;

class ElectroMagnBC2D : public ElectroMagnBC
{
public:

    ElectroMagnBC2D( Params &params, Patch *patch, unsigned int i_boundary );
    virtual ~ElectroMagnBC2D();
    
    virtual void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) = 0;
    
protected:
    
};

#endif
