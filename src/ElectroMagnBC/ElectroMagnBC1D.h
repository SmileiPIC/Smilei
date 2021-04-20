
#ifndef ELECTROMAGNBC1D_H
#define ELECTROMAGNBC1D_H

#include "ElectroMagnBC.h"

class Params;
class Patch;
class ElectroMagn;

class ElectroMagnBC1D : public ElectroMagnBC
{
public:

    ElectroMagnBC1D( Params &params, Patch *patch, unsigned int i_boundary );
    virtual ~ElectroMagnBC1D();
    
    virtual void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) = 0;
    
    void applyBConEdges( ElectroMagn *EMfields, Patch *patch );
    
protected:
    
};

#endif
