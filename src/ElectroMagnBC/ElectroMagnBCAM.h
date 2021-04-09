
#ifndef ELECTROMAGNBCAM_H
#define ELECTROMAGNBCAM_H

#include "ElectroMagnBC.h"

class Params;
class Patch;
class ElectroMagn;

class ElectroMagnBCAM : public ElectroMagnBC
{
public:

    ElectroMagnBCAM( Params &params, Patch *patch, unsigned int i_boundary );
    virtual ~ElectroMagnBCAM();
    
    virtual void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) = 0;
    
    void applyBConEdges( ElectroMagn *EMfields, Patch *patch );
    
protected:

    //! Number of ghost cells of the region in the longitudinal direction 
    int region_oversize_l;

    //! Number of modes
    unsigned int Nmode;
};

#endif
