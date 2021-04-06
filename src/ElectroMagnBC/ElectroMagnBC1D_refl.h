
#ifndef ELECTROMAGNBC1D_REFL_H
#define ELECTROMAGNBC1D_REFL_H

#include "ElectroMagnBC1D.h"

class Params;
class ElectroMagn;

class ElectroMagnBC1D_refl : public ElectroMagnBC1D
{
public:
    ElectroMagnBC1D_refl( Params &param, Patch *patch, unsigned int i_boundary );
    ~ElectroMagnBC1D_refl() {};
    
    void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) override;
    
private:

    //! Oversize
    unsigned int oversize_;
    
};

#endif

