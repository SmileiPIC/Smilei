
#ifndef ELECTROMAGNBC2D_refl_H
#define ELECTROMAGNBC2D_refl_H

#include "ElectroMagnBC2D.h"

class Params;
class ElectroMagn;

class ElectroMagnBC2D_refl : public ElectroMagnBC2D
{
public:
    ElectroMagnBC2D_refl( Params &params, Patch *patch, unsigned int i_boundary );
    ~ElectroMagnBC2D_refl() {};
    
    void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) override;
    
private:

    //! Oversize (nb of ghost cells)
    unsigned int oversize_;
    
};

#endif

