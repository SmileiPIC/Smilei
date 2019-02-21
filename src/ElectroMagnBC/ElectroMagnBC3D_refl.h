
#ifndef ELECTROMAGNBC3D_refl_H
#define ELECTROMAGNBC3D_refl_H

#include "ElectroMagnBC3D.h"

class Params;
class ElectroMagn;

class ElectroMagnBC3D_refl : public ElectroMagnBC3D
{
public:
    ElectroMagnBC3D_refl( Params &params, Patch *patch, unsigned int _min_max );
    ~ElectroMagnBC3D_refl() {};
    
    void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) override;
    
private:

    //! Oversize (nb of ghost cells)
    unsigned int oversize_x, oversize_y, oversize_z;
    
};

#endif

