
#ifndef ELECTROMAGNBCAM_ramp_H
#define ELECTROMAGNBCAM_ramp_H


#include <vector>
#include <complex>
#include "Tools.h"
#include "ElectroMagnBCAM.h"
#include "ElectroMagnAM.h"
#include "cField2D.h"
#include "dcomplex.h"

class Params;
class ElectroMagn;
class Field;

class ElectroMagnBCAM_ramp : public ElectroMagnBCAM
{
public:

    ElectroMagnBCAM_ramp( Params &params, Patch *patch, unsigned int i_boundary, unsigned int ncells );
    ~ElectroMagnBCAM_ramp() {};
    
    virtual void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) override;
    
private:

    //! Number of modes
    unsigned int Nmode;
    
    // number of damping cells
    unsigned int number_of_cells_;
    
};

#endif

