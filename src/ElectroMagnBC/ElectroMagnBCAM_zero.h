
#ifndef ELECTROMAGNBCAM_zero_H
#define ELECTROMAGNBCAM_zero_H


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

class ElectroMagnBCAM_zero : public ElectroMagnBCAM
{
public:

    ElectroMagnBCAM_zero( Params &params, Patch *patch, unsigned int _min_max );
    ~ElectroMagnBCAM_zero() {};
    
    virtual void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) override;
    
private:

    //! Number of modes
    unsigned int Nmode;
    
};

#endif

