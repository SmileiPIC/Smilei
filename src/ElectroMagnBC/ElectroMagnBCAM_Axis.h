
#ifndef ELECTROMAGNBCAM_AXIS_H
#define ELECTROMAGNBCAM_AXIS_H


#include <vector>
#include "Tools.h"
#include "ElectroMagnBCAM.h"
#include "ElectroMagnAM.h"
#include "cField2D.h"


class Params;
class ElectroMagn;
class Field;

class ElectroMagnBCAM_Axis : public ElectroMagnBCAM
{
public:

    ElectroMagnBCAM_Axis( Params &params, Patch *patch, unsigned int _min_max );
    ~ElectroMagnBCAM_Axis() {};
    
    virtual void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) override;
    
    void save_fields( Field *, Patch *patch ) override;
    void disableExternalFields() override;
    
private:


    //! Conversion factor from degree to radian
    double conv_deg2rad;
    
    //! Number of modes
    unsigned int Nmode;
    
    //! Oversize along the radial direction
    //unsigned int oversize_R;
    
};

#endif

