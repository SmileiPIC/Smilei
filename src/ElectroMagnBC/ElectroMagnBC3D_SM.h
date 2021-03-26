
#ifndef ELECTROMAGNBC3D_SM_H
#define ELECTROMAGNBC3D_SM_H


#include <vector>
#include "Tools.h"
#include "ElectroMagnBC3D.h"
#include "ElectroMagn3D.h"
#include "Field3D.h"
#include "Field2D.h"


class Params;
class ElectroMagn;
class Field;

class ElectroMagnBC3D_SM : public ElectroMagnBC3D
{
public:

    ElectroMagnBC3D_SM( Params &params, Patch *patch, unsigned int i_boundary );
    ~ElectroMagnBC3D_SM();
    
    void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) override;
    void save_fields( Field *, Patch *patch ) override;
    void disableExternalFields() override;
    
    //! Save external fields for silver muller EM Boundary condition
    std::vector<Field2D *> B_val;
    
    //! Constants used for the Silver-Mueller boundary conditions 
    double Alpha_, Beta_, Gamma_, Delta_, Epsilon_, Zeta_, Eta_;
    
    unsigned int axis0_, axis1_, axis2_;
    int sign_;
    std::vector<unsigned int> iB_;
};

#endif

