
#ifndef ELECTROMAGNBC2D_SM_H
#define ELECTROMAGNBC2D_SM_H


#include <vector>
#include "Tools.h"
#include "ElectroMagnBC2D.h"
#include "ElectroMagn2D.h"
#include "Field2D.h"


class Params;
class ElectroMagn;
class Field;

class ElectroMagnBC2D_SM : public ElectroMagnBC2D
{
public:

    ElectroMagnBC2D_SM( Params &params, Patch *patch, unsigned int i_boundary );
    ~ElectroMagnBC2D_SM() {};
    
    virtual void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) override;
    
    void save_fields( Field *, Patch *patch ) override;
    void disableExternalFields() override;
    
    //! Save external fields for silver muller EM Boundary condition
    std::vector<std::vector<double> > B_val;
    
private:


    //! Constant used for the Silver-Mueller boundary conditions
    double Alpha_, Beta_, Gamma_, Delta_, Epsilon_;
    
    unsigned int axis0_, axis1_;
    int sign_;
    std::vector<unsigned int> iB_;
};

#endif

