
#ifndef ELECTROMAGNBC1D_SM_H
#define ELECTROMAGNBC1D_SM_H

#include "ElectroMagnBC1D.h"

class Params;
class ElectroMagn;

class ElectroMagnBC1D_SM : public ElectroMagnBC1D
{
public:
    ElectroMagnBC1D_SM( Params &param, Patch *patch, unsigned int i_boundary );
    ~ElectroMagnBC1D_SM();
    
    void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) override;
    
    void save_fields( Field *, Patch *patch ) override;
    
    double By_val, Bz_val;
    
    
private:

    //! Constants used for the Silver-Mueller boundary conditions
    double Alpha, Beta, Gamma;
    
    //! Locations to apply the profile
    unsigned int iE, iB, iB_old;
    int sign_;
    
};

#endif

