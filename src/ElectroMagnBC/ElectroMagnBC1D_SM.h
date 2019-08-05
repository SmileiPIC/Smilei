
#ifndef ELECTROMAGNBC1D_SM_H
#define ELECTROMAGNBC1D_SM_H

#include "ElectroMagnBC1D.h"

class Params;
class ElectroMagn;

class ElectroMagnBC1D_SM : public ElectroMagnBC1D
{
public:
    ElectroMagnBC1D_SM( Params &param, Patch *patch, unsigned int _min_max );
    ~ElectroMagnBC1D_SM();
    
    void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) override;
    
    void save_fields( Field *, Patch *patch ) override;
    
    double By_val, Bz_val;
    
    
private:

    //! \todo Create properties the laser time-profile (MG & TV)
    //! Constant used for the Silver-Mueller boundary conditions
    double Alpha_SM;
    
    //! Constant used for the Silver-Mueller boundary conditions
    double Beta_SM;
    
    //! Constant used for the Silver-Mueller boundary conditions
    double Gamma_SM;
    
};

#endif

