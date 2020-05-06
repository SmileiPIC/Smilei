
#ifndef ELECTROMAGNBCAM_BM_H
#define ELECTROMAGNBCAM_BM_H


#include <vector>
#include <complex>
#include "Tools.h"
#include "ElectroMagnBCAM.h"
#include "ElectroMagnAM.h"
#include "cField2D.h"


class Params;
class ElectroMagn;
class Field;

class ElectroMagnBCAM_BM : public ElectroMagnBCAM
{
public:

    ElectroMagnBCAM_BM( Params &params, Patch *patch, unsigned int _min_max );
    ~ElectroMagnBCAM_BM() {};
    
    virtual void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) override;
    
    void save_fields( Field *, Patch *patch ) override;
    void disableExternalFields() override;
    
    //! Save external fields for Buneman EM Boundary condition
    std::vector< std::complex<double> > Bl_val,  Br_val,  Bt_val;
    
private:

    
    //! Constant used for the Buneman boundary conditions (+R)
    double Alpha_Bl_Rmax, Beta_Bl_Rmax, Gamma_Bl_Rmax ;
    
    //! Constant used for the Buneman boundary conditions (+R)
    double  Alpha_Bt_Rmax, Beta_Bt_Rmax, Gamma_Bt_Rmax, Delta_Bt_Rmax, Epsilon_Bt_Rmax ;
    
    //! Constant used for the Buneman boundary conditions (+R)
    double CB_BM;
    //! Constant used for the Buneman boundary conditions (+R)
    double CE_BM;
};

#endif


