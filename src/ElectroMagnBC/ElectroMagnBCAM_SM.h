
#ifndef ELECTROMAGNBCAM_SM_H
#define ELECTROMAGNBCAM_SM_H


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

class ElectroMagnBCAM_SM : public ElectroMagnBCAM
{
public:

    ElectroMagnBCAM_SM( Params &params, Patch *patch, unsigned int _min_max );
    ~ElectroMagnBCAM_SM() {};
    
    virtual void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) override;
    
    void save_fields( Field *, Patch *patch ) override;
    void disableExternalFields() override;
    
    //! Save external fields for silver muller EM Boundary condition
    std::vector< std::complex<double> > Bl_val,  Br_val,  Bt_val;
    
private:


    //! Conversion factor from degree to radian
    double conv_deg2rad;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmin)
    double Alpha_SM_Xmin;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmin)
    double Beta_SM_Xmin;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmin)
    double Gamma_SM_Xmin;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmin)
    double Delta_SM_Xmin;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmin)
    std::complex<double> Epsilon_SM_Xmin;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmax)
    double Alpha_SM_Xmax;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmax)
    double Beta_SM_Xmax;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmax)
    double Gamma_SM_Xmax;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmax)
    double Delta_SM_Xmax;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmax)
    std::complex<double> Epsilon_SM_Xmax;
    
};

#endif

