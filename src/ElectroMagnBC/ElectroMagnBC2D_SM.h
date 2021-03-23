
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

    ElectroMagnBC2D_SM( Params &params, Patch *patch, unsigned int _min_max );
    ~ElectroMagnBC2D_SM() {};
    
    virtual void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) override;
    
    void save_fields( Field *, Patch *patch ) override;
    void disableExternalFields() override;
    
    //! Save external fields for silver muller EM Boundary condition
    std::vector<double> Bx_val,  By_val,  Bz_val;
    
private:


    //! Constant used for the Silver-Mueller boundary conditions (Xmin)
    double Alpha_Xmin, Beta_Xmin, Gamma_Xmin, Delta_Xmin, Epsilon_Xmin;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmax)
    double Alpha_Xmax, Beta_Xmax, Gamma_Xmax, Delta_Xmax, Epsilon_Xmax;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse)
    double Alpha_Ymin, Beta_Ymin, Gamma_Ymin, Delta_Ymin, Epsilon_Ymin;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse)
    double Alpha_Ymax, Beta_Ymax, Gamma_Ymax, Delta_Ymax, Epsilon_Ymax;
    
};

#endif

