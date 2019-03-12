
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
    double Alpha_SM_W;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmin)
    double Beta_SM_W;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmin)
    double Gamma_SM_W;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmin)
    double Delta_SM_W;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmin)
    double Epsilon_SM_W;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmax)
    double Alpha_SM_E;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmax)
    double Beta_SM_E;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmax)
    double Gamma_SM_E;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmax)
    double Delta_SM_E;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmax)
    double Epsilon_SM_E;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse)
    double Alpha_SM_S;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse)
    double Beta_SM_S;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse)
    double Delta_SM_S;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse)
    double Epsilon_SM_S;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse)
    double Alpha_SM_N;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse)
    double Beta_SM_N;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse)
    double Delta_SM_N;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse)
    double Epsilon_SM_N;
    
};

#endif

