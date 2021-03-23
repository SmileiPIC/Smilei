
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

    ElectroMagnBC3D_SM( Params &params, Patch *patch, unsigned int _min_max );
    ~ElectroMagnBC3D_SM();
    
    void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) override;
    void save_fields( Field *, Patch *patch ) override;
    void disableExternalFields() override;
    
    //! Save external fields for silver muller EM Boundary condition
    Field2D *Bx_val, *By_val, *Bz_val;
    
    //! Constants used for the Silver-Mueller boundary conditions (Xmin)
    double Alpha_Xmin, Beta_Xmin, Gamma_Xmin, Delta_Xmin, Epsilon_Xmin, Zeta_Xmin, Eta_Xmin;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmax)
    double Alpha_Xmax, Beta_Xmax, Gamma_Xmax, Delta_Xmax, Epsilon_Xmax, Zeta_Xmax, Eta_Xmax;
    
    
    //! Constant used for the Silver-Mueller boundary conditions (Ymin)
    double Alpha_Ymin, Beta_Ymin, Gamma_Ymin, Delta_Ymin, Epsilon_Ymin, Zeta_Ymin, Eta_Ymin;
    
    //! Constant used for the Silver-Mueller boundary conditions (Ymax)
    double Alpha_Ymax, Beta_Ymax, Gamma_Ymax, Delta_Ymax, Epsilon_Ymax, Zeta_Ymax, Eta_Ymax;
    
    
    //! Constant used for the Silver-Mueller boundary conditions (Zmin)
    double Alpha_Zmin, Beta_Zmin, Delta_Zmin, Epsilon_Zmin, Zeta_Zmin, Eta_Zmin;
    
    //! Constant used for the Silver-Mueller boundary conditions (Zmax)
    double Alpha_Zmax, Beta_Zmax, Delta_Zmax, Epsilon_Zmax, Zeta_Zmax, Eta_Zmax;
    
};

#endif

