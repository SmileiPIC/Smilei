
#ifndef ELECTROMAGNBC3D_BM_H
#define ELECTROMAGNBC3D_BM_H


#include <vector>
#include "Tools.h"
#include "ElectroMagnBC3D.h"
#include "ElectroMagn3D.h"
#include "Field3D.h"
#include "Field2D.h"


class Params;
class ElectroMagn;
class Field;

class ElectroMagnBC3D_BM : public ElectroMagnBC3D
{
public:

    ElectroMagnBC3D_BM( Params &params, Patch *patch, unsigned int _W_E );
    ~ElectroMagnBC3D_BM();
    
    void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) override;
    void save_fields( Field *, Patch *patch ) override;
    void disableExternalFields() override;
    
    //! Save external fields for Buneman EM Boundary condition
    Field2D *Bx_val, *By_val, *Bz_val;
    
private:

    //! Constants used for the Buneman boundary conditions depending on theta and phi.
    double cb[3][2], ce[3][2];
    
    //! Constant used for the Buneman boundary conditions (Xmin)
    double Alpha_BM_x;
    
    //! Constant used for the Buneman boundary conditions (Xmin)
    double Beta_BM_x;
    
    //! Constant used for the Buneman boundary conditions (Xmin)
    double Gamma_BM_x;
    
    //! Constant used for the Buneman boundary conditions (Transverse Y)
    double Alpha_BM_y;
    
    //! Constant used for the Buneman boundary conditions (Transverse Y)
    double Beta_BM_y;
    
    //! Constant used for the Buneman boundary conditions (Transverse Y)
    double Gamma_BM_y;
    
    //! Constant used for the Buneman boundary conditions (Transverse Z)
    double Alpha_BM_z;
    
    //! Constant used for the Buneman boundary conditions (Transverse Z)
    double Beta_BM_z;
    
    //! Constant used for the Buneman boundary conditions (Trnsverse Z)
    double Gamma_BM_z;
    
};

#endif

