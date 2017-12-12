
#ifndef ELECTROMAGNBC3D_BM_H
#define ELECTROMAGNBC3D_BM_H


#include <vector>
#include "Tools.h"
#include "ElectroMagnBC.h" 
#include "ElectroMagn3D.h" 
#include "Field3D.h"
#include "Field2D.h"


class Params;
class ElectroMagn;
class Field;

class ElectroMagnBC3D_BM : public ElectroMagnBC {
public:

    ElectroMagnBC3D_BM( Params &params, Patch* patch, unsigned int _W_E );
    ~ElectroMagnBC3D_BM();

    void apply(ElectroMagn* EMfields, double time_dual, Patch* patch) override;
    void save_fields(Field*, Patch* patch) override;
    void disableExternalFields() override;

    //! Save external fields for Buneman EM Boundary condition
    Field2D *Bx_val, *By_val, *Bz_val;
    
private:
    
    //! Conversion factor from degree to radian
    double conv_deg2rad;
    
    //! Number of nodes on the primal grid in the x-direction
    unsigned int nx_p;
    
    //! Number of nodes on the dual grid in the x-direction
    unsigned int nx_d;
    
    //! Number of nodes on the primal grid in the y-direction
    unsigned int ny_p;
    
    //! Number of nodes on the dual grid in the y-direction
    unsigned int ny_d;
    
    //! Number of nodes on the primal grid in the y-direction
    unsigned int nz_p;
    
    //! Number of nodes on the dual grid in the z-direction
    unsigned int nz_d;
    
    //! Spatial step dx for 3D3V cartesian simulations
    double dx;
    
    //! Spatial step dy for 3D3V cartesian simulations
    double dy;

    //! Spatial step dz for 3D3V cartesian simulations
    double dz;

    //! Ratio of the time-step by the spatial-step dt/dx for 3D3V cartesian simulations
    double dt_ov_dx;
    
    //! Ratio of the time-step by the spatial-step dt/dy for 3D3V cartesian simulations
    double dt_ov_dy;
  
    //! Ratio of the time-step by the spatial-step dt/dz for 3D3V cartesian simulations
    double dt_ov_dz;
  
    //! Ratio of the spatial-step by the time-step dx/dt for 3D3V cartesian simulations
    double dx_ov_dt;
    
    //! Ratio of the spatial-step by the time-step dy/dt for 3D3V cartesian simulations
    double dy_ov_dt;
    
    //! Ratio of the spatial-step by the time-step dz/dt for 3D3V cartesian simulations
    double dz_ov_dt;



    //! Constant used for the Buneman boundary conditions (Xmin)
    double Alpha_BM_xmin;
    
    //! Constant used for the Buneman boundary conditions (Xmin)
    double Beta_BM_xmin;
    
    //! Constant used for the Buneman boundary conditions (Xmin)
    double Gamma_BM_xmin;
    
    //! Constant used for the Buneman boundary conditions (Xmax)
    double Alpha_BM_xmax;
    
    //! Constant used for the Buneman boundary conditions (Xmax)
    double Beta_BM_xmax;
    
    //! Constant used for the Buneman boundary conditions (Xmax)
    double Gamma_BM_xmax;
    

    //! Constant used for the Buneman boundary conditions (Transverse Y)
    double Alpha_BM_ymin;
    
    //! Constant used for the Buneman boundary conditions (Transverse Y)
    double Beta_BM_ymin;
    
    //! Constant used for the Buneman boundary conditions (Transverse Y)
    double Gamma_BM_ymin;

    //! Constant used for the Buneman boundary conditions (Transverse Y)
    double Alpha_BM_ymax;
    
    //! Constant used for the Buneman boundary conditions (Transverse Y)
    double Beta_BM_ymax;
    
    //! Constant used for the Buneman boundary conditions (Trnsverse Y)
    double Gamma_BM_ymax;

     //! Constant used for the Buneman boundary conditions (Transverse Z)
    double Alpha_BM_zmin;
    
    //! Constant used for the Buneman boundary conditions (Transverse Z)
    double Beta_BM_zmin;
    
    //! Constant used for the Buneman boundary conditions (Trnsverse Z)
    double Gamma_BM_zmin;

    //! Constant used for the Buneman boundary conditions (Transverse Z)
    double Alpha_BM_zmax;
    
    //! Constant used for the Buneman boundary conditions (Transverse Z)
    double Beta_BM_zmax;
    
    //! Constant used for the Buneman boundary conditions (Trnsverse Z)
    double Gamma_BM_zmax;

};

#endif

