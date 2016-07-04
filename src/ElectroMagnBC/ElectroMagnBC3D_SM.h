
#ifndef ELECTROMAGNBC3D_SM_H
#define ELECTROMAGNBC3D_SM_H


#include <vector>
#include "Tools.h"
#include "ElectroMagnBC.h" 
#include "ElectroMagn3D.h" 
#include "Field3D.h"
#include "Field2D.h"


class Params;
class ElectroMagn;
class Field;

class ElectroMagnBC3D_SM : public ElectroMagnBC {
public:

    ElectroMagnBC3D_SM( Params &params, Patch* patch );
    ~ElectroMagnBC3D_SM();

    virtual void apply_xmin(ElectroMagn* EMfields, double time_dual, Patch* patch);
    virtual void apply_xmax(ElectroMagn* EMfields, double time_dual, Patch* patch);
    virtual void apply_ymin(ElectroMagn* EMfields, double time_dual, Patch* patch);
    virtual void apply_ymax(ElectroMagn* EMfields, double time_dual, Patch* patch);
    virtual void apply_zmin(ElectroMagn* EMfields, double time_dual, Patch* patch);
    virtual void apply_zmax(ElectroMagn* EMfields, double time_dual, Patch* patch);

private:
    
    virtual void save_fields_BC3D_Long(Field*);
    virtual void save_fields_BC3D_TransY(Field*);
    virtual void save_fields_BC3D_TransZ(Field*);

 	//! Save external fields for silver muller EM Boundary condition
     Field2D *Bz_xvalmin_Long,  *Bz_xvalmax_Long,  *By_xvalmin_Long,  *By_xvalmax_Long,  *Bx_xvalmin_Long,  *Bx_xvalmax_Long,
         *Bz_yvalmin_TransY, *Bz_yvalmax_TransY, *By_yvalmin_TransY, *By_yvalmax_TransY, *Bx_yvalmin_TransY, *Bx_yvalmax_TransY,
         *Bz_zvalmin_TransZ, *Bz_zvalmax_TransZ, *By_zvalmin_TransZ, *By_zvalmax_TransZ, *Bx_zvalmin_TransZ, *Bx_zvalmax_TransZ;
    
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



    //! Constant used for the Silver-Mueller boundary conditions (West)
    double Alpha_SM_W;
    
    //! Constant used for the Silver-Mueller boundary conditions (West)
    double Beta_SM_W;
    
    //! Constant used for the Silver-Mueller boundary conditions (West)
    double Gamma_SM_W;
    
    //! Constant used for the Silver-Mueller boundary conditions (West)
    double Delta_SM_W;
    
    //! Constant used for the Silver-Mueller boundary conditions (West)
    double Epsilon_SM_W;
    
    //! Constant used for the Silver-Mueller boundary conditions (East)
    double Alpha_SM_E;
    
    //! Constant used for the Silver-Mueller boundary conditions (East)
    double Beta_SM_E;
    
    //! Constant used for the Silver-Mueller boundary conditions (East)
    double Gamma_SM_E;
    
    //! Constant used for the Silver-Mueller boundary conditions (East)
    double Delta_SM_E;
    
    //! Constant used for the Silver-Mueller boundary conditions (East)
    double Epsilon_SM_E;
    

    //! Constant used for the Silver-Mueller boundary conditions (Transverse Y)
    double Alpha_SM_S;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Y)
    double Beta_SM_S;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Y)
    double Delta_SM_S;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Y)
    double Epsilon_SM_S;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Y)
    double Alpha_SM_N;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Y)
    double Beta_SM_N;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Y)
    double Delta_SM_N;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Y)
    double Epsilon_SM_N;


     //! Constant used for the Silver-Mueller boundary conditions (Transverse Z)
    double Alpha_SM_B;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Z)
    double Beta_SM_B;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Z)
    double Delta_SM_B;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Z)
    double Epsilon_SM_B;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Z)
    double Alpha_SM_T;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Z)
    double Beta_SM_T;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Z)
    double Delta_SM_T;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Z)
    double Epsilon_SM_T;
   
};

#endif

