
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

    ElectroMagnBC3D_SM( Params &params, Patch* patch, unsigned int _min_max );
    ~ElectroMagnBC3D_SM();

    void apply(ElectroMagn* EMfields, double time_dual, Patch* patch) override;
    void save_fields(Field*, Patch* patch) override;
    void disableExternalFields() override;

    //! Save external fields for silver muller EM Boundary condition
    Field2D *Bx_val, *By_val, *Bz_val;
    
private:
    
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
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmin)
    double Zeta_SM_W;

    //! Constant used for the Silver-Mueller boundary conditions (Xmin)
    double Eta_SM_W;

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
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmin)
    double Zeta_SM_E;

    //! Constant used for the Silver-Mueller boundary conditions (Xmin)
    double Eta_SM_E;


    //! Constant used for the Silver-Mueller boundary conditions (Transverse Y)
    double Alpha_SM_S;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Y)
    double Beta_SM_S;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Y)
    double Delta_SM_S;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Y)
    double Epsilon_SM_S;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmin)
    double Zeta_SM_S;

    //! Constant used for the Silver-Mueller boundary conditions (Xmin)
    double Eta_SM_S;

    //! Constant used for the Silver-Mueller boundary conditions (Transverse Y)
    double Alpha_SM_N;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Y)
    double Beta_SM_N;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Y)
    double Delta_SM_N;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Y)
    double Epsilon_SM_N;

    //! Constant used for the Silver-Mueller boundary conditions (Xmin)
    double Zeta_SM_N;

    //! Constant used for the Silver-Mueller boundary conditions (Xmin)
    double Eta_SM_N;


     //! Constant used for the Silver-Mueller boundary conditions (Transverse Z)
    double Alpha_SM_B;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Z)
    double Beta_SM_B;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Z)
    double Delta_SM_B;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Z)
    double Epsilon_SM_B;
    
    //! Constant used for the Silver-Mueller boundary conditions (Xmin)
    double Zeta_SM_B;

    //! Constant used for the Silver-Mueller boundary conditions (Xmin)
    double Eta_SM_B;

    //! Constant used for the Silver-Mueller boundary conditions (Transverse Z)
    double Alpha_SM_T;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Z)
    double Beta_SM_T;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Z)
    double Delta_SM_T;
    
    //! Constant used for the Silver-Mueller boundary conditions (Transverse Z)
    double Epsilon_SM_T;
   
    //! Constant used for the Silver-Mueller boundary conditions (Xmin)
    double Zeta_SM_T;

    //! Constant used for the Silver-Mueller boundary conditions (Xmin)
    double Eta_SM_T;

};

#endif

