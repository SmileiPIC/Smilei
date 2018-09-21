
#ifndef ELECTROMAGNBCRZ_SM_H
#define ELECTROMAGNBCRZ_SM_H


#include <vector>
#include <complex>
#include "Tools.h"
#include "ElectroMagnBC.h"
#include "ElectroMagnAM.h"
#include "cField2D.h"
#include "dcomplex.h"

class Params;
class ElectroMagn;
class Field;

class ElectroMagnBCAM_SM : public ElectroMagnBC {
public:
    
    ElectroMagnBCAM_SM( Params &params, Patch* patch, unsigned int _min_max );
    ~ElectroMagnBCAM_SM() {};
    
    virtual void apply(ElectroMagn* EMfields, double time_dual, Patch* patch) override;
    
    void save_fields(Field*, Patch* patch) override;
    void disableExternalFields() override;

    //! Save external fields for silver muller EM Boundary condition
    std::vector< std::complex<double> > Bl_val,  Br_val,  Bt_val;
    
private:
    
    
    //! Conversion factor from degree to radian
    double conv_deg2rad;
    
    //! Number of nodes on the primal grid in the x-direction
    unsigned int nl_p;
    
    //! Number of nodes on the dual grid in the x-direction
    unsigned int nl_d;
    
    //! Number of nodes on the primal grid in the y-direction
    unsigned int nr_p;
    
    //! Number of nodes on the dual grid in the y-direction
    unsigned int nr_d;
    
    //! Spatial step dx for 2D3V cartesian simulations
    double dl;
    
    //! Spatial step dy for 2D3V cartesian simulations
    double dr;
    
    //! Ratio of the time-step by the spatial-step dt/dx for 2D3V cartesian simulations
    double dt_ov_dl;
    
    //! Ratio of the time-step by the spatial-step dt/dy for 2D3V cartesian simulations
    double dt_ov_dr;
    
    //! Ratio of the spatial-step by the time-step dx/dt for 2D3V cartesian simulations
    double dr_ov_dt;
    
    //! Ratio of the spatial-step by the time-step dy/dt for 2D3V cartesian simulations
    double dl_ov_dt;
    
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
	//! Number of modes
	unsigned int Nmode;
    
};

#endif

