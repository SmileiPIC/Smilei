
#ifndef ELECTROMAGNBCAM_BM_H
#define ELECTROMAGNBCAM_BM_H


#include <vector>
#include <complex>
#include "Tools.h"
#include "ElectroMagnBC.h"
#include "ElectroMagnAM.h"
#include "cField2D.h"


class Params;
class ElectroMagn;
class Field;

class ElectroMagnBCAM_BM : public ElectroMagnBC {
public:
    
    ElectroMagnBCAM_BM( Params &params, Patch* patch, unsigned int _min_max );
    ~ElectroMagnBCAM_BM() {};
    
    virtual void apply(ElectroMagn* EMfields, double time_dual, Patch* patch) override;
    
    void save_fields(Field*, Patch* patch) override;
    void disableExternalFields() override;

    //! Save external fields for Buneman EM Boundary condition
    std::vector< std::complex<double> > Bl_val,  Br_val,  Bt_val;
    
private:
    
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
   
    
    //! Constant used for the Buneman boundary conditions (+R)
    double Alpha_Bl_Rmax, Beta_Bl_Rmax, Gamma_Bl_Rmax ;
    
    //! Constant used for the Buneman boundary conditions (+R)
    double  Alpha_Bt_Rmax, Beta_Bt_Rmax, Gamma_Bt_Rmax, Delta_Bt_Rmax, Epsilon_Bt_Rmax ;

    //! Constant used for the Buneman boundary conditions (+R)
    double CB_BM;
    //! Constant used for the Buneman boundary conditions (+R)
    double CE_BM;
    //! Number of modes
    unsigned int Nmode;
    
};

#endif


