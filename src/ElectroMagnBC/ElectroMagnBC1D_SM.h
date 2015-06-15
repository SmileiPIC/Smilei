
#ifndef ELECTROMAGNBC1D_SM_H
#define ELECTROMAGNBC1D_SM_H

#include "ElectroMagnBC.h" 

class PicParams;
class ElectroMagn;

class ElectroMagnBC1D_SM : public ElectroMagnBC {
public:
    ElectroMagnBC1D_SM( PicParams &param, LaserParams &laser_params);
    ~ElectroMagnBC1D_SM();

    virtual void apply(ElectroMagn* EMfields, double time_dual, Patch* patch);

 private:
    //! Number of nodes on the primal grid
    unsigned int nx_p;

    //! Number of nodes on the dual grid
    unsigned int nx_d;

    //! Spatial step dx for 1d3v cartesian simulations
    double dx;

    //! Ratio of the time-step by the spatial-step dt/dx for 1d3v cartesian simulations
    double dt_ov_dx;

    //! Ratio of the spatial-step by the time-step dx/dt for 1d3v cartesian simulations
    double dx_ov_dt;

    //! \todo Create properties the laser time-profile (MG & TV)
    //! Constant used for the Silver-Mueller boundary conditions
    double Alpha_SM;

    //! Constant used for the Silver-Mueller boundary conditions
    double Beta_SM;

    //! Constant used for the Silver-Mueller boundary conditions
    double Gamma_SM;

    
};

#endif

