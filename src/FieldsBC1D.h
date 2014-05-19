
#ifndef FIELDSBC1D_H
#define FIELDSBC1D_H

#include "FieldsBC.h" 

class PicParams;
class ElectroMagn;

class FieldsBC1D : public FieldsBC {
public:
    FieldsBC1D( PicParams *param );
    ~FieldsBC1D();

    virtual void apply(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi);

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

