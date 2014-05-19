
#ifndef FIELDSBC2D_H
#define FIELDSBC2D_H

#include "FieldsBC.h" 

class PicParams;
class ElectroMagn;

class FieldsBC2D : public FieldsBC {
public:
    FieldsBC2D( PicParams *params );
    ~FieldsBC2D();

    virtual void apply(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi);

 private:
    //! Number of nodes on the primal grid in the x-direction
    unsigned int nx_p;

    //! Number of nodes on the dual grid in the x-direction
    unsigned int nx_d;

    //! Number of nodes on the primal grid in the y-direction
    unsigned int ny_p;

    //! Number of nodes on the dual grid in the y-direction
    unsigned int ny_d;

    //! Spatial step dx for 2D3V cartesian simulations
    double dx;

    //! Spatial step dy for 2D3V cartesian simulations
    double dy;

    //! Ratio of the time-step by the spatial-step dt/dx for 2D3V cartesian simulations
    double dt_ov_dx;

    //! Ratio of the time-step by the spatial-step dt/dy for 2D3V cartesian simulations
    double dt_ov_dy;

    //! Ratio of the spatial-step by the time-step dx/dt for 2D3V cartesian simulations
    double dx_ov_dt;

    //! Ratio of the spatial-step by the time-step dy/dt for 2D3V cartesian simulations
    double dy_ov_dt;

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
    
};

#endif

