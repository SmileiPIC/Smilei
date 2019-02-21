
#ifndef ELECTROMAGNBC2D_H
#define ELECTROMAGNBC2D_H

#include "ElectroMagnBC.h"

class Params;
class Patch;
class ElectroMagn;

class ElectroMagnBC2D : public ElectroMagnBC
{
public:

    ElectroMagnBC2D( Params &params, Patch *patch, unsigned int _min_max );
    virtual ~ElectroMagnBC2D();
    
    virtual void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) = 0;
    
    void applyBConEdges( ElectroMagn *EMfields, Patch *patch );
    
protected:

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
    
};

#endif
