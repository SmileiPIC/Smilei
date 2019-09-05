
#ifndef ELECTROMAGNBC1D_H
#define ELECTROMAGNBC1D_H

#include "ElectroMagnBC.h"

class Params;
class Patch;
class ElectroMagn;

class ElectroMagnBC1D : public ElectroMagnBC
{
public:

    ElectroMagnBC1D( Params &params, Patch *patch, unsigned int _min_max );
    virtual ~ElectroMagnBC1D();
    
    virtual void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) = 0;
    
    void applyBConEdges( ElectroMagn *EMfields, Patch *patch );
    
protected:

    //! Number of nodes on the primal grid in the x-direction
    unsigned int nx_p;
    
    //! Number of nodes on the dual grid in the x-direction
    unsigned int nx_d;
    
    //! Spatial step dx for 1D3V cartesian simulations
    double dx;
    
    //! Ratio of the time-step by the spatial-step dt/dx for 1D3V cartesian simulations
    double dt_ov_dx;
    
    //! Ratio of the spatial-step by the time-step dx/dt for 1D3V cartesian simulations
    double dx_ov_dt;
    
};

#endif
