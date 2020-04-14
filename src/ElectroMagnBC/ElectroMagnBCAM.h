
#ifndef ELECTROMAGNBCAM_H
#define ELECTROMAGNBCAM_H

#include "ElectroMagnBC.h"

class Params;
class Patch;
class ElectroMagn;

class ElectroMagnBCAM : public ElectroMagnBC
{
public:

    ElectroMagnBCAM( Params &params, Patch *patch, unsigned int _min_max );
    virtual ~ElectroMagnBCAM();
    
    virtual void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) = 0;
    
    void applyBConEdges( ElectroMagn *EMfields, Patch *patch );
    
protected:

    //! Number of nodes on the primal grid in the l-direction
    unsigned int nl_p;
    
    //! Number of nodes on the dual grid in the l-direction
    unsigned int nl_d;
    
    //! Number of nodes on the primal grid in the r-direction
    unsigned int nr_p;
    
    //! Number of nodes on the dual grid in the r-direction
    unsigned int nr_d;
    
    //! Spatial step dl for AM simulations
    double dl;
    
    //! Spatial step dr for AM simulations
    double dr;
    
    //! Ratio of the time-step by the spatial-step dt/dl for AM simulations
    double dt_ov_dl;
    
    //! Ratio of the time-step by the spatial-step dt/dr for AM simulations
    double dt_ov_dr;
    
    //! Ratio of the spatial-step by the time-step dl/dt for AM simulations
    double dl_ov_dt;
    
    //! Ratio of the spatial-step by the time-step dr/dt for AM simulations
    double dr_ov_dt;
   
    //! Number of ghost cells of the region in the longitudinal direction 
    int region_oversize_l;

    //! Number of modes
    unsigned int Nmode;
};

#endif
