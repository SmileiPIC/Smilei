
#ifndef ELECTROMAGNBC2D_refl_H
#define ELECTROMAGNBC2D_refl_H

#include "ElectroMagnBC.h"

class Params;
class ElectroMagn;

class ElectroMagnBC2D_refl : public ElectroMagnBC {
public:
    ElectroMagnBC2D_refl( Params &params, Patch* patch );
    ~ElectroMagnBC2D_refl();
    
    virtual void apply_xmin(ElectroMagn* EMfields, double time_dual, Patch* patch);
    virtual void apply_xmax(ElectroMagn* EMfields, double time_dual, Patch* patch);
    virtual void apply_ymin(ElectroMagn* EMfields, double time_dual, Patch* patch);
    virtual void apply_ymax(ElectroMagn* EMfields, double time_dual, Patch* patch);
    virtual void apply_zmin(ElectroMagn* EMfields, double time_dual, Patch* patch);
    virtual void apply_zmax(ElectroMagn* EMfields, double time_dual, Patch* patch);
    
private:
    
    //! Oversize (nb of ghost cells)
    unsigned int oversize_;
    
    //! Number of nodes on the primal grid in the x-direction
    unsigned int nx_p;
    
    //! Number of nodes on the dual grid in the x-direction
    unsigned int nx_d;
    
    //! Number of nodes on the primal grid in the y-direction
    unsigned int ny_p;
    
    //! Number of nodes on the dual grid in the y-direction
    unsigned int ny_d;
    
};

#endif

