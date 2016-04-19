
#ifndef ELECTROMAGNBC1D_REFL_H
#define ELECTROMAGNBC1D_REFL_H

#include "ElectroMagnBC.h" 

class Params;
class ElectroMagn;

class ElectroMagnBC1D_refl : public ElectroMagnBC {
public:
    ElectroMagnBC1D_refl( Params &param, Patch* patch );
    ~ElectroMagnBC1D_refl();

    virtual void apply_xmin(ElectroMagn* EMfields, double time_dual, Patch* patch);
    virtual void apply_xmax(ElectroMagn* EMfields, double time_dual, Patch* patch);
    virtual void apply_ymin(ElectroMagn* EMfields, double time_dual, Patch* patch);
    virtual void apply_ymax(ElectroMagn* EMfields, double time_dual, Patch* patch);
    
 private:
    
    //! Oversize
    unsigned int oversize_;
    
    //! Number of nodes on the primal grid
    unsigned int nx_p;

    //! Number of nodes on the dual grid
    unsigned int nx_d;
    
};

#endif

