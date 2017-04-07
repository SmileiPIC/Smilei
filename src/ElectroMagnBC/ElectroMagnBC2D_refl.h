
#ifndef ELECTROMAGNBC2D_refl_H
#define ELECTROMAGNBC2D_refl_H

#include "ElectroMagnBC.h"

class Params;
class ElectroMagn;

class ElectroMagnBC2D_refl : public ElectroMagnBC {
public:
    ElectroMagnBC2D_refl( Params &params, Patch* patch, unsigned int _min_max );
    ~ElectroMagnBC2D_refl() {};
    
    void apply(ElectroMagn* EMfields, double time_dual, Patch* patch) override;

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

