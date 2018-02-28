
#ifndef ELECTROMAGNBC3D_refl_H
#define ELECTROMAGNBC3D_refl_H

#include "ElectroMagnBC.h"

class Params;
class ElectroMagn;

class ElectroMagnBC3D_refl : public ElectroMagnBC {
public:
    ElectroMagnBC3D_refl( Params &params, Patch* patch, unsigned int _min_max );
    ~ElectroMagnBC3D_refl() {};
    
    void apply(ElectroMagn* EMfields, double time_dual, Patch* patch) override;

private:
    
    //! Oversize (nb of ghost cells)
    unsigned int oversize_x, oversize_y, oversize_z;
    
    //! Number of nodes on the primal grid in the x-direction
    unsigned int nx_p;
    
    //! Number of nodes on the dual grid in the x-direction
    unsigned int nx_d;
    
    //! Number of nodes on the primal grid in the y-direction
    unsigned int ny_p;
    
    //! Number of nodes on the dual grid in the y-direction
    unsigned int ny_d;
    
    //! Number of nodes on the primal grid in the z-direction
    unsigned int nz_p;
    
    //! Number of nodes on the dual grid in the z-direction
    unsigned int nz_d;
    
};

#endif

