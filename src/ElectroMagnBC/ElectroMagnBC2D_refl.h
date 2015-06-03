
#ifndef ELECTROMAGNBC2D_refl_H
#define ELECTROMAGNBC2D_refl_H

#include "ElectroMagnBC.h"

class PicParams;
class ElectroMagn;

class ElectroMagnBC2D_refl : public ElectroMagnBC {
public:
    ElectroMagnBC2D_refl( PicParams &params, LaserParams &laser_params );
    ~ElectroMagnBC2D_refl();
    
    virtual void apply_xmin(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi);
    virtual void apply_xmax(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi);
    virtual void apply_ymin(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi);
    virtual void apply_ymax(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi);
    
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

