
#ifndef ELECTROMAGNBC_H
#define ELECTROMAGNBC_H

#include <vector>

#include "Patch.h"

class Params;
class Patch;
class ElectroMagn;
class Laser;
class Field;

class ElectroMagnBC
{
public:
    ElectroMagnBC( Params &params, Patch *patch, unsigned int i_boundary );
    virtual ~ElectroMagnBC();
    void clean();
    
    virtual void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) = 0;
    
    void laserDisabled();
    
    virtual void save_fields( Field *, Patch *patch ) {};
    virtual void disableExternalFields() {};
    
    //! Vector for the various lasers
    std::vector<Laser *> vecLaser;
    
protected:
    
    // side of BC is applied 0:xmin 1:xmax 2:ymin 3:ymax 4:zmin 5:zmax
    unsigned int i_boundary_;
    
    //! Number of nodes on the primal grid
    std::vector<unsigned int> n_p;
    
    //! Number of nodes on the dual grid
    std::vector<unsigned int> n_d;
    
    //! time-step
    double dt;
    
    //! Spatial step
    std::vector<double> d;
    
    //! Ratio of the time-step by the spatial-step
    std::vector<double> dt_ov_d;
};

#endif

