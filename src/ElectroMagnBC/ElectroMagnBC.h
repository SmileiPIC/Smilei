
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
    ElectroMagnBC( Params &params, Patch *patch, unsigned int _min_max );
    virtual ~ElectroMagnBC();
    void clean();
    
    virtual void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) = 0;
    
    void laserDisabled();
    
    virtual void save_fields( Field *, Patch *patch ) {};
    virtual void disableExternalFields() {};
    
    //! Vector for the various lasers
    std::vector<Laser *> vecLaser;
    
protected:

    //! time-step
    double dt;
    
    // side of BC is applied 0:xmin 1:xmax 2:ymin 3:ymax 4:zmin 5:zmax
    unsigned int min_max;

    // number of damping cells
    std::vector<unsigned int> number_of_damping_cells;
    
};

#endif

