
#ifndef ELECTROMAGNBC_H
#define ELECTROMAGNBC_H

#include <vector>

#include "Patch.h"

class Params;
class Patch;
class ElectroMagn;
class Laser;
class Field;
class Solver;

class ElectroMagnBC
{
public:
    ElectroMagnBC( Params &params, Patch *patch, unsigned int i_boundary );
    virtual ~ElectroMagnBC();
    void clean();
    
    virtual void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) = 0;
    
    void laserDisabled();
    
    virtual void save_fields( Field *, Patch * ) {};
    virtual void disableExternalFields() {};
    
    //! Vector for the various lasers
    std::vector<Laser *> vecLaser;

    virtual Field* getExPML() { ERROR("Not using PML");return NULL;}
    virtual Field* getEyPML() { ERROR("Not using PML");return NULL;}
    virtual Field* getEzPML() { ERROR("Not using PML");return NULL;}
    virtual Field* getBxPML() { ERROR("Not using PML");return NULL;}
    virtual Field* getByPML() { ERROR("Not using PML");return NULL;}
    virtual Field* getBzPML() { ERROR("Not using PML");return NULL;}
    virtual Field* getDxPML() { ERROR("Not using PML");return NULL;}
    virtual Field* getDyPML() { ERROR("Not using PML");return NULL;}
    virtual Field* getDzPML() { ERROR("Not using PML");return NULL;}
    virtual Field* getHxPML() { ERROR("Not using PML");return NULL;}
    virtual Field* getHyPML() { ERROR("Not using PML");return NULL;}
    virtual Field* getHzPML() { ERROR("Not using PML");return NULL;}
    
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
    Solver* pml_solver_;
};

#endif

