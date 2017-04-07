
#ifndef ELECTROMAGNBC_H
#define ELECTROMAGNBC_H

#include <vector>

#include "Patch.h"

class Params;
class Patch;
class ElectroMagn;
class Laser;
class Field;

class ElectroMagnBC {
public:
    ElectroMagnBC( Params &params, Patch* patch );
    virtual ~ElectroMagnBC();
    void clean();
    
    virtual void apply_xmin(ElectroMagn* EMfields, double time_dual, Patch* patch) = 0;
    virtual void apply_xmax(ElectroMagn* EMfields, double time_dual, Patch* patch) = 0;
    virtual void apply_ymin(ElectroMagn* EMfields, double time_dual, Patch* patch) = 0;
    virtual void apply_ymax(ElectroMagn* EMfields, double time_dual, Patch* patch) = 0;
    virtual void apply_zmin(ElectroMagn* EMfields, double time_dual, Patch* patch) = 0;
    virtual void apply_zmax(ElectroMagn* EMfields, double time_dual, Patch* patch) = 0;
    void laserDisabled();
    
    virtual void save_fields_BC1D(Field*) {}
    virtual void save_fields_BC2D_Long(Field*) {}
    virtual void save_fields_BC2D_Trans(Field*) {}
    virtual void save_fields_BC3D_Long(Field*) {}
    virtual void save_fields_BC3D_TransY(Field*) {}
    virtual void save_fields_BC3D_TransZ(Field*) {}
    
    //! Vector for the various lasers
    std::vector<Laser*> vecLaser;
    
protected:

    //! time-step
    double dt;

};

#endif

