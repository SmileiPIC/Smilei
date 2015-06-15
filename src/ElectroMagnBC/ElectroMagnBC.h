
#ifndef ELECTROMAGNBC_H
#define ELECTROMAGNBC_H

#include <vector>

#include "Patch.h"

class PicParams;
class LaserParams;
class SmileiMPI;
class Patch;
class ElectroMagn;
class Laser;

class ElectroMagnBC {
public:
    ElectroMagnBC( PicParams &params,  LaserParams &laser_params );
    ~ElectroMagnBC();

    //virtual void apply(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi) = 0;
    virtual void apply(ElectroMagn* EMfields, double time_dual, Patch* patch) = 0;
    void laserDisabled();

 protected:

    //! Vector for the various lasers
    std::vector<Laser*> laser_;
    
    //! time-step
    double dt;

};

#endif

