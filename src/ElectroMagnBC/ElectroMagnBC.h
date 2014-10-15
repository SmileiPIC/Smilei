
#ifndef ELECTROMAGNBC_H
#define ELECTROMAGNBC_H

#include <vector>

class PicParams;
class LaserParams;
class SmileiMPI;
class ElectroMagn;
class Laser;

class ElectroMagnBC {
public:
    ElectroMagnBC( PicParams &params,  LaserParams &laser_params );
    ~ElectroMagnBC();

    virtual void apply(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi) = 0;
    void laserDisabled();

 protected:

    //! Vector for the various lasers
    std::vector<Laser*> laser_;
    
    //! time-step
    double dt;

};

#endif

