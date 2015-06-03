
#ifndef ELECTROMAGNBC_H
#define ELECTROMAGNBC_H

#include <vector>

class PicParams;
class LaserParams;
class SmileiMPI;
class ElectroMagn;
class LaserProfile;

class ElectroMagnBC {
public:
    ElectroMagnBC( PicParams &params,  LaserParams &laser_params );
    ~ElectroMagnBC();

    virtual void apply_xmin(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi) = 0;
    virtual void apply_xmax(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi) = 0;
    virtual void apply_ymin(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi) = 0;
    virtual void apply_ymax(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi) = 0;
    void laserDisabled();

 protected:

    //! Vector for the various lasers
    std::vector<LaserProfile*> laser_;
    
    //! time-step
    double dt;

};

#endif

