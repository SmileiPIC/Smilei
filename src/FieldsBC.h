
#ifndef FIELDSBC_H
#define FIELDSBC_H

#include <vector>

class PicParams;
class LaserParams;
class SmileiMPI;
class ElectroMagn;
class Laser;

class FieldsBC {
public:
    FieldsBC( PicParams &params,  LaserParams &laser_params );
    ~FieldsBC();

    virtual void apply(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi) = 0;
    void laserDisabled();

 protected:

    //! Vector for the various lasers
    std::vector<Laser*> laser_;
    
    //! time-step
    double dt;

};

#endif

