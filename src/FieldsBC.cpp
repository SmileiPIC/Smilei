#include "FieldsBC.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "PicParams.h"
#include "Laser.h"
#include "Tools.h"

using namespace std;

FieldsBC::FieldsBC( PicParams *params )
{

    // check for laser conditions
    laser_.resize(params->n_laser);

    for (unsigned int i=0; i<laser_.size(); i++) {
        DEBUG(5,"Initializing Laser "<<i);
        laser_[i] = new Laser(params->sim_time, params->laser_param[i]);
    }

    double dt = params->timestep;

}

FieldsBC::~FieldsBC()
{
    for (unsigned int i=0; i< laser_.size(); i++) {
        delete laser_[i];
    }


}
