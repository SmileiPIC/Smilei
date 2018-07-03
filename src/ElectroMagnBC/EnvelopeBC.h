
#ifndef ENVELOPEBC_H
#define ENVELOPEBC_H

#include <vector>

#include "Patch.h"

class Params;
class Patch;
class LaserEnvelope;
class Field;

class EnvelopeBC {
public:
    EnvelopeBC( Params &params, Patch* patch, unsigned int _min_max );
    virtual ~EnvelopeBC();
    void clean();
    
    virtual void apply(LaserEnvelope *envelope, double time_dual, Patch* patch) = 0;
    
protected:
  
    //! time-step
    double dt;

    // side of BC is applied 0:xmin 1:xmax 2:ymin 3:ymax 4:zmin 5:zmax
    unsigned int min_max;

};

#endif

