
#ifndef ENVELOPEBC1D_refl_H
#define ENVELOPEBC1D_refl_H


#include <vector>
#include "Tools.h"
#include "EnvelopeBC.h" 
#include "LaserEnvelope.h" 
#include "Field1D.h"
#include "cField1D.h"

class Params;
class LaserEnvelope;
//class Field;

class EnvelopeBC1D_refl : public EnvelopeBC {
public:

    EnvelopeBC1D_refl( Params &params, Patch* patch, unsigned int _min_max );
    ~EnvelopeBC1D_refl(){};

    void apply(LaserEnvelope* envelope, double time_dual, Patch* patch) override;
    
    
private:
    
    //! Oversize (nb of ghost cells)
    unsigned int oversize_;
    
    //! Number of nodes on the primal grid in the x-direction
    unsigned int nx_p;
    
    //! Number of nodes on the dual grid in the x-direction
    unsigned int nx_d;
  
};

#endif

