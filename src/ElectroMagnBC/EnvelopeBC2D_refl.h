
#ifndef ENVELOPEBC2D_refl_H
#define ENVELOPEBC2D_refl_H


#include <vector>
#include "Tools.h"
#include "EnvelopeBC.h"
#include "LaserEnvelope.h"
#include "Field2D.h"
#include "cField2D.h"
#include "ElectroMagn.h"

class Params;
class LaserEnvelope;
//class Field;

class EnvelopeBC2D_refl : public EnvelopeBC
{
public:

    EnvelopeBC2D_refl( Params &params, Patch *patch, unsigned int _min_max );
    ~EnvelopeBC2D_refl() {};
    
    void apply( LaserEnvelope *envelope, ElectroMagn *EMfields, double time_dual, Patch *patch ) override;
    
    
private:

    //! Oversize (nb of ghost cells)
    unsigned int oversize_;
    
    //! Number of nodes on the primal grid in the x-direction
    unsigned int nx_p;
    
    //! Number of nodes on the dual grid in the x-direction
    unsigned int nx_d;
    
    //! Number of nodes on the primal grid in the y-direction
    unsigned int ny_p;
    
    //! Number of nodes on the dual grid in the y-direction
    unsigned int ny_d;
    
    
};

#endif

