
#ifndef ENVELOPEBC3D_refl_H
#define ENVELOPEBC3D_refl_H


#include <vector>
#include "Tools.h"
#include "EnvelopeBC.h"
#include "LaserEnvelope.h"
#include "Field3D.h"
#include "cField3D.h"

class Params;
class LaserEnvelope;
//class Field;

class EnvelopeBC3D_refl : public EnvelopeBC
{
public:

    EnvelopeBC3D_refl( Params &params, Patch *patch, unsigned int _min_max );
    ~EnvelopeBC3D_refl() {};
    
    void apply( LaserEnvelope *envelope, double time_dual, Patch *patch ) override;
    
    
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
    
    //! Number of nodes on the primal grid in the y-direction
    unsigned int nz_p;
    
    //! Number of nodes on the dual grid in the z-direction
    unsigned int nz_d;
    
};

#endif

