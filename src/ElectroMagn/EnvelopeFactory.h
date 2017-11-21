
#ifndef ENVELOPEFACTORY_H
#define ENVELOPEFACTORY_H

#include "LaserEnvelope.h"
#include "Params.h"
#include "Patch.h"

class EnvelopeFactory {
public:
    static LaserEnvelope* create( Params& params, Patch* patch ) {
        if ( params.geometry == "3Dcartesian" ) {
            return new LaserEnvelope3D( params, patch );
        }
        else
            return NULL;
    }
    static LaserEnvelope* clone( LaserEnvelope *envelope, Patch* patch ) {
        if  (patch->envelope == NULL)
            return NULL;
        
        if ( dynamic_cast<LaserEnvelope3D*>( patch->envelope ) ) {
            return new LaserEnvelope3D( envelope, patch );
        }
        else
            return NULL;
    }
};


#endif
