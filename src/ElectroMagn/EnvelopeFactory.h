
#ifndef ENVELOPEFACTORY_H
#define ENVELOPEFACTORY_H

#include "LaserEnvelope.h"
#include "Params.h"
#include "Patch.h"

class EnvelopeFactory
{
public:
    static LaserEnvelope *create( Params &params, Patch *patch )
    {
        if( params.geometry == "1Dcartesian" ) {
            return new LaserEnvelope1D( params, patch );
        } else if( params.geometry == "2Dcartesian" ) {
            return new LaserEnvelope2D( params, patch );
        } else if( params.geometry == "3Dcartesian" ) {
            return new LaserEnvelope3D( params, patch );
        } else if( params.geometry == "AMcylindrical" ) {
            if (params.nmodes!=1){
                WARNING("The envelope will be modeled with only 1 azimuthal mode, \n while the electromagentic fields will use all the modes"); 
            }
            return new LaserEnvelopeAM( params, patch );
        } else {
            return NULL;
        }
    }
    static LaserEnvelope *clone( LaserEnvelope *envelope, Patch *patch, Params &params, unsigned int n_moved )
    {
        if( envelope == NULL ) {
            return NULL;
        }
        
        if( dynamic_cast<LaserEnvelope1D *>( envelope ) ) {
            return new LaserEnvelope1D( envelope, patch, params, n_moved );
        } else if( dynamic_cast<LaserEnvelope2D *>( envelope ) ) {
            return new LaserEnvelope2D( envelope, patch, params, n_moved );
        } else if( dynamic_cast<LaserEnvelope3D *>( envelope ) ) {
            return new LaserEnvelope3D( envelope, patch, params, n_moved );
        } else if( dynamic_cast<LaserEnvelopeAM *>( envelope ) ) {
            return new LaserEnvelopeAM( envelope, patch, params, n_moved );
        } else {
            return NULL;
        }
    }
};

#endif
