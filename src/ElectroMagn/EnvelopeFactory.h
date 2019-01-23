
#ifndef ENVELOPEFACTORY_H
#define ENVELOPEFACTORY_H

#include "LaserEnvelope.h"
#include "Params.h"
#include "Patch.h"

class EnvelopeFactory {
public:
    static LaserEnvelope* create( Params& params, Patch* patch, ElectroMagn* EMfields ) {
        if ( params.geometry == "1Dcartesian" ) {
            return new LaserEnvelope1D( params, patch, EMfields );
        } else if ( params.geometry == "2Dcartesian" ) {
            return new LaserEnvelope2D( params, patch, EMfields );
        } else if ( params.geometry == "3Dcartesian" ) {
            return new LaserEnvelope3D( params, patch, EMfields );
        }
        else
            return NULL;
    }
    static LaserEnvelope* clone( LaserEnvelope *envelope, Patch* patch, ElectroMagn* EMfields, Params& params, unsigned int n_moved ) {
        if  (envelope == NULL)
            return NULL;

        if ( dynamic_cast<LaserEnvelope1D*>( envelope ) ) {
            return new LaserEnvelope1D( envelope, patch , EMfields, params, n_moved );
        } else if ( dynamic_cast<LaserEnvelope2D*>( envelope ) ) {
            return new LaserEnvelope2D( envelope, patch , EMfields, params, n_moved );
        } else if ( dynamic_cast<LaserEnvelope3D*>( envelope ) ) {
            return new LaserEnvelope3D( envelope, patch , EMfields, params, n_moved );
        }
        else
            return NULL;
    }
};

#endif
