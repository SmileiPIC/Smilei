#ifndef PUSHERFACTORY_H
#define PUSHERFACTORY_H

#include "Pusher.h"
#include "PusherBoris.h"

#include "PicParams.h"

#include "Tools.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherFactory

//  --------------------------------------------------------------------------------------------------------------------
class PusherFactory {
public:
    //  --------------------------------------------------------------------------------------------------------------------
    //! Create appropriate pusher for the species ispec
    //! \param ispec SpeciesId
    //! \param params Parameters
    //  --------------------------------------------------------------------------------------------------------------------
    static Pusher* create(PicParams& params, int ispec) {
        Pusher* Push = NULL;

        // assign the correct Pusher to Push
        if ( params.species_param[ispec].dynamics_type == "norm" )
            Push = new PusherBoris( params, ispec );
        else
            ERROR( "Unknown dynamics : " << params.species_param[ispec].dynamics_type );

        return Push;
    }

};

#endif
