#ifndef PUSHERFACTORY_H
#define PUSHERFACTORY_H

#include "Pusher.h"
#include "PusherBoris.h"

#include "Params.h"
#include "Species.h"

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
    static Pusher* create(Params& params, SpeciesStructure& sparams) {
        Pusher* Push = NULL;

        // assign the correct Pusher to Push
        if ( sparams.dynamics_type == "norm" )
            Push = new PusherBoris( params, sparams );
        else
            ERROR( "For species " << sparams.species_type << ": unknown dynamics_type `" << sparams.dynamics_type);

        return Push;
    }

};

#endif
