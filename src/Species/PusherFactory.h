#ifndef PUSHERFACTORY_H
#define PUSHERFACTORY_H

#include "Pusher.h"
#include "PusherBoris.h"
#include "PusherVay.h"
#include "PusherBorisNR.h"
#include "PusherRRLL.h"
#include "PusherHigueraCary.h"

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
    static Pusher* create(Params& params, Species * species) {
        Pusher* Push = NULL;

        // assign the correct Pusher to Push
        if ( species->pusher == "boris" )
        {
            Push = new PusherBoris( params, species );
        }
        else if ( species->pusher == "borisnr" )
        {
            Push = new PusherBorisNR( params, species );
        }
        else if ( species->pusher == "rrll" )
        {
            Push = new PusherRRLL( params, species );
        }
        else if ( species->pusher == "vay" )
        {
            Push = new PusherVay( params, species );
        }
        else if ( species->pusher == "higueracary" )
        {
            Push = new PusherHigueraCary( params, species );
        }
        else {
            ERROR( "For species " << species->name << ": unknown pusher `" << species->pusher << "`");
        }

        return Push;
    }

};

#endif
