// --------------------------------------------------------------------------------------------------------------------
//
//! \file PusherFactory.h
//
//! \brief Class PusherFactory that manage the pusher association to the species.
//
// --------------------------------------------------------------------------------------------------------------------

#ifndef PUSHERFACTORY_H
#define PUSHERFACTORY_H

#include "Pusher.h"
#include "PusherBoris.h"
#include "PusherVay.h"
#include "PusherBorisNR.h"
#include "PusherRRLL.h"
#include "PusherHigueraCary.h"
#include "PusherPhoton.h"

#include "Params.h"
#include "Species.h"

#include "Tools.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherFactory
//
//! \brief Class PusherFactory that manage the pusher association to the species.
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

        // Particle of matter
        if (species->mass > 0) {
            // assign the correct Pusher to Push
            if ( species->dynamics_type == "norm")
            {
                if (species->mass > 0)
                {
                    Push = new PusherBoris( params, species );
                }
                else
                {
                    Push = new PusherPhoton( params, species );
                }
            }
            else if ( species->dynamics_type == "borisnr" )
            {
                Push = new PusherBorisNR( params, species );
            }
            /*else if ( species->dynamics_type == "rrll" )
            {
                Push = new PusherRRLL( params, species );
            }*/
            else if ( species->dynamics_type == "vay" )
            {
                Push = new PusherVay( params, species );
            }
            else if ( species->dynamics_type == "higueracary" )
            {
                Push = new PusherHigueraCary( params, species );
            }
            else {
                ERROR( "For species " << species->species_type
                                      << ": unknown dynamics_type `"
                                      << species->dynamics_type << "`");
            }
        }
        // Photon
        else if (species->mass == 0)
        {
            if ( species->dynamics_type == "norm")
            {
                Push = new PusherPhoton( params, species );
            }
            else {
                ERROR( "For photon species " << species->species_type
                                      << ": unknown dynamics_type `"
                                      << species->dynamics_type << "`");
            }
        }

        return Push;
    }

};

#endif
