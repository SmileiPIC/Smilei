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
#include "PusherBorisBTIS3.h"
#include "PusherPonderomotiveBoris.h"
#include "PusherPonderomotivePositionBoris.h"
#include "PusherPonderomotiveBorisBTIS3.h"
#include "PusherVay.h"
#include "PusherBorisNR.h"
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
class PusherFactory
{
public:
    //  --------------------------------------------------------------------------------------------------------------------
    //! Create appropriate pusher for the species ispec
    //! \param ispec SpeciesId
    //! \param params Parameters
    //  --------------------------------------------------------------------------------------------------------------------
    static Pusher *create( Params &params, Species *species )
    {
        Pusher *Push = NULL;
        
        // Particle of matter
        if( species->mass_ > 0 ) {
            // assign the correct Pusher to Push
            // Pusher of Boris
            if( species->pusher_name_ == "boris" ) {
                    Push = new PusherBoris( params, species );
            } else if( species->pusher_name_ == "ponderomotive_boris" ) {
            
                int n_envlaser = params.Laser_Envelope_model;
                if( n_envlaser <1 ) {
                    ERROR_NAMELIST( "No Laser Envelope present." 
                                <<  " The pusher ponderomotive_boris can be used only in presence of a Laser Envelope.",
                    LINK_NAMELIST + std::string("#pusher"));
                }
                
                Push = new PusherPonderomotiveBoris( params, species );
            // Non-relativistic Boris pusher
            } else if( species->pusher_name_ == "borisnr" ) {
                Push = new PusherBorisNR( params, species );
            }
            // Pusher of J.L. Vay
            else if( species->pusher_name_ == "vay" ) {
                Push = new PusherVay( params, species );
            // Pusher of Higuera Cary
            } else if( species->pusher_name_ == "higueracary" ) {
                Push = new PusherHigueraCary( params, species );
            } else if (species->pusher_name_ == "borisBTIS3"){
                if (!params.use_BTIS3){
                    ERROR("Pusher borisBTIS3 can be used only if use_BTIS3 = True in Main block");
                } else {
                    Push = new PusherBorisBTIS3( params, species );
                }
            } else if (species->pusher_name_ == "ponderomotive_borisBTIS3"){
                if (!params.use_BTIS3){
                    ERROR("Pusher ponderomotive_borisBTIS3 can be used only if use_BTIS3 = True in Main block");
                }
                int n_envlaser = params.Laser_Envelope_model;
                if( n_envlaser <1 ) {
                    ERROR_NAMELIST( "No Laser Envelope present." 
                                <<  " The pusher ponderomotive_borisBTIS3 can be used only in presence of a Laser Envelope.",
                    LINK_NAMELIST + std::string("#pusher"));
                }
                if (params.use_BTIS3 && (n_envlaser >=1)) {
                    Push = new PusherPonderomotiveBorisBTIS3( params, species );
                }
            } else {
                ERROR_NAMELIST( "For species " << species->name_
                       << ": unknown pusher `"
                       << species->pusher_name_ << "`",
                    LINK_NAMELIST + std::string("#pusher") );
            }
        }
        // Photon
        else if( species->mass_ == 0 ) {
            if( species->pusher_name_ == "norm" ) {
                Push = new PusherPhoton( params, species );
            } else {
                ERROR_NAMELIST( "For photon species " << species->name_
                       << ": unknown pusher `"
                       << species->pusher_name_ << "`",
                   LINK_NAMELIST + std::string("#pusher") );
            }
        }
        if (params.use_BTIS3){
            if( !params.Laser_Envelope_model && species->pusher_name_ != "borisBTIS3" ) {
                ERROR_NAMELIST( "For species " << species->name_ 
                << " the only pusher available with BTIS3 interpolation scheme is borisBTIS3",
                LINK_NAMELIST + std::string("#pusher") );
            }
            if( params.Laser_Envelope_model && species->pusher_name_ != "ponderomotive_borisBTIS3" ) {
                ERROR_NAMELIST( "For species " << species->name_ 
                << " the only pusher available with BTIS3 interpolation scheme with envelope is ponderomotive_borisBTIS3",
                LINK_NAMELIST + std::string("#pusher") );
            }
            
        }
        if( params.Laser_Envelope_model && !params.use_BTIS3) {
            if( species->pusher_name_ != "ponderomotive_boris" ) {
                ERROR_NAMELIST( "For species " << species->name_ 
                << " the only pusher available to interact with the envelope without BTIS3 interpolation scheme is ponderomotive_boris",
                LINK_NAMELIST + std::string("#pusher") );
            }
        }
        return Push;
    }
    
    static Pusher *create_ponderomotive_position_updater( Params &params, Species *species )
    {
        Pusher *Push_ponderomotive_position = NULL;
        
        // Particle of matter
        if( species->mass_ > 0 ) {
            // assign the correct Pusher to Push_ponderomotive_position
            if( (species->pusher_name_ == "ponderomotive_boris") || (species->pusher_name_ == "ponderomotive_borisBTIS3") ) {
                    Push_ponderomotive_position = new PusherPonderomotivePositionBoris( params, species );
            }
            
            else {
                ERROR_NAMELIST( "For species " << species->name_
                       << ": unknown pusher `"
                       << species->pusher_name_ << "`",
                   LINK_NAMELIST + std::string("#pusher") );
            }
        } else {
            ERROR_NAMELIST( "Ponderomotive pusher is not a valid choice for photons",
            LINK_NAMELIST + std::string("#pusher"));
        }
        return Push_ponderomotive_position;
    }
    
};

#endif
