// ----------------------------------------------------------------------------
//! \file MultiphotonBreitWheelerFactory.h
//
//! \brief This file contains the class MultiphotonBreitWheelerFactory that
// manages the MultiphotonBreitWheeler process initialization.
//
// ----------------------------------------------------------------------------

#ifndef MULTIPHOTONBREITWHEELERFACTORY_H
#define MULTIPHOTONBREITWHEELERFACTORY_H

#include "MultiphotonBreitWheeler.h"
#include "Params.h"
#include "Tools.h"

//  ----------------------------------------------------------------------------
//! Class MultiphotonBreitWheelerFactory
//
//  ----------------------------------------------------------------------------

class MultiphotonBreitWheelerFactory
{

public:

    //  --------------------------------------------------------------------
    //! Create appropriate multiphoton Breit Wheeler model
    //! for the species ispec
    //! \param params Parameters
    //! \param species species object
    //  --------------------------------------------------------------------
    static MultiphotonBreitWheeler *create( Params &params,
                                            Species *species,
                                            Random * rand)
    {
        MultiphotonBreitWheeler *Multiphoton_Breit_Wheeler_process = NULL;
        
        // Assign the correct Radiation model to Radiate
        if( !species->multiphoton_Breit_Wheeler_[0].empty() ) {
            Multiphoton_Breit_Wheeler_process = new MultiphotonBreitWheeler( params, species, rand );
            if( params.Laser_Envelope_model & species->ponderomotive_dynamics ) {
                ERROR( "Multiphoton Breit-Wheeler model is not yet implemented for species interacting with Laser Envelope model." );
            }
        }
        
        return Multiphoton_Breit_Wheeler_process;
        
    }
    
};

#endif
