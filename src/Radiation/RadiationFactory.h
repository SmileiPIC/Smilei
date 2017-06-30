// ----------------------------------------------------------------------------
//! \file RadiationFactory.h
//
//! \brief This file contains the header for the class RadiationFactory that
// manages the different radiation models.
//
// ----------------------------------------------------------------------------

#ifndef RADIATIONFACTORY_H
#define RADIATIONFACTORY_H

#include "Radiation.h"
#include "RadiationNlicsMC.h"
#include "RadiationNlicsCont.h"
#include "RadiationLL.h"
#include "RadiationNiel.h"

#include "Params.h"
#include "Species.h"

#include "Tools.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class RadiationFactory
//
//  --------------------------------------------------------------------------------------------------------------------

class RadiationFactory {
public:
    //  --------------------------------------------------------------------------------------------------------------------
    //! Create appropriate radiation model for the species ispec
    //! \param ispec Species Id
    //! \param params Parameters
    //  --------------------------------------------------------------------------------------------------------------------
    static Radiation* create(Params& params, Species * species) {
        Radiation* Radiate = NULL;

        // assign the correct Radiation model to Radiate
        if ( species->radiation_model == "Monte-Carlo" )
        {
            Radiate = new RadiationNlicsMC( params, species );
        }
        // Corrected LL + stochastic diffusive operator
        else if ( species->radiation_model == "Niel")
        {
            Radiate = new RadiationNiel( params, species );
        }
        // Corrected continuous radiation loss model
        else if ( species->radiation_model == "continuous" )
        {
            Radiate = new RadiationNlicsCont( params, species );
        }
        // Classical continuous radiation loss model from Landau-Lifshitz (LL)
        else if ( species->radiation_model == "Landau-Lifshitz" )
        {
            Radiate = new RadiationLL( params, species );
        }
        else if ( species->radiation_model != "none" )
        {
            ERROR( "For species " << species->species_type
                                  << ": unknown radiation_model `"
                                  << species->radiation_model << "`");
        }

        return Radiate;
    }

};

#endif
