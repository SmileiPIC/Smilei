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
#include "RadiationNlcisMC.h"
#include "RadiationNlcisCont.h"

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
        if ( species->radiation_type == "Monte-Carlo" )
        {
            //Radiate = new RadiationNlcisMC( params, species );
        }
        // Monte-Carlo nonlinear inverse Compton scattering
        else if ( species->radiation_type == "continuous" )
        {
            //Radiate = new RadiationNlcisCont( params, species );
        }
        else if ( species->radiation_type != "none" )
        {
            ERROR( "For species " << species->species_type << ": unknown radiation_type `" << species->radiation_type << "`");
        }

        return Radiate;
    }

};

#endif
