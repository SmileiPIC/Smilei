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
#include "RadiationMonteCarlo.h"
#include "RadiationCorrLandauLifshitz.h"
#include "RadiationLandauLifshitz.h"
#include "RadiationDiagRadiationSpectrum.h"
#include "RadiationNiel.h"

#include "Params.h"
#include "Species.h"

#include "Tools.h"

//  ----------------------------------------------------------------------------
//! Class RadiationFactory
//
//  ----------------------------------------------------------------------------

class RadiationFactory
{
public:
    //  ------------------------------------------------------------------------
    //! Create appropriate radiation model for the species `species`
    //! \param species Species object
    //! \param params Parameters
    //  --------------------------------------------------------------------------------------------------------------------
    static Radiation *create( Params &params, Species *species, Random * rand   )
    {
        Radiation *Radiate = NULL;

        // assign the correct Radiation model to Radiate
        if( species->radiation_model_ == "mc" ) {
            Radiate = new RadiationMonteCarlo( params, species, rand  );
        }
        // Corrected LL + stochastic diffusive operator
        else if( species->radiation_model_ == "niel" ) {
            Radiate = new RadiationNiel( params, species, rand  );
        }
        // Corrected continuous radiation loss model
        else if( species->radiation_model_ == "cll" ) {
            Radiate = new RadiationCorrLandauLifshitz( params, species, rand  );
        }
        // Classical continuous radiation loss model from Landau-Lifshitz (LL)
        else if( species->radiation_model_ == "ll" ) {
            Radiate = new RadiationLandauLifshitz( params, species, rand );
        }
        // Radiation is only a diagnostic (DiagRadiationSpectrum can be called for this species): only compute
        // the electron quantum parameter
        else if( species->radiation_model_ == "diagradiationspectrum" ) {
            Radiate = new RadiationDiagRadiationSpectrum( params, species, rand );
        } else if( species->radiation_model_ != "none" ) {
            ERROR( "For species " << species->name_
                   << ": unknown radiation_model `"
                   << species->radiation_model_ << "`" );
        }

        int n_envlaser = params.Laser_Envelope_model;
        if( ( n_envlaser >=1 ) & ( species->ponderomotive_dynamics ) & ( species->radiation_model_ != "none" ) ) {
            ERROR( "Radiation model is not yet implemented for species interacting with Laser Envelope model." );
        }


        return Radiate;
    }

};

#endif
