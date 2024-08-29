#ifndef IonizationFactory_H
#define IonizationFactory_H

#include "Ionization.h"
#include "IonizationFromRate.h"
#include "IonizationTunnelEnvelopeAveraged.h"
#include "Params.h"
#include "Species.h"
#include "IonizationTunnel.h"
#include "Tools.h"

//! this class create and associate the right ionization model to species
class IonizationFactory
{
   public:
    static Ionization *create(Params &params, Species *species)
    {
        Ionization *Ionize = NULL;
        std::string model = species->ionization_model_;

        if( model == "tunnel" ) {
            checkMaxCharge(species);
            checkNotLaserEnvelopeModel(params);
            Ionize = new IonizationTunnel<0>( params, species ); // The original model included in Smilei
            
        } else if( model == "tunnel_envelope_averaged" ) {
            checkMaxCharge(species);
            checkTestParticle(species);
            if ( !params.Laser_Envelope_model ) {
                ERROR( "The ionization model tunnel_envelope_averaged needs a laser envelope");
            }

            Ionize = new IonizationTunnelEnvelopeAveraged( params, species );

         } else if( model == "from_rate" ) {
            if ( species->max_charge_ > ( int ) species->maximum_charge_state_ ) {
                ERROR( "For species '" << species->name_ << ": charge > maximum_charge_state" );
            }

            Ionize = new IonizationFromRate( params, species );

        } else if (model == "tunnel_full_PPT") {  
            checkMaxCharge(species);
            checkNotLaserEnvelopeModel(params);
            Ionize = new IonizationTunnel<1>(params, species); // FullPPT
        } else if (model == "tunnel_TL") {  
            checkMaxCharge(species);
            checkNotLaserEnvelopeModel(params);
            Ionize = new IonizationTunnel<2>(params, species); // Tong&Ling
        } else if (model == "tunnel_BSI") {  
            checkMaxCharge(species);
            checkNotLaserEnvelopeModel(params);
            Ionize = new IonizationTunnel<3>(params, species); // BSI
        }

        return Ionize;
    }
           
  private:
    inline static void checkMaxCharge(const Species *species) {
        if ( species->max_charge_ > ( int )species->atomic_number_ ) {
            ERROR( "Charge > atomic_number for species " << species->name_ );
        }
    }

    inline static void checkTestParticle(const Species *species) {
        if( species->particles->is_test ) {
            ERROR( "Cannot ionize test species " << species->name_ );
        }
    }

    inline static void checkNotLaserEnvelopeModel(const Params &params) {
        if ( params.Laser_Envelope_model ) {
            ERROR( "The ionization model for species interacting with envelope is tunnel_envelope_averaged" );
        }
    }
};

#endif
