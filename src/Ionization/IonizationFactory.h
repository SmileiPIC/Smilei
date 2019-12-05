#ifndef IonizationFactory_H
#define IonizationFactory_H

#include "Ionization.h"
#include "IonizationTunnel.h"
#include "IonizationFromRate.h"

#include "Params.h"

#include "Tools.h"

#include "Species.h"

//! this class create and associate the right ionization model to species
class IonizationFactory
{
public:
    static Ionization *create( Params &params, Species *species )
    {
        Ionization *Ionize = NULL;
        std::string model=species->ionization_model;
        
        if( model == "tunnel" ) {
            
            if( species->max_charge_ > ( int )species->atomic_number_ ) {
                ERROR( "Charge > atomic_number for species " << species->name_ );
            }
            
            Ionize = new IonizationTunnel( params, species );
            
        } else if( model == "from_rate" ) {
            
            if( species->max_charge_ > ( int )species->maximum_charge_state_ ) {
                ERROR( "For species '" << species->name_ << ": charge > maximum_charge_state" );
            }
            
            Ionize = new IonizationFromRate( params, species );
            
        }
        
        return Ionize;
    }

};

#endif
