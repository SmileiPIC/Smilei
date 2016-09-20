#ifndef IonizationFactory_H
#define IonizationFactory_H

#include "Ionization.h"
#include "IonizationTunnel.h"

#include "Params.h"

#include "Tools.h"

#include "Species.h"

//! this class create and associate the right ionization model to species
class IonizationFactory {
public:
    static Ionization* create(Params& params, Species * species) {
        Ionization* Ionize = NULL;
        std::string model=species->ionization_model;
        
        if ( model == "tunnel" ) {
            if (species->max_charge > (int)species->atomic_number)
                ERROR( "Charge > atomic_number for species " << species->species_type );
            
            Ionize = new IonizationTunnel( params, species );
            
        } else if ( model != "none" ) {
            WARNING( "For species " << species->species_type << ": unknown ionization model `" << model << "` ... assuming no ionization");
        }
        return Ionize;
    }

};

#endif
