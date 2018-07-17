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
                ERROR( "Charge > atomic_number for species " << species->name );
            if( species->particles->is_test )
                ERROR( "Cannot ionize test species " << species->name );
            
            Ionize = new IonizationTunnel( params, species );

            int n_envlaser = PyTools::nComponents("LaserEnvelope");
            if ( (n_envlaser >=1) & (species->ponderomotive_dynamics) ){
                ERROR( "Ionization is not yet implemented for species interacting with Laser Envelope model.");
                                                                       }
                       
        } else if ( model != "none" ) {
            WARNING( "For species " << species->name << ": unknown ionization model `" << model << "` ... assuming no ionization");
        }
        return Ionize;
    }

};

#endif
