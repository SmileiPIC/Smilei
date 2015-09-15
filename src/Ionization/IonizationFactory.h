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
    static Ionization* create(Params& params, SpeciesStructure& sparams, double max_charge) {
        Ionization* Ionize = NULL;
        std::string model=sparams.ionization_model;

        if ( model == "tunnel" ) {
            if (max_charge > (int)sparams.atomic_number)
                ERROR( "Charge > atomic_number for species " << sparams.species_type );

            Ionize = new IonizationTunnel( params, sparams );

        } else if ( model != "none" ) {
            WARNING( "For species " << sparams.species_type << ": unknown ionization model `" << model << "` ... assuming no ionization");
        }
        return Ionize;
    }

};

#endif
