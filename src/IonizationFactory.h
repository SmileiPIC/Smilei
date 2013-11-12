#ifndef IonizationFactory_H
#define IonizationFactory_H

#include "Ionization.h"
#include "IonizationTunnel.h"

#include "PicParams.h"

#include "Tools.h"

class IonizationFactory {
public:
	static Ionization* create(PicParams* params, int ispec) {
		Ionization* Ionize = NULL;
		std::string model=params->species_param[ispec].ionization_model;
		
		if ( model.empty() || model == "none") {
			Ionize=NULL;
		} else if ( model == "tunnel" ) {
            DEBUG("0 IONIZE IONIZE IONIZE IONIZE IONIZE IONIZE "<< Ionize);
			if (params->species_param[ispec].charge > params->species_param[ispec].atomic_number)
				ERROR( "Charge > atomic_number for specie " << ispec );

		    Ionize = new IonizationTunnel( params, ispec );
            DEBUG("1 IONIZE IONIZE IONIZE IONIZE IONIZE IONIZE "<< Ionize);
            
	    } else {
		    ERROR( "Unknown Ionization Model : " << model );
		}
		
        DEBUG("2 IONIZE IONIZE IONIZE IONIZE IONIZE IONIZE "<< Ionize);
		return Ionize;
	}
	
};

#endif
