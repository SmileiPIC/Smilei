#include "Ionization.h"

Ionization::Ionization(PicParams *params, int ispec) {
	int 	atomic_number   = params->species_param[ispec].atomic_number;

	DEBUG(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> "<< atomic_number);
	
}


Ionization::~Ionization() {
}
