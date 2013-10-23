#include "Ionization.h"

Ionization::Ionization(PicParams *params, int ispec) {

	dt      = params->timestep;
	nDim_   = params->nDim_particle;
	atomic_number_   =	params->species_param[ispec].atomic_number;
	

	Potential.resize(atomic_number_);
	azimuthal_quantum_number.resize(atomic_number_);
	switch (atomic_number_) {
		case 1:
			Potential[0]=13.6;
			azimuthal_quantum_number[0]=0;
			break;
		case 2:
			Potential[0]=24.98;
			Potential[1]=54.33;
			azimuthal_quantum_number[0]=0;
			azimuthal_quantum_number[1]=0;
			break;
		default:
			break;
	}
	
	for (int i=0;i<atomic_number_; i++) {
		DEBUG(5,"i " << i << " potential: " << Potential[i] << " Az.q.num: " << azimuthal_quantum_number[i]);
	}
}


Ionization::~Ionization() {
}
