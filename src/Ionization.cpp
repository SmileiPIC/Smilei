#include "Ionization.h"

Ionization::Ionization(PicParams *params, int ispec) {

    wavelength_SI        = params->wavelength_SI;
    
	dt                   = params->timestep;
    nDim_field           = params->nDim_field;
	nDim_particle        = params->nDim_particle;
	atomic_number_       = params->species_param[ispec].atomic_number;
    ionized_species_mass = params->species_param[ispec].mass;
	
    // Normalization constant from Smilei normalization to/from atomic units
    eV_to_au = 1.0 / 27.2116;
    EC_to_au = 6.24381e-6 / wavelength_SI;
    au_to_w0 = 2.19475e7  * wavelength_SI;  //wavelength_SI / 21.9465e6;

    // Ionization potential & quantum numbers (all in atomic units 1 au = 27.2116 eV)
	Potential.resize(atomic_number_);
	Azimuthal_quantum_number.resize(atomic_number_);
	switch (atomic_number_) {
		case 1:
			Potential[0] = 13.6 * eV_to_au;
			Azimuthal_quantum_number[0]=0;
			break;
		case 2:
			Potential[0] = 24.98*eV_to_au;
			Potential[1] = 54.33*eV_to_au;
			Azimuthal_quantum_number[0] = 0;
			Azimuthal_quantum_number[1] = 0;
			break;
		default:
			break;
	}
	
	for (int i=0;i<atomic_number_; i++) {
		DEBUG(5,"i " << i << " potential: " << Potential[i] << " Az.q.num: " << Azimuthal_quantum_number[i]);
	}
}


Ionization::~Ionization() {
}
