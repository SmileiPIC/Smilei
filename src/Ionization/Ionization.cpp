#include "Ionization.h"

Ionization::Ionization(PicParams& params, int ispec) {

    wavelength_SI        = params.wavelength_SI;

    dt                   = params.timestep;
    nDim_field           = params.nDim_field;
    nDim_particle        = params.nDim_particle;
    atomic_number_       = params.species_param[ispec].atomic_number;
    ionized_species_mass = params.species_param[ispec].mass;

    // Normalization constant from Smilei normalization to/from atomic units
    eV_to_au = 1.0 / 27.2116;
    EC_to_au = 6.24381e-6 / wavelength_SI;
    au_to_w0 = 2.19475e7  * wavelength_SI;  //wavelength_SI / 21.9465e6;

    // Ionization potential & quantum numbers (all in atomic units 1 au = 27.2116 eV)
    Potential.resize(atomic_number_);
    Azimuthal_quantum_number.resize(atomic_number_);
    switch (atomic_number_) {
    case 1:
        Potential[0] = 13.61 * eV_to_au;
        Azimuthal_quantum_number[0]=0;
        break;
    case 2:
        Potential[0] = 24.98*eV_to_au;
        Potential[1] = 54.33*eV_to_au;
        Azimuthal_quantum_number[0] = 0;
        Azimuthal_quantum_number[1] = 0;
        break;
    case 3:
        Potential[0] = 	5.34*eV_to_au;
        Potential[1] = 74.45*eV_to_au;
        Potential[2] = 121.9*eV_to_au;
        Azimuthal_quantum_number[0] = 0;
        Azimuthal_quantum_number[1] = 0;
        Azimuthal_quantum_number[2] = 0;
    case 4:
        Potential[0] = 	8.41*eV_to_au;
        Potential[1] = 18.69*eV_to_au;
        Potential[2] = 149.3*eV_to_au;
        Potential[3] = 214.9*eV_to_au;
        Azimuthal_quantum_number[0] = 0;
        Azimuthal_quantum_number[1] = 0;
        Azimuthal_quantum_number[2] = 0;
        Azimuthal_quantum_number[3] = 0;
    case 5:
        Potential[0] = 	8.43*eV_to_au;
        Potential[1] = 25.81*eV_to_au;
        Potential[2] = 39.58*eV_to_au;
        Potential[3] = 249.3*eV_to_au;
        Potential[4] = 332.9*eV_to_au;
        Azimuthal_quantum_number[0] = 1;
        Azimuthal_quantum_number[1] = 0;
        Azimuthal_quantum_number[2] = 0;
        Azimuthal_quantum_number[3] = 0;
        Azimuthal_quantum_number[4] = 0;
    case 6:
        Potential[0] = 	10.6*eV_to_au;
        Potential[1] = 26.05*eV_to_au;
        Potential[2] = 50.40*eV_to_au;
        Potential[3] = 67.59*eV_to_au;
        Potential[4] = 374.2*eV_to_au;
        Potential[5] = 475.6*eV_to_au;
        Azimuthal_quantum_number[0] = 1;
        Azimuthal_quantum_number[1] = 1;
        Azimuthal_quantum_number[2] = 0;
        Azimuthal_quantum_number[3] = 0;
        Azimuthal_quantum_number[4] = 0;
        Azimuthal_quantum_number[5] = 0;
    case 7:
        Potential[0] = 	13.37*eV_to_au;
        Potential[1] = 32.02*eV_to_au;
        Potential[2] = 51.55*eV_to_au;
        Potential[3] = 82.57*eV_to_au;
        Potential[4] = 103.1*eV_to_au;
        Potential[5] = 524.0*eV_to_au;
        Potential[6] = 643.3*eV_to_au;
        Azimuthal_quantum_number[0] = 1;
        Azimuthal_quantum_number[1] = 1;
        Azimuthal_quantum_number[2] = 1;
        Azimuthal_quantum_number[3] = 0;
        Azimuthal_quantum_number[4] = 0;
        Azimuthal_quantum_number[5] = 0;
        Azimuthal_quantum_number[6] = 0;
    default:
        break;
    }

    for (unsigned int i=0; i<atomic_number_; i++) {
        DEBUG(5,"i " << i << " potential: " << Potential[i] << " Az.q.num: " << Azimuthal_quantum_number[i]);
    }
}


Ionization::~Ionization() {
}
