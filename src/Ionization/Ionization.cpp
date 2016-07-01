#include "Ionization.h"
#include "Species.h"
#include "IonizationTables.h"

Ionization::Ionization(Params& params, Species * species) {

    referenceAngularFrequency_SI = params.referenceAngularFrequency_SI;

    dt                   = params.timestep;
    nDim_field           = params.nDim_field;
    nDim_particle        = params.nDim_particle;
    atomic_number_       = species->atomic_number;
    ionized_species_mass = species->mass;

    // Normalization constant from Smilei normalization to/from atomic units
    eV_to_au = 1.0 / 27.2116;
    EC_to_au = 3.314742578e-15 * referenceAngularFrequency_SI; // hbar omega / (me c^2 alpha^3)
    au_to_w0 = 4.134137172e+16 / referenceAngularFrequency_SI; // alpha^2 me c^2 / (hbar omega)
    
    // Ionization potential & quantum numbers (all in atomic units 1 au = 27.2116 eV)
    Potential.resize(atomic_number_);
    Azimuthal_quantum_number.resize(atomic_number_);
    for( int Zstar=0; Zstar<(int)atomic_number_; Zstar++) {
        Potential               [Zstar] = IonizationTables::ionization_energy      (atomic_number_, Zstar) * eV_to_au;
        Azimuthal_quantum_number[Zstar] = IonizationTables::azimuthal_atomic_number(atomic_number_, Zstar);
    }

    for (unsigned int i=0; i<atomic_number_; i++) {
        DEBUG("ioniz: i " << i << " potential: " << Potential[i] << " Az.q.num: " << Azimuthal_quantum_number[i]);
    }
}


Ionization::~Ionization() {
}
