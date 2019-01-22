#include "Ionization.h"
#include "Species.h"

Ionization::Ionization( Params &params, Species *species )
{

    reference_angular_frequency_SI = params.reference_angular_frequency_SI;
    
    dt                      = params.timestep;
    invdt                   = 1./dt;
    nDim_field              = params.nDim_field;
    nDim_particle           = params.nDim_particle;
    ionized_species_invmass = 1./species->mass;
    
    // Normalization constant from Smilei normalization to/from atomic units
    eV_to_au   = 1.0 / 27.2116;
    au_to_mec2 = 27.2116/510.998e3;
    EC_to_au   = 3.314742578e-15 * reference_angular_frequency_SI; // hbar omega / (me c^2 alpha^3)
    au_to_w0   = 4.134137172e+16 / reference_angular_frequency_SI; // alpha^2 me c^2 / (hbar omega)
    
    
}


Ionization::~Ionization()
{
}
