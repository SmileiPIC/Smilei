#include <vector>

#include "Ionization.h"
#include "Species.h"

using namespace std;

Ionization::Ionization( Params &params, Species *species )
{

    reference_angular_frequency_SI = params.reference_angular_frequency_SI;
    
    dt                      = params.timestep;
    invdt                   = 1./dt;
    nDim_field              = params.nDim_field;
    nDim_particle           = params.nDim_particle;
    ionized_species_invmass = 1./species->mass_;
    
    // Normalization constant from Smilei normalization to/from atomic units
    eV_to_au   = 1.0 / 27.2116;
    au_to_mec2 = 27.2116/510.998e3;
    EC_to_au   = 3.314742578e-15 * reference_angular_frequency_SI; // hbar omega / (me c^2 alpha^3)
    au_to_w0   = 4.134137172e+16 / reference_angular_frequency_SI; // alpha^2 me c^2 / (hbar omega)
    
#ifdef _OMPTASKS
    new_electrons_per_bin = new Particles[species->Nbins];
    ion_charge_per_bin_.resize( species->Nbins );
#endif
}


Ionization::~Ionization()
{
}


void Ionization::joinNewElectrons( unsigned int Nbins )
{
    // if tasks on bins are used for ionization, join the lists of new electrons 
    // created in each bin, to have the list of new electrons for this species and patch
    
    size_t start = new_electrons.size();
    
    // Resize new_electrons
    size_t total_n_new = 0;
    for( size_t ibin = 0 ; ibin < Nbins ; ibin++ ) {
        total_n_new += new_electrons_per_bin[ibin].size();
    }
    new_electrons.createParticles( total_n_new );
    // Also resize ion_charge_ if necessary
    if( save_ion_charge_ ) {
        ion_charge_.resize( start + total_n_new );
    }
    
    // Move each new_electrons_per_bin into new_electrons
    for( size_t ibin = 0 ; ibin < Nbins ; ibin++ ) {
        size_t n_new = new_electrons_per_bin[ibin].size();
        for( size_t i=0; i<new_electrons.dimension(); i++ ) {
            copy( &new_electrons_per_bin[ibin].position( i, 0 ), &new_electrons_per_bin[ibin].position( i, n_new ), &new_electrons.position( i, start ) );
        }
        for( size_t i=0; i<new_electrons.dimension(); i++ ) {
            copy( &new_electrons_per_bin[ibin].momentum( i, 0 ), &new_electrons_per_bin[ibin].momentum( i, n_new ), &new_electrons.momentum( i, start ) );
        }
        copy( &new_electrons_per_bin[ibin].weight( 0 ), &new_electrons_per_bin[ibin].weight( n_new ), &new_electrons.weight( start ) );
        copy( &new_electrons_per_bin[ibin].charge( 0 ), &new_electrons_per_bin[ibin].charge( n_new ), &new_electrons.charge( start ) );
        new_electrons_per_bin[ibin].clear();
        if( save_ion_charge_ ) {
            copy( ion_charge_per_bin_[ibin].begin(), ion_charge_per_bin_[ibin].end(), ion_charge_.begin() + start );
            ion_charge_per_bin_[ibin].clear();
        }
        start += n_new;
    }
}
