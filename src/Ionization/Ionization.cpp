#include "Ionization.h"
#include "Species.h"

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
        new_electrons_per_bin = new Particles[species->particles->first_index.size()];
#endif
}


Ionization::~Ionization()
{
}


void Ionization::joinNewElectrons(unsigned int Nbins)
{

    // if tasks on bins are used for ionization, join the lists of new electrons 
    // created in each bin, to have the list of new electrons for this species and patch
    for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
        // number of particles to add from the bin
        unsigned int n_particles_to_create = new_electrons_per_bin[ibin].size();
        new_electrons.createParticles(n_particles_to_create);

        for (unsigned int ipart = 0; ipart < n_particles_to_create ; ipart++){

            int idNew = (new_electrons.size() - n_particles_to_create) + ipart;

            for( unsigned int i=0; i<new_electrons.dimension(); i++ ) {
                new_electrons.position( i, idNew ) = (new_electrons_per_bin[ibin]).position( i, ipart );
            }
            for( unsigned int i=0; i<3; i++ ) {
                new_electrons.momentum( i, idNew ) = (new_electrons_per_bin[ibin]).momentum( i, ipart );
            }
            new_electrons.weight( idNew ) = (new_electrons_per_bin[ibin]).weight( ipart );
            new_electrons.charge( idNew ) = (new_electrons_per_bin[ibin]).charge( ipart );
       
        } // end ipart
        new_electrons_per_bin[ibin].clear();
    } // end ibin

    

}
