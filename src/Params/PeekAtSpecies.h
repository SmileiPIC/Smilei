#ifndef PEEKATSPECIES_H
#define PEEKATSPECIES_H

#include "PyTools.h"
#include "Params.h"
#include "Profile.h"

//! Class to get a quick view of a species profile, even when the species particles have not been created
class PeekAtSpecies {
public:
    PeekAtSpecies(Params& params, unsigned int species_id);
    ~PeekAtSpecies();
    
    double numberOfParticlesInPatch( unsigned int hindex );
    
    double numberOfParticlesInPatch( std::vector<double> x_cell );
    
    double totalNumberofParticles( );
    
private:
    Profile *densityProfile;
    Profile *ppcProfile;
    Params * params;
};

#endif