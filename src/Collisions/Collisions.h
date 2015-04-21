/*

Collisions class - Frederic Perez - 03/2015

This is based on the work described here
http://dx.doi.org/10.1063/1.4742167

Binary collisions, between macro-particles, are treated according to a scheme
first described by Nanbu (http://dx.doi.org/10.1103/PhysRevE.55.4642).

To include collisions in the simulations, add a block in the input file, 
similar to the following:

# species1    : "type" or "name" of the first  species that collide
#               (can be a list of species)
# species2    : "type" or "name" of the second species that collide
#               (can be a list of species) (can be the same as the first species)
# coulomb_log : value of the Coulomb logarithm. If negative or zero, then automatically computed.
collisions
	species1 = ion1
	species2 = electron1
	coulomb_log = 2.0
end

Several collision types can be defined. For each type, add a group "collisions".

*/

#ifndef COLLISIONS_H
#define COLLISIONS_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "Species.h"

class Collisions
{

public:
    //! Constructor for Collisions between two species
    Collisions(PicParams&,unsigned int,std::vector<unsigned int>,std::vector<unsigned int>,double,bool);
    
    //! Identification number of the Collisions object
    int n_collisions;
    
    //! Group of the species numbers that are associated for Collisions.
    std::vector<unsigned int> species_group1, species_group2;
    
    //! Coulomb logarithm (zero or negative means automatic)
    double coulomb_log;
    
    //! True if collisions inside a group of species, False if collisions between different groups of species
    double intra_collisions;
    
    //! Method to calculate the Debye length in each cluster
    static void calculate_debye_length(PicParams&,std::vector<Species*>&);
    
    //! is true if any of the collisions objects need automatically-computed coulomb log
    static bool debye_length_required;
    
    //! Method called in the main smilei loop to apply collisions at each timestep
    void collide(PicParams&,std::vector<Species*>&);
    
private:
    
    //! Contains the debye length in each cluster, computed each timestep
    static std::vector<double> debye_length_squared; 
    
    static double cos_chi(double);


};


#endif
