#ifndef COLLISIONS_H
#define COLLISIONS_H

#include <vector>

#include "Tools.h"
#include "H5.h"
#include "CollisionalIonization.h"

class Patch;
class Params;
class Species;
class VectorPatch;

class Collisions
{

public:
    //! Constructor for Collisions between two species
    Collisions( Patch* patch, unsigned int n_collisions, std::vector<unsigned int>,
        std::vector<unsigned int>, double coulomb_log, bool intra_collisions,
        int debug_every, int Z, bool ionizing, bool tracked_electrons, int nDim,
        double,std::string);
    //! Cloning Constructor
    Collisions(Collisions*, int);
    //! destructor
    ~Collisions();
    
    //! Method that creates a vector of Collisions objects: one for each group in the input file.
    static std::vector<Collisions*> create(Params&, Patch*, std::vector<Species*>&);
    //! Method that clones a vector of Collisions objects
    static std::vector<Collisions*> clone(std::vector<Collisions*>, Params&);
    
    //! Method to calculate the Debye length in each cluster
    static void calculate_debye_length(Params&, Patch*);
    
    //! is true if any of the collisions objects need automatically-computed coulomb log
    static bool debye_length_required;
    
    //! Method called in the main smilei loop to apply collisions at each timestep
    void collide(Params&, Patch* ,int, std::vector<Diagnostic*>&);
    
    //! Outputs the debug info if requested
    static void debug(Params& params, int itime, unsigned int icoll, VectorPatch& vecPatches);
    
    // Technique given by Nanbu in http://dx.doi.org/10.1103/PhysRevE.55.4642
    //   to pick randomly the deflection angle cosine, in the center-of-mass frame.
    // It involves the "s" parameter (~ collision frequency * deflection expectation)
    //   and a random number "U".
    // Technique slightly modified in http://dx.doi.org/10.1063/1.4742167
    inline double cos_chi(double s)
    {
        double A, invA;
        //!\todo make a faster rand by preallocating ??
        double U = Rand::uniform();
        
        if( s < 0.1 ) {
            if ( U<0.0001 ) U=0.0001; // ensures cos_chi > 0
            return 1. + s*log(U);
        }
        if( s < 3.  ) {
            // the polynomial has been modified from the article in order to have a better form
            invA = 0.00569578 +(0.95602 + (-0.508139 + (0.479139 + ( -0.12789 + 0.0238957*s )*s )*s )*s )*s;
            A = 1./invA;
            return  invA  * log( exp(-A) + 2.*U*sinh(A) );
        }
        if( s < 6.  ) {
            A = 3.*exp(-s);
            return (1./A) * log( exp(-A) + 2.*U*sinh(A) );
        }
        return 2.*U - 1.;
    }
    //! CollisionalIonization object, created if ionization required
    CollisionalIonization * Ionization;
    
private:
    
    //! Identification number of the Collisions object
    int n_collisions;
    
    //! Group of the species numbers that are associated for Collisions.
    std::vector<unsigned int> species_group1, species_group2;
    
    //! Coulomb logarithm (zero or negative means automatic)
    double coulomb_log;
    
    //! True if collisions inside a group of species, False if collisions between different groups of species
    bool intra_collisions;
    
    //! Number of timesteps between each dump of collisions debugging
    int debug_every;
    
    //! Species atomic number, in case of ionization
    int atomic_number;
    
    //! Hdf5 file name
    std::string filename;
    
    //! Temporary variables for the debugging file
    double smean, logLmean, ncol;//, temperature
};


#endif
