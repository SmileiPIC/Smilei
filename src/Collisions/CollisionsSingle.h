#ifndef COLLISIONSSINGLE_H
#define COLLISIONSSINGLE_H

#include <vector>

#include "Tools.h"
#include "H5.h"
#include "Collisions.h"
#include "CollisionalIonization.h"

class Patch;
class Params;
class Species;
class VectorPatch;

class CollisionsSingle : public Collisions
{

public:
    //! Constructor for Collisions between two species
    CollisionsSingle(
        Params &params,
        unsigned int n_collisions,
        std::vector<unsigned int> sg1,
        std::vector<unsigned int> sg2,
        double coulomb_log,
        double coulomb_log_factor,
        bool intra_collisions,
        int every,
        int debug_every,
        CollisionalIonization *ionization,
        CollisionalNuclearReaction *nuclear_reaction,
        std::string fname
    ) : Collisions(
        params,
        n_collisions,
        sg1,
        sg2,
        coulomb_log,
        coulomb_log_factor,
        intra_collisions,
        every,
        debug_every,
        ionization,
        nuclear_reaction,
        fname
    ) {} ;
    
    //! Cloning Constructor
    CollisionsSingle( Collisions *coll ) : Collisions( coll ) {};
    //! destructor
    ~CollisionsSingle() {};
    
    //! Method called in the main smilei loop to apply collisions at each timestep
    void collide( Params &, Patch *, int, std::vector<Diagnostic *> & ) override;
    
};


#endif
