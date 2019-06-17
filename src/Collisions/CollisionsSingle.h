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
    CollisionsSingle( Params &params,
                      unsigned int n_collisions,
                      std::vector<unsigned int> sg1,
                      std::vector<unsigned int> sg2,
                      double coulomb_log,
                      bool intra_collisions,
                      int debug_every,
                      int Z,
                      int ionization_electrons,
                      bool tracked,
                      int nDim,
                      std::string fname
                    ) : Collisions( params,
                                        n_collisions,
                                        sg1,
                                        sg2,
                                        coulomb_log,
                                        intra_collisions,
                                        debug_every,
                                        Z,
                                        ionization_electrons,
                                        tracked,
                                        nDim,
                                        fname
                                      ) {} ;
    //! Cloning Constructor
    CollisionsSingle( Collisions *coll, int ndim ) : Collisions( coll, ndim ) {};
    //! destructor
    ~CollisionsSingle() {};
    
    //! Method called in the main smilei loop to apply collisions at each timestep
    void collide( Params &, Patch *, int, std::vector<Diagnostic *> & ) override;
    
};


#endif
