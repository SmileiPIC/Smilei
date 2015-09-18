#ifndef PartWall_H
#define PartWall_H

#include "Params.h"
#include "Species.h"
#include "Particles.h"
#include "tabulatedFunctions.h"

class SmileiMPI;

//  --------------------------------------------------------------------------------------------------------------------
//! Class PartWall
//  --------------------------------------------------------------------------------------------------------------------
class PartWall {
public:
    //! PartWall creator
    PartWall(double, unsigned short);
    //! PartWall destructor
    ~PartWall(){};
    
    //! Method that creates a vector of PartWall objects: one for each group in the input file.
    static std::vector<PartWall*> create(Params&, SmileiMPI*);
    
    //! Wall boundary condition pointer (same prototypes for all conditions)
    //! @see BoundaryConditionType.h for functions that this pointer will target
    int (*wall) ( Particles &particles, int ipart, int direction, double limit_pos, Species *species, double &nrj_iPart);
    
    //! Method which applies particles wall
    int apply (Particles &particles, int ipart, Species *species, double &nrj_iPart);
    
private:
    //! position of a wall in its direction
    double position;
    unsigned short direction;

};

#endif

