#ifndef PartWall_H
#define PartWall_H

#include "Params.h"
#include "Species.h"
#include "Particles.h"
#include "tabulatedFunctions.h"

class Patch;

//  --------------------------------------------------------------------------------------------------------------------
//! Class PartWall
//  --------------------------------------------------------------------------------------------------------------------
class PartWall {
public:
    //! PartWall constructor
    PartWall(double, unsigned short, std::string);
    //! PartWall destructor
    ~PartWall(){};
    
    //! Method that creates a vector of PartWall objects: one for each group in the input file.
    static std::vector<PartWall*> create(Params&, Patch*);
    //! Method that clones a vector of PartWall objects
    static std::vector<PartWall*> clone(std::vector<PartWall*>);
    
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


class PartWalls {
public:
    //! PartWalls constructor
    PartWalls(Params&, Patch*);
    //! PartWalls cloning constructor
    PartWalls(PartWalls*, Patch*);
    //! PartWalls destructor
    ~PartWalls();
    
    unsigned int size()  { return vecPartWall.size();  };
    void clear() { vecPartWall.clear(); };
    void resize(int i) { vecPartWall.resize(i); };
    void push_back(PartWall* pw) { vecPartWall.push_back(pw); };
    PartWall* operator[](int i) const { return vecPartWall[i]; }
    
private:
    std::vector<PartWall*> vecPartWall;
    
    // save all walls parameters even though the current patch doesn't contain them
    std::vector<short> direction;
    std::vector<double> position;
    std::vector<std::string> kind;
};

#endif

