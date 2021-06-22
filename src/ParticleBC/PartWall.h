#ifndef PartWall_H
#define PartWall_H

#include "Params.h"
#include "tabulatedFunctions.h"

class Patch;
class Species;
class Particles;

//  --------------------------------------------------------------------------------------------------------------------
//! Class PartWall
//  --------------------------------------------------------------------------------------------------------------------
class PartWall
{
public:
    //! PartWall constructor
    PartWall( double, unsigned short, std::string, double dt );
    //! PartWall destructor
    ~PartWall() {};
    
    //! Wall boundary condition pointer (same prototypes for all conditions)
    //! @see BoundaryConditionType.h for functions that this pointer will target
    void ( *wall )( Particles &particles, SmileiMPI *smpi, int imin, int imax, int direction, double limit_pos, double dt, Species *species, int ithread, double &nrj_iPart );
    
    //! Method which applies particles wall
    void apply( Particles &particles, SmileiMPI *smpi, int imin, int imax, Species *species, int ithread, double &nrj_iPart );
    
    double dt_;

private:
    //! position of a wall in its direction
    double position;
    
    //! direction of the partWall (x=0, y=1, z=2)
    unsigned short direction;
    
};


class PartWalls
{
public:
    //! PartWalls constructor
    PartWalls( Params &, Patch * );
    //! PartWalls cloning constructor
    PartWalls( PartWalls *, Patch * );
    //! PartWalls destructor
    ~PartWalls();
    
    //! Returns the number of walls
    unsigned int size()
    {
        return vecPartWall.size();
    };
    
    //! Clears the vector of walls
    void clear()
    {
        vecPartWall.clear();
    };
    
    //! Resizes the vector of walls
    void resize( int i )
    {
        vecPartWall.resize( i );
    };
    
    //! Pushes a pointer to a wall at the end of the vector
    void push_back( PartWall *pw )
    {
        vecPartWall.push_back( pw );
    };
    
    //! Accesses the wall of index i
    PartWall *operator[]( int i ) const
    {
        return vecPartWall[i];
    }
    
private:
    //! The vector of partWall objects
    std::vector<PartWall *> vecPartWall;
    
    // save all walls parameters even though the current patch doesn't contain them
    
    //! direction of all walls
    std::vector<short> direction;
    //! position of all walls
    std::vector<double> position;
    //! kind of all walls (type of boundary condition applied)
    std::vector<std::string> kind;

};

#endif

