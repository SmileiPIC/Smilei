/*! @file Pusher.h

  @brief Pusher.h  generic class for the particle pusher

  @date 2013-02-15
*/

#ifndef PartWall_H
#define PartWall_H

#include "Params.h"
#include "Particles.h"
#include "tabulatedFunctions.h"

class SmileiMPI;


//  --------------------------------------------------------------------------------------------------------------------
//! Class PartWall
//  --------------------------------------------------------------------------------------------------------------------
class PartWall {
public:
    //! PartWall creator, (default no MPI set)
    PartWall( short, std::string, double );
    //! PartWall destructor
    ~PartWall(){};
    
    //! Method that creates a vector of PartWall objects: one for each group in the input file.
    static std::vector<PartWall*> create(Params&, SmileiMPI*);
    
    
    //! Wall boundary condition pointer (same prototypes for all conditions)
    //! @see BoundaryConditionType.h for functions that this pointer will target
    int (*wall) ( Particles &particles, int ipart, int direction, double limit_pos, SpeciesStructure &params, double &nrj_iPart);
    
    //! Method which applies particles wall.
    int (PartWall::*apply) ( Particles &particles, int ipart, SpeciesStructure &params, double &nrj_iPart);
    // Preliminary functions
    inline int apply_x( Particles &particles, int ipart, SpeciesStructure &params, double &nrj_iPart) {
        int keep_part = 1;
        if ((particles.position_old(0, ipart) < pos && particles.position(0, ipart) > pos) ||
            (particles.position_old(0, ipart) > pos && particles.position(0, ipart) < pos) ) {
            keep_part = (*wall)( particles, ipart, 0, 2.*pos, params, nrj_iPart );
        }
        return keep_part;
    };
    inline int apply_y( Particles &particles, int ipart, SpeciesStructure &params, double &nrj_iPart) {
        int keep_part = 1;
        if ((particles.position_old(1, ipart) < pos && particles.position(1, ipart) > pos) ||
            (particles.position_old(1, ipart) > pos && particles.position(1, ipart) < pos) ) {
            keep_part = (*wall)( particles, ipart, 1, 2.*pos, params, nrj_iPart );
        }
        return keep_part;
    };
    inline int apply_z( Particles &particles, int ipart, SpeciesStructure &params, double &nrj_iPart) {
        int keep_part = 1;
        if ((particles.position_old(2, ipart) < pos && particles.position(2, ipart) > pos) ||
            (particles.position_old(2, ipart) > pos && particles.position(2, ipart) < pos) ) {
            keep_part = (*wall)( particles, ipart, 2, 2.*pos, params, nrj_iPart );
        }
        return keep_part;
    };


private:
    //! position of a wall in the direction
    double pos;
    
};

#endif

