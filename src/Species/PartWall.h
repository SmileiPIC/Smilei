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
    PartWall( Params& params, SmileiMPI* smpi );
    //! PartWall destructor
    ~PartWall(){};

    //! West particles boundary conditions pointers (same prototypes for all conditions)
    //! @see BoundaryConditionType.h for functions that this pointers will target
    int (*wall_x) ( Particles &particles, int ipart, int direction, double limit_pos, SpeciesStructure &params, double &nrj_iPart);
    //! South particles boundary conditions pointers
    int (*wall_y) ( Particles &particles, int ipart, int direction, double limit_pos, SpeciesStructure &params, double &nrj_iPart);
    //! North particles boundary conditions pointers
    int (*wall_z) ( Particles &particles, int ipart, int direction, double limit_pos, SpeciesStructure &params, double &nrj_iPart);
    

    //! Method which applies particles wall.
    inline int apply( Particles &particles, int ipart, SpeciesStructure &params, double &nrj_iPart ) {//, bool &contribute ) {
        
        int keep_part = 1;

        if (direction == 0) {
            if ((particles.position_old(0, ipart) < pos && particles.position(0, ipart) > pos) ||
                (particles.position_old(0, ipart) > pos && particles.position(0, ipart) < pos) ) {
                if (wall_x==NULL) {
                    keep_part = 0;
                } else {
                    keep_part = (*wall_x)( particles, ipart, 0, 2.*pos, params,nrj_iPart );
                }
            }
        } else if (direction==1) {
            if ((particles.position_old(1, ipart) < pos && particles.position(1, ipart) > pos) ||
                (particles.position_old(1, ipart) > pos && particles.position(1, ipart) < pos) ) {
                if (wall_y==NULL) {
                    keep_part = 0;
                } else {
                    keep_part = (*wall_y)( particles, ipart, 1, 2.*pos, params,nrj_iPart );
                }
            }
        } else if (direction==2) {
            if ((particles.position_old(2, ipart) < pos && particles.position(2, ipart) > pos) ||
                (particles.position_old(2, ipart) > pos && particles.position(2, ipart) < pos) ) {
                if (wall_z==NULL) {
                    keep_part = 0;
                } else {
                    keep_part = (*wall_z)( particles, ipart, 2, 2.*pos, params,nrj_iPart );
                }
            }
        }

        return keep_part;
    };

    //! true if present on this processor
    bool is_here;
    
private:
    //! direction of the wall
    unsigned int direction;
    
    //! position of a wall in the direction
    double pos;

    //! Space dimension of a particle
    int nDim_particle;
    
};

#endif

