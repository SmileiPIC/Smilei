#ifndef PARTBOUNDCOND_H
#define PARTBOUNDCOND_H

#include "Params.h"
#include "Species.h"
#include "Particles.h"
#include "tabulatedFunctions.h"

class SmileiMPI;

//  --------------------------------------------------------------------------------------------------------------------
//! Class PartBoundCond
//  --------------------------------------------------------------------------------------------------------------------
class PartBoundCond {
public:
    //! partBoundCond creator, (default no MPI set)
    PartBoundCond( Params& params, SpeciesStructure& sparams, SmileiMPI* smpi );
    //! partBoundCond destructor
    ~PartBoundCond();

    //! West particles boundary conditions pointers (same prototypes for all conditions)
    //! @see BoundaryConditionType.h for functions that this pointers will target
    int (*bc_west)  ( Particles &particles, int ipart, int direction, double limit_pos, SpeciesStructure &params, double &nrj_iPart);
    //! East particles boundary conditions pointers
    int (*bc_east)  ( Particles &particles, int ipart, int direction, double limit_pos, SpeciesStructure &params, double &nrj_iPart);
    //! South particles boundary conditions pointers
    int (*bc_south) ( Particles &particles, int ipart, int direction, double limit_pos, SpeciesStructure &params, double &nrj_iPart);
    //! North particles boundary conditions pointers
    int (*bc_north) ( Particles &particles, int ipart, int direction, double limit_pos, SpeciesStructure &params, double &nrj_iPart);
    //! Bottom particles boundary conditions pointers
    int (*bc_bottom)( Particles &particles, int ipart, int direction, double limit_pos, SpeciesStructure &params, double &nrj_iPart);
    //! Up particles boundary conditions pointers
    int (*bc_up)    ( Particles &particles, int ipart, int direction, double limit_pos, SpeciesStructure &params, double &nrj_iPart);

    //! Method which applies particles boundary conditions.
    //! If the MPI process is not a border process, particles will be flagged as an exchange particle returning 0
    //! Conditions along X are applied first, then Y, then Z.
    //! The decision whether the particle is added or not on the Exchange Particle List is defined by the final
    //! value of keep_part. 
    //! Be careful, once an a BC along a given dimension set keep_part to 0, it will remain to 0. 
    inline int apply( Particles &particles, int ipart, SpeciesStructure &params, double &nrj_iPart ) {//, bool &contribute ) {
        
        int keep_part = 1;
        if ( particles.position(0, ipart) <  x_min ) {
            if (bc_west==NULL) keep_part = 0;
            else {
                keep_part = (*bc_west)( particles, ipart, 0, 2.*x_min, params,nrj_iPart );
            }
        }
        else if ( particles.position(0, ipart) >= x_max ) {
            if (bc_east==NULL) keep_part = 0;
            else {
                keep_part = (*bc_east)( particles, ipart, 0, 2.*x_max, params,nrj_iPart );
            }
        }
        if (nDim_particle >= 2) {
            
            if ( particles.position(1, ipart) <  y_min ) {
                if (bc_south==NULL) keep_part = 0;
                else {
                    keep_part *= (*bc_south)( particles, ipart, 1, 2.*y_min, params,nrj_iPart );
                }
            }
            else if ( particles.position(1, ipart) >= y_max ) {
                if (bc_north==NULL) keep_part = 0;
                else {
                    keep_part *= (*bc_north)( particles, ipart, 1, 2.*y_max, params,nrj_iPart );
                }
            }
            
            if (nDim_particle == 3) {
                
                if ( particles.position(2, ipart) <  z_min ) {
                    if (bc_bottom==NULL) keep_part = 0;
                    else {
                        keep_part *= (*bc_bottom)( particles, ipart, 2, 2.*z_min, params,nrj_iPart );
                    }
                }
                else if ( particles.position(2, ipart) >= z_max ) {
                    if (bc_up==NULL) keep_part = 0;
                    else {
                        keep_part *= (*bc_up)( particles, ipart, 2, 2.*z_max, params,nrj_iPart );
                    }
                }
            } // end if (nDim_particle == 3)
        } // end if (nDim_particle >= 2)


        return keep_part;
    };

    //! Move the condition window, not for simulation limits but for MPI exchange conditions
    void moveWindow_x(double shift, SmileiMPI* smpi );

    //! Set the condition window if restart
    inline void updateMvWinLimits( double x_moved ) {
        x_min += x_moved;
        x_max += x_moved;
    }

private:
    //! Min value of the x coordinate of particles on the current processor
    //! Real value, oversize is not considered (same for all)
    double x_min;
    //! Max value of the x coordinate of particles on the current processor
    double x_max;
    //! Min value of the y coordinate of particles on the current processor
    double y_min;
    //! Max value of the y coordinate of particles on the current processor
    double y_max;
    //! Min value of the z coordinate of particles on the current processor
    double z_min;
    //! Max value of the z coordinate of particles on the current processor
    double z_max;

    //! Space dimension of a particle
    int nDim_particle;
    
};

#endif

