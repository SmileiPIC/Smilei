#ifndef PARTBOUNDCOND_H
#define PARTBOUNDCOND_H

#include "Params.h"
#include "Species.h"
#include "Particles.h"
#include "tabulatedFunctions.h"

class Patch;

//  --------------------------------------------------------------------------------------------------------------------
//! Class PartBoundCond
//  --------------------------------------------------------------------------------------------------------------------
class PartBoundCond
{
public:
    //! partBoundCond creator, (default no MPI set)
    PartBoundCond( Params &params, Species *species, Patch *patch );
    //! partBoundCond destructor
    ~PartBoundCond();
    
    //! Xmin particles boundary conditions pointers (same prototypes for all conditions)
    //! @see BoundaryConditionType.h for functions that this pointers will target
    void ( *bc_xmin )( Particles &particles, int ifirst, int ilast, int direction, double limit_pos, Species *species, double &nrj_iPart );
    //! Xmax particles boundary conditions pointers
    void ( *bc_xmax )( Particles &particles, int ifirst, int ilast, int direction, double limit_pos, Species *species, double &nrj_iPart );
    //! Ymin particles boundary conditions pointers
    void ( *bc_ymin )( Particles &particles, int ifirst, int ilast, int direction, double limit_pos, Species *species, double &nrj_iPart );
    //! Ymax particles boundary conditions pointers
    void ( *bc_ymax )( Particles &particles, int ifirst, int ilast, int direction, double limit_pos, Species *species, double &nrj_iPart );
    //! Zmin particles boundary conditions pointers
    void ( *bc_zmin )( Particles &particles, int ifirst, int ilast, int direction, double limit_pos, Species *species, double &nrj_iPart );
    //! Zmax particles boundary conditions pointers
    void ( *bc_zmax )( Particles &particles, int ifirst, int ilast, int direction, double limit_pos, Species *species, double &nrj_iPart );
    
    //! Method which applies particles boundary conditions.
    //! If the MPI process is not a border process, particles will be flagged as an exchange particle returning 0
    //! Conditions along X are applied first, then Y, then Z.
    //! The decision whether the particle is added or not on the Exchange Particle List is defined by the final
    //! value of keep_part.
    //! Be careful, once an a BC along a given dimension set keep_part to 0, it will remain to 0.
    inline void apply( Particles &particles, int imin, int imax, Species *species, double &nrj_tot )
    {
        for (int ipart=imin ; ipart<imax ; ipart++ ) {
            particles.cell_keys[ipart] = 0;
        }

        double nrj_iPart = 0.;
        ( *bc_xmin )( particles, imin, imax, 0, x_min, species, nrj_iPart );
        nrj_tot += nrj_iPart;
        ( *bc_xmax )( particles, imin, imax, 0, x_max, species, nrj_iPart );
        nrj_tot += nrj_iPart;
        if( nDim_particle >= 2 ) {
            ( *bc_ymin )( particles, imin, imax, 1, y_min, species, nrj_iPart );
            nrj_tot += nrj_iPart;
            ( *bc_ymax )( particles, imin, imax, 1, y_max, species, nrj_iPart );
            nrj_tot += nrj_iPart;
            if( nDim_particle == 3 ) {
                ( *bc_zmin )( particles, imin, imax, 2, z_min, species, nrj_iPart );
                nrj_tot += nrj_iPart;
                ( *bc_zmax )( particles, imin, imax, 2, z_max, species, nrj_iPart );
                nrj_tot += nrj_iPart;
            }
        }

    };
    
    ////! Set the condition window if restart (patch position not read)
    //inline void updateMvWinLimits( double x_moved ) {
    //}
    
private:
    //! Min value of the x coordinate of particles on the current processor
    //! Real value, oversize is not considered (same for all)
    double x_min;
    //! Max value of the x coordinate of particles on the current processor
    double x_max;
    //! Min value of the y coordinate of particles on the current processor
    double y_min;
    double y_min2;
    //! Max value of the y coordinate of particles on the current processor
    double y_max;
    double y_max2;
    //! Min value of the z coordinate of particles on the current processor
    double z_min;
    //! Max value of the z coordinate of particles on the current processor
    double z_max;
    
    //! Space dimension of a particle
    int nDim_particle;
    //! Space dimension of field
    int nDim_field;
//<<<<<<< HEAD
    bool isAM;
    
};

#endif
