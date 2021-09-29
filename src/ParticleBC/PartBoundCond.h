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
    void ( *bc_xmin )( Species *species, int imin, int imax, int direction, double limit_inf, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );
    //! Xmax particles boundary conditions pointers
    void ( *bc_xmax )( Species *species, int imin, int imax, int direction, double limit_inf, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );
    //! Ymin particles boundary conditions pointers
    void ( *bc_ymin )( Species *species, int imin, int imax, int direction, double limit_inf, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );
    //! Ymax particles boundary conditions pointers
    void ( *bc_ymax )( Species *species, int imin, int imax, int direction, double limit_inf, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );
    //! Zmin particles boundary conditions pointers
    void ( *bc_zmin )( Species *species, int imin, int imax, int direction, double limit_inf, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );
    //! Zmax particles boundary conditions pointers
    void ( *bc_zmax )( Species *species, int imin, int imax, int direction, double limit_inf, double dt, std::vector<double> &invgf, Random * rand, double &energy_change );
    
    //! Method which applies particles boundary conditions.
    //! If the MPI process is not a border process, particles will be flagged as an exchange particle returning 0
    //! Conditions along X are applied first, then Y, then Z.
    inline void apply( Species *species, int imin, int imax, std::vector<double> &invgf, Random * rand, double &energy_tot )
    {
        for (int ipart=imin ; ipart<imax ; ipart++ ) {
            species->particles->cell_keys[ipart] = 0;
        }

        double energy_change = 0.;
        ( *bc_xmin )( species, imin, imax, 0, x_min, dt_, invgf, rand, energy_change );
        energy_tot += energy_change;
        ( *bc_xmax )( species, imin, imax, 0, x_max, dt_, invgf, rand, energy_change );
        energy_tot += energy_change;
        if( nDim_particle >= 2 ) {
            ( *bc_ymin )( species, imin, imax, 1, y_min, dt_, invgf, rand, energy_change );
            energy_tot += energy_change;
            ( *bc_ymax )( species, imin, imax, 1, y_max, dt_, invgf, rand, energy_change );
            energy_tot += energy_change;
            if( ( nDim_particle == 3 ) && (!isAM) ) {
                ( *bc_zmin )( species, imin, imax, 2, z_min, dt_, invgf, rand, energy_change );
                energy_tot += energy_change;
                ( *bc_zmax )( species, imin, imax, 2, z_max, dt_, invgf, rand, energy_change );
                energy_tot += energy_change;
            }
        }

    };
    
    ////! Set the condition window if restart (patch position not read)
    //inline void updateMvWinLimits( double x_moved ) {
    //}
    
private:
    //! Min x coordinate of particles on the current processor (oversize not considered)
    double x_min;
    //! Max x coordinate of particles on the current processor (oversize not considered)
    double x_max;
    //! Min y coordinate of particles on the current processor (oversize not considered)
    double y_min;
    //! Max y coordinate of particles on the current processor (oversize not considered)
    double y_max;
    //! Min z coordinate of particles on the current processor (oversize not considered)
    double z_min;
    //! Max z coordinate of particles on the current processor (oversize not considered)
    double z_max;
    
    //! Space dimension of a particle
    int nDim_particle;
    //! Space dimension of field
    int nDim_field;
    bool isAM;

    double dt_;

};

#endif
