
#ifndef SIMWINDOW_H
#define SIMWINDOW_H

#include <vector>
#include "Patch.h"

class Params;
class Species;
class ElectroMagn;
class Interpolator;
class Projector;
class SmileiMPI;
class Region;

//  --------------------------------------------------------------------------------------------------------------------
//! Class SimWindow
//  --------------------------------------------------------------------------------------------------------------------
class SimWindow
{

public:
    //! SimWindow creator
    SimWindow( Params &params );
    //! SimWindow destructor
    ~SimWindow();
    //! Move the simulation window (particles, fields, MPI environment & operator related to the grid)

    void shift( VectorPatch &vecPatches, SmileiMPI *smpi, Params &param, unsigned int itime, double time_dual, Region& region );
    
    void operate( Region& region,  VectorPatch& vecPatches, SmileiMPI* smpi, Params& param, double time_dual );
    void operate( Region& region,  VectorPatch& vecPatches, SmileiMPI* smpi, Params& param, double time_dual, unsigned int nmodes );

    //! Tells whether there is a moving window or not
    inline bool isActive()
    {
        return active;
    }
    
    //! Returns a boolean : True if the window should be moved, False if it should not.
    //! Warning : Actually moving the window (function shift) changes the value of x_moved so the returned value of isMoving changes
    //! directly after moving the window.
    //! isMoving is called once in Smilei.cpp, in isProj and in solveMaxwell. Since this is BEFORE shift, it is correct. Take care not to
    //! call isMoving AFTER shift because the returned result might not be the expected one.
    bool isMoving( double time_dual );
    
    //! Return total length the window has moved
    double getXmoved()
    {
        return x_moved;
    }
    //! Return the number of cells the window has moved
    unsigned int getNmoved()
    {
        return n_moved;
    }
    //! Set total length the window has moved (restart case)
    void   setXmoved( double new_val )
    {
        x_moved = new_val;
    }
    //! Set total number of cells the window has moved (restart case)
    void   setNmoved( int new_val )
    {
        n_moved = new_val;
    }
    //! Return number of shifts operations
    unsigned int getNumberOfAdditionalShifts()
    {
        return number_of_additional_shifts;
    }
    //! Return additional shifts iteration
    unsigned int getAdditionalShiftsIteration()
    {
        return additional_shifts_iteration;
    }

    //! Return additional shifts time
    double getAdditionalShiftsTime()
    {
        return additional_shifts_time;
    }
    
    
private:
    //! Tells whether there is a moving window or not
    bool active;
    
    //! Store locally params.cell_length[0], window slides only in x
    double cell_length_x_;
    //! Store locally params.n_space[0], window slides only in x
    double n_space_x_;
    //! Total length the window has moved along x up to now.
    double x_moved;
    //! Total number of cell the window has moved along x up to now.
    unsigned int n_moved;
    //! Velocity of the moving window along x expressed in c.
    double velocity_x;
    //! Time at which the window starts moving.
    double time_start;
    //! Keep track of old patches assignement
    std::vector<Patch *> vecPatches_old;
    //! Keep track of patches to create
    std::vector< std::vector<unsigned int>> patch_to_be_created;
    //! Keep track of patches that receive particles
    std::vector< std::vector<bool>> patch_particle_created;
    //! Max number of threads
    int max_threads;
    //! Time of additional moving window shifts
    double additional_shifts_time;
    //! Iteration of additional moving window shifts
    unsigned int additional_shifts_iteration;
    //! Number of additional moving window shifts
    unsigned int number_of_additional_shifts;
    
    
};

#endif /* SIMWINDOW_H */
