
#ifndef SIMWINDOW_H
#define SIMWINDOW_H

#include <vector>

class Params;
class Species;
class ElectroMagn;
class Interpolator;
class Projector;
class SmileiMPI;

//  --------------------------------------------------------------------------------------------------------------------
//! Class SimWindow
//  --------------------------------------------------------------------------------------------------------------------
class SimWindow {
    
 public:
    //! SimWindow creator
    SimWindow(Params& params);
    //! SimWindow destructor
    ~SimWindow();
    //! Move the simulation window (particles, fields, MPI environment & operator related to the grid)
    void operate(std::vector<Species*> vecSpecies, ElectroMagn* EMfields, Interpolator* Interp, Projector* Proj, SmileiMPI* smpi, Params& param);

    //! Returns a boolean : True if the window should be moved, False if it should not.
    //! Warning : Actually moving the window (function operate) changes the value of x_moved so the returned value of isMoving changes
    //! directly after moving the window.
    //! isMoving is called once in Smilei.cpp, in isProj and in solveMaxwell. Since this is BEFORE operate, it is correct. Take care not to
    //! call isMoving AFTER operate because the returned result might not be the expected one.
    bool isMoving(double time_dual);

    //! Return total length the window has moved
    double getXmoved() {return x_moved;}
    //! Return total number of cell of the window
    double getNspace_win_x() {return nspace_win_x_;}
    //! Set total length the window has moved (restart case)
    void   setXmoved(double new_val) {x_moved = new_val;}

    //! Set the simulation window (particles, fields, MPI environment & operator related to the grid) in restart case
    void setOperators(std::vector<Species*> vecSpecies, Interpolator* Interp, Projector* Proj, SmileiMPI* smpi);
    

 private:
    //! Number of points of the window
    int nspace_win_x_;
    //! Store locally params.cell_length[0], window slides only in x
    double cell_length_x_;
    //! Total length the window has moved along x up to now.
    double x_moved;
    //! Velocity of the moving window along x expressed in c.
    double vx_win_;
    //! Time at which the window starts moving.
    double t_move_win_;

};

#endif /* SIMWINDOW_H */

