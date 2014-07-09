
#ifndef SIMWINDOW_H
#define SIMWINDOW_H

#include <vector>

class PicParams;
class Species;
class ElectroMagn;
class Interpolator;
class Projector;
class SmileiMPI;

class SimWindow {
    
 public:
    SimWindow(PicParams& params);
    ~SimWindow();
    void operate(std::vector<Species*> vecSpecies, ElectroMagn* EMfields, Interpolator* Interp, Projector* Proj, SmileiMPI* smpi);
    bool isMoving(int itime);
    

 private:
    int res_space_win_x_;
    double cell_length_x_;
    double x_moved;    //Total length the window has moved along x up to now.
    double vx_win;     //Velocity of the moving window along x expressed in c.
    double t_move;     //Time at which the window starts moving.

};

#endif /* SIMWINDOW_H */

