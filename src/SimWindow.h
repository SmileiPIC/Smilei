
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
    int cell_length_x_;

};

#endif /* SIMWINDOW_H */

