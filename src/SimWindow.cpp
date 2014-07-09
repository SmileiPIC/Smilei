
#include "SimWindow.h"
#include "PicParams.h"
#include "Species.h"
#include "ElectroMagn.h"
#include "Interpolator.h"
#include "Projector.h"
#include "SmileiMPI.h"
using namespace std;

SimWindow::SimWindow(PicParams& params)
{
    res_space_win_x_ = params.res_space_win_x;
    cell_length_x_   = params.cell_length[0];
    x_moved = 0.;      //The window has not moved at t=0. Warning: not true anymore for restarts.
    vx_win_ = params.vx_win; 
    t_move_win_ = params.t_move_win;      
}

SimWindow::~SimWindow()
{
}

void SimWindow::operate(vector<Species*> vecSpecies, ElectroMagn* EMfields, Interpolator* Interp, Projector* Proj, SmileiMPI* smpi)
{

    unsigned int clrw = vecSpecies[0]->clrw; //clrw must be the same for all species

    smpi->getCellStartingGlobalIndex(0)+= clrw;
    smpi->getDomainLocalMin(0)+= cell_length_x_*clrw;
    smpi->getDomainLocalMax(0)+= cell_length_x_*clrw;

    
    for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) {
	vecSpecies[ispec]->movingWindow_x(clrw,smpi);
    }

    Interp->mv_win(clrw);
    Proj->mv_win(clrw);
    EMfields->movingWindow_x(clrw, smpi);
    x_moved += cell_length_x_*clrw;


}

bool SimWindow::isMoving(int itime)
{
//Warning: isMoving function returns a boolean. True if the window should be moved, False if it should not.
//Actually moving the window (function operate) changes the value of x_moved so the returned value of isMoving changes
//directly after moving the window.
//isMoving is called once in Electromagn. Since this is BEFORE operate, it is correct. Take care not to
//call isMoving AFTER operate because the returned result might not be the expected one.

    return ( (res_space_win_x_) && ((itime - t_move_win_)*vx_win_ > x_moved) );
}
