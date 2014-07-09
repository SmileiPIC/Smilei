
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
    vx_win  = 0.99987; //Should be read from input file.
    t_move = 0.0;      //Should be read from input file.
}

SimWindow::~SimWindow()
{
}

void SimWindow::operate(vector<Species*> vecSpecies, ElectroMagn* EMfields, Interpolator* Interp, Projector* Proj, SmileiMPI* smpi)
{

    smpi->getCellStartingGlobalIndex(0)+= 1;
    smpi->getDomainLocalMin(0)+= cell_length_x_;
    smpi->getDomainLocalMax(0)+= cell_length_x_;

    
    for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) {
	vecSpecies[ispec]->movingWindow_x(1, smpi);
    }

    Interp->mv_win(1);
    Proj->mv_win(1);
    EMfields->movingWindow_x(1, smpi);
    x_moved += cell_length_x_;


}

bool SimWindow::isMoving(int itime)
{
//Warning: isMoving function returns a boolean. True if the window should be moved, False if it should not.
//Actually moving the window (function operate) changes the value of x_moved so the returned value of isMoving changes
//directly after moving the window.
//isMoving is called once in Electromagn. Since this is BEFORE operate, it is correct. Take care not to
//call isMoving AFTER operate because the returned result might not be the expected one.

    return ( (res_space_win_x_) && ((itime - t_move)*vx_win > x_moved) );
}
