
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


}

bool SimWindow::isMoving(int itime)
{
    return ( (res_space_win_x_)&&(itime>res_space_win_x_*0.75) );
}
