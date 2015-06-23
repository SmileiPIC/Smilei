
#include "SimWindow.h"
#include "PicParams.h"
#include "Species.h"
#include "ElectroMagn.h"
#include "Interpolator.h"
#include "Projector.h"
#include "SmileiMPI.h"
#include "Patch.h"
#include "Hilbert_functions.h"

using namespace std;

SimWindow::SimWindow(PicParams& params)
{
    nspace_win_x_ = params.nspace_win_x;
    cell_length_x_   = params.cell_length[0];
    x_moved = 0.;      //The window has not moved at t=0. Warning: not true anymore for restarts.
    vx_win_ = params.vx_win; 
    t_move_win_ = params.t_move_win;      
}

SimWindow::~SimWindow()
{
}

void SimWindow::operate(vector<Species*> vecSpecies, ElectroMagn* EMfields, Interpolator* Interp, Projector* Proj, SmileiMPI* smpi, PicParams& params)
{

    unsigned int clrw;

    if (vecSpecies.size() > 0) {
        clrw = vecSpecies[0]->clrw; //clrw must be the same for all species
    } else {
        clrw = params.clrw; //This is not correct, just an ugly patch. clrw should take the correct global value even when there are no particles.
    }

    smpi->getCellStartingGlobalIndex(0)+= clrw;
    smpi->getDomainLocalMin(0)+= cell_length_x_*clrw;
    smpi->getDomainLocalMax(0)+= cell_length_x_*clrw;

    if (vecSpecies.size() > 0) {
        for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) {
            vecSpecies[ispec]->movingWindow_x(clrw,smpi,params);
        }
        Interp->mv_win(clrw);
        Proj->mv_win(clrw);
    }
    EMfields->movingWindow_x(clrw, smpi);
    x_moved += cell_length_x_*clrw;


}

void SimWindow::operate(std::vector<Patch*> vecPatches, SmileiMPI* smpi, PicParams& params)
{
    int xcall, ycall;
    #pragma omp for
    for (unsigned int ipatch = 0 ; ipatch < vecPatches.size() ; ipatch++) {
        //Si je ne possede pas mon voisin de gauche...
        if (vecPatches[ipatch]->MPI_neighborhood_[3] != vecPatches[ipatch]->MPI_neighborhood_[4]) {
            //...je l'envois si necessaire ... 
            if (vecPatches[ipatch]->MPI_neighborhood_[3] != MPI_PROC_NULL)
                cout << "Sending Patch" << endl;
                //Mpi_Send_Patch(MpiLNeighbour);
            //... et je detruit mes donnees.
            for (unsigned int ispec=0 ; ispec<vecPatches[ipatch]->vecSpecies.size(); ispec++) delete vecPatches[ipatch]->vecSpecies[ispec];
	    vecPatches[ipatch]->vecSpecies.clear();
            delete (vecPatches[ipatch]->EMfields);
            delete (vecPatches[ipatch]->Interp);
            delete (vecPatches[ipatch]->Proj);
        //Sinon, je deviens mon voisin de gauche.
        } else {
            vecPatches[ipatch]->Pcoordinates[0] -= 1;
            vecPatches[ipatch]->min_local[0] -= params.n_space[0]*cell_length_x_;
            vecPatches[ipatch]->max_local[0] -= params.n_space[0]*cell_length_x_;
            vecPatches[ipatch]->cell_starting_global_index[0] -= params.n_space[0];

            //Shift neighborhood tables.
	    for ( int z = 0 ; z < 1+2*(params.nDim_field == 3) ; z++ ) {
	        for ( int y = 0 ; y < 1+2*(params.nDim_field >= 2) ; y++ ) {
	            for ( int x = 2 ; x > 0 ; x-- ) {
	            vecPatches[ipatch]->patch_neighborhood_[z*9+y*3+x] = vecPatches[ipatch]->patch_neighborhood_[z*9+y*3+x-1];
	            vecPatches[ipatch]->MPI_neighborhood_[z*9+y*3+x] = vecPatches[ipatch]->MPI_neighborhood_[z*9+y*3+x-1];
                    }
                }
            }
            //Compute missing part of the new neighborhood tables.
            xcall = vecPatches[ipatch]->Pcoordinates[0]-1;
            ycall = vecPatches[ipatch]->Pcoordinates[1]-1;
            if (params.bc_em_type_long=="periodic") xcall = xcall%((1<<vecPatches[ipatch]->mi[0]));
            if (params.bc_em_type_trans=="periodic") ycall = ycall%((1<<vecPatches[ipatch]->mi[1]));
	    vecPatches[ipatch]->patch_neighborhood_[0] = generalhilbertindex(vecPatches[ipatch]->mi[0] , vecPatches[ipatch]->mi[1], xcall, ycall);
	    vecPatches[ipatch]->patch_neighborhood_[1] = generalhilbertindex(vecPatches[ipatch]->mi[0] , vecPatches[ipatch]->mi[1], xcall, vecPatches[ipatch]->Pcoordinates[1]);
            ycall = vecPatches[ipatch]->Pcoordinates[1]+1;
            if (params.bc_em_type_trans=="periodic") ycall = ycall%((1<<vecPatches[ipatch]->mi[1]));
	    vecPatches[ipatch]->patch_neighborhood_[2] = generalhilbertindex(vecPatches[ipatch]->mi[0] , vecPatches[ipatch]->mi[1], xcall, ycall);
	    for ( int y = 0 ; y < 1+2*(params.nDim_field >= 2) ; y++ ) {
	        for ( int z = 0 ; z < 1+2*(params.nDim_field == 3) ; z++ ) {
                    vecPatches[ipatch]->MPI_neighborhood_[y*3+z] = smpi->hrank(vecPatches[ipatch]->patch_neighborhood_[y*3+z]);
                }
            }
            
             
        }
        //Si je ne possede pas mon voisin de droite...
        if (vecPatches[ipatch]->MPI_neighborhood_[5] != vecPatches[ipatch]->MPI_neighborhood_[4]) {
            //...je reçois ou je cré.
            if (vecPatches[ipatch]->MPI_neighborhood_[3] != MPI_PROC_NULL){
                //Mpi_Receive_Patch(MpiRNeighbour);
                cout << "Receiving Patch" << endl;
            } else {
                //Create_Patch();
                cout << "Patch creation" << endl;
            }
        }
        //Les fonctions SendPatch ou ReceivePatch doivent transformer le Patch comme ci dessus.
        //Ce serait plus simple de mettre hindex, neighbor_ et corner_neighbor_ dans un seul et meme tableau.

    }//End loop on Patches
}



bool SimWindow::isMoving(double time_dual)
{
    return ( (nspace_win_x_) && ((time_dual - t_move_win_)*vx_win_ > x_moved) );
}

void SimWindow::setOperators(vector<Species*> vecSpecies, Interpolator* Interp, Projector* Proj, SmileiMPI* smpi)
{
    smpi->updateMvWinLimits( x_moved, round(x_moved/cell_length_x_) );

    for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) {
	vecSpecies[ispec]->updateMvWinLimits(x_moved);
    }

    Interp->setMvWinLimits( smpi->getCellStartingGlobalIndex(0) );
    Proj->setMvWinLimits  ( smpi->getCellStartingGlobalIndex(0) );

}
