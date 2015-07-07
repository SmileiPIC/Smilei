
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
    n_moved = 0;      //The window has not moved at t=0. Warning: not true anymore for restarts.
    vx_win_ = params.vx_win; 
    t_move_win_ = params.t_move_win;      
}

SimWindow::~SimWindow()
{
}

//void SimWindow::operate(vector<Species*> vecSpecies, ElectroMagn* EMfields, Interpolator* Interp, Projector* Proj, SmileiMPI* smpi, PicParams& params)
//{
//
//    unsigned int clrw;
//
//    if (vecSpecies.size() > 0) {
//        clrw = vecSpecies[0]->clrw; //clrw must be the same for all species
//    } else {
//        clrw = params.clrw; //This is not correct, just an ugly patch. clrw should take the correct global value even when there are no particles.
//    }
//
//    smpi->getCellStartingGlobalIndex(0)+= clrw;
//    smpi->getDomainLocalMin(0)+= cell_length_x_*clrw;
//    smpi->getDomainLocalMax(0)+= cell_length_x_*clrw;
//
//    if (vecSpecies.size() > 0) {
//        for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) {
//            vecSpecies[ispec]->movingWindow_x(clrw,smpi,params);
//        }
//        Interp->mv_win(clrw);
//        Proj->mv_win(clrw);
//    }
//    EMfields->movingWindow_x(clrw, smpi);
//    x_moved += cell_length_x_*clrw;
//
//
//}

void SimWindow::operate(VectorPatch& vecPatches, SmileiMPI* smpi, PicParams& params, DiagParams &diag_params, LaserParams& laser_params)
{
    int xcall, ycall, h0, do_right, receive_flag, old_index;
    Patch* mypatch;
    //std::vector<Patch*> vecPatches_old;

    h0 = vecPatches(0)->hindex;
    #pragma omp single
    {
    vecPatches_old.resize(vecPatches.size());
    }

    cout << " old_size = " << vecPatches_old.size() << endl;

    #pragma omp for schedule(static)
    for (unsigned int ipatch = 0 ; ipatch < vecPatches.size() ; ipatch++) {
        vecPatches_old[ipatch] = vecPatches(ipatch);
    } //Barrier at the end of this omp for is important to prevent an update of x_moved before resolution of isMoving in the main loop.
    #pragma omp single
    {
        x_moved += cell_length_x_*params.n_space[0];
        n_moved += params.n_space[0];
    }

    #pragma omp for private(xcall,ycall,mypatch, do_right, receive_flag, old_index) schedule(runtime)
    for (unsigned int ipatch = 0 ; ipatch < vecPatches.size() ; ipatch++) {
         mypatch = vecPatches_old[ipatch];
         do_right = 0;
         receive_flag = 0;
        //If my right neighbor does not belong to me ...
        if (mypatch->MPI_neighborhood_[2+3*(params.nDim_field >= 2)+9*(params.nDim_field == 3)] != mypatch->MPI_neighborhood_[1+3*(params.nDim_field >= 2)+9*(params.nDim_field == 3)]){
            do_right = 1;
            old_index = mypatch->hindex;
            if (mypatch->MPI_neighborhood_[2+3*(params.nDim_field >= 2)+9*(params.nDim_field == 3)] != MPI_PROC_NULL) receive_flag = 1;
        }
        //If my left neighbor does not belong to me ...
        if (mypatch->MPI_neighborhood_[3*(params.nDim_field >= 2)+9*(params.nDim_field == 3)] != mypatch->MPI_neighborhood_[1+3*(params.nDim_field >= 2)+9*(params.nDim_field == 3)]) {
            //... I might have to MPI send myself to the left...
            if (mypatch->MPI_neighborhood_[3*(params.nDim_field >= 2)+9*(params.nDim_field == 3)] != MPI_PROC_NULL)
                cout << "Sending Patch" << endl;
                //Mpi_Send_Patch(MpiLNeighbour);
            //... and, for sure, destroy my data.
            for (unsigned int ispec=0 ; ispec<mypatch->vecSpecies.size(); ispec++) delete (mypatch->vecSpecies[ispec]);
	    mypatch->vecSpecies.clear();
            delete (mypatch->EMfields);
            delete (mypatch->Interp);
            delete (mypatch->Proj);
            cout << "end of data destruction" << endl;
        } else {
            //Else, if my left neighbor belongs to my Patch vector, I become my left neighbor.
            //Nothing to do on global indexes or min_locals. The Patch structure remains the same and unmoved.
            cout << "Moving Left" << endl;
            mypatch->hindex = mypatch->patch_neighborhood_[0+3*(params.nDim_field >= 2)+9*(params.nDim_field == 3)];
            mypatch->Pcoordinates[0] -= 1;

            //Shift neighborhood tables.
	    for ( int z = 0 ; z < 1+2*(params.nDim_field == 3) ; z++ ) {
	        for ( int y = 0 ; y < 1+2*(params.nDim_field >= 2) ; y++ ) {
	            for ( int x = 2 ; x > 0 ; x-- ) {
	            mypatch->patch_neighborhood_[z*9+y*3+x] = mypatch->patch_neighborhood_[z*9+y*3+x-1];
	            mypatch->MPI_neighborhood_[z*9+y*3+x] = mypatch->MPI_neighborhood_[z*9+y*3+x-1];
                    }
                }
            }
            //Compute missing part of the new neighborhood tables.
            xcall = mypatch->Pcoordinates[0]-1;
            ycall = mypatch->Pcoordinates[1]-1;
            if (params.bc_em_type_long=="periodic" && xcall < 0) xcall += (1<<params.mi[0]);
            if (params.bc_em_type_trans=="periodic" && ycall <0) ycall += (1<<params.mi[1]);
	    mypatch->patch_neighborhood_[0] = generalhilbertindex(params.mi[0] , params.mi[1], xcall, ycall);
	    mypatch->patch_neighborhood_[3] = generalhilbertindex(params.mi[0] , params.mi[1], xcall, mypatch->Pcoordinates[1]);
            ycall = mypatch->Pcoordinates[1]+1;
            if (params.bc_em_type_trans=="periodic" && ycall >= 1<<params.mi[1]) ycall -= (1<<params.mi[1]);
	    mypatch->patch_neighborhood_[6] = generalhilbertindex(params.mi[0] , params.mi[1], xcall, ycall);
	    for ( int y = 0 ; y < 1+2*(params.nDim_field >= 2) ; y++ ) {
	        for ( int z = 0 ; z < 1+2*(params.nDim_field == 3) ; z++ ) {
                    mypatch->MPI_neighborhood_[y*3+z] = smpi->hrank(mypatch->patch_neighborhood_[y*3+z]);
                }
            }
            //For now also need to update neighbor_, corner_neighbor and their MPI counterparts even if these will be obsolete eventually.
           mypatch->corner_neighbor_[0][0]= mypatch->patch_neighborhood_[0] ;
           mypatch->neighbor_[1][0]=        mypatch->patch_neighborhood_[1] ;
           mypatch->corner_neighbor_[1][0]= mypatch->patch_neighborhood_[2] ;
           mypatch->neighbor_[0][0]=        mypatch->patch_neighborhood_[3] ;
           mypatch->neighbor_[0][1]=        mypatch->patch_neighborhood_[5] ;
           mypatch->corner_neighbor_[0][1]= mypatch->patch_neighborhood_[6] ;
           mypatch->neighbor_[1][1]=        mypatch->patch_neighborhood_[7] ;
           mypatch->corner_neighbor_[1][1]= mypatch->patch_neighborhood_[8] ;

	   mypatch->MPI_neighbor_[0][0] = mypatch->MPI_neighborhood_[3];
	   mypatch->MPI_neighbor_[0][1] = mypatch->MPI_neighborhood_[5];
	   mypatch->MPI_neighbor_[1][0] = mypatch->MPI_neighborhood_[1];
	   mypatch->MPI_neighbor_[1][1] = mypatch->MPI_neighborhood_[7];
	   mypatch->MPI_corner_neighbor_[0][0] = mypatch->MPI_neighborhood_[0];
	   mypatch->MPI_corner_neighbor_[0][1] = mypatch->MPI_neighborhood_[6];
	   mypatch->MPI_corner_neighbor_[1][0] = mypatch->MPI_neighborhood_[2];
	   mypatch->MPI_corner_neighbor_[1][1] = mypatch->MPI_neighborhood_[8];

            //And finally put the patch at the correct rank in vecPatches.
            vecPatches.patches_[mypatch->hindex - h0 ] = mypatch ; 
             
        }
        //Warning, from here on the neigbors might have been updated !!
        if (do_right) {
            //...je reçois ou je cré.
            if (receive_flag){
                //Mpi_Receive_Patch(MpiRNeighbour);
                cout << "Receiving Patch" << endl;
            } else {
                //Create_Patch();
                cout << "Patch creation n_moved = " << n_moved << " old_index = " << old_index << endl;
                vecPatches.patches_[ipatch] = new Patch(params, diag_params, laser_params, smpi, old_index, n_moved);

            }
        }
        //Les fonctions SendPatch ou ReceivePatch doivent transformer le Patch comme ci dessus.
        //Ce serait plus simple de mettre hindex, neighbor_ et corner_neighbor_ dans un seul et meme tableau.

    }//End loop on Patches
    cout << "end of operate Simwin" << endl;
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
