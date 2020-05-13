
#include "DoubleGrids.h"

#include <vector>

#include "Region.h"
#include "VectorPatch.h"
#include "Params.h"
#include "SmileiMPI.h"
#include "PatchesFactory.h"

using namespace std;

// ------------------------------------------------------------
// Gather Currents on Region to apply Maxwell solvers on Region
// ------------------------------------------------------------
void DoubleGrids::syncCurrentsOnRegion( VectorPatch &vecPatches, Region &region, Params &params, SmileiMPI *smpi, Timers &timers, int itime )
{
    timers.grids.restart();

    // Loop / additional_patches_ ( = patches included in the local vecPatches but not in local Region, the Region of others MPI will need this data )
    //        additional_patches_ranks stores the MPI rank of the Region which owns additional_patches_
    for ( unsigned int i=0 ; i<region.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = region.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGrids::currentsOnRegionSend( vecPatches(ipatch)->EMfields,
                                           region.additional_patches_[i], region.additional_patches_ranks[i], smpi, vecPatches(ipatch), params );

    }


    // Loop / missing_patches_ ( = patches whose data is needed by the local Region but not own by the local vecPatches )
    //        missing_patches_ranks stores the MPI rank of the vecPacthes which own missing_patches_
    for ( unsigned int i=0 ; i<region.missing_patches_.size() ; i++ ) {

        DoubleGrids::currentsOnRegionRecv( region.patch_->EMfields,
                                           region.missing_patches_[i], region.missing_patches_ranks[i], vecPatches, params, smpi, region );

    }


    // Loop / additional_patches_ to finalize send ( currentsOnRegionSend relies on MPI_Isend )
    for ( unsigned int i=0 ; i<region.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = region.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGrids::currentsOnRegionSendFinalize( vecPatches(ipatch)->EMfields,
                                                   region.additional_patches_[i], region.additional_patches_ranks[i], smpi, vecPatches(ipatch), params );

    }

    // Loop / local_patches_ ( patches own by the local vePatches whose data are used by the local Region )
    for ( unsigned int i=0 ; i<region.local_patches_.size() ; i++ ) {

        unsigned int ipatch = region.local_patches_[i]-vecPatches.refHindex_;
        vecPatches(ipatch)->EMfields->Jx_->add( region.patch_->EMfields->Jx_, params, smpi, vecPatches(ipatch), region.patch_ );
        vecPatches(ipatch)->EMfields->Jy_->add( region.patch_->EMfields->Jy_, params, smpi, vecPatches(ipatch), region.patch_ );
        vecPatches(ipatch)->EMfields->Jz_->add( region.patch_->EMfields->Jz_, params, smpi, vecPatches(ipatch), region.patch_ );

	if(params.is_spectral){
            vecPatches(ipatch)->EMfields->rho_->add( region.patch_->EMfields->rho_, params, smpi, vecPatches(ipatch), region.patch_ );
            // rho_old is save directly on the Region after the resolution of the Maxwell solver
	}

    }
    timers.grids.update();
}

void DoubleGrids::currentsOnRegionSend( ElectroMagn* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params )
{
    // isend( Fields, targeted_mpi_rank, tag, requests );
    //               tag = *5 ? 5 communications are required per patch : 3 currents + rho + rho_old
    smpi->isend( localfields->Jx_, send_to_global_patch_rank, hindex*5  , patch->requests_[0] );
    smpi->isend( localfields->Jy_, send_to_global_patch_rank, hindex*5+1, patch->requests_[1] );
    smpi->isend( localfields->Jz_, send_to_global_patch_rank, hindex*5+2, patch->requests_[2] );

    if(params.is_spectral) {
        smpi->isend( localfields->rho_,    send_to_global_patch_rank, hindex*5+3, patch->requests_[3] );
    }

}

void DoubleGrids::currentsOnRegionSendFinalize( ElectroMagn* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params )
{
    MPI_Status status;
    // Wait for currentsOnRegionSend (isend)
    MPI_Wait( &(patch->requests_[0]), &status );
    MPI_Wait( &(patch->requests_[1]), &status );
    MPI_Wait( &(patch->requests_[2]), &status );

    if(params.is_spectral) {
        MPI_Wait( &(patch->requests_[3]), &status );
    }
}

void DoubleGrids::currentsOnRegionRecv( ElectroMagn* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Region& region )
{
    // fake_patch consists in a piece of the local Region to handle naturally patches communications
    //            need to update its hindex and coordinates to put recv data at the good place in the local Region (add)
    region.fake_patch->hindex = hindex;
    region.fake_patch->Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );

    // recv( Fields, sender_mpi_rank, tag );
    //       tag = *5 ? 5 communications are required per patch : 3 currents + rho + rho_old
    smpi->recv( region.fake_patch->EMfields->Jx_, local_patch_rank, hindex*5 );
    region.fake_patch->EMfields->Jx_->add( globalfields->Jx_, params, smpi, region.fake_patch, region.patch_ );

    smpi->recv( region.fake_patch->EMfields->Jy_, local_patch_rank, hindex*5+1 );
    region.fake_patch->EMfields->Jy_->add( globalfields->Jy_, params, smpi, region.fake_patch, region.patch_ );

    smpi->recv( region.fake_patch->EMfields->Jz_, local_patch_rank, hindex*5+2 );
    region.fake_patch->EMfields->Jz_->add( globalfields->Jz_, params, smpi, region.fake_patch, region.patch_ );

    if(params.is_spectral) {
        smpi->recv( region.fake_patch->EMfields->rho_, local_patch_rank, hindex*5+3 );
        region.fake_patch->EMfields->rho_->add( globalfields->rho_, params, smpi, region.fake_patch, region.patch_ );
    }

}


// ---------------------------------------------------------------------------
// Scatter Fields on Patches for particles interpolation or divergece cleaning
// ---------------------------------------------------------------------------
void DoubleGrids::syncFieldsOnPatches( Region &region, VectorPatch &vecPatches, Params &params, SmileiMPI *smpi, Timers &timers, int itime )
{
    timers.grids.restart();

    // Loop / additional_patches_ ( within local vecPatches but not in local Region )
    //                            get data from Region of others MPI
    for ( unsigned int i=0 ; i<region.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = region.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGrids::fieldsOnPatchesRecv( vecPatches(ipatch)->EMfields,
                                          region.additional_patches_[i], region.additional_patches_ranks[i], smpi,  vecPatches(ipatch) );

    }

    // Loop / missing_patches_ ( within local Region but not in local vecPatches,  )
    //                         send data which do not concern local Region
    for ( unsigned int i=0 ; i<region.missing_patches_.size() ; i++ ) {

        DoubleGrids::fieldsOnPatchesSend( region.patch_->EMfields,
                                          region.missing_patches_[i], region.missing_patches_ranks[i], vecPatches, params, smpi, region );

    }

    // Loop / additional_patches_ to finalize recv ( fieldsOnPatchesRecv relies on MPI_Irecv )
    for ( unsigned int i=0 ; i<region.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = region.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGrids::fieldsOnPatchesRecvFinalize( vecPatches(ipatch)->EMfields,
                                                  region.additional_patches_[i], region.additional_patches_ranks[i], smpi, vecPatches(ipatch) );

    }

    // Loop / local_patches_ ( patches own by the local vePatches whose data are used by the local Region )
    for ( unsigned int i=0 ; i<region.local_patches_.size() ; i++ ) {

        unsigned int ipatch = region.local_patches_[i]-vecPatches.refHindex_;

        vecPatches(ipatch)->EMfields->Ex_->get( region.patch_->EMfields->Ex_, params, smpi, region.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->Ey_->get( region.patch_->EMfields->Ey_, params, smpi, region.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->Ez_->get( region.patch_->EMfields->Ez_, params, smpi, region.patch_, vecPatches(ipatch) );
   
        vecPatches(ipatch)->EMfields->Bx_m->get( region.patch_->EMfields->Bx_m, params, smpi, region.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->By_m->get( region.patch_->EMfields->By_m, params, smpi, region.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->Bz_m->get( region.patch_->EMfields->Bz_m, params, smpi, region.patch_, vecPatches(ipatch) );

    }

    timers.grids.update();
}


void DoubleGrids::fieldsOnPatchesRecv( ElectroMagn* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch )
{
    // irecv( Fields, sender_mpi_rank, tag, requests );
    //        tag = *6 ? 6 communications could be are required per patch
    //        clarify which usage need B, B_m or both
    smpi->irecv( localfields->Ex_, recv_from_global_patch_rank, hindex*6  , patch->requests_[0] );
    smpi->irecv( localfields->Ey_, recv_from_global_patch_rank, hindex*6+1, patch->requests_[1] );
    smpi->irecv( localfields->Ez_, recv_from_global_patch_rank, hindex*6+2, patch->requests_[2] );
   
    smpi->irecv( localfields->Bx_m, recv_from_global_patch_rank, hindex*6+3, patch->requests_[3] );
    smpi->irecv( localfields->By_m, recv_from_global_patch_rank, hindex*6+4, patch->requests_[4] );
    smpi->irecv( localfields->Bz_m, recv_from_global_patch_rank, hindex*6+5, patch->requests_[5] );

}

void DoubleGrids::fieldsOnPatchesRecvFinalize( ElectroMagn* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch )
{
    MPI_Status status;
    // Wait for fieldsOnPatchesRecv (irecv)
    MPI_Wait( &(patch->requests_[0]), &status );
    MPI_Wait( &(patch->requests_[1]), &status );
    MPI_Wait( &(patch->requests_[2]), &status );
    MPI_Wait( &(patch->requests_[3]), &status );
    MPI_Wait( &(patch->requests_[4]), &status );
    MPI_Wait( &(patch->requests_[5]), &status );
}

void DoubleGrids::fieldsOnPatchesSend( ElectroMagn* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Region& region )
{
    // fake_patch consists in a piece of the local Region to handle naturally patches communications
    //            need to update its hindex and coordinates to extract (get) appropriate data from the local Region data before send it
    region.fake_patch->hindex = hindex;
    region.fake_patch->Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );

    // send( Fields, targeted_mpi_rank, tag );
    //       tag = *6 ? 6 communications could be required per patch
    //       clarify which usage need B, B_m or both
    region.fake_patch->EMfields->Ex_->get( globalfields->Ex_, params, smpi, region.patch_, region.fake_patch );
    smpi->send( region.fake_patch->EMfields->Ex_, local_patch_rank, hindex*6 );

    region.fake_patch->EMfields->Ey_->get( globalfields->Ey_, params, smpi, region.patch_, region.fake_patch );
    smpi->send( region.fake_patch->EMfields->Ey_, local_patch_rank, hindex*6+1 );

    region.fake_patch->EMfields->Ez_->get( globalfields->Ez_, params, smpi, region.patch_, region.fake_patch );
    smpi->send( region.fake_patch->EMfields->Ez_, local_patch_rank, hindex*6+2 );

    region.fake_patch->EMfields->Bx_m->get( globalfields->Bx_m, params, smpi, region.patch_, region.fake_patch );
    smpi->send( region.fake_patch->EMfields->Bx_m, local_patch_rank, hindex*6+3 );

    region.fake_patch->EMfields->By_m->get( globalfields->By_m, params, smpi, region.patch_, region.fake_patch );
    smpi->send( region.fake_patch->EMfields->By_m, local_patch_rank, hindex*6+4 );

    region.fake_patch->EMfields->Bz_m->get( globalfields->Bz_m, params, smpi, region.patch_, region.fake_patch );
    smpi->send( region.fake_patch->EMfields->Bz_m, local_patch_rank, hindex*6+5 );

}


// -----------------------------------------------------------
// Gather Fields on Region to compute the divergence cleaning
//                         could be used for the moving window
// -----------------------------------------------------------
void DoubleGrids::syncFieldsOnRegion( VectorPatch& vecPatches, Region& region, Params &params, SmileiMPI* smpi )
{
    // Loop / additional_patches_ ( = patches included in the local vecPatches but not in local Region, the Region of others MPI will need this data )
    //        additional_patches_ranks stores the MPI rank of the Region which owns additional_patches_
    for ( unsigned int i=0 ; i<region.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = region.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGrids::fieldsOnRegionSend( vecPatches(ipatch)->EMfields, region.additional_patches_[i], region.additional_patches_ranks[i], smpi, vecPatches(ipatch), params );
    }

    // Loop / missing_patches_ ( = patches whose data is needed by the local Region but not own by the local vecPatches )
    //        missing_patches_ranks stores the MPI rank of the vecPacthes which own missing_patches_
    for ( unsigned int i=0 ; i<region.missing_patches_.size() ; i++ ) {

        DoubleGrids::fieldsOnRegionRecv( region.patch_->EMfields,
                                         region.missing_patches_[i], region.missing_patches_ranks[i], vecPatches, params, smpi, region );

    }


    // Loop / additional_patches_ to finalize send ( currentsOnRegionSend relies on MPI_Isend )
    for ( unsigned int i=0 ; i<region.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = region.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGrids::fieldsOnRegionSendFinalize( vecPatches(ipatch)->EMfields,
                                                 region.additional_patches_[i], region.additional_patches_ranks[i], smpi, vecPatches(ipatch), params );
    }

    // Loop / local_patches_ ( patches own by the local vePatches whose data are used by the local Region )
    for ( unsigned int i=0 ; i<region.local_patches_.size() ; i++ ) {

        unsigned int ipatch = region.local_patches_[i]-vecPatches.refHindex_;

        vecPatches(ipatch)->EMfields->Ex_->put( region.patch_->EMfields->Ex_, params, smpi, vecPatches(ipatch), region.patch_ );
        vecPatches(ipatch)->EMfields->Ey_->put( region.patch_->EMfields->Ey_, params, smpi, vecPatches(ipatch), region.patch_ );
        vecPatches(ipatch)->EMfields->Ez_->put( region.patch_->EMfields->Ez_, params, smpi, vecPatches(ipatch), region.patch_ );

        vecPatches(ipatch)->EMfields->Bx_->put( region.patch_->EMfields->Bx_, params, smpi, vecPatches(ipatch), region.patch_ );
        vecPatches(ipatch)->EMfields->By_->put( region.patch_->EMfields->By_, params, smpi, vecPatches(ipatch), region.patch_ );
        vecPatches(ipatch)->EMfields->Bz_->put( region.patch_->EMfields->Bz_, params, smpi, vecPatches(ipatch), region.patch_ );

        vecPatches(ipatch)->EMfields->Bx_m->put( region.patch_->EMfields->Bx_m, params, smpi, vecPatches(ipatch), region.patch_ );
        vecPatches(ipatch)->EMfields->By_m->put( region.patch_->EMfields->By_m, params, smpi, vecPatches(ipatch), region.patch_ );
        vecPatches(ipatch)->EMfields->Bz_m->put( region.patch_->EMfields->Bz_m, params, smpi, vecPatches(ipatch), region.patch_ );

    }
}

void DoubleGrids::fieldsOnRegionSend( ElectroMagn* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params )
{
    // isend( Fields, targeted_mpi_rank, tag, requests );
    //        tag = *9 ? 9 communications could be required per patch : B, B_m, both (at least for the moving window)
    smpi->isend( localfields->Ex_, send_to_global_patch_rank, hindex*9  , patch->requests_[0] );
    smpi->isend( localfields->Ey_, send_to_global_patch_rank, hindex*9+1, patch->requests_[1] );
    smpi->isend( localfields->Ez_, send_to_global_patch_rank, hindex*9+2, patch->requests_[2] );

    smpi->isend( localfields->Bx_, send_to_global_patch_rank, hindex*9+3, patch->requests_[3] );
    smpi->isend( localfields->By_, send_to_global_patch_rank, hindex*9+4, patch->requests_[4] );
    smpi->isend( localfields->Bz_, send_to_global_patch_rank, hindex*9+5, patch->requests_[5] );

    smpi->isend( localfields->Bx_m, send_to_global_patch_rank, hindex*9+6, patch->requests_[6] );
    smpi->isend( localfields->By_m, send_to_global_patch_rank, hindex*9+7, patch->requests_[7] );
    smpi->isend( localfields->Bz_m, send_to_global_patch_rank, hindex*9+8, patch->requests_[8] );
    
}

void DoubleGrids::fieldsOnRegionSendFinalize( ElectroMagn* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params )
{
    MPI_Status status;
    // Wait for fieldsOnRegionSend (isend)
    MPI_Wait( &(patch->requests_[0]), &status );
    MPI_Wait( &(patch->requests_[1]), &status );
    MPI_Wait( &(patch->requests_[2]), &status );

    MPI_Wait( &(patch->requests_[3]), &status );
    MPI_Wait( &(patch->requests_[4]), &status );
    MPI_Wait( &(patch->requests_[5]), &status );

    MPI_Wait( &(patch->requests_[6]), &status );
    MPI_Wait( &(patch->requests_[7]), &status );
    MPI_Wait( &(patch->requests_[8]), &status );

}

void DoubleGrids::fieldsOnRegionRecv( ElectroMagn* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Region& region )
{
    // fake_patch consists in a piece of the local Region to handle naturally patches communications
    //            need to update its hindex and coordinates to put recv data at the good place in the local Region (put)
    region.fake_patch->hindex = hindex;
    region.fake_patch->Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );

    // recv( Fields, sender_mpi_rank, tag );
    //       tag = *9 ? 9 communications could be required per patch : B, B_m, both (at least for the moving window)
    smpi->recv( region.fake_patch->EMfields->Ex_, local_patch_rank, hindex*9 );
    region.fake_patch->EMfields->Ex_->put( globalfields->Ex_, params, smpi, region.fake_patch, region.patch_ );

    smpi->recv( region.fake_patch->EMfields->Ey_, local_patch_rank, hindex*9+1 );
    region.fake_patch->EMfields->Ey_->put( globalfields->Ey_, params, smpi, region.fake_patch, region.patch_ );

    smpi->recv( region.fake_patch->EMfields->Ez_, local_patch_rank, hindex*9+2 );
    region.fake_patch->EMfields->Ez_->put( globalfields->Ez_, params, smpi, region.fake_patch, region.patch_ );

    smpi->recv( region.fake_patch->EMfields->Bx_, local_patch_rank, hindex*9+3 );
    region.fake_patch->EMfields->Bx_->put( globalfields->Bx_, params, smpi, region.fake_patch, region.patch_ );

    smpi->recv( region.fake_patch->EMfields->By_, local_patch_rank, hindex*9+4 );
    region.fake_patch->EMfields->By_->put( globalfields->By_, params, smpi, region.fake_patch, region.patch_ );

    smpi->recv( region.fake_patch->EMfields->Bz_, local_patch_rank, hindex*9+5 );
    region.fake_patch->EMfields->Bz_->put( globalfields->Bz_, params, smpi, region.fake_patch, region.patch_ );

    smpi->recv( region.fake_patch->EMfields->Bx_m, local_patch_rank, hindex*9+6 );
    region.fake_patch->EMfields->Bx_m->put( globalfields->Bx_m, params, smpi, region.fake_patch, region.patch_ );

    smpi->recv( region.fake_patch->EMfields->By_m, local_patch_rank, hindex*9+7 );
    region.fake_patch->EMfields->By_m->put( globalfields->By_m, params, smpi, region.fake_patch, region.patch_ );

    smpi->recv( region.fake_patch->EMfields->Bz_m, local_patch_rank, hindex*9+8 );
    region.fake_patch->EMfields->Bz_m->put( globalfields->Bz_m, params, smpi, region.fake_patch, region.patch_ );

}













// ---------------------------------------------------------------------------
// Scatter Fields on Patches for particles interpolation or divergece cleaning
// ---------------------------------------------------------------------------
void DoubleGrids::syncBOnPatches( Region &region, VectorPatch &vecPatches, Params &params, SmileiMPI *smpi, Timers &timers, int itime )
{
    timers.grids.restart();

    // Loop / additional_patches_ ( within local vecPatches but not in local Region )
    //                            get data from Region of others MPI
    for ( unsigned int i=0 ; i<region.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = region.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGrids::bOnPatchesRecv( vecPatches(ipatch)->EMfields,
                                     region.additional_patches_[i], region.additional_patches_ranks[i], smpi,  vecPatches(ipatch) );

    }

    // Loop / missing_patches_ ( within local Region but not in local vecPatches,  )
    //                         send data which do not concern local Region
    for ( unsigned int i=0 ; i<region.missing_patches_.size() ; i++ ) {

        DoubleGrids::bOnPatchesSend( region.patch_->EMfields,
                                     region.missing_patches_[i], region.missing_patches_ranks[i], vecPatches, params, smpi, region );

    }

    // Loop / additional_patches_ to finalize recv ( fieldsOnPatchesRecv relies on MPI_Irecv )
    for ( unsigned int i=0 ; i<region.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = region.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGrids::bOnPatchesRecvFinalize( vecPatches(ipatch)->EMfields,
                                             region.additional_patches_[i], region.additional_patches_ranks[i], smpi, vecPatches(ipatch) );

    }

    // Loop / local_patches_ ( patches own by the local vePatches whose data are used by the local Region )
    for ( unsigned int i=0 ; i<region.local_patches_.size() ; i++ ) {

        unsigned int ipatch = region.local_patches_[i]-vecPatches.refHindex_;

        vecPatches(ipatch)->EMfields->Bx_->get( region.patch_->EMfields->Bx_, params, smpi, region.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->By_->get( region.patch_->EMfields->By_, params, smpi, region.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->Bz_->get( region.patch_->EMfields->Bz_, params, smpi, region.patch_, vecPatches(ipatch) );

    }

    timers.grids.update();
}


void DoubleGrids::bOnPatchesRecv( ElectroMagn* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch )
{
    // irecv( Fields, sender_mpi_rank, tag, requests );
    //        tag = *6 ? 6 communications could be are required per patch
    //        clarify which usage need B, B_m or both
    smpi->irecv( localfields->Bx_, recv_from_global_patch_rank, hindex*9+6, patch->requests_[6] );
    smpi->irecv( localfields->By_, recv_from_global_patch_rank, hindex*9+7, patch->requests_[7] );
    smpi->irecv( localfields->Bz_, recv_from_global_patch_rank, hindex*9+8, patch->requests_[8] );

}

void DoubleGrids::bOnPatchesRecvFinalize( ElectroMagn* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch )
{
    MPI_Status status;
    // Wait for fieldsOnPatchesRecv (irecv)
    MPI_Wait( &(patch->requests_[6]), &status );
    MPI_Wait( &(patch->requests_[7]), &status );
    MPI_Wait( &(patch->requests_[8]), &status );
}

void DoubleGrids::bOnPatchesSend( ElectroMagn* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Region& region )
{
    // fake_patch consists in a piece of the local Region to handle naturally patches communications
    //            need to update its hindex and coordinates to extract (get) appropriate data from the local Region data before send it
    region.fake_patch->hindex = hindex;
    region.fake_patch->Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );

    // send( Fields, targeted_mpi_rank, tag );
    //       tag = *9 ? 6 communications could be required per patch
    //       clarify which usage need B, B_m or both
    region.fake_patch->EMfields->Bx_->get( globalfields->Bx_, params, smpi, region.patch_, region.fake_patch );
    smpi->send( region.fake_patch->EMfields->Bx_, local_patch_rank, hindex*9+6 );

    region.fake_patch->EMfields->By_->get( globalfields->By_, params, smpi, region.patch_, region.fake_patch );
    smpi->send( region.fake_patch->EMfields->By_, local_patch_rank, hindex*9+7 );

    region.fake_patch->EMfields->Bz_->get( globalfields->Bz_, params, smpi, region.patch_, region.fake_patch );
    smpi->send( region.fake_patch->EMfields->Bz_, local_patch_rank, hindex*9+8 );

}

// ---------------------------------------------------------------------------
// Scatter Currents on Patches for diags after filtering
// ---------------------------------------------------------------------------
void DoubleGrids::syncCurrentsOnPatches( Region &region, VectorPatch &vecPatches, Params &params, SmileiMPI *smpi, Timers &timers, int itime )
{
    timers.grids.restart();

    // Loop / additional_patches_ ( within local vecPatches but not in local Region )
    //                            get data from Region of others MPI
    for ( unsigned int i=0 ; i<region.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = region.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGrids::currentsOnPatchesRecv( vecPatches(ipatch)->EMfields,
                                          region.additional_patches_[i], region.additional_patches_ranks[i], smpi,  vecPatches(ipatch) );

    }

    // Loop / missing_patches_ ( within local Region but not in local vecPatches,  )
    //                         send data which do not concern local Region
    for ( unsigned int i=0 ; i<region.missing_patches_.size() ; i++ ) {

        DoubleGrids::currentsOnPatchesSend( region.patch_->EMfields,
                                          region.missing_patches_[i], region.missing_patches_ranks[i], vecPatches, params, smpi, region );

    }

    // Loop / additional_patches_ to finalize recv ( fieldsOnPatchesRecv relies on MPI_Irecv )
    for ( unsigned int i=0 ; i<region.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = region.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGrids::currentsOnPatchesRecvFinalize( vecPatches(ipatch)->EMfields,
                                                  region.additional_patches_[i], region.additional_patches_ranks[i], smpi, vecPatches(ipatch) );

    }

    // Loop / local_patches_ ( patches own by the local vePatches whose data are used by the local Region )
    for ( unsigned int i=0 ; i<region.local_patches_.size() ; i++ ) {

        unsigned int ipatch = region.local_patches_[i]-vecPatches.refHindex_;

        vecPatches(ipatch)->EMfields->Jx_->get( region.patch_->EMfields->Jx_, params, smpi, region.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->Jy_->get( region.patch_->EMfields->Jy_, params, smpi, region.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->Jz_->get( region.patch_->EMfields->Jz_, params, smpi, region.patch_, vecPatches(ipatch) );

        vecPatches(ipatch)->EMfields->rho_->get( region.patch_->EMfields->rho_, params, smpi, region.patch_, vecPatches(ipatch) );

    }

    timers.grids.update();
}


void DoubleGrids::currentsOnPatchesRecv( ElectroMagn* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch )
{
    // irecv( Fields, sender_mpi_rank, tag, requests );
    //        tag = *6 ? 6 communications could be are required per patch
    //        clarify which usage need B, B_m or both
    smpi->irecv( localfields->Jx_, recv_from_global_patch_rank, hindex*6  , patch->requests_[0] );
    smpi->irecv( localfields->Jy_, recv_from_global_patch_rank, hindex*6+1, patch->requests_[1] );
    smpi->irecv( localfields->Jz_, recv_from_global_patch_rank, hindex*6+2, patch->requests_[2] );

    smpi->irecv( localfields->rho_, recv_from_global_patch_rank, hindex*6+3, patch->requests_[3] );
}

void DoubleGrids::currentsOnPatchesRecvFinalize( ElectroMagn* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch )
{
    MPI_Status status;
    // Wait for fieldsOnPatchesRecv (irecv)
    MPI_Wait( &(patch->requests_[0]), &status );
    MPI_Wait( &(patch->requests_[1]), &status );
    MPI_Wait( &(patch->requests_[2]), &status );


    MPI_Wait( &(patch->requests_[3]), &status );
}

void DoubleGrids::currentsOnPatchesSend( ElectroMagn* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Region& region )
{
    // fake_patch consists in a piece of the local Region to handle naturally patches communications
    //            need to update its hindex and coordinates to extract (get) appropriate data from the local Region data before send it
    region.fake_patch->hindex = hindex;
    region.fake_patch->Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );

    // send( Fields, targeted_mpi_rank, tag );
    //       tag = *6 ? 6 communications could be required per patch
    //       clarify which usage need B, B_m or both
    region.fake_patch->EMfields->Jx_->get( globalfields->Jx_, params, smpi, region.patch_, region.fake_patch );
    smpi->send( region.fake_patch->EMfields->Jx_, local_patch_rank, hindex*6 );

    region.fake_patch->EMfields->Jy_->get( globalfields->Jy_, params, smpi, region.patch_, region.fake_patch );
    smpi->send( region.fake_patch->EMfields->Jy_, local_patch_rank, hindex*6+1 );

    region.fake_patch->EMfields->Jz_->get( globalfields->Jz_, params, smpi, region.patch_, region.fake_patch );
    smpi->send( region.fake_patch->EMfields->Jz_, local_patch_rank, hindex*6+2 );

    region.fake_patch->EMfields->rho_->get( globalfields->rho_, params, smpi, region.patch_, region.fake_patch );
    smpi->send( region.fake_patch->EMfields->rho_, local_patch_rank, hindex*6+3 );

}
