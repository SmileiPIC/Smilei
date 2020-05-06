
#include "DoubleGridsAM.h"

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
void DoubleGridsAM::syncCurrentsOnRegion( VectorPatch &vecPatches, Region &region, Params &params, SmileiMPI *smpi, Timers &timers, int itime, unsigned int imode )
{
    timers.grids.restart();

    // Loop / additional_patches_ ( = patches included in the local vecPatches but not in local Region, the Region of others MPI will need this data )
    //        additional_patches_ranks stores the MPI rank of the Region which owns additional_patches_
    for ( unsigned int i=0 ; i<region.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = region.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGridsAM::currentsOnRegionSend( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields),
                                             region.additional_patches_[i], region.additional_patches_ranks[i], smpi, vecPatches(ipatch), params, imode );

    }

    // Loop / missing_patches_ ( = patches whose data is needed by the local Region but not own by the local vecPatches )
    //        missing_patches_ranks stores the MPI rank of the vecPacthes which own missing_patches_
    for ( unsigned int i=0 ; i<region.missing_patches_.size() ; i++ ) {

        DoubleGridsAM::currentsOnRegionRecv( static_cast<ElectroMagnAM *>(region.patch_->EMfields)
                                             , region.missing_patches_[i], region.missing_patches_ranks[i], vecPatches, params, smpi, region, imode );

    }

    // Loop / additional_patches_ to finalize send ( currentsOnRegionSend relies on MPI_Isend )
    for ( unsigned int i=0 ; i<region.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = region.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGridsAM::currentsOnRegionSendFinalize( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields)
                                                     , region.additional_patches_[i], region.additional_patches_ranks[i], smpi, vecPatches(ipatch), params, imode );

    }


    ElectroMagnAM * region_fields = static_cast<ElectroMagnAM *>( region.patch_->EMfields );
    // Loop / local_patches_ ( patches own by the local vePatches whose data are used by the local Region )
    for ( unsigned int i=0 ; i<region.local_patches_.size() ; i++ ) {

        unsigned int ipatch = region.local_patches_[i]-vecPatches.refHindex_;        
        ElectroMagnAM * patch_fields = static_cast<ElectroMagnAM *>( vecPatches(ipatch)->EMfields );

        patch_fields->Jl_[imode]->add( region_fields->Jl_[imode], params, smpi, vecPatches(ipatch), region.patch_ );
        patch_fields->Jr_[imode]->add( region_fields->Jr_[imode], params, smpi, vecPatches(ipatch), region.patch_ );
        patch_fields->Jt_[imode]->add( region_fields->Jt_[imode], params, smpi, vecPatches(ipatch), region.patch_ );
        if(params.is_spectral){
            patch_fields->rho_AM_[imode]->add( region_fields->rho_AM_[imode], params, smpi, vecPatches(ipatch), region.patch_ );
            // rho_old is save directly on the Region after the resolution of the Maxwell solver
        }

    }

    timers.grids.update();
}

void DoubleGridsAM::currentsOnRegionSend( ElectroMagnAM* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params, unsigned int imode )
{
    // isendComplex( cFields, targeted_mpi_rank, tag, requests );
    //               tag = *5 ? 5 communications are required per patch : 3 currents + rho + rho_old
    smpi->isendComplex( localfields->Jl_[imode], send_to_global_patch_rank, hindex*5  , patch->requests_[0] );
    smpi->isendComplex( localfields->Jr_[imode], send_to_global_patch_rank, hindex*5+1, patch->requests_[1] );
    smpi->isendComplex( localfields->Jt_[imode], send_to_global_patch_rank, hindex*5+2, patch->requests_[2] );    

    if(params.is_spectral) {
        smpi->isendComplex( localfields->rho_AM_[imode], send_to_global_patch_rank, hindex*5+3, patch->requests_[3] );
    }

}

void DoubleGridsAM::currentsOnRegionSendFinalize( ElectroMagnAM* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params, unsigned int imode )
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

void DoubleGridsAM::currentsOnRegionRecv( ElectroMagnAM* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Region& region, unsigned int imode )
{
    ElectroMagnAM * fake_fields = static_cast<ElectroMagnAM *>( region.fake_patch->EMfields );

    // fake_patch consists in a piece of the local Region to handle naturally patches communications
    //            need to update its hindex and coordinates to put recv data at the good place in the local Region (add)
    region.fake_patch->hindex = hindex;
    region.fake_patch->Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );

    // recvComplex( cFields, sender_mpi_rank, tag );
    //              tag = *5 ? 5 communications are required per patch : 3 currents + rho + rho_old
    smpi->recvComplex( fake_fields->Jl_[imode], local_patch_rank, hindex*5 );
    fake_fields->Jl_[imode]->add( globalfields->Jl_[imode], params, smpi, region.fake_patch, region.patch_ );

    smpi->recvComplex( fake_fields->Jr_[imode], local_patch_rank, hindex*5+1 );
    fake_fields->Jr_[imode]->add( globalfields->Jr_[imode], params, smpi, region.fake_patch, region.patch_ );

    smpi->recvComplex( fake_fields->Jt_[imode], local_patch_rank, hindex*5+2 );
    fake_fields->Jt_[imode]->add( globalfields->Jt_[imode], params, smpi, region.fake_patch, region.patch_ );

    if(params.is_spectral) {
        smpi->recvComplex( fake_fields->rho_AM_[imode], local_patch_rank, hindex*5+3 );
        fake_fields->rho_AM_[imode]->add( globalfields->rho_AM_[imode], params, smpi, region.fake_patch, region.patch_ );
    }

}

// ---------------------------------------------------------------------------
// Scatter Fields on Patches for particles interpolation or divergece cleaning
// ---------------------------------------------------------------------------
void DoubleGridsAM::syncFieldsOnPatches( Region &region, VectorPatch &vecPatches, Params &params, SmileiMPI *smpi, Timers &timers, int itime, unsigned int imode )
{
    timers.grids.restart();

    // Loop / additional_patches_ ( within local vecPatches but not in local Region )
    //                            get data from Region of others MPI
    for ( unsigned int i=0 ; i<region.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = region.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGridsAM::fieldsOnPatchesRecv( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields),
                                            region.additional_patches_[i], region.additional_patches_ranks[i], smpi,  vecPatches(ipatch), params, imode );

    }

    // Loop / missing_patches_ ( within local Region but not in local vecPatches,  )
    //                         send data which do not concern local Region
    for ( unsigned int i=0 ; i<region.missing_patches_.size() ; i++ ) {

        DoubleGridsAM::fieldsOnPatchesSend( static_cast<ElectroMagnAM *>(region.patch_->EMfields),
                                            region.missing_patches_[i], region.missing_patches_ranks[i], vecPatches, params, smpi, region, imode );

    }

    // Loop / additional_patches_ to finalize recv ( fieldsOnPatchesRecv relies on MPI_Irecv )
    for ( unsigned int i=0 ; i<region.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = region.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGridsAM::fieldsOnPatchesRecvFinalize( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields),
                                                    region.additional_patches_[i], region.additional_patches_ranks[i], smpi, vecPatches(ipatch), imode );

    }

    ElectroMagnAM * region_fields = NULL;
    if ( region.local_patches_.size() ) // could be empty for region_global
        region_fields = static_cast<ElectroMagnAM *>( region.patch_->EMfields );

    // Loop / local_patches_ ( patches own by the local vePatches whose data are used by the local Region )
    for ( unsigned int i=0 ; i<region.local_patches_.size() ; i++ ) {

        unsigned int ipatch = region.local_patches_[i]-vecPatches.refHindex_;
        ElectroMagnAM * patch_fields = static_cast<ElectroMagnAM *>( vecPatches(ipatch)->EMfields );

        patch_fields->El_[imode]->get( region_fields->El_[imode], params, smpi, region.patch_, vecPatches(ipatch) );
        patch_fields->Er_[imode]->get( region_fields->Er_[imode], params, smpi, region.patch_, vecPatches(ipatch) );
        patch_fields->Et_[imode]->get( region_fields->Et_[imode], params, smpi, region.patch_, vecPatches(ipatch) );
       
        //Temporary synchronize B_m even if it is not necessary since Bm = B 
        patch_fields->Bl_m[imode]->get( region_fields->Bl_m[imode], params, smpi, region.patch_, vecPatches(ipatch) );
        patch_fields->Br_m[imode]->get( region_fields->Br_m[imode], params, smpi, region.patch_, vecPatches(ipatch) );
        patch_fields->Bt_m[imode]->get( region_fields->Bt_m[imode], params, smpi, region.patch_, vecPatches(ipatch) );

    }

    timers.grids.update();
}


void DoubleGridsAM::fieldsOnPatchesRecv( ElectroMagnAM* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params, unsigned int imode )
{
    // irecvComplex( cFields, sender_mpi_rank, tag, requests );
    //               tag = *9 ? 9 communications could be are required per patch
    //               clarify which usage need B, B_m or both
    smpi->irecvComplex( localfields->El_[imode], recv_from_global_patch_rank, hindex*9  , patch->requests_[0] );
    smpi->irecvComplex( localfields->Er_[imode], recv_from_global_patch_rank, hindex*9+1, patch->requests_[1] );
    smpi->irecvComplex( localfields->Et_[imode], recv_from_global_patch_rank, hindex*9+2, patch->requests_[2] );
   
    smpi->irecvComplex( localfields->Bl_m[imode], recv_from_global_patch_rank, hindex*9+6, patch->requests_[6] );
    smpi->irecvComplex( localfields->Br_m[imode], recv_from_global_patch_rank, hindex*9+7, patch->requests_[7] );
    smpi->irecvComplex( localfields->Bt_m[imode], recv_from_global_patch_rank, hindex*9+8, patch->requests_[8] );

}

void DoubleGridsAM::fieldsOnPatchesRecvFinalize( ElectroMagnAM* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch, unsigned int imode )
{
    MPI_Status status;
    // Wait for fieldsOnPatchesRecv (irecv)
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

void DoubleGridsAM::fieldsOnPatchesSend( ElectroMagnAM* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Region& region, unsigned int imode )
{
    ElectroMagnAM * fake_fields = static_cast<ElectroMagnAM *>( region.fake_patch->EMfields );

    // fake_patch consists in a piece of the local Region to handle naturally patches communications
    //            need to update its hindex and coordinates to extract (get) appropriate data from the local Region data before send it
    region.fake_patch->hindex = hindex;
    region.fake_patch->Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );

    // sendComplex( cFields, targeted_mpi_rank, tag );
    //               tag = *9 ? 9 communications could be required per patch
    //               clarify which usage need B, B_m or both
    fake_fields->El_[imode]->get( globalfields->El_[imode], params, smpi, region.patch_, region.fake_patch );
    smpi->sendComplex( fake_fields->El_[imode], local_patch_rank, hindex*9 );

    fake_fields->Er_[imode]->get( globalfields->Er_[imode], params, smpi, region.patch_, region.fake_patch );
    smpi->sendComplex( fake_fields->Er_[imode], local_patch_rank, hindex*9+1 );

    fake_fields->Et_[imode]->get( globalfields->Et_[imode], params, smpi, region.patch_, region.fake_patch );
    smpi->sendComplex( fake_fields->Et_[imode], local_patch_rank, hindex*9+2 );

    fake_fields->Bl_m[imode]->get( globalfields->Bl_m[imode], params, smpi, region.patch_, region.fake_patch );
    smpi->sendComplex( fake_fields->Bl_m[imode], local_patch_rank, hindex*9+6 );
       
    fake_fields->Br_m[imode]->get( globalfields->Br_m[imode], params, smpi, region.patch_, region.fake_patch );
    smpi->sendComplex( fake_fields->Br_m[imode], local_patch_rank, hindex*9+7 );
       
    fake_fields->Bt_m[imode]->get( globalfields->Bt_m[imode], params, smpi, region.patch_, region.fake_patch );
    smpi->sendComplex( fake_fields->Bt_m[imode], local_patch_rank, hindex*9+8 );

}


// -----------------------------------------------------------
// Gather Fields on Region to compute the divergence cleaning
//                         could be used for the moving window
// -----------------------------------------------------------
void DoubleGridsAM::syncFieldsOnRegion( VectorPatch& vecPatches, Region& region, Params &params, SmileiMPI* smpi, unsigned int imode )
{
    // Loop / additional_patches_ ( = patches included in the local vecPatches but not in local Region, the Region of others MPI will need this data )
    //        additional_patches_ranks stores the MPI rank of the Region which owns additional_patches_
    for ( unsigned int i=0 ; i<region.additional_patches_.size() ; i++ ) {
        
        unsigned int ipatch = region.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGridsAM::fieldsOnRegionSend( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields),
                                           region.additional_patches_[i], region.additional_patches_ranks[i], smpi, vecPatches(ipatch), params, imode );
    }

    // Loop / missing_patches_ ( = patches whose data is needed by the local Region but not own by the local vecPatches )
    //        missing_patches_ranks stores the MPI rank of the vecPacthes which own missing_patches_
    for ( unsigned int i=0 ; i<region.missing_patches_.size() ; i++ ) {
        
        DoubleGridsAM::fieldsOnRegionRecv( static_cast<ElectroMagnAM *>(region.patch_->EMfields),
                                           region.missing_patches_[i], region.missing_patches_ranks[i], vecPatches, params, smpi, region, imode );

    }

    // Loop / additional_patches_ to finalize send ( currentsOnRegionSend relies on MPI_Isend )
    for ( unsigned int i=0 ; i<region.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = region.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGridsAM::fieldsOnRegionSendFinalize( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields),
                                                   region.additional_patches_[i], region.additional_patches_ranks[i], smpi, vecPatches(ipatch), params, imode );
    }

    ElectroMagnAM * region_fields = NULL;
    if ( region.local_patches_.size() )
        region_fields = static_cast<ElectroMagnAM *>( region.patch_->EMfields );

    // Loop / local_patches_ ( patches own by the local vePatches whose data are used by the local Region )
    for ( unsigned int i=0 ; i<region.local_patches_.size() ; i++ ) {

        unsigned int ipatch = region.local_patches_[i]-vecPatches.refHindex_;
        ElectroMagnAM * patch_fields = static_cast<ElectroMagnAM *>( vecPatches(ipatch)->EMfields );

        patch_fields->El_[imode]->put( region_fields->El_[imode], params, smpi, vecPatches(ipatch), region.patch_ );
        patch_fields->Er_[imode]->put( region_fields->Er_[imode], params, smpi, vecPatches(ipatch), region.patch_ );
        patch_fields->Et_[imode]->put( region_fields->Et_[imode], params, smpi, vecPatches(ipatch), region.patch_ );
        
        patch_fields->Bl_[imode]->put( region_fields->Bl_[imode], params, smpi, vecPatches(ipatch), region.patch_ );
        patch_fields->Br_[imode]->put( region_fields->Br_[imode], params, smpi, vecPatches(ipatch), region.patch_ );
        patch_fields->Bt_[imode]->put( region_fields->Bt_[imode], params, smpi, vecPatches(ipatch), region.patch_ );
        
        if (!params.is_spectral) {
            patch_fields->Bl_m[imode]->put( region_fields->Bl_m[imode], params, smpi, vecPatches(ipatch), region.patch_ );
            patch_fields->Br_m[imode]->put( region_fields->Br_m[imode], params, smpi, vecPatches(ipatch), region.patch_ );
            patch_fields->Bt_m[imode]->put( region_fields->Bt_m[imode], params, smpi, vecPatches(ipatch), region.patch_ );
        }
    }
}

void DoubleGridsAM::fieldsOnRegionSend( ElectroMagnAM* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params, unsigned int imode )
{
    // isendComplex( cFields, targeted_mpi_rank, tag, requests );
    //               tag = *9 ? 9 communications could be required per patch : B, B_m, both (at least for the moving window)
    smpi->isendComplex( localfields->El_[imode], send_to_global_patch_rank, hindex*9  , patch->requests_[0] );
    smpi->isendComplex( localfields->Er_[imode], send_to_global_patch_rank, hindex*9+1, patch->requests_[1] );
    smpi->isendComplex( localfields->Et_[imode], send_to_global_patch_rank, hindex*9+2, patch->requests_[2] );

    smpi->isendComplex( localfields->Bl_[imode], send_to_global_patch_rank, hindex*9+3, patch->requests_[3] );
    smpi->isendComplex( localfields->Br_[imode], send_to_global_patch_rank, hindex*9+4, patch->requests_[4] );
    smpi->isendComplex( localfields->Bt_[imode], send_to_global_patch_rank, hindex*9+5, patch->requests_[5] );

    if (!params.is_spectral) {
        smpi->isendComplex( localfields->Bl_m[imode], send_to_global_patch_rank, hindex*9+6, patch->requests_[6] );
        smpi->isendComplex( localfields->Br_m[imode], send_to_global_patch_rank, hindex*9+7, patch->requests_[7] );
        smpi->isendComplex( localfields->Bt_m[imode], send_to_global_patch_rank, hindex*9+8, patch->requests_[8] );
    }
    
}

void DoubleGridsAM::fieldsOnRegionSendFinalize( ElectroMagnAM* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params, unsigned int imode )
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

void DoubleGridsAM::fieldsOnRegionRecv( ElectroMagnAM* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Region& region, unsigned int imode )
{
    ElectroMagnAM * fake_fields = static_cast<ElectroMagnAM *>( region.fake_patch->EMfields );

    // fake_patch consists in a piece of the local Region to handle naturally patches communications
    //            need to update its hindex and coordinates to put recv data at the good place in the local Region (put)
    region.fake_patch->hindex = hindex;
    region.fake_patch->Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );

    // recvComplex( cFields, sender_mpi_rank, tag );
    //              tag = *9 ? 9 communications could be required per patch : B, B_m, both (at least for the moving window)
    smpi->recvComplex( fake_fields->El_[imode], local_patch_rank, hindex*9 );
    fake_fields->El_[imode]->put( globalfields->El_[imode], params, smpi, region.fake_patch, region.patch_ );

    smpi->recvComplex( fake_fields->Er_[imode], local_patch_rank, hindex*9+1 );
    fake_fields->Er_[imode]->put( globalfields->Er_[imode], params, smpi, region.fake_patch, region.patch_ );

    smpi->recvComplex( fake_fields->Et_[imode], local_patch_rank, hindex*9+2 );
    fake_fields->Et_[imode]->put( globalfields->Et_[imode], params, smpi, region.fake_patch, region.patch_ );

    smpi->recvComplex( fake_fields->Bl_[imode], local_patch_rank, hindex*9+3 );
    fake_fields->Bl_[imode]->put( globalfields->Bl_[imode], params, smpi, region.fake_patch, region.patch_ );

    smpi->recvComplex( fake_fields->Br_[imode], local_patch_rank, hindex*9+4 );
    fake_fields->Br_[imode]->put( globalfields->Br_[imode], params, smpi, region.fake_patch, region.patch_ );

    smpi->recvComplex( fake_fields->Bt_[imode], local_patch_rank, hindex*9+5 );
    fake_fields->Bt_[imode]->put( globalfields->Bt_[imode], params, smpi, region.fake_patch, region.patch_ );

    if (!params.is_spectral) {
        smpi->recvComplex( fake_fields->Bl_m[imode], local_patch_rank, hindex*9+6 );
        fake_fields->Bl_m[imode]->put( globalfields->Bl_m[imode], params, smpi, region.fake_patch, region.patch_ );
        
        smpi->recvComplex( fake_fields->Br_m[imode], local_patch_rank, hindex*9+7 );
        fake_fields->Br_m[imode]->put( globalfields->Br_m[imode], params, smpi, region.fake_patch, region.patch_ );
        
        smpi->recvComplex( fake_fields->Bt_m[imode], local_patch_rank, hindex*9+8 );
        fake_fields->Bt_m[imode]->put( globalfields->Bt_m[imode], params, smpi, region.fake_patch, region.patch_ );
    }

}



// ---------------------------------------------------------------------------
// Scatter Fields on Patches for particles interpolation or divergece cleaning
// ---------------------------------------------------------------------------
void DoubleGridsAM::syncBOnPatches( Region &region, VectorPatch &vecPatches, Params &params, SmileiMPI *smpi, Timers &timers, int itime, unsigned int imode )
{
    timers.grids.restart();

    // Loop / additional_patches_ ( within local vecPatches but not in local Region )
    //                            get data from Region of others MPI
    for ( unsigned int i=0 ; i<region.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = region.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGridsAM::bOnPatchesRecv( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields),
                                            region.additional_patches_[i], region.additional_patches_ranks[i], smpi,  vecPatches(ipatch), params, imode );

    }

    // Loop / missing_patches_ ( within local Region but not in local vecPatches,  )
    //                         send data which do not concern local Region
    for ( unsigned int i=0 ; i<region.missing_patches_.size() ; i++ ) {

        DoubleGridsAM::bOnPatchesSend( static_cast<ElectroMagnAM *>(region.patch_->EMfields),
                                            region.missing_patches_[i], region.missing_patches_ranks[i], vecPatches, params, smpi, region, imode );

    }

    // Loop / additional_patches_ to finalize recv ( fieldsOnPatchesRecv relies on MPI_Irecv )
    for ( unsigned int i=0 ; i<region.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = region.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGridsAM::bOnPatchesRecvFinalize( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields),
                                                    region.additional_patches_[i], region.additional_patches_ranks[i], smpi, vecPatches(ipatch), imode );

    }

    ElectroMagnAM * region_fields = NULL;
    if ( region.local_patches_.size() ) // could be empty for region_global
        region_fields = static_cast<ElectroMagnAM *>( region.patch_->EMfields );

    // Loop / local_patches_ ( patches own by the local vePatches whose data are used by the local Region )
    for ( unsigned int i=0 ; i<region.local_patches_.size() ; i++ ) {

        unsigned int ipatch = region.local_patches_[i]-vecPatches.refHindex_;
        ElectroMagnAM * patch_fields = static_cast<ElectroMagnAM *>( vecPatches(ipatch)->EMfields );

        patch_fields->Bl_[imode]->get( region_fields->Bl_[imode], params, smpi, region.patch_, vecPatches(ipatch) );
        patch_fields->Br_[imode]->get( region_fields->Br_[imode], params, smpi, region.patch_, vecPatches(ipatch) );
        patch_fields->Bt_[imode]->get( region_fields->Bt_[imode], params, smpi, region.patch_, vecPatches(ipatch) );

    }

    timers.grids.update();
}


void DoubleGridsAM::bOnPatchesRecv( ElectroMagnAM* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params, unsigned int imode )
{
    // irecvComplex( cFields, sender_mpi_rank, tag, requests );
    //               tag = *9 ? 9 communications could be are required per patch
    //               clarify which usage need B, B_m or both
    smpi->irecvComplex( localfields->Bl_[imode], recv_from_global_patch_rank, hindex*9+3, patch->requests_[3] );
    smpi->irecvComplex( localfields->Br_[imode], recv_from_global_patch_rank, hindex*9+4, patch->requests_[4] );
    smpi->irecvComplex( localfields->Bt_[imode], recv_from_global_patch_rank, hindex*9+5, patch->requests_[5] );

}

void DoubleGridsAM::bOnPatchesRecvFinalize( ElectroMagnAM* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch, unsigned int imode )
{
    MPI_Status status;
    // Wait for fieldsOnPatchesRecv (irecv)
    MPI_Wait( &(patch->requests_[3]), &status );
    MPI_Wait( &(patch->requests_[4]), &status );
    MPI_Wait( &(patch->requests_[5]), &status );
}

void DoubleGridsAM::bOnPatchesSend( ElectroMagnAM* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Region& region, unsigned int imode )
{
    ElectroMagnAM * fake_fields = static_cast<ElectroMagnAM *>( region.fake_patch->EMfields );

    // fake_patch consists in a piece of the local Region to handle naturally patches communications
    //            need to update its hindex and coordinates to extract (get) appropriate data from the local Region data before send it
    region.fake_patch->hindex = hindex;
    region.fake_patch->Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );

    // sendComplex( cFields, targeted_mpi_rank, tag );
    //               tag = *9 ? 9 communications could be required per patch
    //               clarify which usage need B, B_m or both
    fake_fields->Bl_[imode]->get( globalfields->Bl_[imode], params, smpi, region.patch_, region.fake_patch );
    smpi->sendComplex( fake_fields->Bl_[imode], local_patch_rank, hindex*9+3 );

    fake_fields->Br_[imode]->get( globalfields->Br_[imode], params, smpi, region.patch_, region.fake_patch );
    smpi->sendComplex( fake_fields->Br_[imode], local_patch_rank, hindex*9+4 );

    fake_fields->Bt_[imode]->get( globalfields->Bt_[imode], params, smpi, region.patch_, region.fake_patch );
    smpi->sendComplex( fake_fields->Bt_[imode], local_patch_rank, hindex*9+5 );

}


// ---------------------------------------------------------------------------
// Scatter Currents on Patches for diags after filtering
// ---------------------------------------------------------------------------
void DoubleGridsAM::syncCurrentsOnPatches( Region &region, VectorPatch &vecPatches, Params &params, SmileiMPI *smpi, Timers &timers, int itime, unsigned int imode )
{
    timers.grids.restart();

    // Loop / additional_patches_ ( within local vecPatches but not in local Region )
    //                            get data from Region of others MPI
    for ( unsigned int i=0 ; i<region.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = region.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGridsAM::currentsOnPatchesRecv( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields),
                                            region.additional_patches_[i], region.additional_patches_ranks[i], smpi,  vecPatches(ipatch), params, imode );

    }

    // Loop / missing_patches_ ( within local Region but not in local vecPatches,  )
    //                         send data which do not concern local Region
    for ( unsigned int i=0 ; i<region.missing_patches_.size() ; i++ ) {

        DoubleGridsAM::currentsOnPatchesSend( static_cast<ElectroMagnAM *>(region.patch_->EMfields),
                                            region.missing_patches_[i], region.missing_patches_ranks[i], vecPatches, params, smpi, region, imode );

    }

    // Loop / additional_patches_ to finalize recv ( fieldsOnPatchesRecv relies on MPI_Irecv )
    for ( unsigned int i=0 ; i<region.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = region.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGridsAM::currentsOnPatchesRecvFinalize( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields),
                                                    region.additional_patches_[i], region.additional_patches_ranks[i], smpi, vecPatches(ipatch), imode );

    }

    ElectroMagnAM * region_fields = NULL;
    if ( region.local_patches_.size() ) // could be empty for region_global
        region_fields = static_cast<ElectroMagnAM *>( region.patch_->EMfields );

    // Loop / local_patches_ ( patches own by the local vePatches whose data are used by the local Region )
    for ( unsigned int i=0 ; i<region.local_patches_.size() ; i++ ) {

        unsigned int ipatch = region.local_patches_[i]-vecPatches.refHindex_;
        ElectroMagnAM * patch_fields = static_cast<ElectroMagnAM *>( vecPatches(ipatch)->EMfields );

        patch_fields->Jl_[imode]->get( region_fields->Jl_[imode], params, smpi, region.patch_, vecPatches(ipatch) );
        patch_fields->Jr_[imode]->get( region_fields->Jr_[imode], params, smpi, region.patch_, vecPatches(ipatch) );
        patch_fields->Jt_[imode]->get( region_fields->Jt_[imode], params, smpi, region.patch_, vecPatches(ipatch) );
       
        patch_fields->rho_AM_[imode]->get( region_fields->rho_AM_[imode], params, smpi, region.patch_, vecPatches(ipatch) );
    }

    timers.grids.update();
}


void DoubleGridsAM::currentsOnPatchesRecv( ElectroMagnAM* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params, unsigned int imode )
{
    // irecvComplex( cFields, sender_mpi_rank, tag, requests );
    //               tag = *9 ? 9 communications could be are required per patch
    //               clarify which usage need B, B_m or both
    smpi->irecvComplex( localfields->Jl_[imode], recv_from_global_patch_rank, hindex*9  , patch->requests_[0] );
    smpi->irecvComplex( localfields->Jr_[imode], recv_from_global_patch_rank, hindex*9+1, patch->requests_[1] );
    smpi->irecvComplex( localfields->Jt_[imode], recv_from_global_patch_rank, hindex*9+2, patch->requests_[2] );

    smpi->irecvComplex( localfields->rho_AM_[imode], recv_from_global_patch_rank, hindex*9+3, patch->requests_[3] );

}

void DoubleGridsAM::currentsOnPatchesRecvFinalize( ElectroMagnAM* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch, unsigned int imode )
{
    MPI_Status status;
    // Wait for fieldsOnPatchesRecv (irecv)
    MPI_Wait( &(patch->requests_[0]), &status );
    MPI_Wait( &(patch->requests_[1]), &status );
    MPI_Wait( &(patch->requests_[2]), &status );
    MPI_Wait( &(patch->requests_[3]), &status );
}

void DoubleGridsAM::currentsOnPatchesSend( ElectroMagnAM* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Region& region, unsigned int imode )
{
    ElectroMagnAM * fake_fields = static_cast<ElectroMagnAM *>( region.fake_patch->EMfields );

    // fake_patch consists in a piece of the local Region to handle naturally patches communications
    //            need to update its hindex and coordinates to extract (get) appropriate data from the local Region data before send it
    region.fake_patch->hindex = hindex;
    region.fake_patch->Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );

    // sendComplex( cFields, targeted_mpi_rank, tag );
    //               tag = *9 ? 9 communications could be required per patch
    //               clarify which usage need B, B_m or both
    fake_fields->Jl_[imode]->get( globalfields->Jl_[imode], params, smpi, region.patch_, region.fake_patch );
    smpi->sendComplex( fake_fields->Jl_[imode], local_patch_rank, hindex*9 );

    fake_fields->Jr_[imode]->get( globalfields->Jr_[imode], params, smpi, region.patch_, region.fake_patch );
    smpi->sendComplex( fake_fields->Jr_[imode], local_patch_rank, hindex*9+1 );

    fake_fields->Jt_[imode]->get( globalfields->Jt_[imode], params, smpi, region.patch_, region.fake_patch );
    smpi->sendComplex( fake_fields->Jt_[imode], local_patch_rank, hindex*9+2 );

    fake_fields->rho_AM_[imode]->get( globalfields->rho_AM_[imode], params, smpi, region.patch_, region.fake_patch );
    smpi->sendComplex( fake_fields->rho_AM_[imode], local_patch_rank, hindex*9+3 );

}
