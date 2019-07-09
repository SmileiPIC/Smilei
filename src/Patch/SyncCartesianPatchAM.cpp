
#include "SyncCartesianPatchAM.h"

#include <vector>

#include "Domain.h"
#include "VectorPatch.h"
#include "Params.h"
#include "SmileiMPI.h"
#include "PatchesFactory.h"

using namespace std;

void SyncCartesianPatchAM::sendCurrentsToDomain( VectorPatch &vecPatches, Domain &domain, Params &params, SmileiMPI *smpi, Timers &timers, int itime, unsigned int imode )
{
    timers.grids.restart();

    // Loop / local_patches_ -> OK
    //     put()

    // Loop / additional_patches_ (identify where goes the additionnal patches regarding patch_count)
    //     patch().send()
    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {
        //cout << smpi->getRank() << " will send " << domain.additional_patches_[i] << " to " << domain.additional_patches_ranks[i] << endl;
        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;

        SyncCartesianPatchAM::sendPatchedToCartesian( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields), domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi, vecPatches(ipatch), params, imode );

        // who = smpi->hrank_THEORIQUE( domain.additional_patches_[i] );
        //     MANAGED WITH --> domain.additional_patches_ranks !!!
        // Send( vecPatches(ipatch)->EMfields->Jx_, who )
    }


    // Loop / missing_patches_ (identify where are the missing patches)
    //     domain.recv( domain.missing_patches_[i] from smpi->hrank( domain.missing_patches_[i] ) )
    //     put to domain
    for ( unsigned int i=0 ; i<domain.missing_patches_.size() ; i++ ) {
        unsigned int ipatch = domain.missing_patches_[i]-vecPatches.refHindex_;
        //cout << smpi->getRank() << " will recv " << domain.missing_patches_[i] << " from " << domain.missing_patches_ranks[i] << endl;

        SyncCartesianPatchAM::recvPatchedToCartesian( static_cast<ElectroMagnAM *>(domain.patch_->EMfields), domain.missing_patches_[i], domain.missing_patches_ranks[i], vecPatches, params, smpi, domain, imode );

    }


    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {
        //cout << smpi->getRank() << " will send " << domain.additional_patches_[i] << endl;
        //cout << smpi->getRank() << " finaliser send " << domain.additional_patches_[i] << " to " << domain.additional_patches_ranks[i] << endl;
        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;

        SyncCartesianPatchAM::finalize_sendPatchedToCartesian( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields), domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi, vecPatches(ipatch), params, imode );

        // who = smpi->hrank_THEORIQUE( domain.additional_patches_[i] );
        //     MANAGED WITH --> domain.additional_patches_ranks !!!
        // Send( vecPatches(ipatch)->EMfields->Jx_, who )
    }

    //cout << smpi->getRank() << " - "; 
    //for ( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {

    ElectroMagnAM * domain_fields = static_cast<ElectroMagnAM *>( domain.patch_->EMfields );

    for ( unsigned int i=0 ; i<domain.local_patches_.size() ; i++ ) {
        unsigned int ipatch = domain.local_patches_[i]-vecPatches.refHindex_;
        
        ElectroMagnAM * patch_fields = static_cast<ElectroMagnAM *>( vecPatches(ipatch)->EMfields );

        patch_fields->Jl_[imode]->put( domain_fields->Jl_[imode], params, smpi, vecPatches(ipatch), domain.patch_ );
        patch_fields->Jr_[imode]->put( domain_fields->Jr_[imode], params, smpi, vecPatches(ipatch), domain.patch_ );
        patch_fields->Jt_[imode]->put( domain_fields->Jt_[imode], params, smpi, vecPatches(ipatch), domain.patch_ );
        if(params.is_spectral){
            patch_fields->rho_AM_[imode]->put( domain_fields->rho_AM_[imode],       params, smpi, vecPatches(ipatch), domain.patch_ );

    // useless rho_old is save directly on vecPatches concerned by the Maxwell soler see VectorPatches::solveMaxwell()
    //vecPatches(ipatch)->EMfields->rhoold_->put( domain.patch_->EMfields->rhoold_, params, smpi, vecPatches(ipatch), domain.patch_ );
 }

    }
    //cout << endl;
    timers.grids.update();
}

void SyncCartesianPatchAM::sendPatchedToCartesian( ElectroMagnAM* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params, unsigned int imode )
{
    //smpi->send( localfields->Jx_, hindex, send_to_global_patch_rank );
    //    isend( EM->Bz_m, to, mpi_tag+tag, requests[tag]); tag++;
    smpi->isendComplex( localfields->Jl_[imode], send_to_global_patch_rank, hindex*5  , patch->requests_[0] );
    smpi->isendComplex( localfields->Jr_[imode], send_to_global_patch_rank, hindex*5+1, patch->requests_[1] );
    smpi->isendComplex( localfields->Jt_[imode], send_to_global_patch_rank, hindex*5+2, patch->requests_[2] );
    

    if(params.is_spectral) {
        smpi->isendComplex( localfields->rho_AM_[imode],    send_to_global_patch_rank, hindex*5+3, patch->requests_[3] );
        //smpi->isend( localfields->rhoold_, send_to_global_patch_rank, hindex*5+4, patch->requests_[4] );
    }

}

void SyncCartesianPatchAM::finalize_sendPatchedToCartesian( ElectroMagnAM* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params, unsigned int imode )
{
    MPI_Status status;
    MPI_Wait( &(patch->requests_[0]), &status );
    MPI_Wait( &(patch->requests_[1]), &status );
    MPI_Wait( &(patch->requests_[2]), &status );
    if(params.is_spectral) {
        MPI_Wait( &(patch->requests_[3]), &status );
        //MPI_Wait( &(patch->requests_[4]), &status ); 
    }
}

void SyncCartesianPatchAM::recvPatchedToCartesian( ElectroMagnAM* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Domain& domain, unsigned int imode )
{
    // Jx_
    // define fake_patch
    unsigned int n_moved = 0;
    //Patch* fake_patch = PatchesFactory::clone(vecPatches(0), params, smpi, vecPatches.domain_decomposition_, hindex, n_moved, false);
    //smpi->recv( fake_patch->Jx_, hindex, local_patch_rank );
    //    recv(  EM->Bz_m, from, tag ); tag++;

    //smpi->recv( fake_patch->EMfields->Jx_, local_patch_rank, hindex );
    ElectroMagnAM * fake_fields = static_cast<ElectroMagnAM *>( domain.fake_patch->EMfields );

    //  vecPatches(ipatch) -> need sender patch coordinates : vecPatches.getDomainCoordinates( hindex )
    // Buffer will be resized for each component, fake local patch, which wil have Jxyz, and coordinates to update ?
    domain.fake_patch->hindex = hindex;
    domain.fake_patch->Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );

    smpi->recvComplex( fake_fields->Jl_[imode], local_patch_rank, hindex*5 );
    fake_fields->Jl_[imode]->put( globalfields->Jl_[imode], params, smpi, domain.fake_patch, domain.patch_ );

    smpi->recvComplex( fake_fields->Jr_[imode], local_patch_rank, hindex*5+1 );
    fake_fields->Jr_[imode]->put( globalfields->Jr_[imode], params, smpi, domain.fake_patch, domain.patch_ );

    smpi->recvComplex( fake_fields->Jt_[imode], local_patch_rank, hindex*5+2 );
    fake_fields->Jt_[imode]->put( globalfields->Jt_[imode], params, smpi, domain.fake_patch, domain.patch_ );

    if(params.is_spectral) {
        smpi->recvComplex( fake_fields->rho_AM_[imode], local_patch_rank, hindex*5+3 );
        fake_fields->rho_AM_[imode]->put( globalfields->rho_AM_[imode], params, smpi, domain.fake_patch, domain.patch_ );
        //smpi->recv( domain.fake_patch->EMfields->rhoold_, local_patch_rank, hindex*5+4 );
        //domain.fake_patch->EMfields->rhoold_->put( globalfields->rhoold_, params, smpi, domain.fake_patch, domain.patch_ );
    }

    //delete fake_patch;

}


void SyncCartesianPatchAM::recvFieldsFromDomain( Domain &domain, VectorPatch &vecPatches, Params &params, SmileiMPI *smpi, Timers &timers, int itime, unsigned int imode )
{
    timers.grids.restart();

    // Loop / additional_patches_ (regarding cartesian), get data from cartesian
    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {
        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;

        SyncCartesianPatchAM::recvCartesianToPatches( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields), domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi,  vecPatches(ipatch), imode );

    }

    // Loop / missing_patches_ (regarding cartesian), send data which do not concern myself
    //     patch().send()
    for ( unsigned int i=0 ; i<domain.missing_patches_.size() ; i++ ) {
        unsigned int ipatch = domain.missing_patches_[i]-vecPatches.refHindex_;

        SyncCartesianPatchAM::sendCartesianToPatches( static_cast<ElectroMagnAM *>(domain.patch_->EMfields), domain.missing_patches_[i], domain.missing_patches_ranks[i], vecPatches, params, smpi, domain, imode );

    }

    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {
        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;

        SyncCartesianPatchAM::finalize_recvCartesianToPatches( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields), domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi, vecPatches(ipatch), imode );

    }

    ElectroMagnAM * domain_fields = NULL;
    if ( domain.local_patches_.size() )
        domain_fields = static_cast<ElectroMagnAM *>( domain.patch_->EMfields );

    //for ( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
    for ( unsigned int i=0 ; i<domain.local_patches_.size() ; i++ ) {
        unsigned int ipatch = domain.local_patches_[i]-vecPatches.refHindex_;

        ElectroMagnAM * patch_fields = static_cast<ElectroMagnAM *>( vecPatches(ipatch)->EMfields );

        patch_fields->El_[imode]->get( domain_fields->El_[imode], params, smpi, domain.patch_, vecPatches(ipatch) );
        patch_fields->Er_[imode]->get( domain_fields->Er_[imode], params, smpi, domain.patch_, vecPatches(ipatch) );
        patch_fields->Et_[imode]->get( domain_fields->Et_[imode], params, smpi, domain.patch_, vecPatches(ipatch) );
        
        if (!params.is_spectral) {
            patch_fields->Bl_m[imode]->get( domain_fields->Bl_m[imode], params, smpi, domain.patch_, vecPatches(ipatch) );
            patch_fields->Br_m[imode]->get( domain_fields->Br_m[imode], params, smpi, domain.patch_, vecPatches(ipatch) );
            patch_fields->Bt_m[imode]->get( domain_fields->Bt_m[imode], params, smpi, domain.patch_, vecPatches(ipatch) );
        }

        patch_fields->Bl_[imode]->get( domain_fields->Bl_[imode], params, smpi, domain.patch_, vecPatches(ipatch) );
        patch_fields->Br_[imode]->get( domain_fields->Br_[imode], params, smpi, domain.patch_, vecPatches(ipatch) );
        patch_fields->Bt_[imode]->get( domain_fields->Bt_[imode], params, smpi, domain.patch_, vecPatches(ipatch) );

//        vecPatches(ipatch)->EMfields->Bx_->get( domain.patch_->EMfields->Bx_, params, smpi, domain.patch_, vecPatches(ipatch) );
//        vecPatches(ipatch)->EMfields->By_->get( domain.patch_->EMfields->By_, params, smpi, domain.patch_, vecPatches(ipatch) );
//        vecPatches(ipatch)->EMfields->Bz_->get( domain.patch_->EMfields->Bz_, params, smpi, domain.patch_, vecPatches(ipatch) );

//        vecPatches(ipatch)->EMfields->Jx_->get( domain.patch_->EMfields->Jx_, params, smpi, domain.patch_, vecPatches(ipatch) );
//        vecPatches(ipatch)->EMfields->Jy_->get( domain.patch_->EMfields->Jy_, params, smpi, domain.patch_, vecPatches(ipatch) );
//        vecPatches(ipatch)->EMfields->Jz_->get( domain.patch_->EMfields->Jz_, params, smpi, domain.patch_, vecPatches(ipatch) );
//	if(params.is_spectral) {
//          vecPatches(ipatch)->EMfields->rho_->get( domain.patch_->EMfields->rho_, params, smpi, domain.patch_, vecPatches(ipatch) );
//          vecPatches(ipatch)->EMfields->rhoold_->get( domain.patch_->EMfields->rhoold_, params, smpi, domain.patch_,
//	  vecPatches(ipatch) );
//	}

    }

    timers.grids.update();
}


void SyncCartesianPatchAM::recvCartesianToPatches( ElectroMagnAM* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch, unsigned int imode )
{
    smpi->irecvComplex( localfields->El_[imode], recv_from_global_patch_rank, hindex*9  , patch->requests_[0] );
    smpi->irecvComplex( localfields->Er_[imode], recv_from_global_patch_rank, hindex*9+1, patch->requests_[1] );
    smpi->irecvComplex( localfields->Et_[imode], recv_from_global_patch_rank, hindex*9+2, patch->requests_[2] );
   
    smpi->irecvComplex( localfields->Bl_[imode], recv_from_global_patch_rank, hindex*9+3, patch->requests_[3] );
    smpi->irecvComplex( localfields->Br_[imode], recv_from_global_patch_rank, hindex*9+4, patch->requests_[4] );
    smpi->irecvComplex( localfields->Bt_[imode], recv_from_global_patch_rank, hindex*9+5, patch->requests_[5] );

    if (localfields->Bl_m[imode]->cdata_!=localfields->Bl_m[imode]->cdata_) {
        smpi->irecvComplex( localfields->Bl_m[imode], recv_from_global_patch_rank, hindex*9+6, patch->requests_[6] );
        smpi->irecvComplex( localfields->Br_m[imode], recv_from_global_patch_rank, hindex*9+7, patch->requests_[7] );
        smpi->irecvComplex( localfields->Bt_m[imode], recv_from_global_patch_rank, hindex*9+8, patch->requests_[8] );
    }

}

void SyncCartesianPatchAM::finalize_recvCartesianToPatches( ElectroMagnAM* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch, unsigned int imode )
{
    MPI_Status status;
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

void SyncCartesianPatchAM::sendCartesianToPatches( ElectroMagnAM* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Domain& domain, unsigned int imode )
{
    // Jx_
    // define fake_patch
    unsigned int n_moved = 0;
    //Patch* fake_patch = PatchesFactory::clone(vecPatches(0), params, smpi, vecPatches.domain_decomposition_, hindex, n_moved, false);
    //smpi->recv( fake_patch->Jx_, hindex, local_patch_rank );
    //    recv(  EM->Bz_m, from, tag ); tag++;

    //smpi->recv( fake_patch->EMfields->Jx_, local_patch_rank, hindex );
    ElectroMagnAM * fake_fields = static_cast<ElectroMagnAM *>( domain.fake_patch->EMfields );

    //  vecPatches(ipatch) -> need sender patch coordinates : vecPatches.getDomainCoordinates( hindex )
    // Buffer will be resized for each component, fake local patch, which wil have Jxyz, and coordinates to update ?


    domain.fake_patch->hindex = hindex;
    domain.fake_patch->Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );

    fake_fields->El_[imode]->get( globalfields->El_[imode], params, smpi, domain.patch_, domain.fake_patch );
    smpi->sendComplex( fake_fields->El_[imode], local_patch_rank, hindex*9 );

    fake_fields->Er_[imode]->get( globalfields->Er_[imode], params, smpi, domain.patch_, domain.fake_patch );
    smpi->sendComplex( fake_fields->Er_[imode], local_patch_rank, hindex*9+1 );

    fake_fields->Et_[imode]->get( globalfields->Et_[imode], params, smpi, domain.patch_, domain.fake_patch );
    smpi->sendComplex( fake_fields->Et_[imode], local_patch_rank, hindex*9+2 );

    fake_fields->Bl_[imode]->get( globalfields->Bl_[imode], params, smpi, domain.patch_, domain.fake_patch );
    smpi->sendComplex( fake_fields->Bl_[imode], local_patch_rank, hindex*9+3 );

    fake_fields->Br_[imode]->get( globalfields->Br_[imode], params, smpi, domain.patch_, domain.fake_patch );
    smpi->sendComplex( fake_fields->Br_[imode], local_patch_rank, hindex*9+4 );

    fake_fields->Bt_[imode]->get( globalfields->Bt_[imode], params, smpi, domain.patch_, domain.fake_patch );
    smpi->sendComplex( fake_fields->Bt_[imode], local_patch_rank, hindex*9+5 );

    if (globalfields->Bl_m[imode]->cdata_!=globalfields->Bl_[imode]->cdata_) {
       fake_fields->Bl_m[imode]->get( globalfields->Bl_m[imode], params, smpi, domain.patch_, domain.fake_patch );
       smpi->sendComplex( fake_fields->Bl_m[imode], local_patch_rank, hindex*9+6 );
       
       fake_fields->Br_m[imode]->get( globalfields->Br_m[imode], params, smpi, domain.patch_, domain.fake_patch );
       smpi->sendComplex( fake_fields->Br_m[imode], local_patch_rank, hindex*9+7 );
       
       fake_fields->Bt_m[imode]->get( globalfields->Bt_m[imode], params, smpi, domain.patch_, domain.fake_patch );
       smpi->sendComplex( fake_fields->Bt_m[imode], local_patch_rank, hindex*9+8 );
    }

    //if(params.is_spectral){
    //    smpi->send( localfields->rho_, hindex, global_patch_rank, params, smpi );
    //
    //}

    //delete fake_patch;
}




// ---------------------------------
// --------------- MW --------------
// ---------------------------------


void SyncCartesianPatchAM::sendFieldsToDomain( VectorPatch& vecPatches, Domain& domain, Params &params, SmileiMPI* smpi, unsigned int imode )
{

    // Loop / local_patches_ -> OK
    //     put()

    // Loop / additional_patches_ (identify where goes the additionnal patches regarding patch_count)
    //     patch().send()
    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {
        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;
        SyncCartesianPatchAM::sendPatchedToCartesian_MW( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields), domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi, vecPatches(ipatch), params, imode );
    }


    // Loop / missing_patches_ (identify where are the missing patches)
    //     domain.recv( domain.missing_patches_[i] from smpi->hrank( domain.missing_patches_[i] ) )
    //     put to domain
    for ( unsigned int i=0 ; i<domain.missing_patches_.size() ; i++ ) {
        unsigned int ipatch = domain.missing_patches_[i]-vecPatches.refHindex_;
        SyncCartesianPatchAM::recvPatchedToCartesian_MW( static_cast<ElectroMagnAM *>(domain.patch_->EMfields), domain.missing_patches_[i], domain.missing_patches_ranks[i], vecPatches, params, smpi, domain, imode );

    }


    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {
        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;
        SyncCartesianPatchAM::finalize_sendPatchedToCartesian_MW( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields), domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi, vecPatches(ipatch), params, imode );
    }

    ElectroMagnAM * domain_fields = NULL;
    if ( domain.local_patches_.size() )
        domain_fields = static_cast<ElectroMagnAM *>( domain.patch_->EMfields );

    for ( unsigned int i=0 ; i<domain.local_patches_.size() ; i++ ) {
        unsigned int ipatch = domain.local_patches_[i]-vecPatches.refHindex_;

        ElectroMagnAM * patch_fields = static_cast<ElectroMagnAM *>( vecPatches(ipatch)->EMfields );

        patch_fields->El_[imode]->put( domain_fields->El_[imode], params, smpi, vecPatches(ipatch), domain.patch_ );
        patch_fields->Er_[imode]->put( domain_fields->Er_[imode], params, smpi, vecPatches(ipatch), domain.patch_ );
        patch_fields->Et_[imode]->put( domain_fields->Et_[imode], params, smpi, vecPatches(ipatch), domain.patch_ );
        
        patch_fields->Bl_[imode]->put( domain_fields->Bl_[imode], params, smpi, vecPatches(ipatch), domain.patch_ );
        patch_fields->Br_[imode]->put( domain_fields->Br_[imode], params, smpi, vecPatches(ipatch), domain.patch_ );
        patch_fields->Bt_[imode]->put( domain_fields->Bt_[imode], params, smpi, vecPatches(ipatch), domain.patch_ );
        
        if (!params.is_spectral) {
            patch_fields->Bl_m[imode]->put( domain_fields->Bl_m[imode], params, smpi, vecPatches(ipatch), domain.patch_ );
            patch_fields->Br_m[imode]->put( domain_fields->Br_m[imode], params, smpi, vecPatches(ipatch), domain.patch_ );
            patch_fields->Bt_m[imode]->put( domain_fields->Bt_m[imode], params, smpi, vecPatches(ipatch), domain.patch_ );
        }
    }
}

void SyncCartesianPatchAM::sendPatchedToCartesian_MW( ElectroMagnAM* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params, unsigned int imode )
{
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

void SyncCartesianPatchAM::finalize_sendPatchedToCartesian_MW( ElectroMagnAM* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params, unsigned int imode )
{
    MPI_Status status;
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

void SyncCartesianPatchAM::recvPatchedToCartesian_MW( ElectroMagnAM* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Domain& domain, unsigned int imode )
{
    // Jx_
    unsigned int n_moved = 0;

    ElectroMagnAM * fake_fields = static_cast<ElectroMagnAM *>( domain.fake_patch->EMfields );

    //  vecPatches(ipatch) -> need sender patch coordinates : vecPatches.getDomainCoordinates( hindex )
    // Buffer will be resized for each component, fake local patch, which wil have Jxyz, and coordinates to update ?
    domain.fake_patch->hindex = hindex;
    domain.fake_patch->Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );

    smpi->recvComplex( fake_fields->El_[imode], local_patch_rank, hindex*9 );
    fake_fields->El_[imode]->put( globalfields->El_[imode], params, smpi, domain.fake_patch, domain.patch_ );

    smpi->recvComplex( fake_fields->Er_[imode], local_patch_rank, hindex*9+1 );
    fake_fields->Er_[imode]->put( globalfields->Er_[imode], params, smpi, domain.fake_patch, domain.patch_ );

    smpi->recvComplex( fake_fields->Et_[imode], local_patch_rank, hindex*9+2 );
    fake_fields->Et_[imode]->put( globalfields->Et_[imode], params, smpi, domain.fake_patch, domain.patch_ );

    smpi->recvComplex( fake_fields->Bl_[imode], local_patch_rank, hindex*9+3 );
    fake_fields->Bl_[imode]->put( globalfields->Bl_[imode], params, smpi, domain.fake_patch, domain.patch_ );

    smpi->recvComplex( fake_fields->Br_[imode], local_patch_rank, hindex*9+4 );
    fake_fields->Br_[imode]->put( globalfields->Br_[imode], params, smpi, domain.fake_patch, domain.patch_ );

    smpi->recvComplex( fake_fields->Bt_[imode], local_patch_rank, hindex*9+5 );
    fake_fields->Bt_[imode]->put( globalfields->Bt_[imode], params, smpi, domain.fake_patch, domain.patch_ );

    if (!params.is_spectral) {
        smpi->recvComplex( fake_fields->Bl_m[imode], local_patch_rank, hindex*9+6 );
        fake_fields->Bl_m[imode]->put( globalfields->Bl_m[imode], params, smpi, domain.fake_patch, domain.patch_ );
        
        smpi->recvComplex( fake_fields->Br_m[imode], local_patch_rank, hindex*9+7 );
        fake_fields->Br_m[imode]->put( globalfields->Br_m[imode], params, smpi, domain.fake_patch, domain.patch_ );
        
        smpi->recvComplex( fake_fields->Bt_m[imode], local_patch_rank, hindex*9+8 );
        fake_fields->Bt_m[imode]->put( globalfields->Bt_m[imode], params, smpi, domain.fake_patch, domain.patch_ );
    }

}
