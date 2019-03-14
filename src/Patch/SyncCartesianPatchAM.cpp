
#include "SyncCartesianPatchAM.h"

#include <vector>

#include "Domain.h"
#include "VectorPatch.h"
#include "Params.h"
#include "SmileiMPI.h"
#include "PatchesFactory.h"

using namespace std;

void SyncCartesianPatchAM::patchedToCartesian( VectorPatch &vecPatches, Domain &domain, Params &params, SmileiMPI *smpi, Timers &timers, int itime, unsigned int imode )
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
    //smpi->isend( localfields->Jl_[imode], send_to_global_patch_rank, hindex*5  , patch->requests_[0] );
    //smpi->isend( localfields->Jr_[imode], send_to_global_patch_rank, hindex*5+1, patch->requests_[1] );
    //smpi->isend( localfields->Jt_[imode], send_to_global_patch_rank, hindex*5+2, patch->requests_[2] );
    

    if(params.is_spectral) {
        //smpi->isend( localfields->rho_AM_[imode],    send_to_global_patch_rank, hindex*5+3, patch->requests_[3] );
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

    //smpi->recv( fake_fields->Jl_[imode], local_patch_rank, hindex*5 );
    fake_fields->Jl_[imode]->put( globalfields->Jl_[imode], params, smpi, domain.fake_patch, domain.patch_ );

    //smpi->recv( fake_fields->Jr_[imode], local_patch_rank, hindex*5+1 );
    fake_fields->Jr_[imode]->put( globalfields->Jr_[imode], params, smpi, domain.fake_patch, domain.patch_ );

    //smpi->recv( fake_fields->Jt_[imode], local_patch_rank, hindex*5+2 );
    fake_fields->Jt_[imode]->put( globalfields->Jt_[imode], params, smpi, domain.fake_patch, domain.patch_ );

    if(params.is_spectral) {
        //smpi->recv( fake_fields->rho_AM_[imode], local_patch_rank, hindex*5+3 );
        fake_fields->rho_AM_[imode]->put( globalfields->rho_AM_[imode], params, smpi, domain.fake_patch, domain.patch_ );
        //smpi->recv( domain.fake_patch->EMfields->rhoold_, local_patch_rank, hindex*5+4 );
        //domain.fake_patch->EMfields->rhoold_->put( globalfields->rhoold_, params, smpi, domain.fake_patch, domain.patch_ );
    }

    //delete fake_patch;

}


void SyncCartesianPatchAM::cartesianToPatches( Domain &domain, VectorPatch &vecPatches, Params &params, SmileiMPI *smpi, Timers &timers, int itime, unsigned int imode )
{
    timers.grids.restart();

    // Loop / additional_patches_ (regarding cartesian), get data from cartesian
    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {
        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;

        SyncCartesianPatchAM::recvCartesianToPatches( vecPatches(ipatch)->EMfields, domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi,  vecPatches(ipatch), imode );

    }


    // Loop / missing_patches_ (regarding cartesian), send data which do not concern myself
    //     patch().send()
    for ( unsigned int i=0 ; i<domain.missing_patches_.size() ; i++ ) {
        unsigned int ipatch = domain.missing_patches_[i]-vecPatches.refHindex_;

        SyncCartesianPatchAM::sendCartesianToPatches( domain.patch_->EMfields, domain.missing_patches_[i], domain.missing_patches_ranks[i], vecPatches, params, smpi, domain, imode );

    }


    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {
        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;

        SyncCartesianPatchAM::finalize_recvCartesianToPatches( vecPatches(ipatch)->EMfields, domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi, vecPatches(ipatch), imode );

    }


    ElectroMagnAM * domain_fields = static_cast<ElectroMagnAM *>( domain.patch_->EMfields );

    //for ( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
    for ( unsigned int i=0 ; i<domain.local_patches_.size() ; i++ ) {
        unsigned int ipatch = domain.local_patches_[i]-vecPatches.refHindex_;

        ElectroMagnAM * patch_fields = static_cast<ElectroMagnAM *>( vecPatches(ipatch)->EMfields );

        patch_fields->El_[imode]->get( domain_fields->El_[imode], params, smpi, domain.patch_, vecPatches(ipatch) );
        patch_fields->Er_[imode]->get( domain_fields->Er_[imode], params, smpi, domain.patch_, vecPatches(ipatch) );
        patch_fields->Et_[imode]->get( domain_fields->Et_[imode], params, smpi, domain.patch_, vecPatches(ipatch) );
        
        patch_fields->Bl_m[imode]->get( domain_fields->Bl_m[imode], params, smpi, domain.patch_, vecPatches(ipatch) );
        patch_fields->Br_m[imode]->get( domain_fields->Br_m[imode], params, smpi, domain.patch_, vecPatches(ipatch) );
        patch_fields->Bt_m[imode]->get( domain_fields->Bt_m[imode], params, smpi, domain.patch_, vecPatches(ipatch) );

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


void SyncCartesianPatchAM::recvCartesianToPatches( ElectroMagn* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch, unsigned int imode )
{
    smpi->irecv( localfields->Ex_, recv_from_global_patch_rank, hindex*6  , patch->requests_[0] );
    smpi->irecv( localfields->Ey_, recv_from_global_patch_rank, hindex*6+1, patch->requests_[1] );
    smpi->irecv( localfields->Ez_, recv_from_global_patch_rank, hindex*6+2, patch->requests_[2] );
   

    smpi->irecv( localfields->Bx_m, recv_from_global_patch_rank, hindex*6+3, patch->requests_[3] );
    smpi->irecv( localfields->By_m, recv_from_global_patch_rank, hindex*6+4, patch->requests_[4] );
    smpi->irecv( localfields->Bz_m, recv_from_global_patch_rank, hindex*6+5, patch->requests_[5] );

}

void SyncCartesianPatchAM::finalize_recvCartesianToPatches( ElectroMagn* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch, unsigned int imode )
{
    MPI_Status status;
    MPI_Wait( &(patch->requests_[0]), &status );
    MPI_Wait( &(patch->requests_[1]), &status );
    MPI_Wait( &(patch->requests_[2]), &status );
    MPI_Wait( &(patch->requests_[3]), &status );
    MPI_Wait( &(patch->requests_[4]), &status );
    MPI_Wait( &(patch->requests_[5]), &status );
}

void SyncCartesianPatchAM::sendCartesianToPatches( ElectroMagn* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Domain& domain, unsigned int imode )
{
    // Jx_
    // define fake_patch
    unsigned int n_moved = 0;
    //Patch* fake_patch = PatchesFactory::clone(vecPatches(0), params, smpi, vecPatches.domain_decomposition_, hindex, n_moved, false);
    //smpi->recv( fake_patch->Jx_, hindex, local_patch_rank );
    //    recv(  EM->Bz_m, from, tag ); tag++;

    //smpi->recv( fake_patch->EMfields->Jx_, local_patch_rank, hindex );

    //  vecPatches(ipatch) -> need sender patch coordinates : vecPatches.getDomainCoordinates( hindex )
    // Buffer will be resized for each component, fake local patch, which wil have Jxyz, and coordinates to update ?


    domain.fake_patch->hindex = hindex;
    domain.fake_patch->Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );

    domain.fake_patch->EMfields->Ex_->get( globalfields->Ex_, params, smpi, domain.patch_, domain.fake_patch );
    smpi->send( domain.fake_patch->EMfields->Ex_, local_patch_rank, hindex*6 );

    domain.fake_patch->EMfields->Ey_->get( globalfields->Ey_, params, smpi, domain.patch_, domain.fake_patch );
    smpi->send( domain.fake_patch->EMfields->Ey_, local_patch_rank, hindex*6+1 );

    domain.fake_patch->EMfields->Ez_->get( globalfields->Ez_, params, smpi, domain.patch_, domain.fake_patch );
    smpi->send( domain.fake_patch->EMfields->Ez_, local_patch_rank, hindex*6+2 );

    domain.fake_patch->EMfields->Bx_m->get( globalfields->Bx_m, params, smpi, domain.patch_, domain.fake_patch );
    smpi->send( domain.fake_patch->EMfields->Bx_m, local_patch_rank, hindex*6+3 );

    domain.fake_patch->EMfields->By_m->get( globalfields->By_m, params, smpi, domain.patch_, domain.fake_patch );
    smpi->send( domain.fake_patch->EMfields->By_m, local_patch_rank, hindex*6+4 );

    domain.fake_patch->EMfields->Bz_m->get( globalfields->Bz_m, params, smpi, domain.patch_, domain.fake_patch );
    smpi->send( domain.fake_patch->EMfields->Bz_m, local_patch_rank, hindex*6+5 );


    //if(params.is_spectral){
    //    smpi->send( localfields->rho_, hindex, global_patch_rank, params, smpi );
    //
    //}

    //delete fake_patch;
}




// ---------------------------------
// --------------- MW --------------
// ---------------------------------


void SyncCartesianPatchAM::patchedToCartesian_MW( VectorPatch& vecPatches, Domain& domain, Params &params, SmileiMPI* smpi )
{

    // Loop / local_patches_ -> OK
    //     put()

    // Loop / additional_patches_ (identify where goes the additionnal patches regarding patch_count)
    //     patch().send()
    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {
        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;
        SyncCartesianPatchAM::sendPatchedToCartesian_MW( vecPatches(ipatch)->EMfields, domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi, vecPatches(ipatch), params );
    }


    // Loop / missing_patches_ (identify where are the missing patches)
    //     domain.recv( domain.missing_patches_[i] from smpi->hrank( domain.missing_patches_[i] ) )
    //     put to domain
    for ( unsigned int i=0 ; i<domain.missing_patches_.size() ; i++ ) {
        unsigned int ipatch = domain.missing_patches_[i]-vecPatches.refHindex_;
        SyncCartesianPatchAM::recvPatchedToCartesian_MW( domain.patch_->EMfields, domain.missing_patches_[i], domain.missing_patches_ranks[i], vecPatches, params, smpi, domain );

    }


    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {
        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;
        SyncCartesianPatchAM::finalize_sendPatchedToCartesian_MW( vecPatches(ipatch)->EMfields, domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi, vecPatches(ipatch), params );
    }

    for ( unsigned int i=0 ; i<domain.local_patches_.size() ; i++ ) {
        unsigned int ipatch = domain.local_patches_[i]-vecPatches.refHindex_;

        vecPatches(ipatch)->EMfields->Ex_->put( domain.patch_->EMfields->Ex_, params, smpi, vecPatches(ipatch), domain.patch_ );
        vecPatches(ipatch)->EMfields->Ey_->put( domain.patch_->EMfields->Ey_, params, smpi, vecPatches(ipatch), domain.patch_ );
        vecPatches(ipatch)->EMfields->Ez_->put( domain.patch_->EMfields->Ez_, params, smpi, vecPatches(ipatch), domain.patch_ );

        vecPatches(ipatch)->EMfields->Bx_->put( domain.patch_->EMfields->Bx_, params, smpi, vecPatches(ipatch), domain.patch_ );
        vecPatches(ipatch)->EMfields->By_->put( domain.patch_->EMfields->By_, params, smpi, vecPatches(ipatch), domain.patch_ );
        vecPatches(ipatch)->EMfields->Bz_->put( domain.patch_->EMfields->Bz_, params, smpi, vecPatches(ipatch), domain.patch_ );

        vecPatches(ipatch)->EMfields->Bx_m->put( domain.patch_->EMfields->Bx_m, params, smpi, vecPatches(ipatch), domain.patch_ );
        vecPatches(ipatch)->EMfields->By_m->put( domain.patch_->EMfields->By_m, params, smpi, vecPatches(ipatch), domain.patch_ );
        vecPatches(ipatch)->EMfields->Bz_m->put( domain.patch_->EMfields->Bz_m, params, smpi, vecPatches(ipatch), domain.patch_ );


    }
}

void SyncCartesianPatchAM::sendPatchedToCartesian_MW( ElectroMagn* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params )
{
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

void SyncCartesianPatchAM::finalize_sendPatchedToCartesian_MW( ElectroMagn* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params )
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

void SyncCartesianPatchAM::recvPatchedToCartesian_MW( ElectroMagn* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Domain& domain )
{
    // Jx_
    unsigned int n_moved = 0;

    //  vecPatches(ipatch) -> need sender patch coordinates : vecPatches.getDomainCoordinates( hindex )
    // Buffer will be resized for each component, fake local patch, which wil have Jxyz, and coordinates to update ?
    domain.fake_patch->hindex = hindex;
    domain.fake_patch->Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );

    smpi->recv( domain.fake_patch->EMfields->Ex_, local_patch_rank, hindex*9 );
    domain.fake_patch->EMfields->Ex_->put( globalfields->Ex_, params, smpi, domain.fake_patch, domain.patch_ );

    smpi->recv( domain.fake_patch->EMfields->Ey_, local_patch_rank, hindex*9+1 );
    domain.fake_patch->EMfields->Ey_->put( globalfields->Ey_, params, smpi, domain.fake_patch, domain.patch_ );

    smpi->recv( domain.fake_patch->EMfields->Ez_, local_patch_rank, hindex*9+2 );
    domain.fake_patch->EMfields->Ez_->put( globalfields->Ez_, params, smpi, domain.fake_patch, domain.patch_ );

    smpi->recv( domain.fake_patch->EMfields->Bx_, local_patch_rank, hindex*9+3 );
    domain.fake_patch->EMfields->Bx_->put( globalfields->Bx_, params, smpi, domain.fake_patch, domain.patch_ );

    smpi->recv( domain.fake_patch->EMfields->By_, local_patch_rank, hindex*9+4 );
    domain.fake_patch->EMfields->By_->put( globalfields->By_, params, smpi, domain.fake_patch, domain.patch_ );

    smpi->recv( domain.fake_patch->EMfields->Bz_, local_patch_rank, hindex*9+5 );
    domain.fake_patch->EMfields->Bz_->put( globalfields->Bz_, params, smpi, domain.fake_patch, domain.patch_ );

    smpi->recv( domain.fake_patch->EMfields->Bx_m, local_patch_rank, hindex*9+6 );
    domain.fake_patch->EMfields->Bx_m->put( globalfields->Bx_m, params, smpi, domain.fake_patch, domain.patch_ );

    smpi->recv( domain.fake_patch->EMfields->By_m, local_patch_rank, hindex*9+7 );
    domain.fake_patch->EMfields->By_m->put( globalfields->By_m, params, smpi, domain.fake_patch, domain.patch_ );

    smpi->recv( domain.fake_patch->EMfields->Bz_m, local_patch_rank, hindex*9+8 );
    domain.fake_patch->EMfields->Bz_m->put( globalfields->Bz_m, params, smpi, domain.fake_patch, domain.patch_ );

}
