
#include "DoubleGrids.h"

#include <vector>

#include "Domain.h"
#include "VectorPatch.h"
#include "Params.h"
#include "SmileiMPI.h"
#include "PatchesFactory.h"

using namespace std;

// ------------------------------------------------------------
// Gather Currents on Domain to apply Maxwell solvers on Domain
// ------------------------------------------------------------
void DoubleGrids::syncCurrentsOnDomain( VectorPatch &vecPatches, Domain &domain, Params &params, SmileiMPI *smpi, Timers &timers, int itime )
{
    timers.grids.restart();

    // Loop / additional_patches_ ( = patches included in the local vecPatches but not in local Domain, the Domain of others MPI will need this data )
    //        additional_patches_ranks stores the MPI rank of the Domain which owns additional_patches_
    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGrids::currentsOnDomainSend( vecPatches(ipatch)->EMfields,
                                           domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi, vecPatches(ipatch), params );

    }


    // Loop / missing_patches_ ( = patches whose data is needed by the local Domain but not own by the local vecPatches )
    //        missing_patches_ranks stores the MPI rank of the vecPacthes which own missing_patches_
    for ( unsigned int i=0 ; i<domain.missing_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.missing_patches_[i]-vecPatches.refHindex_;
        DoubleGrids::currentsOnDomainRecv( domain.patch_->EMfields,
                                           domain.missing_patches_[i], domain.missing_patches_ranks[i], vecPatches, params, smpi, domain );

    }


    // Loop / additional_patches_ to finalize send ( currentsOnDomainSend relies on MPI_Isend )
    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGrids::currentsOnDomainSendFinalize( vecPatches(ipatch)->EMfields,
                                                   domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi, vecPatches(ipatch), params );

    }

    // Loop / local_patches_ ( patches own by the local vePatches whose data are used by the local Domain )
    for ( unsigned int i=0 ; i<domain.local_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.local_patches_[i]-vecPatches.refHindex_;
        vecPatches(ipatch)->EMfields->Jx_->put( domain.patch_->EMfields->Jx_, params, smpi, vecPatches(ipatch), domain.patch_ );
        vecPatches(ipatch)->EMfields->Jy_->put( domain.patch_->EMfields->Jy_, params, smpi, vecPatches(ipatch), domain.patch_ );
        vecPatches(ipatch)->EMfields->Jz_->put( domain.patch_->EMfields->Jz_, params, smpi, vecPatches(ipatch), domain.patch_ );

	if(params.is_spectral){
            vecPatches(ipatch)->EMfields->rho_->put( domain.patch_->EMfields->rho_, params, smpi, vecPatches(ipatch), domain.patch_ );
            // rho_old is save directly on the Domain after the resolution of the Maxwell solver
	}

    }
    timers.grids.update();
}

void DoubleGrids::currentsOnDomainSend( ElectroMagn* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params )
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

void DoubleGrids::currentsOnDomainSendFinalize( ElectroMagn* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params )
{
    MPI_Status status;
    // Wait for currentsOnDomainSend (isend)
    MPI_Wait( &(patch->requests_[0]), &status );
    MPI_Wait( &(patch->requests_[1]), &status );
    MPI_Wait( &(patch->requests_[2]), &status );

    if(params.is_spectral) {
        MPI_Wait( &(patch->requests_[3]), &status );
    }
}

void DoubleGrids::currentsOnDomainRecv( ElectroMagn* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Domain& domain )
{
    // fake_patch consists in a piece of the local Domain to handle naturally patches communications
    //            need to update its hindex and coordinates to put recv data at the good place in the local Domain (put)
    domain.fake_patch->hindex = hindex;
    domain.fake_patch->Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );

    // recv( Fields, sender_mpi_rank, tag );
    //       tag = *5 ? 5 communications are required per patch : 3 currents + rho + rho_old
    smpi->recv( domain.fake_patch->EMfields->Jx_, local_patch_rank, hindex*5 );
    domain.fake_patch->EMfields->Jx_->put( globalfields->Jx_, params, smpi, domain.fake_patch, domain.patch_ );

    smpi->recv( domain.fake_patch->EMfields->Jy_, local_patch_rank, hindex*5+1 );
    domain.fake_patch->EMfields->Jy_->put( globalfields->Jy_, params, smpi, domain.fake_patch, domain.patch_ );

    smpi->recv( domain.fake_patch->EMfields->Jz_, local_patch_rank, hindex*5+2 );
    domain.fake_patch->EMfields->Jz_->put( globalfields->Jz_, params, smpi, domain.fake_patch, domain.patch_ );

    if(params.is_spectral) {
        smpi->recv( domain.fake_patch->EMfields->rho_, local_patch_rank, hindex*5+3 );
        domain.fake_patch->EMfields->rho_->put( globalfields->rho_, params, smpi, domain.fake_patch, domain.patch_ );
    }

}


// ---------------------------------------------------------------------------
// Scatter Fields on Patches for particles interpolation or divergece cleaning
// ---------------------------------------------------------------------------
void DoubleGrids::syncFieldsOnPatches( Domain &domain, VectorPatch &vecPatches, Params &params, SmileiMPI *smpi, Timers &timers, int itime )
{
    timers.grids.restart();

    // Loop / additional_patches_ ( within local vecPatches but not in local Domain )
    //                            get data from Domain of others MPI
    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGrids::fieldsOnPatchesRecv( vecPatches(ipatch)->EMfields,
                                          domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi,  vecPatches(ipatch) );

    }

    // Loop / missing_patches_ ( within local Domain but not in local vecPatches,  )
    //                         send data which do not concern local Domain
    for ( unsigned int i=0 ; i<domain.missing_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.missing_patches_[i]-vecPatches.refHindex_;
        DoubleGrids::fieldsOnPatchesSend( domain.patch_->EMfields,
                                          domain.missing_patches_[i], domain.missing_patches_ranks[i], vecPatches, params, smpi, domain );

    }

    // Loop / additional_patches_ to finalize recv ( fieldsOnPatchesRecv relies on MPI_Irecv )
    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGrids::fieldsOnPatchesRecvFinalize( vecPatches(ipatch)->EMfields,
                                                  domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi, vecPatches(ipatch) );

    }

    // Loop / local_patches_ ( patches own by the local vePatches whose data are used by the local Domain )
    for ( unsigned int i=0 ; i<domain.local_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.local_patches_[i]-vecPatches.refHindex_;

        vecPatches(ipatch)->EMfields->Ex_->get( domain.patch_->EMfields->Ex_, params, smpi, domain.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->Ey_->get( domain.patch_->EMfields->Ey_, params, smpi, domain.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->Ez_->get( domain.patch_->EMfields->Ez_, params, smpi, domain.patch_, vecPatches(ipatch) );
   
        vecPatches(ipatch)->EMfields->Bx_m->get( domain.patch_->EMfields->Bx_m, params, smpi, domain.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->By_m->get( domain.patch_->EMfields->By_m, params, smpi, domain.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->Bz_m->get( domain.patch_->EMfields->Bz_m, params, smpi, domain.patch_, vecPatches(ipatch) );

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

void DoubleGrids::fieldsOnPatchesSend( ElectroMagn* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Domain& domain )
{
    // fake_patch consists in a piece of the local Domain to handle naturally patches communications
    //            need to update its hindex and coordinates to extract (get) appropriate data from the local Domain data before send it
    domain.fake_patch->hindex = hindex;
    domain.fake_patch->Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );

    // send( Fields, targeted_mpi_rank, tag );
    //       tag = *6 ? 6 communications could be required per patch
    //       clarify which usage need B, B_m or both
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

}


// -----------------------------------------------------------
// Gather Fields on Domain to compute the divergence cleaning
//                         could be used for the moving window
// -----------------------------------------------------------
void DoubleGrids::syncFieldsOnDomain( VectorPatch& vecPatches, Domain& domain, Params &params, SmileiMPI* smpi )
{
    // Loop / additional_patches_ ( = patches included in the local vecPatches but not in local Domain, the Domain of others MPI will need this data )
    //        additional_patches_ranks stores the MPI rank of the Domain which owns additional_patches_
    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGrids::fieldsOnDomainSend( vecPatches(ipatch)->EMfields, domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi, vecPatches(ipatch), params );
    }

    // Loop / missing_patches_ ( = patches whose data is needed by the local Domain but not own by the local vecPatches )
    //        missing_patches_ranks stores the MPI rank of the vecPacthes which own missing_patches_
    for ( unsigned int i=0 ; i<domain.missing_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.missing_patches_[i]-vecPatches.refHindex_;
        DoubleGrids::fieldsOnDomainRecv( domain.patch_->EMfields,
                                         domain.missing_patches_[i], domain.missing_patches_ranks[i], vecPatches, params, smpi, domain );

    }


    // Loop / additional_patches_ to finalize send ( currentsOnDomainSend relies on MPI_Isend )
    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGrids::fieldsOnDomainSendFinalize( vecPatches(ipatch)->EMfields,
                                                 domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi, vecPatches(ipatch), params );
    }

    // Loop / local_patches_ ( patches own by the local vePatches whose data are used by the local Domain )
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

void DoubleGrids::fieldsOnDomainSend( ElectroMagn* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params )
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

void DoubleGrids::fieldsOnDomainSendFinalize( ElectroMagn* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params )
{
    MPI_Status status;
    // Wait for fieldsOnDomainSend (isend)
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

void DoubleGrids::fieldsOnDomainRecv( ElectroMagn* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Domain& domain )
{
    // fake_patch consists in a piece of the local Domain to handle naturally patches communications
    //            need to update its hindex and coordinates to put recv data at the good place in the local Domain (put)
    domain.fake_patch->hindex = hindex;
    domain.fake_patch->Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );

    // recv( Fields, sender_mpi_rank, tag );
    //       tag = *9 ? 9 communications could be required per patch : B, B_m, both (at least for the moving window)
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
