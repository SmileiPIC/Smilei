
#include "DoubleGridsAM.h"

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
void DoubleGridsAM::syncCurrentsOnDomain( VectorPatch &vecPatches, Domain &domain, Params &params, SmileiMPI *smpi, Timers &timers, int itime, unsigned int imode )
{
    timers.grids.restart();

    // Loop / additional_patches_ ( = patches included in the local vecPatches but not in local Domain, the Domain of others MPI will need this data )
    //        additional_patches_ranks stores the MPI rank of the Domain which owns additional_patches_
    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGridsAM::currentsOnDomainSend( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields),
                                             domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi, vecPatches(ipatch), params, imode );

    }

    // Loop / missing_patches_ ( = patches whose data is needed by the local Domain but not own by the local vecPatches )
    //        missing_patches_ranks stores the MPI rank of the vecPacthes which own missing_patches_
    for ( unsigned int i=0 ; i<domain.missing_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.missing_patches_[i]-vecPatches.refHindex_;
        DoubleGridsAM::currentsOnDomainRecv( static_cast<ElectroMagnAM *>(domain.patch_->EMfields)
                                             , domain.missing_patches_[i], domain.missing_patches_ranks[i], vecPatches, params, smpi, domain, imode );

    }

    // Loop / additional_patches_ to finalize send ( currentsOnDomainSend relies on MPI_Isend )
    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGridsAM::currentsOnDomainSendFinalize( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields)
                                                     , domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi, vecPatches(ipatch), params, imode );

    }


    ElectroMagnAM * domain_fields = static_cast<ElectroMagnAM *>( domain.patch_->EMfields );
    // Loop / local_patches_ ( patches own by the local vePatches whose data are used by the local Domain )
    for ( unsigned int i=0 ; i<domain.local_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.local_patches_[i]-vecPatches.refHindex_;        
        ElectroMagnAM * patch_fields = static_cast<ElectroMagnAM *>( vecPatches(ipatch)->EMfields );

        patch_fields->Jl_[imode]->add( domain_fields->Jl_[imode], params, smpi, vecPatches(ipatch), domain.patch_ );
        patch_fields->Jr_[imode]->add( domain_fields->Jr_[imode], params, smpi, vecPatches(ipatch), domain.patch_ );
        patch_fields->Jt_[imode]->add( domain_fields->Jt_[imode], params, smpi, vecPatches(ipatch), domain.patch_ );
        if(params.is_spectral){
            patch_fields->rho_AM_[imode]->add( domain_fields->rho_AM_[imode], params, smpi, vecPatches(ipatch), domain.patch_ );
            // rho_old is save directly on the Domain after the resolution of the Maxwell solver
        }

    }

    timers.grids.update();
}

void DoubleGridsAM::currentsOnDomainSend( ElectroMagnAM* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params, unsigned int imode )
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

void DoubleGridsAM::currentsOnDomainSendFinalize( ElectroMagnAM* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params, unsigned int imode )
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

void DoubleGridsAM::currentsOnDomainRecv( ElectroMagnAM* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Domain& domain, unsigned int imode )
{
    ElectroMagnAM * fake_fields = static_cast<ElectroMagnAM *>( domain.fake_patch->EMfields );

    // fake_patch consists in a piece of the local Domain to handle naturally patches communications
    //            need to update its hindex and coordinates to put recv data at the good place in the local Domain (add)
    domain.fake_patch->hindex = hindex;
    domain.fake_patch->Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );

    // recvComplex( cFields, sender_mpi_rank, tag );
    //              tag = *5 ? 5 communications are required per patch : 3 currents + rho + rho_old
    smpi->recvComplex( fake_fields->Jl_[imode], local_patch_rank, hindex*5 );
    fake_fields->Jl_[imode]->add( globalfields->Jl_[imode], params, smpi, domain.fake_patch, domain.patch_ );

    smpi->recvComplex( fake_fields->Jr_[imode], local_patch_rank, hindex*5+1 );
    fake_fields->Jr_[imode]->add( globalfields->Jr_[imode], params, smpi, domain.fake_patch, domain.patch_ );

    smpi->recvComplex( fake_fields->Jt_[imode], local_patch_rank, hindex*5+2 );
    fake_fields->Jt_[imode]->add( globalfields->Jt_[imode], params, smpi, domain.fake_patch, domain.patch_ );

    if(params.is_spectral) {
        smpi->recvComplex( fake_fields->rho_AM_[imode], local_patch_rank, hindex*5+3 );
        fake_fields->rho_AM_[imode]->add( globalfields->rho_AM_[imode], params, smpi, domain.fake_patch, domain.patch_ );
    }

}

// ---------------------------------------------------------------------------
// Scatter Fields on Patches for particles interpolation or divergece cleaning
// ---------------------------------------------------------------------------
void DoubleGridsAM::syncFieldsOnPatches( Domain &domain, VectorPatch &vecPatches, Params &params, SmileiMPI *smpi, Timers &timers, int itime, unsigned int imode )
{
    timers.grids.restart();

    // Loop / additional_patches_ ( within local vecPatches but not in local Domain )
    //                            get data from Domain of others MPI
    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGridsAM::fieldsOnPatchesRecv( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields),
                                            domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi,  vecPatches(ipatch), params, imode );

    }

    // Loop / missing_patches_ ( within local Domain but not in local vecPatches,  )
    //                         send data which do not concern local Domain
    for ( unsigned int i=0 ; i<domain.missing_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.missing_patches_[i]-vecPatches.refHindex_;
        DoubleGridsAM::fieldsOnPatchesSend( static_cast<ElectroMagnAM *>(domain.patch_->EMfields),
                                            domain.missing_patches_[i], domain.missing_patches_ranks[i], vecPatches, params, smpi, domain, imode );

    }

    // Loop / additional_patches_ to finalize recv ( fieldsOnPatchesRecv relies on MPI_Irecv )
    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGridsAM::fieldsOnPatchesRecvFinalize( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields),
                                                    domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi, vecPatches(ipatch), imode );

    }

    ElectroMagnAM * domain_fields = NULL;
    if ( domain.local_patches_.size() ) // could be empty for domain_global
        domain_fields = static_cast<ElectroMagnAM *>( domain.patch_->EMfields );

    // Loop / local_patches_ ( patches own by the local vePatches whose data are used by the local Domain )
    for ( unsigned int i=0 ; i<domain.local_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.local_patches_[i]-vecPatches.refHindex_;
        ElectroMagnAM * patch_fields = static_cast<ElectroMagnAM *>( vecPatches(ipatch)->EMfields );

        patch_fields->El_[imode]->get( domain_fields->El_[imode], params, smpi, domain.patch_, vecPatches(ipatch) );
        patch_fields->Er_[imode]->get( domain_fields->Er_[imode], params, smpi, domain.patch_, vecPatches(ipatch) );
        patch_fields->Et_[imode]->get( domain_fields->Et_[imode], params, smpi, domain.patch_, vecPatches(ipatch) );
       
        //Temporary synchronize B_m even if it is not necessary since Bm = B 
        patch_fields->Bl_m[imode]->get( domain_fields->Bl_m[imode], params, smpi, domain.patch_, vecPatches(ipatch) );
        patch_fields->Br_m[imode]->get( domain_fields->Br_m[imode], params, smpi, domain.patch_, vecPatches(ipatch) );
        patch_fields->Bt_m[imode]->get( domain_fields->Bt_m[imode], params, smpi, domain.patch_, vecPatches(ipatch) );

        patch_fields->Bl_[imode]->get( domain_fields->Bl_[imode], params, smpi, domain.patch_, vecPatches(ipatch) );
        patch_fields->Br_[imode]->get( domain_fields->Br_[imode], params, smpi, domain.patch_, vecPatches(ipatch) );
        patch_fields->Bt_[imode]->get( domain_fields->Bt_[imode], params, smpi, domain.patch_, vecPatches(ipatch) );

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
   
    smpi->irecvComplex( localfields->Bl_[imode], recv_from_global_patch_rank, hindex*9+3, patch->requests_[3] );
    smpi->irecvComplex( localfields->Br_[imode], recv_from_global_patch_rank, hindex*9+4, patch->requests_[4] );
    smpi->irecvComplex( localfields->Bt_[imode], recv_from_global_patch_rank, hindex*9+5, patch->requests_[5] );

    if (!params.is_spectral) {
        smpi->irecvComplex( localfields->Bl_m[imode], recv_from_global_patch_rank, hindex*9+6, patch->requests_[6] );
        smpi->irecvComplex( localfields->Br_m[imode], recv_from_global_patch_rank, hindex*9+7, patch->requests_[7] );
        smpi->irecvComplex( localfields->Bt_m[imode], recv_from_global_patch_rank, hindex*9+8, patch->requests_[8] );
    }

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

void DoubleGridsAM::fieldsOnPatchesSend( ElectroMagnAM* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Domain& domain, unsigned int imode )
{
    ElectroMagnAM * fake_fields = static_cast<ElectroMagnAM *>( domain.fake_patch->EMfields );

    // fake_patch consists in a piece of the local Domain to handle naturally patches communications
    //            need to update its hindex and coordinates to extract (get) appropriate data from the local Domain data before send it
    domain.fake_patch->hindex = hindex;
    domain.fake_patch->Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );

    // sendComplex( cFields, targeted_mpi_rank, tag );
    //               tag = *9 ? 9 communications could be required per patch
    //               clarify which usage need B, B_m or both
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

    if (!params.is_spectral) {
       fake_fields->Bl_m[imode]->get( globalfields->Bl_m[imode], params, smpi, domain.patch_, domain.fake_patch );
       smpi->sendComplex( fake_fields->Bl_m[imode], local_patch_rank, hindex*9+6 );
       
       fake_fields->Br_m[imode]->get( globalfields->Br_m[imode], params, smpi, domain.patch_, domain.fake_patch );
       smpi->sendComplex( fake_fields->Br_m[imode], local_patch_rank, hindex*9+7 );
       
       fake_fields->Bt_m[imode]->get( globalfields->Bt_m[imode], params, smpi, domain.patch_, domain.fake_patch );
       smpi->sendComplex( fake_fields->Bt_m[imode], local_patch_rank, hindex*9+8 );
    }

}


// -----------------------------------------------------------
// Gather Fields on Domain to compute the divergence cleaning
//                         could be used for the moving window
// -----------------------------------------------------------
void DoubleGridsAM::syncFieldsOnDomain( VectorPatch& vecPatches, Domain& domain, Params &params, SmileiMPI* smpi, unsigned int imode )
{
    // Loop / additional_patches_ ( = patches included in the local vecPatches but not in local Domain, the Domain of others MPI will need this data )
    //        additional_patches_ranks stores the MPI rank of the Domain which owns additional_patches_
    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {
        
        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGridsAM::fieldsOnDomainSend( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields),
                                           domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi, vecPatches(ipatch), params, imode );
    }

    // Loop / missing_patches_ ( = patches whose data is needed by the local Domain but not own by the local vecPatches )
    //        missing_patches_ranks stores the MPI rank of the vecPacthes which own missing_patches_
    for ( unsigned int i=0 ; i<domain.missing_patches_.size() ; i++ ) {
        
        unsigned int ipatch = domain.missing_patches_[i]-vecPatches.refHindex_;
        DoubleGridsAM::fieldsOnDomainRecv( static_cast<ElectroMagnAM *>(domain.patch_->EMfields),
                                           domain.missing_patches_[i], domain.missing_patches_ranks[i], vecPatches, params, smpi, domain, imode );

    }

    // Loop / additional_patches_ to finalize send ( currentsOnDomainSend relies on MPI_Isend )
    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGridsAM::fieldsOnDomainSendFinalize( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields),
                                                   domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi, vecPatches(ipatch), params, imode );
    }

    ElectroMagnAM * domain_fields = NULL;
    if ( domain.local_patches_.size() )
        domain_fields = static_cast<ElectroMagnAM *>( domain.patch_->EMfields );

    // Loop / local_patches_ ( patches own by the local vePatches whose data are used by the local Domain )
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

void DoubleGridsAM::fieldsOnDomainSend( ElectroMagnAM* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params, unsigned int imode )
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

void DoubleGridsAM::fieldsOnDomainSendFinalize( ElectroMagnAM* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params, unsigned int imode )
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

void DoubleGridsAM::fieldsOnDomainRecv( ElectroMagnAM* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Domain& domain, unsigned int imode )
{
    ElectroMagnAM * fake_fields = static_cast<ElectroMagnAM *>( domain.fake_patch->EMfields );

    // fake_patch consists in a piece of the local Domain to handle naturally patches communications
    //            need to update its hindex and coordinates to put recv data at the good place in the local Domain (put)
    domain.fake_patch->hindex = hindex;
    domain.fake_patch->Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );

    // recvComplex( cFields, sender_mpi_rank, tag );
    //              tag = *9 ? 9 communications could be required per patch : B, B_m, both (at least for the moving window)
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



// ---------------------------------------------------------------------------
// Scatter Fields on Patches for particles interpolation or divergece cleaning
// ---------------------------------------------------------------------------
void DoubleGridsAM::syncBOnPatches( Domain &domain, VectorPatch &vecPatches, Params &params, SmileiMPI *smpi, Timers &timers, int itime, unsigned int imode )
{
    timers.grids.restart();

    // Loop / additional_patches_ ( within local vecPatches but not in local Domain )
    //                            get data from Domain of others MPI
    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGridsAM::bOnPatchesRecv( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields),
                                            domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi,  vecPatches(ipatch), params, imode );

    }

    // Loop / missing_patches_ ( within local Domain but not in local vecPatches,  )
    //                         send data which do not concern local Domain
    for ( unsigned int i=0 ; i<domain.missing_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.missing_patches_[i]-vecPatches.refHindex_;
        DoubleGridsAM::bOnPatchesSend( static_cast<ElectroMagnAM *>(domain.patch_->EMfields),
                                            domain.missing_patches_[i], domain.missing_patches_ranks[i], vecPatches, params, smpi, domain, imode );

    }

    // Loop / additional_patches_ to finalize recv ( fieldsOnPatchesRecv relies on MPI_Irecv )
    for ( unsigned int i=0 ; i<domain.additional_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.additional_patches_[i]-vecPatches.refHindex_;
        DoubleGridsAM::bOnPatchesRecvFinalize( static_cast<ElectroMagnAM *>(vecPatches(ipatch)->EMfields),
                                                    domain.additional_patches_[i], domain.additional_patches_ranks[i], smpi, vecPatches(ipatch), imode );

    }

    ElectroMagnAM * domain_fields = NULL;
    if ( domain.local_patches_.size() ) // could be empty for domain_global
        domain_fields = static_cast<ElectroMagnAM *>( domain.patch_->EMfields );

    // Loop / local_patches_ ( patches own by the local vePatches whose data are used by the local Domain )
    for ( unsigned int i=0 ; i<domain.local_patches_.size() ; i++ ) {

        unsigned int ipatch = domain.local_patches_[i]-vecPatches.refHindex_;
        ElectroMagnAM * patch_fields = static_cast<ElectroMagnAM *>( vecPatches(ipatch)->EMfields );

        patch_fields->Bl_[imode]->get( domain_fields->Bl_[imode], params, smpi, domain.patch_, vecPatches(ipatch) );
        patch_fields->Br_[imode]->get( domain_fields->Br_[imode], params, smpi, domain.patch_, vecPatches(ipatch) );
        patch_fields->Bt_[imode]->get( domain_fields->Bt_[imode], params, smpi, domain.patch_, vecPatches(ipatch) );

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

void DoubleGridsAM::bOnPatchesSend( ElectroMagnAM* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Domain& domain, unsigned int imode )
{
    ElectroMagnAM * fake_fields = static_cast<ElectroMagnAM *>( domain.fake_patch->EMfields );

    // fake_patch consists in a piece of the local Domain to handle naturally patches communications
    //            need to update its hindex and coordinates to extract (get) appropriate data from the local Domain data before send it
    domain.fake_patch->hindex = hindex;
    domain.fake_patch->Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );

    // sendComplex( cFields, targeted_mpi_rank, tag );
    //               tag = *9 ? 9 communications could be required per patch
    //               clarify which usage need B, B_m or both
    fake_fields->Bl_[imode]->get( globalfields->Bl_[imode], params, smpi, domain.patch_, domain.fake_patch );
    smpi->sendComplex( fake_fields->Bl_[imode], local_patch_rank, hindex*9+3 );

    fake_fields->Br_[imode]->get( globalfields->Br_[imode], params, smpi, domain.patch_, domain.fake_patch );
    smpi->sendComplex( fake_fields->Br_[imode], local_patch_rank, hindex*9+4 );

    fake_fields->Bt_[imode]->get( globalfields->Bt_[imode], params, smpi, domain.patch_, domain.fake_patch );
    smpi->sendComplex( fake_fields->Bt_[imode], local_patch_rank, hindex*9+5 );

}

