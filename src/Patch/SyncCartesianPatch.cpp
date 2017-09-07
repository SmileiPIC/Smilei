
#include "SyncCartesianPatch.h"

#include <vector>

#include "VectorPatch.h"
#include "Params.h"
#include "SmileiMPI.h"

using namespace std;

void SyncCartesianPatch::patchedToCartesian( VectorPatch& vecPatches, Patch* patch, Params &params, SmileiMPI* smpi, Timers &timers, int itime )
{
    //patch->EMfields->rho_->put_to( 0. );
    for ( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
        //SyncCartesianPatch::sync( vecPatches(ipatch)->EMfields->rho_, patch->EMfields->rho_, params, smpi, vecPatches(ipatch), patch );
        //SyncCartesianPatch::sync( vecPatches(ipatch)->EMfields->Ex_, patch->EMfields->Ex_, params, smpi, vecPatches(ipatch), patch );
        //SyncCartesianPatch::sync( vecPatches(ipatch)->EMfields->Ey_, patch->EMfields->Ey_, params, smpi, vecPatches(ipatch), patch );
        //SyncCartesianPatch::sync( vecPatches(ipatch)->EMfields->Jx_, patch->EMfields->Jx_, params, smpi, vecPatches(ipatch), patch );

        SyncCartesianPatch::sync( vecPatches(ipatch)->EMfields->Ex_, patch->EMfields->Ex_, params, smpi, vecPatches(ipatch), patch );
        SyncCartesianPatch::sync( vecPatches(ipatch)->EMfields->Ey_, patch->EMfields->Ey_, params, smpi, vecPatches(ipatch), patch );
        SyncCartesianPatch::sync( vecPatches(ipatch)->EMfields->Ez_, patch->EMfields->Ez_, params, smpi, vecPatches(ipatch), patch );
        SyncCartesianPatch::sync( vecPatches(ipatch)->EMfields->Bx_, patch->EMfields->Bx_, params, smpi, vecPatches(ipatch), patch );
        SyncCartesianPatch::sync( vecPatches(ipatch)->EMfields->By_, patch->EMfields->By_, params, smpi, vecPatches(ipatch), patch );
        SyncCartesianPatch::sync( vecPatches(ipatch)->EMfields->Bz_, patch->EMfields->Bz_, params, smpi, vecPatches(ipatch), patch );
    }

}


void SyncCartesianPatch::cartesianToPatches( Patch* patch, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Timers &timers, int itime )
{
    for ( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
        //SyncCartesianPatch::syncBack( patch->EMfields->rho_, vecPatches(ipatch)->EMfields->rho_, params, smpi, patch, vecPatches(ipatch) );
        //SyncCartesianPatch::syncBack( patch->EMfields->Ex_, vecPatches(ipatch)->EMfields->Ex_, params, smpi, patch, vecPatches(ipatch) );
        //SyncCartesianPatch::syncBack( patch->EMfields->Ey_, vecPatches(ipatch)->EMfields->Ey_, params, smpi, patch, vecPatches(ipatch) );
        //SyncCartesianPatch::syncBack( patch->EMfields->Jx_, vecPatches(ipatch)->EMfields->Jx_, params, smpi, patch, vecPatches(ipatch) );

        SyncCartesianPatch::syncBack( patch->EMfields->Ex_, vecPatches(ipatch)->EMfields->Ex_, params, smpi, patch, vecPatches(ipatch) );
        SyncCartesianPatch::syncBack( patch->EMfields->Ey_, vecPatches(ipatch)->EMfields->Ey_, params, smpi, patch, vecPatches(ipatch) );
        SyncCartesianPatch::syncBack( patch->EMfields->Ez_, vecPatches(ipatch)->EMfields->Ez_, params, smpi, patch, vecPatches(ipatch) );
        SyncCartesianPatch::syncBack( patch->EMfields->Bx_, vecPatches(ipatch)->EMfields->Bx_, params, smpi, patch, vecPatches(ipatch) );
        SyncCartesianPatch::syncBack( patch->EMfields->By_, vecPatches(ipatch)->EMfields->By_, params, smpi, patch, vecPatches(ipatch) );
        SyncCartesianPatch::syncBack( patch->EMfields->Bz_, vecPatches(ipatch)->EMfields->Bz_, params, smpi, patch, vecPatches(ipatch) );
        
    }

}


void SyncCartesianPatch::sync( Field* inField, Field* outField, Params &params, SmileiMPI* smpi, Patch* inPatch, Patch* outPatch )
{
    Field2D* in2D  = static_cast<Field2D*>( inField  );
    Field2D* out2D = static_cast<Field2D*>( outField );

    std::vector<unsigned int> dual =  in2D->isDual_;

    int iout = inPatch->Pcoordinates[0]*params.n_space[0] - outPatch->Pcoordinates[0]*params.n_space[0]*params.global_factor[0] ;
    int jout = inPatch->Pcoordinates[1]*params.n_space[1] - outPatch->Pcoordinates[1]*params.n_space[1]*params.global_factor[1] ;
 
    //for ( unsigned int i = params.oversize[0] ; i < in2D->dims_[0]-params.oversize[0] ; i++ ) {
    //    for ( unsigned int j = params.oversize[1] ; j < in2D->dims_[1]-params.oversize[1] ; j++ ) {
    for ( unsigned int i = 0 ; i < in2D->dims_[0] ; i++ ) {
        for ( unsigned int j = 0 ; j < in2D->dims_[1] ; j++ ) {
            ( *out2D )( iout+i, jout+j ) = ( *in2D )( i, j );
            //( *out2D )( iout+i, jout+j ) += (inPatch->hindex+1);
        }
    }    
    inField->put_to( 0. );
    
}


void SyncCartesianPatch::syncBack( Field* inField, Field* outField, Params &params, SmileiMPI* smpi, Patch* inPatch, Patch* outPatch )
{
    Field2D* in2D  = static_cast<Field2D*>( inField  );
    Field2D* out2D = static_cast<Field2D*>( outField );

    std::vector<unsigned int> dual =  in2D->isDual_;

    int iin = outPatch->Pcoordinates[0]*params.n_space[0] - inPatch->Pcoordinates[0]*params.n_space[0]*params.global_factor[0] ;
    int jin = outPatch->Pcoordinates[1]*params.n_space[1] - inPatch->Pcoordinates[1]*params.n_space[1]*params.global_factor[1] ;
 
    //for ( unsigned int i = params.oversize[0] ; i < out2D->dims_[0]-params.oversize[0] ; i++ ) {
    //    for ( unsigned int j = params.oversize[1] ; j < out2D->dims_[1]-params.oversize[1] ; j++ ) {
    for ( unsigned int i = 0 ; i < out2D->dims_[0] ; i++ ) {
        for ( unsigned int j = 0 ; j < out2D->dims_[1] ; j++ ) {
            ( *out2D )( i, j ) = ( *in2D )( iin+i, jin+j );
            //( *out2D )( i, j ) = in2D->hindex;
        }
    }    

    
}
