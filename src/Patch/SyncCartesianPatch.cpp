
#include "SyncCartesianPatch.h"

#include <vector>

#include "VectorPatch.h"
#include "Params.h"
#include "SmileiMPI.h"

using namespace std;

void SyncCartesianPatch::patchedToCartesian( VectorPatch& vecPatches, Patch* patch, Params &params, SmileiMPI* smpi, Timers &timers, int itime )
{
    for ( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
        SyncCartesianPatch::sync( vecPatches(ipatch)->EMfields->rho_, patch->EMfields->rho_, params, smpi, vecPatches(ipatch), patch );
        SyncCartesianPatch::sync( vecPatches(ipatch)->EMfields->Ex_, patch->EMfields->Ex_, params, smpi, vecPatches(ipatch), patch );
    }

}

void SyncCartesianPatch::sync( Field* inField, Field* outField, Params &params, SmileiMPI* smpi, Patch* inPatch, Patch* outPatch )
{
    Field2D* in2D  = static_cast<Field2D*>( inField  );
    Field2D* out2D = static_cast<Field2D*>( outField );

    std::vector<unsigned int> dual =  in2D->isDual_;

    int iout = inPatch->Pcoordinates[0]*params.n_space[0] - outPatch->Pcoordinates[0]*params.n_space[0]*params.global_factor[0] ;
    int jout = inPatch->Pcoordinates[1]*params.n_space[1] - outPatch->Pcoordinates[1]*params.n_space[1]*params.global_factor[1] ;
 
    for ( unsigned int i = 0 ; i < in2D->dims_[0] ; i++ ) {
        for ( unsigned int j = 0 ; j < in2D->dims_[1] ; j++ ) {
            ( *out2D )( iout+i, jout+j ) = ( *in2D )( i, j );
        }
    }    

    
}

