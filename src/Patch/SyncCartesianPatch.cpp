
#include "SyncCartesianPatch.h"

#include <vector>

#include "Domain.h"
#include "VectorPatch.h"
#include "Params.h"
#include "SmileiMPI.h"

using namespace std;

void SyncCartesianPatch::patchedToCartesian( VectorPatch& vecPatches, Domain& domain, Params &params, SmileiMPI* smpi, Timers &timers, int itime )
{
    for ( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {

        vecPatches(ipatch)->EMfields->Jx_->push( domain.patch_->EMfields->Jx_, params, smpi, vecPatches(ipatch), domain.patch_ );
        vecPatches(ipatch)->EMfields->Jy_->push( domain.patch_->EMfields->Jy_, params, smpi, vecPatches(ipatch), domain.patch_ );
        vecPatches(ipatch)->EMfields->Jz_->push( domain.patch_->EMfields->Jz_, params, smpi, vecPatches(ipatch), domain.patch_ );
    }

}


void SyncCartesianPatch::cartesianToPatches( Domain& domain, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Timers &timers, int itime )
{
    for ( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {

        vecPatches(ipatch)->EMfields->Ex_->pull( domain.patch_->EMfields->Ex_, params, smpi, domain.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->Ey_->pull( domain.patch_->EMfields->Ey_, params, smpi, domain.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->Ez_->pull( domain.patch_->EMfields->Ez_, params, smpi, domain.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->Bx_m->pull( domain.patch_->EMfields->Bx_m, params, smpi, domain.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->By_m->pull( domain.patch_->EMfields->By_m, params, smpi, domain.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->Bz_m->pull( domain.patch_->EMfields->Bz_m, params, smpi, domain.patch_, vecPatches(ipatch) );

        // Diags only
        vecPatches(ipatch)->EMfields->Bx_->pull( domain.patch_->EMfields->Bx_, params, smpi, domain.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->By_->pull( domain.patch_->EMfields->By_, params, smpi, domain.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->Bz_->pull( domain.patch_->EMfields->Bz_, params, smpi, domain.patch_, vecPatches(ipatch) );
        
    }

}

