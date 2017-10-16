
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
        vecPatches(ipatch)->EMfields->Ex_->put( domain.patch_->EMfields->Ex_, params, smpi, vecPatches(ipatch), domain.patch_ );
        vecPatches(ipatch)->EMfields->Ey_->put( domain.patch_->EMfields->Ey_, params, smpi, vecPatches(ipatch), domain.patch_ );
        vecPatches(ipatch)->EMfields->Ez_->put( domain.patch_->EMfields->Ez_, params, smpi, vecPatches(ipatch), domain.patch_ );

        vecPatches(ipatch)->EMfields->Bx_->put( domain.patch_->EMfields->Bx_, params, smpi, vecPatches(ipatch), domain.patch_ );
        vecPatches(ipatch)->EMfields->By_->put( domain.patch_->EMfields->By_, params, smpi, vecPatches(ipatch), domain.patch_ );
        vecPatches(ipatch)->EMfields->Bz_->put( domain.patch_->EMfields->Bz_, params, smpi, vecPatches(ipatch), domain.patch_ );

        vecPatches(ipatch)->EMfields->Jx_->put( domain.patch_->EMfields->Jx_, params, smpi, vecPatches(ipatch), domain.patch_ );
        vecPatches(ipatch)->EMfields->Jy_->put( domain.patch_->EMfields->Jy_, params, smpi, vecPatches(ipatch), domain.patch_ );
        vecPatches(ipatch)->EMfields->Jz_->put( domain.patch_->EMfields->Jz_, params, smpi, vecPatches(ipatch), domain.patch_ );

        vecPatches(ipatch)->EMfields->rho_->put( domain.patch_->EMfields->rho_, params, smpi, vecPatches(ipatch), domain.patch_ );
        vecPatches(ipatch)->EMfields->rhoold_->put( domain.patch_->EMfields->rhoold_, params, smpi, vecPatches(ipatch), domain.patch_ );

    }

}


void SyncCartesianPatch::cartesianToPatches( Domain& domain, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Timers &timers, int itime )
{
    for ( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {

        vecPatches(ipatch)->EMfields->Ex_->get( domain.patch_->EMfields->Ex_, params, smpi, domain.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->Ey_->get( domain.patch_->EMfields->Ey_, params, smpi, domain.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->Ez_->get( domain.patch_->EMfields->Ez_, params, smpi, domain.patch_, vecPatches(ipatch) );

        vecPatches(ipatch)->EMfields->Bx_->get( domain.patch_->EMfields->Bx_, params, smpi, domain.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->By_->get( domain.patch_->EMfields->By_, params, smpi, domain.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->Bz_->get( domain.patch_->EMfields->Bz_, params, smpi, domain.patch_, vecPatches(ipatch) );

        vecPatches(ipatch)->EMfields->Jx_->get( domain.patch_->EMfields->Jx_, params, smpi, domain.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->Jy_->get( domain.patch_->EMfields->Jy_, params, smpi, domain.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->Jz_->get( domain.patch_->EMfields->Jz_, params, smpi, domain.patch_, vecPatches(ipatch) );

        vecPatches(ipatch)->EMfields->rho_->get( domain.patch_->EMfields->rho_, params, smpi, domain.patch_, vecPatches(ipatch) );
        vecPatches(ipatch)->EMfields->rhoold_->get( domain.patch_->EMfields->rhoold_, params, smpi, domain.patch_, vecPatches(ipatch) );


    }

}

