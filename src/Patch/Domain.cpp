
#include "Domain.h"

#include "PatchesFactory.h"
#include "GeometryFactory.h"
#include "DiagnosticCartFields2D.h"

#include "Params.h"
#include "OpenPMDparams.h"
#include "SimWindow.h"
#include "Timers.h"
#include "SyncCartesianPatch.h"

Domain::Domain( Params &params ) :
    VecPatchCart_( params ) 
{
}

void Domain::build( Params &params, SmileiMPI* smpi, VectorPatch& vecPatches, OpenPMDparams& openPMD )
{
    cartGeom_ = GeometryFactory::createGlobal( params );
    cartPatch_ = PatchesFactory::create( params, smpi, cartGeom_, vecPatches.refHindex_ / vecPatches.size() );
    cartPatch_->set( params, cartGeom_, vecPatches );
    VecPatchCart_.patches_.push_back( cartPatch_ );

    VecPatchCart_.refHindex_ = vecPatches.refHindex_ / vecPatches.size();
    VecPatchCart_.update_field_list();

    //VecPatchCart_.update_field_list(0);
    VecPatchCart_.patches_[0]->finalizeMPIenvironment();
    VecPatchCart_.nrequests = vecPatches(0)->requests_.size();

    diagCart_ = new DiagnosticCartFields2D( params, smpi, VecPatchCart_, 0, openPMD ); 

    for (unsigned int ifield=0 ; ifield<VecPatchCart_(0)->EMfields->Jx_s.size(); ifield++) {
        if( VecPatchCart_(0)->EMfields->Jx_s[ifield]->data_ == NULL ){
            delete VecPatchCart_(0)->EMfields->Jx_s[ifield];
            VecPatchCart_(0)->EMfields->Jx_s[ifield]=NULL;
        }
    }
    for (unsigned int ifield=0 ; ifield<VecPatchCart_(0)->EMfields->Jy_s.size(); ifield++) {
        if( VecPatchCart_(0)->EMfields->Jy_s[ifield]->data_ == NULL ){
            delete VecPatchCart_(0)->EMfields->Jy_s[ifield];
            VecPatchCart_(0)->EMfields->Jy_s[ifield]=NULL;
        }
    }
    for (unsigned int ifield=0 ; ifield<VecPatchCart_(0)->EMfields->Jz_s.size(); ifield++) {
        if( VecPatchCart_(0)->EMfields->Jz_s[ifield]->data_ == NULL ){
            delete VecPatchCart_(0)->EMfields->Jz_s[ifield];
            VecPatchCart_(0)->EMfields->Jz_s[ifield]=NULL;
        }
    }
    for (unsigned int ifield=0 ; ifield<VecPatchCart_(0)->EMfields->rho_s.size(); ifield++) {
        if( VecPatchCart_(0)->EMfields->rho_s[ifield]->data_ == NULL ){
            delete VecPatchCart_(0)->EMfields->rho_s[ifield];
            VecPatchCart_(0)->EMfields->rho_s[ifield]=NULL;
        }
    }

    diagCart_->init( params, smpi, VecPatchCart_ );
    diagCart_->theTimeIsNow = diagCart_->prepare( 0 );
    //if ( diagCart_->theTimeIsNow )
    //    diagCart_->run( smpi, VecPatchCart_, 0, simWindow );


}

Domain::~Domain()
{
}

void Domain::clean()
{
    if (diagCart_ !=NULL) {
        diagCart_->closeFile();
        delete diagCart_;
    }
    if (cartPatch_!=NULL) delete cartPatch_;
    if (cartGeom_ !=NULL) delete cartGeom_;

}

void Domain::solveMaxwell( Params& params, SimWindow* simWindow, int itime, double time_dual, Timers& timers )
{
    //if ( diagCart_!=NULL ) {
    //timers.diagsNEW.restart();
    //SyncVectorPatch::exchangeE( vecPatches );
    //SyncVectorPatch::finalizeexchangeE( vecPatches );
    //SyncVectorPatch::exchangeB( vecPatches );
    //SyncVectorPatch::finalizeexchangeB( vecPatches ); 

    //SyncCartesianPatch::patchedToCartesian( vecPatches, cartPatch_, params, smpi, timers, itime );

    //SyncVectorPatch::exchangeE( VecPatchCart_ );
    //SyncVectorPatch::finalizeexchangeE( VecPatchCart_ );
    //SyncVectorPatch::exchangeB( VecPatchCart_ );
    //SyncVectorPatch::finalizeexchangeB( VecPatchCart_ );

    VecPatchCart_.solveMaxwell( params, simWindow, itime, time_dual, timers );

    //SyncCartesianPatch::cartesianToPatches( cartPatch_, vecPatches, params, smpi, timers, itime );
                    
    //SyncVectorPatch::exchangeE( vecPatches );
    //SyncVectorPatch::finalizeexchangeE( vecPatches );
    //SyncVectorPatch::exchangeB( vecPatches );
    //SyncVectorPatch::finalizeexchangeB( vecPatches ); 

    //diagCart_->theTimeIsNow = diagCart_->prepare( itime );
    //if ( diagCart_->theTimeIsNow ) {
    //    diagCart_->run( smpi, VecPatchCart_, itime, simWindow );
    //}
    //timers.diagsNEW.update();
    //}

}

