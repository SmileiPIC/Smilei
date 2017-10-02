
#include "Domain.h"

#include "PatchesFactory.h"
#include "DomainDecompositionFactory.h"
#include "DiagnosticCartFields2D.h"

#include "Params.h"
#include "OpenPMDparams.h"
#include "SimWindow.h"
#include "Timers.h"
#include "SyncCartesianPatch.h"

Domain::Domain( Params &params ) :
    vecPatch_( params ),
    decomposition_(NULL),
    patch_(NULL),
    diag_(NULL)
{
}

void Domain::build( Params &params, SmileiMPI* smpi, VectorPatch& vecPatches, OpenPMDparams& openPMD )
{
    decomposition_ = DomainDecompositionFactory::createGlobal( params );
    patch_ = PatchesFactory::create( params, smpi, decomposition_, vecPatches.refHindex_ / vecPatches.size() );
    patch_->set( params, decomposition_, vecPatches );
    vecPatch_.patches_.push_back( patch_ );

    vecPatch_.refHindex_ = vecPatches.refHindex_ / vecPatches.size();
    vecPatch_.update_field_list();

    //vecPatch_.update_field_list(0);
    vecPatch_.patches_[0]->finalizeMPIenvironment();
    vecPatch_.nrequests = vecPatches(0)->requests_.size();
    vecPatch_.initExternals( params );

    diag_ = new DiagnosticCartFields2D( params, smpi, vecPatch_, 0, openPMD ); 

    for (unsigned int ifield=0 ; ifield<vecPatch_(0)->EMfields->Jx_s.size(); ifield++) {
        if( vecPatch_(0)->EMfields->Jx_s[ifield]->data_ == NULL ){
            delete vecPatch_(0)->EMfields->Jx_s[ifield];
            vecPatch_(0)->EMfields->Jx_s[ifield]=NULL;
        }
    }
    for (unsigned int ifield=0 ; ifield<vecPatch_(0)->EMfields->Jy_s.size(); ifield++) {
        if( vecPatch_(0)->EMfields->Jy_s[ifield]->data_ == NULL ){
            delete vecPatch_(0)->EMfields->Jy_s[ifield];
            vecPatch_(0)->EMfields->Jy_s[ifield]=NULL;
        }
    }
    for (unsigned int ifield=0 ; ifield<vecPatch_(0)->EMfields->Jz_s.size(); ifield++) {
        if( vecPatch_(0)->EMfields->Jz_s[ifield]->data_ == NULL ){
            delete vecPatch_(0)->EMfields->Jz_s[ifield];
            vecPatch_(0)->EMfields->Jz_s[ifield]=NULL;
        }
    }
    for (unsigned int ifield=0 ; ifield<vecPatch_(0)->EMfields->rho_s.size(); ifield++) {
        if( vecPatch_(0)->EMfields->rho_s[ifield]->data_ == NULL ){
            delete vecPatch_(0)->EMfields->rho_s[ifield];
            vecPatch_(0)->EMfields->rho_s[ifield]=NULL;
        }
    }

    diag_->init( params, smpi, vecPatch_ );
    diag_->theTimeIsNow = diag_->prepare( 0 );
    //if ( diag_->theTimeIsNow )
    //    diag_->run( smpi, vecPatch_, 0, simWindow );


}

Domain::~Domain()
{
}

void Domain::clean()
{
    if (diag_ !=NULL) {
        diag_->closeFile();
        delete diag_;
    }
    if (patch_!=NULL) delete patch_;
    if (decomposition_ !=NULL) delete decomposition_;

}

void Domain::solveMaxwell( Params& params, SimWindow* simWindow, int itime, double time_dual, Timers& timers )
{
    //if ( diag_!=NULL ) {
    //timers.diagsNEW.restart();
    //SyncVectorPatch::exchangeE( vecPatches );
    //SyncVectorPatch::finalizeexchangeE( vecPatches );
    //SyncVectorPatch::exchangeB( vecPatches );
    //SyncVectorPatch::finalizeexchangeB( vecPatches ); 

    //SyncCartesianPatch::patchedToCartesian( vecPatches, patch_, params, smpi, timers, itime );

    //SyncVectorPatch::exchangeE( vecPatch_ );
    //SyncVectorPatch::finalizeexchangeE( vecPatch_ );
    //SyncVectorPatch::exchangeB( vecPatch_ );
    //SyncVectorPatch::finalizeexchangeB( vecPatch_ );

    vecPatch_.solveMaxwell( params, simWindow, itime, time_dual, timers );

    //SyncCartesianPatch::cartesianToPatches( patch_, vecPatches, params, smpi, timers, itime );
                    
    //SyncVectorPatch::exchangeE( vecPatches );
    //SyncVectorPatch::finalizeexchangeE( vecPatches );
    //SyncVectorPatch::exchangeB( vecPatches );
    //SyncVectorPatch::finalizeexchangeB( vecPatches ); 

    //diag_->theTimeIsNow = diag_->prepare( itime );
    //if ( diag_->theTimeIsNow ) {
    //    diag_->run( smpi, vecPatch_, itime, simWindow );
    //}
    //timers.diagsNEW.update();
    //}

}

