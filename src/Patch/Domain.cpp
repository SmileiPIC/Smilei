
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
    VecPatchCart( params ) 
{
}

void Domain::build( Params &params, SmileiMPI* smpi, VectorPatch& vecPatches, OpenPMDparams& openPMD )
{
    cartGeom = GeometryFactory::createGlobal( params );
    cartPatch = PatchesFactory::create( params, smpi, cartGeom, vecPatches.refHindex_ / vecPatches.size() );
    cartPatch->set( params, cartGeom, vecPatches );
    VecPatchCart.patches_.push_back( cartPatch );

    VecPatchCart.refHindex_ = vecPatches.refHindex_ / vecPatches.size();
    VecPatchCart.update_field_list();

    //VecPatchCart.update_field_list(0);
    VecPatchCart.patches_[0]->finalizeMPIenvironment();
    VecPatchCart.nrequests = vecPatches(0)->requests_.size();

    diagCart = new DiagnosticCartFields2D( params, smpi, VecPatchCart, 0, openPMD ); 

    for (unsigned int ifield=0 ; ifield<VecPatchCart(0)->EMfields->Jx_s.size(); ifield++) {
        if( VecPatchCart(0)->EMfields->Jx_s[ifield]->data_ == NULL ){
            delete VecPatchCart(0)->EMfields->Jx_s[ifield];
            VecPatchCart(0)->EMfields->Jx_s[ifield]=NULL;
        }
    }
    for (unsigned int ifield=0 ; ifield<VecPatchCart(0)->EMfields->Jy_s.size(); ifield++) {
        if( VecPatchCart(0)->EMfields->Jy_s[ifield]->data_ == NULL ){
            delete VecPatchCart(0)->EMfields->Jy_s[ifield];
            VecPatchCart(0)->EMfields->Jy_s[ifield]=NULL;
        }
    }
    for (unsigned int ifield=0 ; ifield<VecPatchCart(0)->EMfields->Jz_s.size(); ifield++) {
        if( VecPatchCart(0)->EMfields->Jz_s[ifield]->data_ == NULL ){
            delete VecPatchCart(0)->EMfields->Jz_s[ifield];
            VecPatchCart(0)->EMfields->Jz_s[ifield]=NULL;
        }
    }
    for (unsigned int ifield=0 ; ifield<VecPatchCart(0)->EMfields->rho_s.size(); ifield++) {
        if( VecPatchCart(0)->EMfields->rho_s[ifield]->data_ == NULL ){
            delete VecPatchCart(0)->EMfields->rho_s[ifield];
            VecPatchCart(0)->EMfields->rho_s[ifield]=NULL;
        }
    }

    diagCart->init( params, smpi, VecPatchCart );
    diagCart->theTimeIsNow = diagCart->prepare( 0 );
    //if ( diagCart->theTimeIsNow )
    //    diagCart->run( smpi, VecPatchCart, 0, simWindow );


}

Domain::~Domain()
{
}

void Domain::clean()
{
    if (diagCart !=NULL) {
        diagCart->closeFile();
        delete diagCart;
    }
    if (cartPatch!=NULL) delete cartPatch;
    if (cartGeom !=NULL) delete cartGeom;

}

void Domain::solveMaxwell( Params& params, SimWindow* simWindow, int itime, double time_dual, Timers& timers )
{
    //if ( diagCart!=NULL ) {
    //timers.diagsNEW.restart();
    //SyncVectorPatch::exchangeE( vecPatches );
    //SyncVectorPatch::finalizeexchangeE( vecPatches );
    //SyncVectorPatch::exchangeB( vecPatches );
    //SyncVectorPatch::finalizeexchangeB( vecPatches ); 

    //SyncCartesianPatch::patchedToCartesian( vecPatches, cartPatch, params, smpi, timers, itime );

    //SyncVectorPatch::exchangeE( VecPatchCart );
    //SyncVectorPatch::finalizeexchangeE( VecPatchCart );
    //SyncVectorPatch::exchangeB( VecPatchCart );
    //SyncVectorPatch::finalizeexchangeB( VecPatchCart );

    VecPatchCart.solveMaxwell( params, simWindow, itime, time_dual, timers );

    //SyncCartesianPatch::cartesianToPatches( cartPatch, vecPatches, params, smpi, timers, itime );
                    
    //SyncVectorPatch::exchangeE( vecPatches );
    //SyncVectorPatch::finalizeexchangeE( vecPatches );
    //SyncVectorPatch::exchangeB( vecPatches );
    //SyncVectorPatch::finalizeexchangeB( vecPatches ); 

    //diagCart->theTimeIsNow = diagCart->prepare( itime );
    //if ( diagCart->theTimeIsNow ) {
    //    diagCart->run( smpi, VecPatchCart, itime, simWindow );
    //}
    //timers.diagsNEW.update();
    //}

}

