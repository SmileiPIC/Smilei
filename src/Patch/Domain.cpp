
#include "Domain.h"

#include "PatchesFactory.h"
#include "DomainDecompositionFactory.h"
#include "DiagnosticFields1D.h"
#include "DiagnosticCartFields2D.h"
#include "DiagnosticCartFields3D.h"

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
    vecPatch_.nAntennas = vecPatch_(0)->EMfields->antennas.size();
    vecPatch_.initExternals( params );

/*    if ( params.nDim_field == 1 )
        diag_ = new DiagnosticFields1D( params, smpi, vecPatch_, 0, openPMD ); 
    else if ( params.nDim_field == 2 )
        diag_ = new DiagnosticCartFields2D( params, smpi, vecPatch_, 0, openPMD ); 
    else if ( params.nDim_field == 3 )
        diag_ = new DiagnosticCartFields3D( params, smpi, vecPatch_, 0, openPMD );

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

*/
    if(params.is_spectral == true){
        init_pxr(params);
    }

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
    if(params.is_spectral == false) {
    vecPatch_.solveMaxwell( params, simWindow, itime, time_dual, timers );
    }
    else{
    vecPatch_.solveMaxwell_Spectral(params,simWindow,itime,time_dual,timers);
    }

}

void Domain::init_pxr(Params &params)
{
int n0,n1,n2;
int ov0,ov1,ov2;
// unable to convert unsigned int to an iso_c_binding supported type 

n0=(int) ( 1 + params.n_space[2]*params.global_factor[2]);
n1=(int) ( 1 + params.n_space[1]*params.global_factor[1]);
n2=(int) ( 1 + params.n_space[0]*params.global_factor[0]);

ov0=(int) params.oversize[2];
ov1=(int) params.oversize[1];
ov2=(int) params.oversize[0];

Field3D* Ex3D_pxr = static_cast<Field3D*>(this->vecPatch_(0)->EMfields->Ex_pxr);
Field3D* Ey3D_pxr = static_cast<Field3D*>(this->vecPatch_(0)->EMfields->Ey_pxr);
Field3D* Ez3D_pxr = static_cast<Field3D*>(this->vecPatch_(0)->EMfields->Ez_pxr);
Field3D* Bx3D_pxr = static_cast<Field3D*>(this->vecPatch_(0)->EMfields->Bx_pxr);
Field3D* By3D_pxr = static_cast<Field3D*>(this->vecPatch_(0)->EMfields->By_pxr);
Field3D* Bz3D_pxr = static_cast<Field3D*>(this->vecPatch_(0)->EMfields->Bz_pxr);
Field3D* Jx3D_pxr = static_cast<Field3D*>(this->vecPatch_(0)->EMfields->Jx_pxr);
Field3D* Jy3D_pxr = static_cast<Field3D*>(this->vecPatch_(0)->EMfields->Jy_pxr);
Field3D* Jz3D_pxr = static_cast<Field3D*>(this->vecPatch_(0)->EMfields->Jz_pxr);
Field3D* rho3D_pxr = static_cast<Field3D*>(this->vecPatch_(0)->EMfields->rho_pxr);
Field3D* rhoold3D_pxr = static_cast<Field3D*>(this->vecPatch_(0)->EMfields->rhoold_pxr);

//call of extern init routine (defined in picsar)

  init_params_picsar(&n0,&n1,&n2,
    &params.cell_length[2],&params.cell_length[1],&params.cell_length[0],&params.timestep,
    &ov0,&ov1,&ov2,
    &params.norderz,&params.nordery,&params.norderx,
    &params.is_spectral,
    &(Ex3D_pxr->data_[0]),
    &(Ey3D_pxr->data_[0]),
    &(Ez3D_pxr->data_[0]),
    &(Bx3D_pxr->data_[0]),
    &(By3D_pxr->data_[0]),
    &(Bz3D_pxr->data_[0]),
    &(Jx3D_pxr->data_[0]),
    &(Jy3D_pxr->data_[0]),
    &(Jz3D_pxr->data_[0]),
    &(rho3D_pxr->data_[0]),
    &(rhoold3D_pxr->data_[0]));


}

