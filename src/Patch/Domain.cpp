
#include "Domain.h"

#include "PatchesFactory.h"
#include "DomainDecompositionFactory.h"
#include "DiagnosticFields1D.h"
#include "DiagnosticCartFields2D.h"
#include "DiagnosticCartFields3D.h"
#include "ElectroMagn.h"
#include "Solver.h"

#include "Params.h"
#include "OpenPMDparams.h"
#include "SimWindow.h"
#include "Timers.h"
#include "SyncCartesianPatch.h"

using namespace std;

Domain::Domain( Params &params ) :
    vecPatch_( params ),
    decomposition_(NULL),
    patch_(NULL),
    diag_(NULL)
{
}

void Domain::build( Params &params, SmileiMPI* smpi, VectorPatch& vecPatches, OpenPMDparams& openPMD )
{
    int rk(0);
    MPI_Comm_rank( MPI_COMM_WORLD, &rk );
    vecPatch_.refHindex_ = rk;


    decomposition_ = DomainDecompositionFactory::createGlobal( params );
    //patch_ = PatchesFactory::create( params, smpi, decomposition_, vecPatches.refHindex_ / vecPatches.size() );
    patch_ = PatchesFactory::create( params, smpi, decomposition_, rk );
    patch_->set( params, decomposition_, vecPatches );
    vecPatch_.patches_.push_back( patch_ );

    //vecPatch_.refHindex_ = vecPatches.refHindex_ / vecPatches.size();
    vecPatch_.update_field_list();

    //vecPatch_.update_field_list(0);
    vecPatch_.patches_[0]->finalizeMPIenvironment(params);
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

    if (params.is_pxr)
        vecPatch_(0)->EMfields->MaxwellAmpereSolver_->coupling( params, vecPatch_(0)->EMfields ); 
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
    MPI_Barrier( MPI_COMM_WORLD );
    MESSAGE( "in domain solvemaxwell" );
    vecPatch_.solveMaxwell( params, simWindow, itime, time_dual, timers );

}


void Domain::identify_additional_patches(SmileiMPI* smpi, VectorPatch& vecPatches)
{
    for ( int ipatch = 0 ; ipatch < vecPatches.size() ; ipatch++ ) {
        
        bool patch_is_in( true );
        for ( int iDim = 0 ; iDim < patch_->getDomainLocalMin().size() ; iDim++ ) {
            double center = ( vecPatches(ipatch)->getDomainLocalMin(iDim) + vecPatches(ipatch)->getDomainLocalMax(iDim) ) /2.;
            if ( ( center < patch_->getDomainLocalMin(iDim) ) || ( center > patch_->getDomainLocalMax(iDim) ) )
                patch_is_in = false;
        }
        if (!patch_is_in)
            additional_patches_.push_back( vecPatches(ipatch)->hindex );
    }
    
    //cout << smpi->getRank() << " - additional : ";
    //for (int i = 0 ; i < additional_patches_.size() ; i++)
    //    cout << additional_patches_[i] << " " ;
    //cout << endl;

}


void Domain::identify_missing_patches(SmileiMPI* smpi, VectorPatch& vecPatches, Params& params)
{
    //missing_patches_.push_back()

    // Loop on theroritical patches

    std::vector<int> patch_min_coord(patch_->getDomainLocalMin().size());
    std::vector<int> patch_max_coord(patch_->getDomainLocalMin().size());
    int npatch_domain(1);
    for ( int iDim = 0 ; iDim < patch_->getDomainLocalMin().size() ; iDim++ ) {
        patch_min_coord[iDim] = (int)( patch_->getDomainLocalMin(iDim) / params.cell_length[iDim] / (double)params.n_space[iDim] );
        patch_max_coord[iDim] = (int)( patch_->getDomainLocalMax(iDim) / params.cell_length[iDim] / (double)params.n_space[iDim] ) - 1;
        npatch_domain *= params.n_space_domain[iDim] / params.n_space[iDim];
    }
    //cout << patch_min_coord[0] << " " << patch_min_coord[1] << endl;
    //cout << patch_max_coord[0] << " " << patch_max_coord[1] << endl;

    int patch_min_id = vecPatches.domain_decomposition_->getDomainId( patch_min_coord );
    int patch_max_id = vecPatches.domain_decomposition_->getDomainId( patch_max_coord );
    //cout << patch_min_id << " " << patch_max_id << "\t npatche theoritical : " << npatch_domain << endl; 

    int hmin = min (patch_min_id,patch_max_id); 
    int hmax = max (patch_min_id,patch_max_id); 
    for ( int ix = patch_min_coord[0] ; ix <= patch_max_coord[0] ; ix++ ) {
        for ( int iy = patch_min_coord[1] ; iy <= patch_max_coord[1] ; iy++ ) {
            std::vector<int> coords(patch_->getDomainLocalMin().size());
            coords[0] = ix;
            coords[1] = iy;
            int hindex = vecPatches.domain_decomposition_->getDomainId( coords );
            if (hindex < hmin)
                 hmin = hindex;
            if (hindex > hmax)
                hmax = hindex;

        }
    }

    for ( int idx = hmin ; idx <=hmax ; idx++ ) {
        if ( (idx < vecPatches(0)->hindex) || (idx > vecPatches(vecPatches.size()-1)->hindex) )
            missing_patches_.push_back( idx );
        else
            local_patches_.push_back( idx );
    }

    //cout << "Theoritical : " << hmin << " " << hmax << endl;
    //cout << smpi->getRank() << " - missing : " ;
    //for (int i = 0 ; i < missing_patches_.size() ; i++)
    //    cout << missing_patches_[i] << " " ;
    //cout << endl;

    //cout << smpi->getRank() << " - local : " ;
    //for (int i = 0 ; i < local_patches_.size() ; i++)
    //    cout << local_patches_[i] << " " ;
    //cout << endl;

}
