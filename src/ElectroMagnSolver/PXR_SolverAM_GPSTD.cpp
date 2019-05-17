
#include "PXR_SolverAM_GPSTD.h"

#include "ElectroMagnAM.h"
#include "cField2D.h"
#include "cField3D.h"
#include "interface.h"

#include <iostream>

using namespace std;


PXR_SolverAM_GPSTD::PXR_SolverAM_GPSTD( Params &params )
    : SolverAM( params )
{
    El_pxr = NULL;
    Er_pxr = NULL;
    Et_pxr = NULL;
    Bl_pxr = NULL;
    Br_pxr = NULL;
    Bt_pxr = NULL;
    Jl_pxr = NULL;
    Jr_pxr = NULL;
    Jt_pxr = NULL;
    rho_pxr = NULL;
    rhoold_pxr = NULL;
}

PXR_SolverAM_GPSTD::~PXR_SolverAM_GPSTD()
{
    if (El_pxr) delete El_pxr  ;
    if (Er_pxr) delete Er_pxr  ;
    if (Et_pxr) delete Et_pxr  ;
    if (Bl_pxr) delete Bl_pxr  ;
    if (Br_pxr) delete Br_pxr  ;
    if (Bt_pxr) delete Bt_pxr  ;
    if (Jl_pxr) delete Jl_pxr  ;
    if (Jr_pxr) delete Jr_pxr  ;
    if (Jt_pxr) delete Jt_pxr  ;
    if (rho_pxr) delete rho_pxr;
    if (rhoold_pxr) delete rhoold_pxr;
}

void PXR_SolverAM_GPSTD::coupling( Params &params, ElectroMagn *EMfields )
{

    int nl, nr;
    int ovl, ovr;
    // unable to convert unsigned int to an iso_c_binding supported type
    
    std::vector<unsigned int> n_space(params.n_space);
    if (params.uncoupled_grids)
        n_space = params.n_space_domain;
    
    nl=( int ) (0 + n_space[0]);
    nr=( int ) (0 + n_space[1]);
    
    ovl=( int ) params.oversize[0];
    ovr=( int ) params.oversize[1];
    



    std::vector<unsigned int> dimPrim( 3 );
    // Dimension of the primal and dual grids
    for( size_t i=0 ; i<params.nDim_field ; i++ ) {
        // Standard scheme
        dimPrim[i+1] = n_space[i]+1;
        dimPrim[i+1] += 2*params.oversize[i];
    }
    dimPrim[0] = params.nmodes;

    El_pxr = new cField3D( dimPrim );
    Er_pxr = new cField3D( dimPrim );
    Et_pxr = new cField3D( dimPrim );
    Bl_pxr = new cField3D( dimPrim );
    Br_pxr = new cField3D( dimPrim );
    Bt_pxr = new cField3D( dimPrim );
    Jl_pxr = new cField3D( dimPrim );
    Jr_pxr = new cField3D( dimPrim );
    Jt_pxr = new cField3D( dimPrim );
    rho_pxr = new cField3D( dimPrim );
    rhoold_pxr = new cField3D( dimPrim );

    _2Dvectors_to_3D(EMfields);


    double pxr_dl = params.cell_length[0];
    double pxr_dr = params.cell_length[1];

//    for ( int imode=0 ; imode<params.nmodes ; imode++ ) {
//
//        El = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
//        Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
//        Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
//        Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_[imode];
//        Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_[imode];
//        Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_[imode];
//        Jl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Jl_[imode];
//        Jr = ( static_cast<ElectroMagnAM *>( EMfields ) )->Jr_[imode];
//        Jt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Jt_[imode];
//        rho = ( static_cast<ElectroMagnAM *>( EMfields ) )->rho_AM_[imode];
//        rhoold = ( static_cast<ElectroMagnAM *>( EMfields ) )->rho_old_AM_[imode];
//
//        double sum = 0.;
//        for ( int il=0 ; il<Et->dims_[0] ; il++ )
//            for ( int ir=0 ; ir<Et->dims_[1] ; ir++ )
//                sum += ((*Et)(il,ir)).real()*((*Et)(il,ir)).real() + ((*Et)(il,ir)).imag()*((*Et)(il,ir)).imag() ;
//        cout << "Et " << sum << endl;
//        sum = 0.;
//        for ( int il=0 ; il<Et->dims_[0] ; il++ )
//            for ( int ir=0 ; ir<Et->dims_[1] ; ir++ )
//                sum += ((*Er)(il,ir)).real()*((*Er)(il,ir)).real() + ((*Er)(il,ir)).imag()*((*Er)(il,ir)).imag() ;
//        cout << "Er " << sum << endl;
//        sum = 0.;
//        for ( int il=0 ; il<Et->dims_[0] ; il++ )
//            for ( int ir=0 ; ir<Et->dims_[1] ; ir++ )
//                sum += ((*El)(il,ir)).real()*((*El)(il,ir)).real() + ((*El)(il,ir)).imag()*((*El)(il,ir)).imag() ;
//        cout << "El " << sum << endl;
//        sum = 0.;
//        for ( int il=0 ; il<Et->dims_[0] ; il++ )
//            for ( int ir=0 ; ir<Et->dims_[1] ; ir++ )
//                sum += ((*Bt)(il,ir)).real()*((*Bt)(il,ir)).real() + ((*Bt)(il,ir)).imag()*((*Bt)(il,ir)).imag() ;
//        cout << "Bt " << sum << endl;
//        sum = 0.;
//        for ( int il=0 ; il<Et->dims_[0] ; il++ )
//            for ( int ir=0 ; ir<Et->dims_[1] ; ir++ )
//                sum += ((*Br)(il,ir)).real()*((*Br)(il,ir)).real() + ((*Br)(il,ir)).imag()*((*Br)(il,ir)).imag() ;
//        cout << "Br " << sum << endl;
//        sum = 0.;
//        for ( int il=0 ; il<Et->dims_[0] ; il++ )
//            for ( int ir=0 ; ir<Et->dims_[1] ; ir++ )
//                sum += ((*Bl)(il,ir)).real()*((*Bl)(il,ir)).real() + ((*Bl)(il,ir)).imag()*((*Bl)(il,ir)).imag() ;
//        cout << "Bl " << sum << endl;
//    }
#ifdef _PICSAR    

        int nmodes ( Nmode );
        //call of extern init routine (defined in picsar)
        picsar::init_params_picsar_AM( &nr, &nl, &nmodes, &nmodes, 
                                    &pxr_dr, &pxr_dl, &params.timestep,
                                    &ovr, &ovl,
                                    &params.norder[1], &params.norder[0],
                                    &( Et_pxr->cdata_[0] ),
                                    &( Er_pxr->cdata_[0] ),
                                    &( El_pxr->cdata_[0] ),
                                    &( Bt_pxr->cdata_[0] ),
                                    &( Br_pxr->cdata_[0] ),
                                    &( Bl_pxr->cdata_[0] ),
                                    &( Jt_pxr->cdata_[0] ),
                                    &( Jr_pxr->cdata_[0] ),
                                    &( Jl_pxr->cdata_[0] ),
                                    &( rho_pxr->cdata_[0] ),
                                    &( rhoold_pxr->cdata_[0] ) );
#else
    ERROR( "Smilei not linked with picsar" );
#endif
                                
                                
}

void PXR_SolverAM_GPSTD::divergence_cleaning( ElectroMagn *fields )
{
    _2Dvectors_to_3D(fields);
    
#ifdef _PICSAR
    //picsar::divergence_cleaning();
#else
    ERROR( "Smilei not linked with picsar" );
#endif
    
    //duplicate_field_into_smilei( fields );
    _3D_to_2Dvectors(fields);

}


void PXR_SolverAM_GPSTD::operator()( ElectroMagn *fields )
{
    //duplicate_field_into_pxr( fields );
    _2Dvectors_to_3D(fields);
    
#ifdef _PICSAR
    picsar::push_psatd_ebfield_();
#else
    ERROR( "Smilei not linked with picsar" );
#endif
    
    //duplicate_field_into_smilei( fields );
    _3D_to_2Dvectors(fields);
    
}

void PXR_SolverAM_GPSTD::_2Dvectors_to_3D( ElectroMagn *fields )
{
    cField2D* El;
    cField2D* Er;
    cField2D* Et;
    cField2D* Bl;
    cField2D* Br;
    cField2D* Bt;
    cField2D* Jl;
    cField2D* Jr;
    cField2D* Jt;
    cField2D* rho;
    cField2D* rhoold;

    for ( int imode=0 ; imode<( static_cast<ElectroMagnAM *>( fields ) )->El_.size() ; imode++ ) {

        El = ( static_cast<ElectroMagnAM *>( fields ) )->El_[imode];
        Er = ( static_cast<ElectroMagnAM *>( fields ) )->Er_[imode];
        Et = ( static_cast<ElectroMagnAM *>( fields ) )->Et_[imode];
        Bl = ( static_cast<ElectroMagnAM *>( fields ) )->Bl_[imode];
        Br = ( static_cast<ElectroMagnAM *>( fields ) )->Br_[imode];
        Bt = ( static_cast<ElectroMagnAM *>( fields ) )->Bt_[imode];
        Jl = ( static_cast<ElectroMagnAM *>( fields ) )->Jl_[imode];
        Jr = ( static_cast<ElectroMagnAM *>( fields ) )->Jr_[imode];
        Jt = ( static_cast<ElectroMagnAM *>( fields ) )->Jt_[imode];
        rho = ( static_cast<ElectroMagnAM *>( fields ) )->rho_AM_[imode];
        rhoold = ( static_cast<ElectroMagnAM *>( fields ) )->rho_old_AM_[imode];

        for (int il=0;il<El->dims_[0];il++) {
            for (int ir=0;ir<El->dims_[1];ir++) {
                (*El_pxr)(imode,il,ir) = (*El)(il,ir); 
                (*Er_pxr)(imode,il,ir) = (*Er)(il,ir); 
                (*Et_pxr)(imode,il,ir) = (*Et)(il,ir); 
                (*Bl_pxr)(imode,il,ir) = (*Bl)(il,ir); 
                (*Br_pxr)(imode,il,ir) = (*Br)(il,ir); 
                (*Bt_pxr)(imode,il,ir) = (*Bt)(il,ir); 
                (*Jl_pxr)(imode,il,ir) = (*Jl)(il,ir); 
                (*Jr_pxr)(imode,il,ir) = (*Jr)(il,ir); 
                (*Jt_pxr)(imode,il,ir) = (*Jt)(il,ir); 
                (*rho_pxr)(imode,il,ir) = (*rho)(il,ir); 
                (*rhoold_pxr)(imode,il,ir) = (*rhoold)(il,ir); 
            }
        }
    }
}

void PXR_SolverAM_GPSTD::_3D_to_2Dvectors( ElectroMagn *fields )
{
    cField2D* El;
    cField2D* Er;
    cField2D* Et;
    cField2D* Bl;
    cField2D* Br;
    cField2D* Bt;
    cField2D* Jl;
    cField2D* Jr;
    cField2D* Jt;
    cField2D* rho;
    cField2D* rhoold;

    for ( int imode=0 ; imode<( static_cast<ElectroMagnAM *>( fields ) )->El_.size() ; imode++ ) {

        El = ( static_cast<ElectroMagnAM *>(  fields ) )->El_[imode];
        Er = ( static_cast<ElectroMagnAM *>(  fields ) )->Er_[imode];
        Et = ( static_cast<ElectroMagnAM *>(  fields ) )->Et_[imode];
        Bl = ( static_cast<ElectroMagnAM *>(  fields ) )->Bl_[imode];
        Br = ( static_cast<ElectroMagnAM *>(  fields ) )->Br_[imode];
        Bt = ( static_cast<ElectroMagnAM *>(  fields ) )->Bt_[imode];
        Jl = ( static_cast<ElectroMagnAM *>(  fields ) )->Jl_[imode];
        Jr = ( static_cast<ElectroMagnAM *>(  fields ) )->Jr_[imode];
        Jt = ( static_cast<ElectroMagnAM *>(  fields ) )->Jt_[imode];
        rho = ( static_cast<ElectroMagnAM *>( fields ) )->rho_AM_[imode];
        rhoold = ( static_cast<ElectroMagnAM *>( fields ) )->rho_old_AM_[imode];

        for (int il=0;il<El->dims_[0];il++) {
            for (int ir=0;ir<El->dims_[1];ir++) {
                (*El)(il,ir) = (*El_pxr)(imode,il,ir); 
                (*Er)(il,ir) = (*Er_pxr)(imode,il,ir); 
                (*Et)(il,ir) = (*Et_pxr)(imode,il,ir); 
                (*Bl)(il,ir) = (*Bl_pxr)(imode,il,ir); 
                (*Br)(il,ir) = (*Br_pxr)(imode,il,ir); 
                (*Bt)(il,ir) = (*Bt_pxr)(imode,il,ir); 
                (*Jl)(il,ir) = (*Jl_pxr)(imode,il,ir); 
                (*Jr)(il,ir) = (*Jr_pxr)(imode,il,ir); 
                (*Jt)(il,ir) = (*Jt_pxr)(imode,il,ir); 
                (*rho)(il,ir) = (*rho_pxr)(imode,il,ir); 
                (*rhoold)(il,ir) = (*rhoold_pxr)(imode,il,ir); 
            }
        }
    }
}

