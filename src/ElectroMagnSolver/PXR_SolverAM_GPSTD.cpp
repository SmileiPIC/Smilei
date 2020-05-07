
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

void PXR_SolverAM_GPSTD::coupling( Params &params, ElectroMagn *EMfields, bool full_domain )
{
#ifdef _PICSAR    

    int nl, nr;
    int ovl, ovr;
    // unable to convert unsigned int to an iso_c_binding supported type
    
    std::vector<unsigned int> n_space(params.n_space);
    if (full_domain)
        n_space = params.n_space_global;
    else if (params.uncoupled_grids)
        n_space = params.n_space_region;
    
    nl=( int ) (0 + n_space[0]);
    nr=( int ) (0 + n_space[1]);
    
    if (params.uncoupled_grids) {
        ovl=( int ) params.region_oversize[0];
        ovr=( int ) params.region_oversize[1];
    }
    else {
        ovl=( int ) params.oversize[0];
        ovr=( int ) params.oversize[1];
    }
    ovr=0;

    std::vector<unsigned int> dimPrim( 3 );
    // Dimension of the primal and dual grids
    for( size_t i=0 ; i<params.nDim_field ; i++ ) {
        // Standard scheme
        dimPrim[i+1] = n_space[i]+1;
        if (params.uncoupled_grids)
            dimPrim[i+1] += 2*params.region_oversize[i];
        else
            dimPrim[i+1] += 2*params.oversize[i];
    }
    dimPrim[0] = params.nmodes;
    dimPrim[2] = n_space[1]+1;

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
    ERROR( "Smilei not linked with picsar, use make config=picsar" );
#endif
                                
                                
}

void PXR_SolverAM_GPSTD::uncoupling()
{
#ifdef _PICSAR
    picsar::free_params_picsar_AM();
#else
    ERROR( "Smilei not linked with picsar, use make config=picsar" );
#endif
}

//rotational_cleaning is called over a reconstruction of the full domain in a single region in order to correct the initial laser field
void PXR_SolverAM_GPSTD::rotational_cleaning( ElectroMagn *fields )
{
    _2Dvectors_to_3D(fields);
    
#ifdef _PICSAR
    picsar::rotational_cleaning();
#else
    ERROR( "Smilei not linked with picsar, use make config=picsar" );
#endif
    
    //duplicate_field_into_smilei( fields );
    _3D_to_2Dvectors(fields);

}

//densities_correction takes care of densities high frequency filtering and divergence cleaning before Maxwell solver is called.
void PXR_SolverAM_GPSTD::densities_correction(ElectroMagn *fields)
{
//0) Transform toward spectral space
//1) Filter + divergence cleaning (current correction)
//2) Back to intermediate space
#ifdef _PICSAR
    picsar::densities_correction();
#else
    ERROR( "Smilei not linked with picsar, use make config=picsar" );
#endif
//3) Communicate J, rho, rho_old and set them to zero in boundary cells



}

void PXR_SolverAM_GPSTD::operator()( ElectroMagn *fields )
{
    //duplicate_field_into_pxr( fields );
    _2Dvectors_to_3D(fields);
    
#ifdef _PICSAR
    picsar::push_psatd_ebfield_();
#else
    ERROR( "Smilei not linked with picsar, use make config=picsar" );
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

    for ( unsigned int imode=0 ; imode<( static_cast<ElectroMagnAM *>( fields ) )->El_.size() ; imode++ ) {

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

        for (unsigned int il=0;il<El->dims_[0];il++) {
            for (unsigned int ir=fields->oversize[1];ir<El->dims_[1]-fields->oversize[1];ir++) {
                (*El_pxr)(imode,il,ir-fields->oversize[1]) = (*El)(il,ir); 
                (*Er_pxr)(imode,il,ir-fields->oversize[1]) = (*Er)(il,ir); 
                (*Et_pxr)(imode,il,ir-fields->oversize[1]) = (*Et)(il,ir); 
                (*Bl_pxr)(imode,il,ir-fields->oversize[1]) = (*Bl)(il,ir); 
                (*Br_pxr)(imode,il,ir-fields->oversize[1]) = (*Br)(il,ir); 
                (*Bt_pxr)(imode,il,ir-fields->oversize[1]) = (*Bt)(il,ir); 
                (*Jl_pxr)(imode,il,ir-fields->oversize[1]) = (*Jl)(il,ir); 
                (*Jr_pxr)(imode,il,ir-fields->oversize[1]) = (*Jr)(il,ir); 
                (*Jt_pxr)(imode,il,ir-fields->oversize[1]) = (*Jt)(il,ir); 
                (*rho_pxr)(imode,il,ir-fields->oversize[1]) = (*rho)(il,ir); 
                (*rhoold_pxr)(imode,il,ir-fields->oversize[1]) = (*rhoold)(il,ir); 
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

    for ( unsigned int imode=0 ; imode<( static_cast<ElectroMagnAM *>( fields ) )->El_.size() ; imode++ ) {

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

        for (unsigned int il=0;il<El->dims_[0];il++) {
            for (unsigned int ir=fields->oversize[1];ir<El->dims_[1]-fields->oversize[1];ir++) {
                (*El)(il,ir) = (*El_pxr)(imode,il,ir-fields->oversize[1]); 
                (*Er)(il,ir) = (*Er_pxr)(imode,il,ir-fields->oversize[1]); 
                (*Et)(il,ir) = (*Et_pxr)(imode,il,ir-fields->oversize[1]); 
                (*Bl)(il,ir) = (*Bl_pxr)(imode,il,ir-fields->oversize[1]); 
                (*Br)(il,ir) = (*Br_pxr)(imode,il,ir-fields->oversize[1]); 
                (*Bt)(il,ir) = (*Bt_pxr)(imode,il,ir-fields->oversize[1]); 
                (*Jl)(il,ir) = (*Jl_pxr)(imode,il,ir-fields->oversize[1]); 
                (*Jr)(il,ir) = (*Jr_pxr)(imode,il,ir-fields->oversize[1]); 
                (*Jt)(il,ir) = (*Jt_pxr)(imode,il,ir-fields->oversize[1]); 
                (*rho)(il,ir) = (*rho_pxr)(imode,il,ir-fields->oversize[1]); 
                (*rhoold)(il,ir) = (*rhoold_pxr)(imode,il,ir-fields->oversize[1]); 
            }
        }
    }
}

