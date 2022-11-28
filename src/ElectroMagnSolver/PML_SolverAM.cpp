#include "PML_SolverAM.h"
#include "ElectroMagnAM.h"
#include "ElectroMagnBCAM_PML.h"
#include "cField2D.h"
#include "Tools.h"
#include <complex>
#include "dcomplex.h"
#include "Patch.h"

PML_SolverAM::PML_SolverAM( Params &params ):
    SolverAM( params ),
    pml_sigma_( 3, NULL ),
    pml_kappa_( 3, NULL )
{
    std::vector<PyObject *> prof;
    if( PyTools::extract_pyProfiles( "pml_sigma", "Main", 0, prof )){
        if( prof.size() < 2 ){
            ERROR(" in pml_sigma, expecting a list of at least 2 profiles.");
        }
    // extracted profile // number of variables of the function // name of the profile extracted // params // try numpy ?? // try file ?? // time variable ??
        pml_sigma_[0] = new Profile( prof[0], 1, "pml_sigma_l_profile", params, true, false, false );
        if( prof.size() == 2){ 
            pml_sigma_[1] = new Profile( prof[0], 1, "pml_sigma_r_profile", params, true, false, false );
        } else {
            pml_sigma_[1] = new Profile( prof[1], 1, "pml_sigma_r_profile", params, true, false, false );
        }
        // The last profile of the least is taken as the integrate_sigma. Profiles after the 2nd and before the last are neglected.
        pml_sigma_[2] = new Profile( prof[prof.size()-1], 1, "pml_integrate_sigma_profile", params, true, false, false );
    }
    if( PyTools::extract_pyProfiles( "pml_kappa", "Main", 0, prof )){
        if(prof.size() < 2){
            ERROR(" in pml_kappa, expecting a list of at least 2 profiles.");
        }
        pml_kappa_[0] = new Profile( prof[0], 1, "pml_kappa_l_profile", params, true, false, false );
        if( prof.size() == 2){ 
            pml_kappa_[1] = new Profile( prof[0], 1, "pml_kappa_r_profile", params, true, false, false );
        } else{
            pml_kappa_[1] = new Profile( prof[1], 1, "pml_kappa_r_profile", params, true, false, false );
        }
        pml_kappa_[2] = new Profile( prof[prof.size()-1], 1, "pml_integrate_kappa_profile", params, true, false, false );
    }
}

PML_SolverAM::~PML_SolverAM()
{
    for( unsigned int i=0; i<pml_sigma_.size(); i++ ) {
        delete pml_sigma_[i];
    }
    for( unsigned int i=0; i<pml_kappa_.size(); i++ ) {
        delete pml_kappa_[i];
    }
}

void PML_SolverAM::operator()( ElectroMagn * )
{
    ERROR( "This is not a solver for the main domain" );
}

void PML_SolverAM::setDomainSizeAndCoefficients( int iDim, int min_or_max, std::vector<unsigned int> dimPrim, int ncells_pml_domain, int startpml, int* ncells_pml_min, int* ncells_pml_max, Patch* patch )
{
    const unsigned int nl_p = dimPrim[0];
    const unsigned int nl_d = dimPrim[0] + 1;
    const unsigned int nr_p = dimPrim[1];
    const unsigned int nr_d = dimPrim[1] + 1;
    
    // Set global radial index where the PML domain begins j_glob_pml
    //  and set PML sizes nl_p, nr_p ... which are first initialized as patches size.
    if ( iDim == 0 ) {
        j_glob_pml = patch->getCellStartingGlobalIndex( 1 );
        // Size in R is the same as a normal patch => no need to adjust.
    }
    else if ( iDim == 1 ) {
        if (min_or_max==0) {
            j_glob_pml = patch->getCellStartingGlobalIndex( 1 )-ncells_pml_domain + patch->oversize[iDim]+1;
        }
        else if (min_or_max==1) {
            j_glob_pml = patch->getCellStartingGlobalIndex( 1 )+patch->size_[1] + patch->oversize[1];
        }
    }

    //PML Coeffs Kappa,Sigma ...
    //Primal
    kappa_r_p.resize( nr_p );
    sigma_r_p.resize( nr_p );
    kappa_l_p.resize( nl_p );
    sigma_l_p.resize( nl_p );
    integrate_kappa_r_p.resize( nr_p ) ;
    integrate_sigma_r_p.resize( nr_p ) ;
    //Dual
    kappa_r_d.resize( nr_d );
    sigma_r_d.resize( nr_d );
    kappa_l_d.resize( nl_d );
    sigma_l_d.resize( nl_d );
    integrate_kappa_r_d.resize( nr_d ) ;
    integrate_sigma_r_d.resize( nr_d ) ;

    //Coeffs for El,Dl,Hl,Bl
    //Primal
    c1_p_lfield.resize( nr_p, 1. ); // j-dependent
    c2_p_lfield.resize( nr_p, dt ); // j-dependent
    c3_p_lfield.resize( nr_p, 0. ); // j-dependent
    c4_p_lfield.resize( nr_p, 1. ); // j-dependent
    c5_p_lfield.resize( nl_p, 1. ); // i-dependent
    c6_p_lfield.resize( nl_p, 0. ); // i-dependent
    //Dual
    c1_d_lfield.resize( nr_d, 1. ); // j-dependent
    c2_d_lfield.resize( nr_d, dt ); // j-dependent
    c3_d_lfield.resize( nr_d, 0 ); // j-dependent
    c4_d_lfield.resize( nr_d, 1. ); // j-dependent
    c5_d_lfield.resize( nl_d, 1. ); // i-dependent
    c6_d_lfield.resize( nl_d, 0. ); // i-dependent
    //Coefs for Er,Dr,Hr,Br
    //Primal
    c1_p_rfield.resize( nl_p, 1. ); // i-dependent
    c2_p_rfield.resize( nl_p, dt ); // i-dependent
    c3_p_rfield.resize( nr_p, 0. ); // j-dependent
    c4_p_rfield.resize( nr_p, 1. ); // j-dependent
    c5_p_rfield.resize( nr_p, 1. ); // j-dependent
    c6_p_rfield.resize( nr_p, 0. ); // j-dependent
    //Dual
    c1_d_rfield.resize( nl_d, 1. ); // i-dependent
    c2_d_rfield.resize( nl_d, dt ); // i-dependent
    c3_d_rfield.resize( nr_d, 0. ); // j-dependent
    c4_d_rfield.resize( nr_d, 1. ); // j-dependent
    c5_d_rfield.resize( nr_d, 1. ); // j-dependent
    c6_d_rfield.resize( nr_d, 0. ); // j-dependent
    //Coefs for Et,Dt,Ht,Bt
    //Primal
    c1_p_tfield.resize( nr_p, 1. ); // j-dependent
    c2_p_tfield.resize( nr_p, dt ); // j-dependent
    c3_p_tfield.resize( nl_p, 0. ); // i-dependent
    c4_p_tfield.resize( nl_p, 1. ); // i-dependent
    c5_p_tfield.resize( nr_p, 1. ); // j-dependent
    c6_p_tfield.resize( nr_p, 0. ); // j-dependent
    //Dual
    c1_d_tfield.resize( nr_d, 1. ); // j-dependent
    c2_d_tfield.resize( nr_d, dt ); // j-dependent
    c3_d_tfield.resize( nl_d, 0. ); // i-dependent
    c4_d_tfield.resize( nl_d, 1. ); // i-dependent
    c5_d_tfield.resize( nr_d, 1. ); // j-dependent
    c6_d_tfield.resize( nr_d, 0. ); // j-dependent

    double r0 = patch->getDomainLocalMax( 1 ) + patch->oversize[1]*dr  ; // Radius at which pmls start absorbing.
    double rmin = r0 -startpml*dr; // radius of the pml at j=0.

    if ( iDim == 0 ) {
        // 2 cells are vaccum so the PML media begin at r0 which is :
        // Eventually the size of PML media is :
        length_l_pml =  (ncells_pml_domain-startpml+0.5)*dl ;
        // Primal grid
        // Longitudinal
        // Params for first cell of PML-patch (vacuum) i = 0,1,2
        for ( int i=0 ; i<startpml ; i++ ) {
            // Coeffs for the first cell
            kappa_l_p[i] = 1. ;
            sigma_l_p[i] = 0. ;
        }
        // Params for other cells (PML Media) when i>=3
        for( int i = startpml; i < (int) nl_p ; i++ ) {
            kappa_l_p[i] = pml_kappa_[0]->valueAt((i-startpml)*dl/length_l_pml);
            sigma_l_p[i] = pml_sigma_[0]->valueAt((i-startpml)*dl/length_l_pml);
        }
        // Radial
        for( int j = 0; j < (int) nr_p ; j++ ) {
            kappa_r_p[j] = 1. ;
            sigma_r_p[j] = 0. ;
            integrate_kappa_r_p[j] = (j-startpml)*dr;
            integrate_sigma_r_p[j] = 0. ;
        }
        // Dual grid
        // Longitudinal
        // Params for first cell of PML-patch (vacuum) i = 0,1,2,3
        for ( int i=0 ; i<startpml+1 ; i++ ) {
            // Coeffs for the first cell
            kappa_l_d[i] = 1. ;
            sigma_l_d[i] = 0 ;
        }
        // Params for other cells (PML Media) when j>=4
        for( int i = startpml+1 ; i< (int) nl_d ; i++ ) {
            kappa_l_d[i] = pml_kappa_[0]->valueAt((i-startpml-0.5)*dl/length_l_pml);
            sigma_l_d[i] = pml_sigma_[0]->valueAt((i-startpml-0.5)*dl/length_l_pml);
        }
        // Radial
        for( int j = 0; j < (int) nr_d ; j++ ) {
            kappa_r_d[j] = 1. ;
            sigma_r_d[j] = 0. ;
            integrate_kappa_r_d[j] = (j-startpml-0.5)*dr;
            integrate_sigma_r_d[j] = 0. ;
        }
    }

    if ( iDim == 1 ) {
        // 3 cells are vaccum so the PML media begin at r0 which is :
        // Eventually the size of PML media is :
        length_r_pml = (ncells_pml_domain-startpml+0.5)*dr ;
        length_l_pml_lmax = (ncells_pml_max[0]+0.5)*dl;
        length_l_pml_lmin = (ncells_pml_min[0]+0.5)*dl;
        // Primal grid
        // Longitudinal
        for( unsigned int i = 0 ; i<nl_p ; i++ ) {
            kappa_l_p[i] = 1. ;
            sigma_l_p[i] = 0. ;
        }
        if (ncells_pml_min[0] != 0 ){
            for ( int i=0 ; i<ncells_pml_min[0] ; i++ ) {
                kappa_l_p[i] = pml_kappa_[0]->valueAt((ncells_pml_min[0] - 1 - i)*dl/length_l_pml_lmin);
                sigma_l_p[i] = pml_sigma_[0]->valueAt((ncells_pml_min[0] - 1 - i)*dl/length_l_pml_lmin);
            }
        }
        if (ncells_pml_max[0] != 0 ){
            for( int i = (nl_p-1)-(ncells_pml_max[0]-1) ; i< (int) nl_p ; i++ ) {
                kappa_l_p[i] = pml_kappa_[0]->valueAt((i - nl_p + ncells_pml_max[0])*dl/length_l_pml_lmax);
                sigma_l_p[i] = pml_sigma_[0]->valueAt((i - nl_p + ncells_pml_max[0])*dl/length_l_pml_lmax);
            }
        }
        // Radial
        // Params for first cell of PML-patch (vacuum) j = 0,1,2
        for ( int j=0 ; j<startpml ; j++ ) {
            // Coeffs for the first cell
            kappa_r_p[j] = 1. ;
            sigma_r_p[j] = 0. ;
            integrate_kappa_r_p[j] = (j-startpml)*dr;
            integrate_sigma_r_p[j] = 0. ;
        }
        // Params for other cells (PML Media) when j>=3
        for( int j = startpml; j < (int) nr_p ; j++) {
            kappa_r_p[j] = pml_kappa_[1]->valueAt((j-startpml)*dr/length_r_pml);
            sigma_r_p[j] = pml_sigma_[1]->valueAt((j-startpml)*dr/length_r_pml);
            integrate_kappa_r_p[j] = length_r_pml * pml_kappa_[2]->valueAt((j-startpml)*dr/length_r_pml);
            integrate_sigma_r_p[j] = length_r_pml * pml_sigma_[2]->valueAt((j-startpml)*dr/length_r_pml);
        }
        // Dual grid
        // Longitudinal
        for( unsigned int i = 0 ; i<nl_d ; i++ ) {
            kappa_l_d[i] = 1. ;
            sigma_l_d[i] = 0. ;
        }
        if (ncells_pml_min[0] != 0 ){
            for ( int i=0 ; i<ncells_pml_min[0] ; i++ ) {
                kappa_l_d[i] = pml_kappa_[0]->valueAt((( 0.5 + ncells_pml_min[0] - 1 - i ))*dl/length_l_pml_lmin);
                sigma_l_d[i] = pml_sigma_[0]->valueAt((( 0.5 + ncells_pml_min[0] - 1 - i ))*dl/length_l_pml_lmin);
            }
        }
        if (ncells_pml_max[0] != 0 ){
            for( int i = (nl_p-1)-(ncells_pml_max[0]-1)+1 ; i< (int) nl_d ; i++ ) {
                kappa_l_d[i] = pml_kappa_[0]->valueAt((i -  nl_p + ncells_pml_max[0] - 0.5 )*dl/length_l_pml_lmax);
                sigma_l_d[i] = pml_sigma_[0]->valueAt((i -  nl_p + ncells_pml_max[0] - 0.5 )*dl/length_l_pml_lmax);
            }
        }
        // Radial
        // Params for first cell of PML-patch (vacuum) j = 0,1,2,3
        for ( int j=0 ; j<startpml+1 ; j++ ) {
            // Coeffs for the first cell
            kappa_r_d[j] = 1. ;
            sigma_r_d[j] = 0. ;
            integrate_kappa_r_d[j] = (j-startpml-0.5)*dr;
            integrate_sigma_r_d[j] = 0 ;
        }
        // Params for other cells (PML Media) when j>=4
        for( int j = startpml+1 ; j< (int) nr_d ; j++) {
            kappa_r_d[j] = pml_kappa_[1]->valueAt((j-startpml-0.5)*dr/length_r_pml);
            sigma_r_d[j] = pml_sigma_[1]->valueAt((j-startpml-0.5)*dr/length_r_pml);
            integrate_kappa_r_d[j] = length_r_pml * pml_sigma_[2]->valueAt((j-startpml-0.5)*dr/length_r_pml);
            integrate_sigma_r_d[j] = length_r_pml * pml_sigma_[2]->valueAt((j-startpml-0.5)*dr/length_r_pml);
        }
    }

    if ((min_or_max==0)&&(iDim==0)){
        for( int i = 0; i < (int) nl_p ; i++ ) {
            //longitudinal-field-coeff
            c5_p_lfield[i] = 2.*kappa_l_p[(nl_p-1)-i] + dt*sigma_l_p[(nl_p-1)-i] ;
            c6_p_lfield[i] = 2.*kappa_l_p[(nl_p-1)-i] - dt*sigma_l_p[(nl_p-1)-i] ;
            //radial-field-coeff
            c1_p_rfield[i] = ( 2.*kappa_l_p[(nl_p-1)-i] - dt*sigma_l_p[(nl_p-1)-i] ) / ( 2.*kappa_l_p[(nl_p-1)-i] + dt*sigma_l_p[(nl_p-1)-i] ) ;
            c2_p_rfield[i] = ( 2*dt ) / ( 2.*kappa_l_p[(nl_p-1)-i] + dt*sigma_l_p[(nl_p-1)-i] ) ;
            //theta-field-coeff
            c3_p_tfield[i] = ( 2.*kappa_l_p[(nl_p-1)-i] - dt*sigma_l_p[(nl_p-1)-i] ) / ( 2.*kappa_l_p[(nl_p-1)-i] + dt*sigma_l_p[(nl_p-1)-i] ) ;
            c4_p_tfield[i] = (1.) / ( 2.*kappa_l_p[(nl_p-1)-i] + dt*sigma_l_p[(nl_p-1)-i] ) ;
        }

        for( int i = 0; i < (int) nl_d ; i++ ) {
            //longitudinal-field-coeff
            c5_d_lfield[i] = 2.*kappa_l_d[(nl_d-1)-i] + dt*sigma_l_d[(nl_d-1)-i] ;
            c6_d_lfield[i] = 2.*kappa_l_d[(nl_d-1)-i] - dt*sigma_l_d[(nl_d-1)-i] ;
            //radial-field-coeff
            c1_d_rfield[i] = ( 2.*kappa_l_d[(nl_d-1)-i] - dt*sigma_l_d[(nl_d-1)-i] ) / ( 2.*kappa_l_d[(nl_d-1)-i] + dt*sigma_l_d[(nl_d-1)-i] ) ;
            c2_d_rfield[i] = ( 2*dt ) / ( 2.*kappa_l_d[(nl_d-1)-i] + dt*sigma_l_d[(nl_d-1)-i] ) ;
            //theta-field-coeff
            c3_d_tfield[i] = ( 2.*kappa_l_d[(nl_d-1)-i] - dt*sigma_l_d[(nl_d-1)-i] ) / ( 2.*kappa_l_d[(nl_d-1)-i] + dt*sigma_l_d[(nl_d-1)-i] ) ;
            c4_d_tfield[i] = (1.) / ( 2.*kappa_l_d[(nl_d-1)-i] + dt*sigma_l_d[(nl_d-1)-i] ) ;
        }
    }
    else {
        for( int i = 0; i < (int) nl_p ; i++ ) {
            //longitudinal-field-coeff
            c5_p_lfield[i] = 2.*kappa_l_p[i] + dt*sigma_l_p[i] ;
            c6_p_lfield[i] = 2.*kappa_l_p[i] - dt*sigma_l_p[i] ;
            //radial-field-coeff
            c1_p_rfield[i] = ( 2.*kappa_l_p[i] - dt*sigma_l_p[i] ) / ( 2.*kappa_l_p[i] + dt*sigma_l_p[i] ) ;
            c2_p_rfield[i] = ( 2*dt ) / ( 2.*kappa_l_p[i] + dt*sigma_l_p[i] ) ;
            //theta-field-coeff
            c3_p_tfield[i] = ( 2.*kappa_l_p[i] - dt*sigma_l_p[i] ) / ( 2.*kappa_l_p[i] + dt*sigma_l_p[i] ) ;
            c4_p_tfield[i] = (1.) / ( 2.*kappa_l_p[i] + dt*sigma_l_p[i] ) ;
        }

        for( int i = 0; i < (int) nl_d ; i++ ) {
            //longitudinal-field-coeff
            c5_d_lfield[i] = 2.*kappa_l_d[i] + dt*sigma_l_d[i] ;
            c6_d_lfield[i] = 2.*kappa_l_d[i] - dt*sigma_l_d[i] ;
            //radial-field-coeff
            c1_d_rfield[i] = ( 2.*kappa_l_d[i] - dt*sigma_l_d[i] ) / ( 2.*kappa_l_d[i] + dt*sigma_l_d[i] ) ;
            c2_d_rfield[i] = ( 2*dt ) / ( 2.*kappa_l_d[i] + dt*sigma_l_d[i] ) ;
            //theta-field-coeff
            c3_d_tfield[i] = ( 2.*kappa_l_d[i] - dt*sigma_l_d[i] ) / ( 2.*kappa_l_d[i] + dt*sigma_l_d[i] ) ;
            c4_d_tfield[i] = (1.) / ( 2.*kappa_l_d[i] + dt*sigma_l_d[i] ) ;
        }
    } // End X

    if ((min_or_max==0)&&(iDim==0)){
        for( int j = 0; j < (int) nr_p ; j++ ) {
            //longitudinal-field-coeff
            c1_p_lfield[j] = ( 2.*kappa_r_p[(nr_p-1)-j] - dt*sigma_r_p[(nr_p-1)-j] ) / ( 2.*kappa_r_p[(nr_p-1)-j] + dt*sigma_r_p[(nr_p-1)-j] ) ;
            c2_p_lfield[j] = ( 2*dt ) / ( 2.*kappa_r_p[(nr_p-1)-j] + dt*sigma_r_p[(nr_p-1)-j] ) ;
            c3_p_lfield[j] = ( 2.*( r0 + integrate_kappa_r_p[(nr_p-1)-j] ) - dt*integrate_sigma_r_p[(nr_p-1)-j] ) / ( 2.*( r0 + integrate_kappa_r_p[(nr_p-1)-j] ) + dt*integrate_sigma_r_p[(nr_p-1)-j] ) ;
            c4_p_lfield[j] = ( rmin + ((nr_p-1)-j)*dr ) / ( 2.*( r0 + integrate_kappa_r_p[(nr_p-1)-j] ) + dt*integrate_sigma_r_p[(nr_p-1)-j] ) ;
            //radial-field-coeff
            c3_p_rfield[j] = ( 2.*( r0 + integrate_kappa_r_p[(nr_p-1)-j] ) - dt*integrate_sigma_r_p[(nr_p-1)-j] ) / ( 2.*( r0 + integrate_kappa_r_p[(nr_p-1)-j] ) + dt*integrate_sigma_r_p[(nr_p-1)-j] ) ;
            c4_p_rfield[j] = ( rmin + ((nr_p-1)-j)*dr ) / ( 2.*( r0 + integrate_kappa_r_p[(nr_p-1)-j] ) + dt*integrate_sigma_r_p[(nr_p-1)-j] ) ;
            c5_p_rfield[j] = 2.*kappa_r_p[(nr_p-1)-j] + dt*sigma_r_p[(nr_p-1)-j] ;
            c6_p_rfield[j] = 2.*kappa_r_p[(nr_p-1)-j] - dt*sigma_r_p[(nr_p-1)-j] ;
            //theta-field-coeff
            c1_p_tfield[j] = ( 2.*kappa_r_p[(nr_p-1)-j] - dt*sigma_r_p[(nr_p-1)-j] ) / ( 2.*kappa_r_p[(nr_p-1)-j] + dt*sigma_r_p[(nr_p-1)-j] ) ;
            c2_p_tfield[j] = ( 2*dt ) / ( 2.*kappa_r_p[(nr_p-1)-j] + dt*sigma_r_p[(nr_p-1)-j] ) ;
            c5_p_tfield[j] = ( 2.*( r0 + integrate_kappa_r_p[(nr_p-1)-j] ) + dt*integrate_sigma_r_p[(nr_p-1)-j] ) / ( rmin + ((nr_p-1)-j)*dr ) ;
            c6_p_tfield[j] = ( 2.*( r0 + integrate_kappa_r_p[(nr_p-1)-j] ) - dt*integrate_sigma_r_p[(nr_p-1)-j] ) / ( rmin + ((nr_p-1)-j)*dr ) ;
        }

        for( int j = 0; j < (int) nr_d ; j++ ) {
            //longitudinal-field-coeff
            c1_d_lfield[j] = ( 2.*kappa_r_d[(nr_d-1)-j] - dt*sigma_r_d[(nr_d-1)-j] ) / ( 2.*kappa_r_d[(nr_d-1)-j] + dt*sigma_r_d[(nr_d-1)-j] ) ;
            c2_d_lfield[j] = ( 2*dt ) / ( 2.*kappa_r_d[(nr_d-1)-j] + dt*sigma_r_d[(nr_d-1)-j] ) ;
            c3_d_lfield[j] = ( 2.*( r0 + integrate_kappa_r_d[(nr_d-1)-j] ) - dt*integrate_sigma_r_d[(nr_d-1)-j] ) / ( 2.*( r0 + integrate_kappa_r_d[(nr_d-1)-j] ) + dt*integrate_sigma_r_d[(nr_d-1)-j] ) ;
            c4_d_lfield[j] = ( rmin + (((nr_d-1)-j)-0.5)*dr ) / ( 2.*( r0 + integrate_kappa_r_d[(nr_d-1)-j] ) + dt*integrate_sigma_r_d[(nr_d-1)-j] ) ;
            //radial-field-coeff
            c3_d_rfield[j] = ( 2.*( r0 + integrate_kappa_r_d[(nr_d-1)-j] ) - dt*integrate_sigma_r_d[(nr_d-1)-j] ) / ( 2.*( r0 + integrate_kappa_r_d[(nr_d-1)-j] ) + dt*integrate_sigma_r_d[(nr_d-1)-j] ) ;
            c4_d_rfield[j] = ( rmin + (((nr_d-1)-j)-0.5)*dr ) / ( 2.*( r0 + integrate_kappa_r_d[(nr_d-1)-j] ) + dt*integrate_sigma_r_d[(nr_d-1)-j] ) ;
            c5_d_rfield[j] = 2.*kappa_r_d[(nr_d-1)-j] + dt*sigma_r_d[(nr_d-1)-j] ;
            c6_d_rfield[j] = 2.*kappa_r_d[(nr_d-1)-j] - dt*sigma_r_d[(nr_d-1)-j] ;
            //theta-field-coeff
            c1_d_tfield[j] = ( 2.*kappa_r_d[(nr_d-1)-j] - dt*sigma_r_d[(nr_d-1)-j] ) / ( 2.*kappa_r_d[(nr_d-1)-j] + dt*sigma_r_d[(nr_d-1)-j] ) ;
            c2_d_tfield[j] = ( 2*dt ) / ( 2.*kappa_r_d[(nr_d-1)-j] + dt*sigma_r_d[(nr_d-1)-j] ) ;
            c5_d_tfield[j] = ( 2.*( r0 + integrate_kappa_r_d[(nr_d-1)-j] ) + dt*integrate_sigma_r_d[(nr_d-1)-j] ) / ( rmin + (((nr_d-1)-j)-0.5)*dr ) ;
            c6_d_tfield[j] = ( 2.*( r0 + integrate_kappa_r_d[(nr_d-1)-j] ) - dt*integrate_sigma_r_d[(nr_d-1)-j] ) / ( rmin + (((nr_d-1)-j)-0.5)*dr ) ;
        }
    }
    else {
        for( int j = 0; j < (int) nr_p ; j++ ) {
            //longitudinal-field-coeff
            c1_p_lfield[j] = ( 2.*kappa_r_p[j] - dt*sigma_r_p[j] ) / ( 2.*kappa_r_p[j] + dt*sigma_r_p[j] ) ;
            c2_p_lfield[j] = ( 2*dt ) / ( 2.*kappa_r_p[j] + dt*sigma_r_p[j] ) ;
            c3_p_lfield[j] = ( 2.*( r0 + integrate_kappa_r_p[j] ) - dt*integrate_sigma_r_p[j] ) / ( 2.*( r0 + integrate_kappa_r_p[j] ) + dt*integrate_sigma_r_p[j] ) ;
            c4_p_lfield[j] = ( rmin + j*dr ) / ( 2.*( r0 + integrate_kappa_r_p[j] ) + dt*integrate_sigma_r_p[j] ) ;
            //radial-field-coeff
            c3_p_rfield[j] = ( 2.*( r0 + integrate_kappa_r_p[j] ) - dt*integrate_sigma_r_p[j] ) / ( 2.*( r0 + integrate_kappa_r_p[j] ) + dt*integrate_sigma_r_p[j] ) ;
            c4_p_rfield[j] = ( rmin + j*dr ) / ( 2.*( r0 + integrate_kappa_r_p[j] ) + dt*integrate_sigma_r_p[j] ) ;
            c5_p_rfield[j] = 2.*kappa_r_p[j] + dt*sigma_r_p[j] ;
            c6_p_rfield[j] = 2.*kappa_r_p[j] - dt*sigma_r_p[j] ;
            //theta-field-coeff
            c1_p_tfield[j] = ( 2.*kappa_r_p[j] - dt*sigma_r_p[j] ) / ( 2.*kappa_r_p[j] + dt*sigma_r_p[j] ) ;
            c2_p_tfield[j] = ( 2*dt ) / ( 2.*kappa_r_p[j] + dt*sigma_r_p[j] ) ;
            c5_p_tfield[j] = ( 2.*( r0 + integrate_kappa_r_p[j] ) + dt*integrate_sigma_r_p[j] ) / ( rmin + j*dr ) ;
            c6_p_tfield[j] = ( 2.*( r0 + integrate_kappa_r_p[j] ) - dt*integrate_sigma_r_p[j] ) / ( rmin + j*dr ) ;
        }

        for( int j = 0; j < (int) nr_d ; j++ ) {
            //longitudinal-field-coeff
            c1_d_lfield[j] = ( 2.*kappa_r_d[j] - dt*sigma_r_d[j] ) / ( 2.*kappa_r_d[j] + dt*sigma_r_d[j] ) ;
            c2_d_lfield[j] = ( 2*dt ) / ( 2.*kappa_r_d[j] + dt*sigma_r_d[j] ) ;
            c3_d_lfield[j] = ( 2.*( r0 + integrate_kappa_r_d[j] ) - dt*integrate_sigma_r_d[j] ) / ( 2.*( r0 + integrate_kappa_r_d[j] ) + dt*integrate_sigma_r_d[j] ) ;
            c4_d_lfield[j] = ( rmin + (j-0.5)*dr ) / ( 2.*( r0 + integrate_kappa_r_d[j] ) + dt*integrate_sigma_r_d[j] ) ;
            //radial-field-coeff
            c3_d_rfield[j] = ( 2.*( r0 + integrate_kappa_r_d[j] ) - dt*integrate_sigma_r_d[j] ) / ( 2.*( r0 + integrate_kappa_r_d[j] ) + dt*integrate_sigma_r_d[j] ) ;
            c4_d_rfield[j] = ( rmin + (j-0.5)*dr ) / ( 2.*( r0 + integrate_kappa_r_d[j] ) + dt*integrate_sigma_r_d[j] ) ;
            c5_d_rfield[j] = 2.*kappa_r_d[j] + dt*sigma_r_d[j] ;
            c6_d_rfield[j] = 2.*kappa_r_d[j] - dt*sigma_r_d[j] ;
            //theta-field-coeff
            c1_d_tfield[j] = ( 2.*kappa_r_d[j] - dt*sigma_r_d[j] ) / ( 2.*kappa_r_d[j] + dt*sigma_r_d[j] ) ;
            c2_d_tfield[j] = ( 2*dt ) / ( 2.*kappa_r_d[j] + dt*sigma_r_d[j] ) ;
            c5_d_tfield[j] = ( 2.*( r0 + integrate_kappa_r_d[j] ) + dt*integrate_sigma_r_d[j] ) / ( rmin + (j-0.5)*dr ) ;
            c6_d_tfield[j] = ( 2.*( r0 + integrate_kappa_r_d[j] ) - dt*integrate_sigma_r_d[j] ) / ( rmin + (j-0.5)*dr ) ;
        }
    } //  End Y
}

void PML_SolverAM::compute_E_from_D( ElectroMagn *fields, int iDim, int min_or_max, std::vector<unsigned int> dimPrim, unsigned int solvermin, unsigned int solvermax )
{
    const unsigned int nl_p = dimPrim[0];
    const unsigned int nl_d = dimPrim[0] + 1;
    const unsigned int nr_p = dimPrim[1];
    const unsigned int nr_d = dimPrim[1] + 1;
    
    ElectroMagnBCAM_PML* pml_fields = static_cast<ElectroMagnBCAM_PML*>( fields->emBoundCond[iDim*2+min_or_max] );
    cField2D* El_pml = NULL;
    cField2D* Er_pml = NULL;
    cField2D* Et_pml = NULL;
    cField2D* Hl_pml = NULL;
    cField2D* Hr_pml = NULL;
    cField2D* Ht_pml = NULL;
    cField2D* Dl_pml = NULL;
    cField2D* Dr_pml = NULL;
    cField2D* Dt_pml = NULL;

    bool isYmin = ( static_cast<ElectroMagnAM *>( fields ) )->isYmin;

    if (min_or_max==0){
        isMin=true;
        isMax=false;
    }
    else if (min_or_max==1){
        isMin=false;
        isMax=true;
    }

    if (iDim==0){
        for( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
            El_pml = pml_fields->El_[imode];
            Er_pml = pml_fields->Er_[imode];
            Et_pml = pml_fields->Et_[imode];
            Hl_pml = pml_fields->Hl_[imode];
            Hr_pml = pml_fields->Hr_[imode];
            Ht_pml = pml_fields->Ht_[imode];
            Dl_pml = pml_fields->Dl_[imode];
            Dr_pml = pml_fields->Dr_[imode];
            Dt_pml = pml_fields->Dt_[imode];
            //Electric field El^(d,p) Remind that in PML, there no current
            for( unsigned int i=solvermin ; i<(unsigned int)solvermax ; i++ ) {
                for( unsigned int j=isYmin*3 ; j<nr_p ; j++ ) {
                    // Standard FDTD
                    // ( *Dl_pml )( i, j ) = + ( *Dl_pml )( i, j )
                    //                     + dt / ( ( j_glob_pml+j )*dr ) * ( ( j+j_glob_pml+0.5 ) * ( *Ht_pml )( i, j+1 ) - ( j+j_glob_pml-0.5 )*( *Ht_pml )( i, j ) )
                    //                     + dt / ( ( j_glob_pml+j )*dr ) * Icpx*( double )imode*( *Hr_pml )( i, j ) ;
                    // ( *El_pml )( i, j ) = 1*( *Dl_pml )( i, j );
                    // PML FDTD
                    Dl_pml_old = std::complex<double>(1,0)*( *Dl_pml )( i, j ) ;
                    ( *Dl_pml )( i, j ) = + c1_p_lfield[j] * ( *Dl_pml )( i, j )
                                          + c2_p_lfield[j] / ( ( j_glob_pml+j )*dr ) * ( ( j+j_glob_pml+0.5 ) * ( *Ht_pml )( i, j+1 ) - ( j+j_glob_pml-0.5 )*( *Ht_pml )( i, j ) )
                                          + c2_p_lfield[j] / ( ( j_glob_pml+j )*dr ) * Icpx*( double )imode*( *Hr_pml )( i, j ) ;
                    ( *El_pml )( i, j ) = + c3_p_lfield[j] * ( *El_pml )( i, j )
                                          + c4_p_lfield[j] * (c5_d_lfield[i]*( *Dl_pml )( i, j ) - c6_d_lfield[i]*Dl_pml_old ) ;
                }
            }
            //Electric field Er^(p,d) Remind that in PML, there no current
            for( unsigned int i=solvermin ; i<(unsigned int)solvermax ; i++ ) {
                for( unsigned int j=isYmin*3 ; j<nr_d ; j++ ) {
                    // Standard FDTD
                    // ( *Dr_pml )( i, j ) = + ( *Dr_pml )( i, j )
                    //                     - dt/dl * ( ( *Ht_pml )( i+1, j ) - ( *Ht_pml )( i, j ) )
                    //                     - dt*Icpx*( double )imode/( ( j+j_glob_pml-0.5 )*dr )*( *Hl_pml )( i, j ) ;
                    // ( *Er_pml )( i, j ) = 1*( *Dr_pml )( i, j );
                    // PML FDTD
                    Dr_pml_old = std::complex<double>(1,0)*( *Dr_pml )( i, j ) ;
                    ( *Dr_pml )( i, j ) = + c1_p_rfield[i] * ( *Dr_pml )( i, j )
                                          - c2_p_rfield[i]/dl * ( ( *Ht_pml )( i+1, j ) - ( *Ht_pml )( i, j ) )
                                          - c2_p_rfield[i]*Icpx*( double )imode/( ( j+j_glob_pml-0.5 )*dr)*( *Hl_pml )( i, j ) ;
                    ( *Er_pml )( i, j ) = + c3_d_rfield[j] * ( *Er_pml )( i, j )
                                          + c4_d_rfield[j] * (c5_d_rfield[j]*( *Dr_pml )( i, j ) - c6_d_rfield[j]*Dr_pml_old );
                }
            }
            //Electric field Et^(p,p) Remind that in PML, there no current
            for( unsigned int i=solvermin ; i<(unsigned int)solvermax ; i++ ) {
                for( unsigned int j=isYmin*3 ; j<nr_p ; j++ ) {
                    // Standard FDTD
                    // ( *Dt_pml )( i, j ) = + 1 * ( *Dt_pml )( i, j )
                    //                     - dt/dr * ( ( *Hl_pml )( i, j+1 ) - ( *Hl_pml )( i, j ) )
                    //                     + dt/dl * ( ( *Hr_pml )( i+1, j ) - ( *Hr_pml )( i, j ) ) ;
                    // ( *Et_pml )( i, j ) = 1*( *Dt_pml )( i, j );
                    // PML FDTD
                    Dt_pml_old = std::complex<double>(1,0)*( *Dt_pml )( i, j ) ;
                    ( *Dt_pml )( i, j ) = + c1_p_tfield[j] * ( *Dt_pml )( i, j )
                                          - c2_p_tfield[j]/dr * ( ( *Hl_pml )( i, j+1 ) - ( *Hl_pml )( i, j ) )
                                          + c2_p_tfield[j]/dl * ( ( *Hr_pml )( i+1, j ) - ( *Hr_pml )( i, j ) ) ;
                    ( *Et_pml )( i, j ) = + c3_p_tfield[i] * ( *Et_pml )( i, j )
                                          + c4_p_tfield[i] * (c5_p_tfield[j]*( *Dt_pml )( i, j ) - c6_p_tfield[j]*Dt_pml_old );
                }
            }
            if( isYmin ) {
                // Conditions on axis
                unsigned int j=2;
                if( imode==0 ) {
                    for( unsigned int i=0 ; i<nl_p  ; i++ ) {
                        ( *Et_pml )( i, j )  = 0;
                        ( *Et_pml )( i, j-1 )= -( *Et_pml )( i, j+1 );
                        //
                        ( *Dt_pml )( i, j ) = ( *Et_pml )( i, j ) ;
                        ( *Dt_pml )( i, j ) = ( *Et_pml )( i, j-1 ) ;
                    }
                    for( unsigned int i=0 ; i<nl_p  ; i++ ) {
                        ( *Er_pml )( i, j ) = -( *Er_pml )( i, j+1 );
                        //
                        ( *Dr_pml )( i, j ) = ( *Er_pml )( i, j ) ;
                    }
                    for( unsigned int i=0 ; i<nl_d ; i++ ) {
                        ( *El_pml )( i, j )  += 4.*dt/dr*( *Ht_pml )( i, j+1 );
                        ( *El_pml )( i, j-1 ) =( *El_pml )( i, j+1 );
                        //
                        ( *Dl_pml )( i, j ) = ( *El_pml )( i, j ) ;
                        ( *Dl_pml )( i, j-1 ) = ( *El_pml )( i, j-1 ) ;
                    }
                } else if( imode==1 ) {
                    for( unsigned int i=0 ; i<nl_d  ; i++ ) {
                        ( *El_pml )( i, j )= 0;
                        ( *El_pml )( i, j-1 )=-( *El_pml )( i, j+1 );
                        //
                        ( *Dl_pml )( i, j ) = ( *El_pml )( i, j ) ;
                        ( *Dl_pml )( i, j-1 ) = ( *El_pml )( i, j-1 ) ;
                    }
                    for( unsigned int i=0 ; i<nl_p  ; i++ ) {
                        ( *Et_pml )( i, j )= -Icpx/8.*( 9.*( *Er_pml )( i, j+1 )-( *Er_pml )( i, j+2 ) );
                        ( *Et_pml )( i, j-1 )=( *Et_pml )( i, j+1 );
                        //
                        ( *Dt_pml )( i, j ) = ( *Et_pml )( i, j ) ; 
                        ( *Dt_pml )( i, j-1 ) = ( *Et_pml )( i, j-1 ) ;
                    }
                    for( unsigned int i=0 ; i<nl_p ; i++ ) {
                        ( *Er_pml )( i, j )=2.*Icpx*( *Et_pml )( i, j )-( *Er_pml )( i, j+1 );
                        //
                        ( *Dr_pml )( i, j ) = ( *Er_pml )( i, j ) ;
                    }
                } else { // mode > 1
                    for( unsigned int  i=0 ; i<nl_d; i++ ) {
                        ( *El_pml )( i, j )= 0;
                        ( *El_pml )( i, j-1 )=-( *El_pml )( i, j+1 );
                        //
                        ( *Dl_pml )( i, j ) = ( *El_pml )( i, j ) ;
                        ( *Dl_pml )( i, j-1 ) = ( *El_pml )( i, j-1 ) ;
                    }
                    for( unsigned int  i=0 ; i<nl_p; i++ ) {
                        ( *Er_pml )( i, j+1 )= ( *Er_pml )( i, j+2 ) / 9.;
                        ( *Er_pml )( i, j )= -( *Er_pml )( i, j+1 );
                        //
                        ( *Dr_pml )( i, j+1 ) = ( *Er_pml )( i, j+1 ) ;
                        ( *Dr_pml )( i, j ) = ( *Er_pml )( i, j ) ;
                    }
                    for( unsigned int i=0 ; i<nl_p; i++ ) {
                        ( *Et_pml )( i, j )= 0;
                        ( *Et_pml )( i, j-1 )=-( *Et_pml )( i, j+1 );
                        //
                        ( *Dt_pml )( i, j ) = ( *Et_pml )( i, j ) ;
                        ( *Dt_pml )( i, j-1 ) = ( *Et_pml )( i, j-1 );
                    }
                }
            }
            // if( isYmin ) {
            //     // Conditions on axis
            //     unsigned int j=2;
            //     if( imode==0 ) {
            //         for( unsigned int i=3*isMax ; i<nl_p-3*isMin  ; i++ ) {
            //             ( *Dt_pml )( i, j )   = 0 ;
            //             ( *Et_pml )( i, j )   = 0 ;

            //             Dt_pml_old = ( *Dt_pml )( i, j-1 ) ;
            //             ( *Dt_pml )( i, j-1 ) = -( *Dt_pml )( i, j+1 ) ;
            //             ( *Et_pml )( i, j-1 ) = + c3_p_tfield[i] * ( *Et_pml )( i, j-1 )
            //                                     + c4_p_tfield[i] * (c5_p_tfield[j]*( *Dt_pml )( i, j-1 ) - c6_p_tfield[j]*Dt_pml_old ) ;
            //         }
            //         for( unsigned int i=3*isMax ; i<nl_p-3*isMin ; i++ ) {
            //             Dr_pml_old = ( *Dr_pml )( i, j ) ;
            //             ( *Dr_pml )( i, j )   = -( *Dr_pml )( i, j+1 ) ;
            //             ( *Er_pml )( i, j )   = + c3_d_rfield[j] * ( *Er_pml )( i, j )
            //                                     + c4_d_rfield[j] * (c5_d_rfield[j]*( *Dr_pml )( i, j ) - c6_d_rfield[j]*Dr_pml_old );
            //         }
            //         for( unsigned int i=3*isMax ; i<nl_d-3*isMin ; i++ ) {
            //             Dl_pml_old            = ( *Dl_pml )( i, j ) ;
            //             ( *Dl_pml )( i, j )   = + c1_p_lfield[j] * ( *Dl_pml )( i, j )
            //                                     + 4*c2_p_lfield[j]/dr * ( *Ht_pml )( i, j+1 ) ;
            //             ( *El_pml )( i, j )   = + c3_p_lfield[j] * ( *El_pml )( i, j )
            //                                     + c4_p_lfield[j] * (c5_d_lfield[i]*( *Dl_pml )( i, j ) - c6_d_lfield[i]*Dl_pml_old ) ;

            //             Dl_pml_old            = ( *Dl_pml )( i, j-1 ) ;
            //             ( *Dl_pml )( i, j-1 ) = ( *Dl_pml )( i, j+1 ) ;
            //             ( *El_pml )( i, j-1 ) = + c3_p_lfield[j] * ( *El_pml )( i, j-1 )
            //                                     + c4_p_lfield[j] * (c5_d_lfield[i]*( *Dl_pml )( i, j-1 ) - c6_d_lfield[i]*Dl_pml_old ) ;
            //         }
            //     } else if( imode==1 ) {
            //         for( unsigned int i=3*isMax ; i<nl_d-3*isMin  ; i++ ) {
            //             ( *Dl_pml )( i, j )   = 0 ;
            //             ( *El_pml )( i, j )   = 0 ;

            //             Dl_pml_old            = ( *Dl_pml )( i, j-1 ) ;
            //             ( *Dl_pml )( i, j-1 ) = -( *Dl_pml )( i, j+1 ) ;
            //             ( *El_pml )( i, j-1 ) = + c3_p_lfield[j] * ( *El_pml )( i, j-1 )
            //                                     + c4_p_lfield[j] * (c5_d_lfield[i]*( *Dl_pml )( i, j-1 ) - c6_d_lfield[i]*Dl_pml_old ) ;
            //         }
            //         for( unsigned int i=3*isMax ; i<nl_p-3*isMin  ; i++ ) {
            //             Dt_pml_old = ( *Dt_pml )( i, j ) ;
            //             ( *Dt_pml )( i, j )   = -Icpx/8.*( 9.*( *Dr_pml )( i, j+1 )-( *Dr_pml )( i, j+2 ) );
            //             ( *Et_pml )( i, j )   = + c3_p_tfield[i] * ( *Et_pml )( i, j )
            //                                     + c4_p_tfield[i] * (c5_p_tfield[j]*( *Dt_pml )( i, j ) - c6_p_tfield[j]*Dt_pml_old ) ;


            //             Dt_pml_old = ( *Dt_pml )( i, j-1 ) ;
            //             ( *Dt_pml )( i, j-1 ) = ( *Dt_pml )( i, j+1 );
            //             ( *Et_pml )( i, j-1 ) = + c3_p_tfield[i] * ( *Et_pml )( i, j-1 )
            //                                     + c4_p_tfield[i] * (c5_p_tfield[j]*( *Dt_pml )( i, j-1 ) - c6_p_tfield[j]*Dt_pml_old ) ;
            //         }
            //         for( unsigned int i=3*isMax ; i<nl_p-3*isMin ; i++ ) {
            //             Dr_pml_old = ( *Dr_pml )( i, j ) ;
            //             ( *Dr_pml )( i, j )   = 2.*Icpx*( *Dt_pml )( i, j )-( *Dr_pml )( i, j+1 );
            //             ( *Er_pml )( i, j )   = + c3_d_rfield[j] * ( *Er_pml )( i, j )
            //                                     + c4_d_rfield[j] * (c5_d_rfield[j]*( *Dr_pml )( i, j ) - c6_d_rfield[j]*Dr_pml_old );
            //         }
            //     } else { // mode > 1
            //         for( unsigned int  i=3*isMax ; i<nl_d-3*isMin; i++ ) {
            //             ( *Dl_pml )( i, j )   = 0 ;
            //             ( *El_pml )( i, j )   = 0 ;

            //             Dl_pml_old = ( *Dl_pml )( i, j-1 ) ;
            //             ( *Dl_pml )( i, j-1 ) = -( *Dl_pml )( i, j+1 ) ;
            //             ( *El_pml )( i, j-1 ) = + c3_p_lfield[j] * ( *El_pml )( i, j-1 )
            //                                     + c4_p_lfield[j] * (c5_d_lfield[i]*( *Dl_pml )( i, j-1 ) - c6_d_lfield[i]*Dl_pml_old ) ;
            //         }
            //         for( unsigned int  i=3*isMax ; i<nl_p-3*isMin; i++ ) {
            //             Dr_pml_old = ( *Dr_pml )( i, j+1 ) ;
            //             ( *Dr_pml )( i, j+1 ) = ( *Dr_pml )( i, j+2 ) / 9.;
            //             ( *Er_pml )( i, j+1 ) = + c3_d_rfield[j] * ( *Er_pml )( i, j+1 )
            //                                     + c4_d_rfield[j] * (c5_d_rfield[j]*( *Dr_pml )( i, j+1 ) - c6_d_rfield[j]*Dr_pml_old );

            //             Dr_pml_old = ( *Dr_pml )( i, j ) ;
            //             ( *Dr_pml )( i, j )   = -( *Dr_pml )( i, j+1 ) ;
            //             ( *Er_pml )( i, j )   = + c3_d_rfield[j] * ( *Er_pml )( i, j )
            //                                     + c4_d_rfield[j] * (c5_d_rfield[j]*( *Dr_pml )( i, j ) - c6_d_rfield[j]*Dr_pml_old );
            //         }
            //         for( unsigned int i=3*isMax ; i<nl_p-3*isMin; i++ ) {
            //             ( *Dt_pml )( i, j )    = 0 ;
            //             ( *Et_pml )( i, j )    = 0 ;

            //             Dt_pml_old = ( *Dt_pml )( i, j-1 );
            //             ( *Dt_pml )( i, j-1 )  = -( *Dt_pml )( i, j+1 );
            //             ( *Et_pml )( i, j-1 ) = + c3_p_tfield[i] * ( *Et_pml )( i, j-1 )
            //                                     + c4_p_tfield[i] * (c5_p_tfield[j]*( *Dt_pml )( i, j-1 ) - c6_p_tfield[j]*Dt_pml_old ) ;
            //         }
            //     }
            // }
        }
    }
    else if (iDim==1){
        for( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
            El_pml = pml_fields->El_[imode];
            Er_pml = pml_fields->Er_[imode];
            Et_pml = pml_fields->Et_[imode];
            Hl_pml = pml_fields->Hl_[imode];
            Hr_pml = pml_fields->Hr_[imode];
            Ht_pml = pml_fields->Ht_[imode];
            Dl_pml = pml_fields->Dl_[imode];
            Dr_pml = pml_fields->Dr_[imode];
            Dt_pml = pml_fields->Dt_[imode];
            //Electric field El^(d,p) Remind that in PML, there no current
            for( unsigned int i=0 ; i<nl_d ; i++ ) {
                for( unsigned int j=solvermin ; j<(unsigned int)solvermax ; j++ ) {
                    // Standard FDTD
                    // ( *Dl_pml )( i, j ) = + ( *Dl_pml )( i, j )
                    //                     + dt / ( ( j_glob_pml+j )*dr ) * ( ( j+j_glob_pml+0.5 ) * ( *Ht_pml )( i, j+1 ) - ( j+j_glob_pml-0.5 )*( *Ht_pml )( i, j ) )
                    //                     + dt / ( ( j_glob_pml+j )*dr ) * Icpx*( double )imode*( *Hr_pml )( i, j ) ;
                    // ( *El_pml )( i, j ) = 1*( *Dl_pml )( i, j );
                    // PML FDTD
                    Dl_pml_old = std::complex<double>(1,0)*( *Dl_pml )( i, j ) ;
                    ( *Dl_pml )( i, j ) = + c1_p_lfield[j] * ( *Dl_pml )( i, j )
                                          + c2_p_lfield[j] / ( ( j_glob_pml+j )*dr ) * ( ( j+j_glob_pml+0.5 ) * ( *Ht_pml )( i, j+1 ) - ( j+j_glob_pml-0.5 )*( *Ht_pml )( i, j ) )
                                          + c2_p_lfield[j] / ( ( j_glob_pml+j )*dr ) * Icpx*( double )imode*( *Hr_pml )( i, j ) ;
                    ( *El_pml )( i, j ) = + c3_p_lfield[j] * ( *El_pml )( i, j )
                                          + c4_p_lfield[j] * (c5_d_lfield[i]*( *Dl_pml )( i, j ) - c6_d_lfield[i]*Dl_pml_old ) ;
                }
            }
            //Electric field Er^(p,d) Remind that in PML, there no current
            for( unsigned int i=0 ; i<nl_p ; i++ ) {
                for( unsigned int j=solvermin ; j<(unsigned int)solvermax ; j++ ) {
                    // Standard FDTD
                    // ( *Dr_pml )( i, j ) = + ( *Dr_pml )( i, j )
                    //                     - dt/dl * ( ( *Ht_pml )( i+1, j ) - ( *Ht_pml )( i, j ) )
                    //                     - dt*Icpx*( double )imode/( ( j+j_glob_pml-0.5 )*dr )*( *Hl_pml )( i, j ) ;
                    // ( *Er_pml )( i, j ) = 1*( *Dr_pml )( i, j );
                    // PML FDTD
                    Dr_pml_old = std::complex<double>(1,0)*( *Dr_pml )( i, j ) ;
                    ( *Dr_pml )( i, j ) = + c1_p_rfield[i] * ( *Dr_pml )( i, j )
                                          - c2_p_rfield[i]/dl * ( ( *Ht_pml )( i+1, j ) - ( *Ht_pml )( i, j ) )
                                          - c2_p_rfield[i]*Icpx*( double )imode/( ( j+j_glob_pml-0.5 )*dr)*( *Hl_pml )( i, j ) ;
                    ( *Er_pml )( i, j ) = + c3_d_rfield[j] * ( *Er_pml )( i, j )
                                          + c4_d_rfield[j] * (c5_d_rfield[j]*( *Dr_pml )( i, j ) - c6_d_rfield[j]*Dr_pml_old );
                }
            }
            //Electric field Et^(p,p) Remind that in PML, there no current
            for( unsigned int i=0 ; i<nl_p ; i++ ) {
                for( unsigned int j=solvermin ; j<(unsigned int)solvermax ; j++ ) {
                    // Standard FDTD
                    // ( *Dt_pml )( i, j ) = + 1 * ( *Dt_pml )( i, j )
                    //                     - dt/dr * ( ( *Hl_pml )( i, j+1 ) - ( *Hl_pml )( i, j ) )
                    //                     + dt/dl * ( ( *Hr_pml )( i+1, j ) - ( *Hr_pml )( i, j ) ) ;
                    // ( *Et_pml )( i, j ) = 1*( *Dt_pml )( i, j );
                    // PML FDTD
                    Dt_pml_old = std::complex<double>(1,0)*( *Dt_pml )( i, j ) ;
                    ( *Dt_pml )( i, j ) = + c1_p_tfield[j] * ( *Dt_pml )( i, j )
                                          - c2_p_tfield[j]/dr * ( ( *Hl_pml )( i, j+1 ) - ( *Hl_pml )( i, j ) )
                                          + c2_p_tfield[j]/dl * ( ( *Hr_pml )( i+1, j ) - ( *Hr_pml )( i, j ) ) ;
                    ( *Et_pml )( i, j ) = + c3_p_tfield[i] * ( *Et_pml )( i, j )
                                          + c4_p_tfield[i] * (c5_p_tfield[j]*( *Dt_pml )( i, j ) - c6_p_tfield[j]*Dt_pml_old );
                }
            }
        }
    }
}

void PML_SolverAM::compute_H_from_B( ElectroMagn *fields, int iDim, int min_or_max, std::vector<unsigned int> dimPrim, unsigned int solvermin, unsigned int solvermax )
{
    const unsigned int nl_p = dimPrim[0];
    const unsigned int nl_d = dimPrim[0] + 1;
    const unsigned int nr_p = dimPrim[1];
    const unsigned int nr_d = dimPrim[1] + 1;
    
    ElectroMagnBCAM_PML* pml_fields = static_cast<ElectroMagnBCAM_PML*>( fields->emBoundCond[iDim*2+min_or_max] );
    cField2D* El_pml = NULL;
    cField2D* Er_pml = NULL;
    cField2D* Et_pml = NULL;
    cField2D* Hl_pml = NULL;
    cField2D* Hr_pml = NULL;
    cField2D* Ht_pml = NULL;
    cField2D* Bl_pml = NULL;
    cField2D* Br_pml = NULL;
    cField2D* Bt_pml = NULL;

    bool isYmin = ( static_cast<ElectroMagnAM *>( fields ) )->isYmin;

    if (min_or_max==0){
        isMin=true;
        isMax=false;
    }
    else if (min_or_max==1){
        isMin=false;
        isMax=true;
    }

    if (iDim==0){
        for( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
            El_pml = pml_fields->El_[imode];
            Er_pml = pml_fields->Er_[imode];
            Et_pml = pml_fields->Et_[imode];
            Hl_pml = pml_fields->Hl_[imode];
            Hr_pml = pml_fields->Hr_[imode];
            Ht_pml = pml_fields->Ht_[imode];
            Bl_pml = pml_fields->Bl_[imode];
            Br_pml = pml_fields->Br_[imode];
            Bt_pml = pml_fields->Bt_[imode];
            //Magnetic field Bl^(p,d)
            for( unsigned int i=solvermin ; i<(unsigned int)solvermax;  i++ ) {
                for( unsigned int j=1+isYmin*2 ; j<nr_d-1 ; j++ ) {
                    // Standard FDTD
                    // ( *Bl_pml )( i, j ) = + 1 * ( *Bl_pml )( i, j )
                    //                       - dt / ( ( j_glob_pml+j-0.5 )*dr ) * ( ( double )( j+j_glob_pml )*( *Et_pml )( i, j ) - ( double )( j+j_glob_pml-1. )*( *Et_pml )( i, j-1 ) )
                    //                       - dt / ( ( j_glob_pml+j-0.5 )*dr ) * ( Icpx*( double )imode*( *Er_pml )( i, j ) ) ;
                    // ( *Hl_pml )( i, j ) = 1*( *Bl_pml )( i, j );
                    // PML FDTD
                    Bl_pml_old = std::complex<double>(1,0)*( *Bl_pml )( i, j ) ;
                    ( *Bl_pml )( i, j ) = + c1_d_lfield[j] * ( *Bl_pml )( i, j )
                                          - c2_d_lfield[j] / ( ( j_glob_pml+j-0.5 )*dr ) * ( ( double )( j+j_glob_pml )*( *Et_pml )( i, j ) - ( double )( j+j_glob_pml-1. )*( *Et_pml )( i, j-1 ) )
                                          - c2_d_lfield[j] / ( ( j_glob_pml+j-0.5 )*dr ) * ( Icpx*( double )imode*( *Er_pml )( i, j ) ) ;
                    ( *Hl_pml )( i, j ) = + c3_d_lfield[j] * ( *Hl_pml )( i, j )
                                          + c4_d_lfield[j] * (c5_p_lfield[i]*( *Bl_pml )( i, j ) - c6_p_lfield[i]*Bl_pml_old ) ;
                }
            }
            //Magnetic field Br^(d,p)
            for( unsigned int i=solvermin ; i<(unsigned int)solvermax ; i++ ) {
                for( unsigned int j=isYmin*3 ; j<nr_p ; j++ ) {
                    //Standard FDTD
                    // ( *Br_pml )( i, j ) = + 1 * ( *Br_pml )( i, j )
                    //                       + dt/dl * ( ( *Et_pml )( i, j ) - ( *Et_pml )( i-1, j ) )
                    //                       + dt*Icpx*( double )imode/( ( double )( j_glob_pml+j )*dr )*( *El_pml )( i, j ) ;
                    // ( *Hr_pml )( i, j ) = 1*( *Br_pml )( i, j );
                    // PML FDTD
                    Br_pml_old = std::complex<double>(1,0)*( *Br_pml )( i, j ) ;
                    ( *Br_pml )( i, j ) = + c1_d_rfield[i] * ( *Br_pml )( i, j )
                                          + c2_d_rfield[i]/dl * ( ( *Et_pml )( i, j ) - ( *Et_pml )( i-1, j ) )
                                          + c2_d_rfield[i]*Icpx*( double )imode/( ( double )( j_glob_pml+j )*dr )*( *El_pml )( i, j ) ;
                    ( *Hr_pml )( i, j ) = + c3_p_rfield[j] * ( *Hr_pml )( i, j )
                                          + c4_p_rfield[j] * (c5_p_rfield[j]*( *Br_pml )( i, j ) - c6_p_rfield[j]*Br_pml_old );
                }
            }
            //Magnetic field Bt^(d,d)
            for( unsigned int i=solvermin ; i<(unsigned int)solvermax ; i++ ) {
                for( unsigned int j=1+isYmin*2 ; j<nr_d-1 ; j++ ) {
                    // Standard FDTD
                    // ( *Bt_pml )( i, j ) = + 1 * ( *Bt_pml )( i, j )
                    //                       + dt/dr * ( ( *El_pml )( i, j ) - ( *El_pml )( i, j-1 ) )
                    //                       - dt/dl * ( ( *Er_pml )( i, j ) - ( *Er_pml )( i-1, j ) ) ;
                    // ( *Ht_pml )( i, j ) = 1*( *Bt_pml )( i, j );
                    // PML FDTD
                    Bt_pml_old = std::complex<double>(1,0)*( *Bt_pml )( i, j ) ;
                    ( *Bt_pml )( i, j ) = + c1_d_tfield[j] * ( *Bt_pml )( i, j )
                                          + c2_d_tfield[j]/dr * ( ( *El_pml )( i, j ) - ( *El_pml )( i, j-1 ) )
                                          - c2_d_tfield[j]/dl * ( ( *Er_pml )( i, j ) - ( *Er_pml )( i-1, j ) ) ;
                    ( *Ht_pml )( i, j ) = + c3_d_tfield[i] * ( *Ht_pml )( i, j )
                                          + c4_d_tfield[i] * (c5_d_tfield[j]*( *Bt_pml )( i, j ) - c6_d_tfield[j]*Bt_pml_old );
                }
            }
            if( isYmin ) {
                unsigned int j=2;
                if( imode==0 ) {
                    for( unsigned int i=0 ; i<nl_d ; i++ ) {
                        ( *Br_pml )( i, j )=0;
                        ( *Br_pml )( i, 1 )=-( *Br_pml )( i, 3 );
                        //
                        ( *Hr_pml )( i, j ) = ( *Br_pml )( i, j ) ;
                        ( *Hr_pml )( i, 1 ) = ( *Br_pml )( i, 1 ) ;
                    }
                    for( unsigned int i=0 ; i<nl_d ; i++ ) {
                        ( *Bt_pml )( i, j )= -( *Bt_pml )( i, j+1 );
                        //
                        ( *Ht_pml )( i, j ) = ( *Bt_pml )( i, j ) ;
                    }
                    for( unsigned int i=0 ; i<nl_p ; i++ ) {
                        ( *Bl_pml )( i, j )= ( *Bl_pml )( i, j+1 );
                        //
                        ( *Hl_pml )( i, j ) = ( *Bl_pml )( i, j ) ;
                    }
                }
                else if( imode==1 ) {
                    for( unsigned int i=0 ; i<nl_p  ; i++ ) {
                        ( *Bl_pml )( i, j )= -( *Bl_pml )( i, j+1 );
                        //
                        ( *Hl_pml )( i, j ) = ( *Bl_pml )( i, j ) ;
                    }

                    for( unsigned int i=1 ; i<nl_d-1 ; i++ ) {
                        ( *Br_pml )( i, j )+=  Icpx*dt_ov_dr*( *El_pml )( i, j+1 )
                                        +                       dt_ov_dl*( ( *Et_pml )( i, j )-( *Et_pml )( i-1, j ) );
                        ( *Br_pml )( i, 1 )=( *Br_pml )( i, 3 );
                        //
                        ( *Hr_pml )( i, j ) = ( *Br_pml )( i, j ) ;
                        ( *Hr_pml )( i, 1 ) = ( *Br_pml )( i, 1 ) ;
                    }
                    for( unsigned int i=0; i<nl_d ; i++ ) {
                        ( *Bt_pml )( i, j )= -2.*Icpx*( *Br_pml )( i, j )-( *Bt_pml )( i, j+1 );
                        //
                        ( *Ht_pml )( i, j ) = ( *Bt_pml )( i, j ) ;
                    }
                }
                else { // modes > 1
                    for( unsigned int  i=0 ; i<nl_p; i++ ) {
                        ( *Bl_pml )( i, j )= -( *Bl_pml )( i, j+1 );
                        //
                        ( *Hl_pml )( i, j ) = ( *Bl_pml )( i, j ) ;
                    }
                    for( unsigned int i=0 ; i<nl_d; i++ ) {
                        ( *Br_pml )( i, j )= 0;
                        ( *Br_pml )( i, 1 )=-( *Br_pml )( i, 3 );
                        //
                        ( *Hr_pml )( i, j ) = ( *Br_pml )( i, j ) ;
                        ( *Hr_pml )( i, 1 ) = ( *Br_pml )( i, 1 ) ;
                    }
                    for( unsigned int  i=0 ; i<nl_d ; i++ ) {
                        ( *Bt_pml )( i, j )= - ( *Bt_pml )( i, j+1 );
                        //
                        ( *Ht_pml )( i, j ) = ( *Bt_pml )( i, j ) ;
                    }
                }
            }
            // if( isYmin ) {
            //     unsigned int j=2;
            //     if( imode==0 ) {
            //         for( unsigned int i=3*isMax ; i<nl_d-3*isMin ; i++ ) {
            //             ( *Br_pml )( i, j ) = 0;
            //             ( *Hr_pml )( i, j ) = 0;

            //             Br_pml_old = ( *Br_pml )( i, 1 );
            //             ( *Br_pml )( i, 1 ) = -( *Br_pml )( i, 3 );
            //             ( *Hr_pml )( i, 1 ) = + c3_p_rfield[j] * ( *Hr_pml )( i, 1 )
            //                                   + c4_p_rfield[j] * (c5_p_rfield[j]*( *Br_pml )( i, 1 ) - c6_p_rfield[j]*Br_pml_old );
            //         }
            //         for( unsigned int i=3*isMax ; i<nl_d-3*isMin ; i++ ) {
            //             Bt_pml_old = ( *Bt_pml )( i, j );
            //             ( *Bt_pml )( i, j ) = -( *Bt_pml )( i, j+1 );
            //             ( *Ht_pml )( i, j ) = + c3_d_tfield[i] * ( *Ht_pml )( i, j )
            //                                   + c4_d_tfield[i] * (c5_d_tfield[j]*( *Bt_pml )( i, j ) - c6_d_tfield[j]*Bt_pml_old );
            //         }
            //         for( unsigned int i=3*isMax ; i<nl_p-3*isMin ; i++ ) {
            //             Bl_pml_old = ( *Bl_pml )( i, j );
            //             ( *Bl_pml )( i, j ) = ( *Bl_pml )( i, j+1 );
            //             ( *Hl_pml )( i, j ) = + c3_d_lfield[j] * ( *Hl_pml )( i, j )
            //                                   + c4_d_lfield[j] * (c5_p_lfield[i]*( *Bl_pml )( i, j ) - c6_p_lfield[i]*Bl_pml_old ) ;
            //         }
            //     }
            //     else if( imode==1 ) {
            //         for( unsigned int i=3*isMax ; i<nl_p-3*isMin  ; i++ ) {
            //             Bl_pml_old = ( *Bl_pml )( i, j );
            //             ( *Bl_pml )( i, j ) = -( *Bl_pml )( i, j+1 );
            //             ( *Hl_pml )( i, j ) = + c3_d_lfield[j] * ( *Hl_pml )( i, j )
            //                                   + c4_d_lfield[j] * (c5_p_lfield[i]*( *Bl_pml )( i, j ) - c6_p_lfield[i]*Bl_pml_old ) ;
            //         }
            //         for( unsigned int i=3*isMax+1*isMin ; i<nl_d-1*isMax-3*isMin ; i++ ) {
            //             Br_pml_old          = ( *Br_pml )( i, j ) ;
            //             ( *Br_pml )( i, j ) = + c1_d_rfield[i] * ( *Br_pml )( i, j )
            //                                   + c2_d_rfield[i]/dl * ( ( *Et_pml )( i, j ) - ( *Et_pml )( i-1, j ) )
            //                                   + c2_d_rfield[i]/dr * Icpx * ( *El_pml )( i, j+1 ) ;
            //             ( *Hr_pml )( i, j ) = + c3_p_rfield[j] * ( *Hr_pml )( i, j )
            //                                   + c4_p_rfield[j] * (c5_p_rfield[j]*( *Br_pml )( i, j ) - c6_p_rfield[j]*Br_pml_old );

            //             Br_pml_old          = ( *Br_pml )( i, 1 );
            //             ( *Br_pml )( i, 1 ) = ( *Br_pml )( i, 3 );
            //             ( *Hr_pml )( i, 1 ) = + c3_p_rfield[j] * ( *Hr_pml )( i, 1 )
            //                                   + c4_p_rfield[j] * (c5_p_rfield[j]*( *Br_pml )( i, 1 ) - c6_p_rfield[j]*Br_pml_old );

            //         }
            //         for( unsigned int i=3*isMax ; i<nl_d-3*isMin ; i++ ) {
            //             Bt_pml_old = ( *Bt_pml )( i, j );
            //             ( *Bt_pml )( i, j ) = -2.*Icpx*( *Br_pml )( i, j )-( *Bt_pml )( i, j+1 );
            //             ( *Ht_pml )( i, j ) = + c3_d_tfield[i] * ( *Ht_pml )( i, j )
            //                                   + c4_d_tfield[i] * (c5_d_tfield[j]*( *Bt_pml )( i, j ) - c6_d_tfield[j]*Bt_pml_old );
            //         }
            //     } else { // modes > 1
            //         for( unsigned int i=3*isMax ; i<nl_p-3*isMin; i++ ) {
            //             Bl_pml_old = ( *Bl_pml )( i, j );
            //             ( *Bl_pml )( i, j ) = -( *Bl_pml )( i, j+1 );
            //             ( *Hl_pml )( i, j ) = + c3_d_lfield[j] * ( *Hl_pml )( i, j )
            //                                   + c4_d_lfield[j] * (c5_p_lfield[i]*( *Bl_pml )( i, j ) - c6_p_lfield[i]*Bl_pml_old ) ;
            //         }
            //         for( unsigned int i=3*isMax ; i<nl_d-3*isMin; i++ ) {
            //             ( *Br_pml )( i, j ) = 0;
            //             ( *Hr_pml )( i, j ) = 0;

            //             Br_pml_old = ( *Br_pml )( i, 1 );
            //             ( *Br_pml )( i, 1 ) = -( *Br_pml )( i, 3 );
            //             ( *Hr_pml )( i, 1 ) = + c3_p_rfield[j] * ( *Hr_pml )( i, 1 )
            //                                   + c4_p_rfield[j] * (c5_p_rfield[j]*( *Br_pml )( i, 1 ) - c6_p_rfield[j]*Br_pml_old );
            //         }
            //         for( unsigned int  i=3*isMax ; i<nl_d-3*isMin ; i++ ) {
            //             Bt_pml_old = ( *Bt_pml )( i, j ) ;
            //             ( *Bt_pml )( i, j ) = - ( *Bt_pml )( i, j+1 );
            //             ( *Ht_pml )( i, j ) = + c3_d_tfield[i] * ( *Ht_pml )( i, j )
            //                                   + c4_d_tfield[i] * (c5_d_tfield[j]*( *Bt_pml )( i, j ) - c6_d_tfield[j]*Bt_pml_old );
            //         }
            //     }
            // }
        }
    }
    else if (iDim==1){
        for( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
            El_pml = pml_fields->El_[imode];
            Er_pml = pml_fields->Er_[imode];
            Et_pml = pml_fields->Et_[imode];
            Hl_pml = pml_fields->Hl_[imode];
            Hr_pml = pml_fields->Hr_[imode];
            Ht_pml = pml_fields->Ht_[imode];
            Bl_pml = pml_fields->Bl_[imode];
            Br_pml = pml_fields->Br_[imode];
            Bt_pml = pml_fields->Bt_[imode];
            //Magnetic field Bl^(p,d)
            for( unsigned int i=0 ; i<nl_p;  i++ ) {
                for( unsigned int j=solvermin ; j<(unsigned int)solvermax ; j++ ) {
                    // Standard FDTD
                    // ( *Bl_pml )( i, j ) = + 1 * ( *Bl_pml )( i, j )
                    //                       - dt / ( ( j_glob_pml+j-0.5 )*dr ) * ( ( double )( j+j_glob_pml )*( *Et_pml )( i, j ) - ( double )( j+j_glob_pml-1. )*( *Et_pml )( i, j-1 ) )
                    //                       - dt / ( ( j_glob_pml+j-0.5 )*dr ) * ( Icpx*( double )imode*( *Er_pml )( i, j ) ) ;
                    // ( *Hl_pml )( i, j ) = 1*( *Bl_pml )( i, j );
                    // PML FDTD
                    Bl_pml_old = std::complex<double>(1,0)*( *Bl_pml )( i, j ) ;
                    ( *Bl_pml )( i, j ) = + c1_d_lfield[j] * ( *Bl_pml )( i, j )
                                          - c2_d_lfield[j] / ( ( j_glob_pml+j-0.5 )*dr ) * ( ( double )( j+j_glob_pml )*( *Et_pml )( i, j ) - ( double )( j+j_glob_pml-1. )*( *Et_pml )( i, j-1 ) )
                                          - c2_d_lfield[j] / ( ( j_glob_pml+j-0.5 )*dr ) * ( Icpx*( double )imode*( *Er_pml )( i, j ) ) ;
                    ( *Hl_pml )( i, j ) = + c3_d_lfield[j] * ( *Hl_pml )( i, j )
                                          + c4_d_lfield[j] * (c5_p_lfield[i]*( *Bl_pml )( i, j ) - c6_p_lfield[i]*Bl_pml_old ) ;
                }
            }
            //Magnetic field Br^(d,p)
            for( unsigned int i=1 ; i<nl_d-1 ; i++ ) {
                for( unsigned int j=solvermin ; j<(unsigned int)solvermax ; j++ ) {
                    //Standard FDTD
                    // ( *Br_pml )( i, j ) = + 1 * ( *Br_pml )( i, j )
                    //                       + dt/dl * ( ( *Et_pml )( i, j ) - ( *Et_pml )( i-1, j ) )
                    //                       + dt*Icpx*( double )imode/( ( double )( j_glob_pml+j )*dr )*( *El_pml )( i, j ) ;
                    // ( *Hr_pml )( i, j ) = 1*( *Br_pml )( i, j );
                    // PML FDTD
                    Br_pml_old = std::complex<double>(1,0)*( *Br_pml )( i, j ) ;
                    ( *Br_pml )( i, j ) = + c1_d_rfield[i] * ( *Br_pml )( i, j )
                                          + c2_d_rfield[i]/dl * ( ( *Et_pml )( i, j ) - ( *Et_pml )( i-1, j ) )
                                          + c2_d_rfield[i]*Icpx*( double )imode/( ( double )( j_glob_pml+j )*dr )*( *El_pml )( i, j ) ;
                    ( *Hr_pml )( i, j ) = + c3_p_rfield[j] * ( *Hr_pml )( i, j )
                                          + c4_p_rfield[j] * (c5_p_rfield[j]*( *Br_pml )( i, j ) - c6_p_rfield[j]*Br_pml_old );
                }
            }
            //Magnetic field Bt^(d,d)
            for( unsigned int i=1 ; i<nl_d-1 ; i++ ) {
                for( unsigned int j=solvermin ; j<(unsigned int)solvermax ; j++ ) {
                    // Standard FDTD
                    // ( *Bt_pml )( i, j ) = + 1 * ( *Bt_pml )( i, j )
                    //                       + dt/dr * ( ( *El_pml )( i, j ) - ( *El_pml )( i, j-1 ) )
                    //                       - dt/dl * ( ( *Er_pml )( i, j ) - ( *Er_pml )( i-1, j ) ) ;
                    // ( *Ht_pml )( i, j ) = 1*( *Bt_pml )( i, j );
                    // PML FDTD
                    Bt_pml_old = std::complex<double>(1,0)*( *Bt_pml )( i, j ) ;
                    ( *Bt_pml )( i, j ) = + c1_d_tfield[j] * ( *Bt_pml )( i, j )
                                          + c2_d_tfield[j]/dr * ( ( *El_pml )( i, j ) - ( *El_pml )( i, j-1 ) )
                                          - c2_d_tfield[j]/dl * ( ( *Er_pml )( i, j ) - ( *Er_pml )( i-1, j ) ) ;
                    ( *Ht_pml )( i, j ) = + c3_d_tfield[i] * ( *Ht_pml )( i, j )
                                          + c4_d_tfield[i] * (c5_d_tfield[j]*( *Bt_pml )( i, j ) - c6_d_tfield[j]*Bt_pml_old );
                }
            }
        }
    }
}
