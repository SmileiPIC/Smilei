#include "PML_Solver2D_Yee.h"
#include "ElectroMagn.h"
#include "ElectroMagnBC2D_PML.h"
#include "Field2D.h"
#include "Patch.h"

PML_Solver2D_Yee::PML_Solver2D_Yee( Params &params ):
    Solver2D( params ),
    pml_sigma_( 2, NULL ),
    pml_kappa_( 2, NULL )
{

    std::vector<PyObject *> prof;
    if( PyTools::extract_pyProfiles( "pml_sigma", "Main", 0, prof )){
        if( prof.size() == 0 ){
            ERROR(" in pml_sigma, expecting a list of at least 1 profile.");
        }
    // extracted profile // number of variables of the function // name of the profile extracted // params // try numpy ?? // try file ?? // time variable ??
        pml_sigma_[0] = new Profile( prof[0], 1, "pml_sigma_x_profile", params, true, false, false );
        if( prof.size() == 1){ 
            pml_sigma_[1] = new Profile( prof[0], 1, "pml_sigma_y_profile", params, true, false, false );
        } else {
            pml_sigma_[1] = new Profile( prof[1], 1, "pml_sigma_y_profile", params, true, false, false );
        }
    }
    if( PyTools::extract_pyProfiles( "pml_kappa", "Main", 0, prof )){
        if(prof.size() == 0){
            ERROR(" in pml_kappa, expecting a list of at least 1 profile.");
        }
        pml_kappa_[0] = new Profile( prof[0], 1, "pml_kappa_x_profile", params, true, false, false );
        if( prof.size() == 1){ 
            pml_kappa_[1] = new Profile( prof[0], 1, "pml_kappa_y_profile", params, true, false, false );
        } else {
            pml_kappa_[1] = new Profile( prof[1], 1, "pml_kappa_y_profile", params, true, false, false );
        }
    }
}

PML_Solver2D_Yee::~PML_Solver2D_Yee()
{
    for( unsigned int i=0; i<pml_sigma_.size(); i++ ) {
        delete pml_sigma_[i];
    }
    for( unsigned int i=0; i<pml_kappa_.size(); i++ ) {
        delete pml_kappa_[i];
    }
}

void PML_Solver2D_Yee::operator()( ElectroMagn * )
{
    ERROR( "This is not a solver for the main domain" );

}

void PML_Solver2D_Yee::setDomainSizeAndCoefficients( int iDim, int min_or_max, std::vector<unsigned int> dimPrim, int ncells_pml_domain, int startpml, int* ncells_pml_min, int* ncells_pml_max, Patch* )
{
    const unsigned int nx_p = dimPrim[0];
    const unsigned int nx_d = dimPrim[0] + 1;
    const unsigned int ny_p = dimPrim[1];
    const unsigned int ny_d = dimPrim[1] + 1;
    
    //PML Coeffs Kappa,Sigma ...
    //Primal
    kappa_x_p.resize( nx_p );
    sigma_x_p.resize( nx_p );
    kappa_y_p.resize( ny_p );
    sigma_y_p.resize( ny_p );
    kappa_z_p.resize( 1 );
    sigma_z_p.resize( 1 );
    //Dual
    kappa_x_d.resize( nx_d );
    sigma_x_d.resize( nx_d );
    kappa_y_d.resize( ny_d );
    sigma_y_d.resize( ny_d );
    kappa_z_d.resize( 1 );
    sigma_z_d.resize( 1 );

    // While min and max of each dim have got the same ncells_pml, coefficient are the same
    // for min and max PML

    //Coeffs for El,Dl,Hl,Bl
    //Primal
    c1_p_xfield.resize( ny_p ); // j-dependent
    c2_p_xfield.resize( ny_p ); // j-dependent
    c3_p_xfield.resize( 1 ); // k-dependent
    c4_p_xfield.resize( 1 ); // k-dependent
    c5_p_xfield.resize( nx_p ); // i-dependent
    c6_p_xfield.resize( nx_p ); // i-dependent
    //Dual
    c1_d_xfield.resize( ny_d ); // j-dependent
    c2_d_xfield.resize( ny_d ); // j-dependent
    c3_d_xfield.resize( 1 ); // j-dependent
    c4_d_xfield.resize( 1 ); // j-dependent
    c5_d_xfield.resize( nx_d ); // i-dependent
    c6_d_xfield.resize( nx_d ); // i-dependent
    //Coefs for Er,Dr,Hr,Br
    //Primal
    c1_p_yfield.resize( 1 ); // k-dependent
    c2_p_yfield.resize( 1 ); // k-dependent
    c3_p_yfield.resize( nx_p ); // i-dependent
    c4_p_yfield.resize( nx_p ); // i-dependent
    c5_p_yfield.resize( ny_p ); // j-dependent
    c6_p_yfield.resize( ny_p ); // j-dependent
    //Dual
    c1_d_yfield.resize( 1 ); // k-dependent
    c2_d_yfield.resize( 1 ); // k-dependent
    c3_d_yfield.resize( nx_d ); // i-dependent
    c4_d_yfield.resize( nx_d ); // i-dependent
    c5_d_yfield.resize( ny_d ); // k-dependent
    c6_d_yfield.resize( ny_d ); // k-dependent
    //Coefs for Et,Dt,Ht,Bt
    //Primal
    c1_p_zfield.resize( nx_p ); // j-dependent
    c2_p_zfield.resize( nx_p ); // j-dependent
    c3_p_zfield.resize( ny_p ); // i-dependent
    c4_p_zfield.resize( ny_p ); // i-dependent
    c5_p_zfield.resize( 1 ); // j-dependent
    c6_p_zfield.resize( 1 ); // j-dependent
    //Dual
    c1_d_zfield.resize( nx_d ); // i-dependent
    c2_d_zfield.resize( nx_d ); // i-dependent
    c3_d_zfield.resize( ny_d ); // j-dependent
    c4_d_zfield.resize( ny_d ); // j-dependent
    c5_d_zfield.resize( 1 ); // k-dependent
    c6_d_zfield.resize( 1 ); // k-dependent

    if ( iDim == 0 ) {
        // 3 cells (oversize) are vaccum so the PML media begin at r0 which is :
        // Eventually the size of PML media is :
        length_x_pml = (ncells_pml_domain-startpml+0.5)*dx ;
        // Primal grid
        // X-direction
        // Params for first cell of PML-patch (vacuum) i = 0,1,2
        for ( int i=0 ; i<startpml ; i++ ) {
            // Coeffs for the first cell
            kappa_x_p[i] = 1. ;
            sigma_x_p[i] = 0. ;
        }
        // Params for other cells (PML Media) when i>=3
        for( int i = startpml; i< (int) nx_p ; i++ ) {
            kappa_x_p[i] = pml_kappa_[0]->valueAt((i-startpml)*dx/length_x_pml);
            sigma_x_p[i] = pml_sigma_[0]->valueAt((i-startpml)*dx/length_x_pml);
        }
        // Y-direction
        for( unsigned int j = 0 ; j<ny_p ; j++ ) {
            kappa_y_p[j] = 1. ;
            sigma_y_p[j] = 0. ;
        }
        // Z-direction
        for ( int k=0 ;k<1 ; k++ ) {
            kappa_z_p[k] = 1. ;
            sigma_z_p[k] = 0. ;
        }
        // Dual grid
        // X-direction
        // Params for first cell of PML-patch (vacuum) i = 0,1,2,3
        for ( int i=0 ; i<startpml+1 ; i++ ) {
            // Coeffs for the first cell
            kappa_x_d[i] = 1. ;
            sigma_x_d[i] = 0. ;
        }
        // Params for other cells (PML Media) when j>=4
        for( int i = startpml+1 ; i< (int) nx_d ; i++ ) {
            kappa_x_d[i] = pml_kappa_[0]->valueAt((i-startpml-0.5)*dx/length_x_pml);
            sigma_x_d[i] = pml_sigma_[0]->valueAt((i-startpml-0.5)*dx/length_x_pml);
        }
        // Y-direction
        for( unsigned int j = 0 ; j<ny_d ; j++ ) {
            kappa_y_d[j] = 1. ;
            sigma_y_d[j] = 0. ;
        }
        // Z-direction
        for ( int k=0 ; k<1 ; k++ ) {
            kappa_z_d[k] = 1. ;
            sigma_z_d[k] = 0. ;
        }
    }

    if ( iDim == 1 ) {
        // 3 cells are vaccum so the PML media begin at r0 which is :
        // Eventually the size of PML media is :
        length_y_pml = (ncells_pml_domain-startpml+0.5)*dy ;
        length_x_pml_xmax = (ncells_pml_max[0]+0.5)*dx ;
        length_x_pml_xmin = (ncells_pml_min[0]+0.5)*dx ;
        // Primal grid
        // X-direction
        for( unsigned int i = 0 ; i<nx_p ; i++ ) {
            kappa_x_p[i] = 1. ;
            sigma_x_p[i] = 0. ;
        }
        if (ncells_pml_min[0] != 0 ){
            for ( int i=0 ; i<ncells_pml_min[0] ; i++ ) {
                kappa_x_p[i] = pml_kappa_[0]->valueAt((ncells_pml_min[0] - 1 - i)*dx/length_x_pml_xmin);
                sigma_x_p[i] = pml_sigma_[0]->valueAt((ncells_pml_min[0] - 1 - i)*dx/length_x_pml_xmin);
            }
        }
        if (ncells_pml_max[0] != 0 ){
            for( int i = (nx_p-1)-(ncells_pml_max[0]-1) ; i< (int) nx_p ; i++ ) {
                kappa_x_p[i] = pml_kappa_[0]->valueAt((i - nx_p  + ncells_pml_max[0])*dx/length_x_pml_xmax);
                sigma_x_p[i] = pml_sigma_[0]->valueAt((i - nx_p  + ncells_pml_max[0])*dx/length_x_pml_xmax);
            }
        }
        // Y-direction
        // Params for first cell of PML-patch (vacuum) i = 0,1,2
        for ( int j=0 ; j<startpml ; j++ ) {
            // Coeffs for the first cell
            kappa_y_p[j] = 1. ;
            sigma_y_p[j] = 0. ;
        }
        // Params for other cells (PML Media) when j>=3
        for( int j = startpml ; j< (int) ny_p ; j++ ) {
            kappa_y_p[j] = pml_kappa_[1]->valueAt((j-startpml)*dy/length_y_pml);
            sigma_y_p[j] = pml_sigma_[1]->valueAt((j-startpml)*dy/length_y_pml);
        }
        // Z-direction
        for ( int k=0 ;k<1 ; k++ ) {
            kappa_z_p[k] = 1. ;
            sigma_z_p[k] = 0. ;
        }
        // Dual grid
        // X-direction
        for( unsigned int i = 0 ; i<nx_d ; i++ ) {
            kappa_x_d[i] = 1. ;
            sigma_x_d[i] = 0. ;
        }
        if (ncells_pml_min[0] != 0 ){
            for ( int i=0 ; i<ncells_pml_min[0] ; i++ ) {
                kappa_x_d[i] = pml_kappa_[0]->valueAt(( ncells_pml_min[0] - 0.5 - i )*dx/length_x_pml_xmin);
                sigma_x_d[i] = pml_sigma_[0]->valueAt(( ncells_pml_min[0] - 0.5 - i )*dx/length_x_pml_xmin);
            }
        }
        if (ncells_pml_max[0] != 0 ){
            for( int i = (nx_p-1)-(ncells_pml_max[0]-1)+1 ; i< (int) nx_d ; i++ ) {
                kappa_x_d[i] = pml_kappa_[0]->valueAt((i - nx_p + ncells_pml_max[0] - 0.5 )*dx/length_x_pml_xmax);
                sigma_x_d[i] = pml_sigma_[0]->valueAt((i - nx_p  + ncells_pml_max[0] - 0.5)*dx/length_x_pml_xmax);
            }
        }
        // Y-direction
        // Params for first cell of PML-patch (vacuum) i = 0,1,2,3
        for ( int j=0 ; j<startpml+1 ; j++ ) {
            // Coeffs for the first cell
            kappa_y_d[j] = 1. ;
            sigma_y_d[j] = 0. ;
        }
        // Params for other cells (PML Media) when j>=4
        for( int j = startpml+1 ; j< (int) ny_d ; j++ ) {
            kappa_y_d[j] = pml_kappa_[1]->valueAt((j-startpml-0.5)*dy/length_y_pml);
            sigma_y_d[j] = pml_sigma_[1]->valueAt((j-startpml-0.5)*dy/length_y_pml);
        }
        // Z-direction
        for ( int k=0 ; k<1 ; k++ ) {
            kappa_z_d[k] = 1. ;
            sigma_z_d[k] = 0. ;
        }
    }

    if ((min_or_max==0)&&(iDim==0)){
        for( int i = 0 ; i< (int) nx_p ; i++ ) {
            c1_p_zfield[i] = ( 2.*kappa_x_p[(nx_p-1)-i] - dt*sigma_x_p[(nx_p-1)-i] ) / ( 2.*kappa_x_p[(nx_p-1)-i] + dt*sigma_x_p[(nx_p-1)-i] ) ;
            c2_p_zfield[i] = ( 2*dt ) / ( 2.*kappa_x_p[(nx_p-1)-i] + dt*sigma_x_p[(nx_p-1)-i] ) ;
            c3_p_yfield[i] = ( 2.*kappa_x_p[(nx_p-1)-i] - dt*sigma_x_p[(nx_p-1)-i] ) / ( 2.*kappa_x_p[(nx_p-1)-i] + dt*sigma_x_p[(nx_p-1)-i] ) ;
            c4_p_yfield[i] = ( 1. ) / ( 2.*kappa_x_p[(nx_p-1)-i] + dt*sigma_x_p[(nx_p-1)-i] ) ;
            c5_p_xfield[i] = ( 2.*kappa_x_p[(nx_p-1)-i] + dt*sigma_x_p[(nx_p-1)-i] ) ;
            c6_p_xfield[i] = ( 2.*kappa_x_p[(nx_p-1)-i] - dt*sigma_x_p[(nx_p-1)-i] ) ;
        }

        for( int i = 0 ; i< (int) nx_d ; i++ ) {
            c1_d_zfield[i] = ( 2.*kappa_x_d[(nx_d-1)-i] - dt*sigma_x_d[(nx_d-1)-i] ) / ( 2.*kappa_x_d[(nx_d-1)-i] + dt*sigma_x_d[(nx_d-1)-i] ) ;
            c2_d_zfield[i] = ( 2*dt ) / ( 2.*kappa_x_d[(nx_d-1)-i] + dt*sigma_x_d[(nx_d-1)-i] ) ;
            c3_d_yfield[i] = ( 2.*kappa_x_d[(nx_d-1)-i] - dt*sigma_x_d[(nx_d-1)-i] ) / ( 2.*kappa_x_d[(nx_d-1)-i] + dt*sigma_x_d[(nx_d-1)-i] ) ;
            c4_d_yfield[i] = ( 1. ) / ( 2.*kappa_x_d[(nx_d-1)-i] + dt*sigma_x_d[(nx_d-1)-i] ) ;
            c5_d_xfield[i] = ( 2.*kappa_x_d[(nx_d-1)-i] + dt*sigma_x_d[(nx_d-1)-i] ) ;
            c6_d_xfield[i] = ( 2.*kappa_x_d[(nx_d-1)-i] - dt*sigma_x_d[(nx_d-1)-i] ) ;
        }
    }
    else {
        for( int i = 0 ; i< (int) nx_p ; i++ ) {
            c1_p_zfield[i] = ( 2.*kappa_x_p[i] - dt*sigma_x_p[i] ) / ( 2.*kappa_x_p[i] + dt*sigma_x_p[i] ) ;
            c2_p_zfield[i] = ( 2*dt ) / ( 2.*kappa_x_p[i] + dt*sigma_x_p[i] ) ;
            c3_p_yfield[i] = ( 2.*kappa_x_p[i] - dt*sigma_x_p[i] ) / ( 2.*kappa_x_p[i] + dt*sigma_x_p[i] ) ;
            c4_p_yfield[i] = ( 1. ) / ( 2.*kappa_x_p[i] + dt*sigma_x_p[i] ) ;
            c5_p_xfield[i] = ( 2.*kappa_x_p[i] + dt*sigma_x_p[i] ) ;
            c6_p_xfield[i] = ( 2.*kappa_x_p[i] - dt*sigma_x_p[i] ) ;
        }

        for( int i = 0 ; i< (int) nx_d ; i++ ) {
            c1_d_zfield[i] = ( 2.*kappa_x_d[i] - dt*sigma_x_d[i] ) / ( 2.*kappa_x_d[i] + dt*sigma_x_d[i] ) ;
            c2_d_zfield[i] = ( 2*dt ) / ( 2.*kappa_x_d[i] + dt*sigma_x_d[i] ) ;
            c3_d_yfield[i] = ( 2.*kappa_x_d[i] - dt*sigma_x_d[i] ) / ( 2.*kappa_x_d[i] + dt*sigma_x_d[i] ) ;
            c4_d_yfield[i] = ( 1. ) / ( 2.*kappa_x_d[i] + dt*sigma_x_d[i] ) ;
            c5_d_xfield[i] = ( 2.*kappa_x_d[i] + dt*sigma_x_d[i] ) ;
            c6_d_xfield[i] = ( 2.*kappa_x_d[i] - dt*sigma_x_d[i] ) ;
        }
    } // End X

    if (min_or_max==0){
        for( int j = 0 ; j< (int) ny_p ; j++ ) {
            c1_p_xfield[j] = ( 2.*kappa_y_p[(ny_p-1)-j] - dt*sigma_y_p[(ny_p-1)-j] ) / ( 2.*kappa_y_p[(ny_p-1)-j] + dt*sigma_y_p[(ny_p-1)-j] ) ;
            c2_p_xfield[j] = ( 2*dt ) / ( 2.*kappa_y_p[(ny_p-1)-j] + dt*sigma_y_p[(ny_p-1)-j] ) ;
            c3_p_zfield[j] = ( 2.*kappa_y_p[(ny_p-1)-j] - dt*sigma_y_p[(ny_p-1)-j] ) / ( 2.*kappa_y_p[(ny_p-1)-j] + dt*sigma_y_p[(ny_p-1)-j] ) ;
            c4_p_zfield[j] = ( 1. ) / ( 2.*kappa_y_p[(ny_p-1)-j] + dt*sigma_y_p[(ny_p-1)-j] ) ;
            c5_p_yfield[j] = ( 2.*kappa_y_p[(ny_p-1)-j] + dt*sigma_y_p[(ny_p-1)-j] ) ;
            c6_p_yfield[j] = ( 2.*kappa_y_p[(ny_p-1)-j] - dt*sigma_y_p[(ny_p-1)-j] ) ;
        }

        for( int j = 0 ; j< (int) ny_d ; j++ ) {
            c1_d_xfield[j] = ( 2.*kappa_y_d[(ny_d-1)-j] - dt*sigma_y_d[(ny_d-1)-j] ) / ( 2.*kappa_y_d[(ny_d-1)-j] + dt*sigma_y_d[(ny_d-1)-j] ) ;
            c2_d_xfield[j] = ( 2*dt ) / ( 2.*kappa_y_d[(ny_d-1)-j] + dt*sigma_y_d[(ny_d-1)-j] ) ;
            c3_d_zfield[j] = ( 2.*kappa_y_d[(ny_d-1)-j] - dt*sigma_y_d[(ny_d-1)-j] ) / ( 2.*kappa_y_d[(ny_d-1)-j] + dt*sigma_y_d[(ny_d-1)-j] ) ;
            c4_d_zfield[j] = ( 1. ) / ( 2.*kappa_y_d[(ny_d-1)-j] + dt*sigma_y_d[(ny_d-1)-j] ) ;
            c5_d_yfield[j] = ( 2.*kappa_y_d[(ny_d-1)-j] + dt*sigma_y_d[(ny_d-1)-j] ) ;
            c6_d_yfield[j] = ( 2.*kappa_y_d[(ny_d-1)-j] - dt*sigma_y_d[(ny_d-1)-j] ) ;
        }
    }
    else if (min_or_max==1){
        for( int j = 0 ; j< (int) ny_p ; j++ ) {
            c1_p_xfield[j] = ( 2.*kappa_y_p[j] - dt*sigma_y_p[j] ) / ( 2.*kappa_y_p[j] + dt*sigma_y_p[j] ) ;
            c2_p_xfield[j] = ( 2*dt ) / ( 2.*kappa_y_p[j] + dt*sigma_y_p[j] ) ;
            c3_p_zfield[j] = ( 2.*kappa_y_p[j] - dt*sigma_y_p[j] ) / ( 2.*kappa_y_p[j] + dt*sigma_y_p[j] ) ;
            c4_p_zfield[j] = ( 1. ) / ( 2.*kappa_y_p[j] + dt*sigma_y_p[j] ) ;
            c5_p_yfield[j] = ( 2.*kappa_y_p[j] + dt*sigma_y_p[j] ) ;
            c6_p_yfield[j] = ( 2.*kappa_y_p[j] - dt*sigma_y_p[j] ) ;
        }

        for( int j = 0 ; j< (int) ny_d ; j++ ) {
            c1_d_xfield[j] = ( 2.*kappa_y_d[j] - dt*sigma_y_d[j] ) / ( 2.*kappa_y_d[j] + dt*sigma_y_d[j] ) ;
            c2_d_xfield[j] = ( 2*dt ) / ( 2.*kappa_y_d[j] + dt*sigma_y_d[j] ) ;
            c3_d_zfield[j] = ( 2.*kappa_y_d[j] - dt*sigma_y_d[j] ) / ( 2.*kappa_y_d[j] + dt*sigma_y_d[j] ) ;
            c4_d_zfield[j] = ( 1. ) / ( 2.*kappa_y_d[j] + dt*sigma_y_d[j] ) ;
            c5_d_yfield[j] = ( 2.*kappa_y_d[j] + dt*sigma_y_d[j] ) ;
            c6_d_yfield[j] = ( 2.*kappa_y_d[j] - dt*sigma_y_d[j] ) ;
        }
    } // End Y

    if (min_or_max==0){
        for ( int k=0 ; k<1 ; k++ ) {
            c1_p_yfield[k] = ( 2.*kappa_z_p[k] - dt*sigma_z_p[k] ) / ( 2.*kappa_z_p[k] + dt*sigma_z_p[k] ) ;
            c2_p_yfield[k] = ( 2*dt ) / ( 2.*kappa_z_p[k] + dt*sigma_z_p[k] ) ;
            c3_p_xfield[k] = ( 2.*kappa_z_p[k] - dt*sigma_z_p[k] ) / ( 2.*kappa_z_p[k] + dt*sigma_z_p[k] ) ;
            c4_p_xfield[k] = ( 1. ) / ( 2.*kappa_z_p[k] + dt*sigma_z_p[k] ) ;
            c5_p_zfield[k] = ( 2.*kappa_z_p[k] + dt*sigma_z_p[k] ) ;
            c6_p_zfield[k] = ( 2.*kappa_z_p[k] - dt*sigma_z_p[k] ) ;
        }

        for ( int k=0 ; k<1 ; k++ ) {
            c1_d_yfield[k] = ( 2.*kappa_z_d[k] - dt*sigma_z_d[k] ) / ( 2.*kappa_z_d[k] + dt*sigma_z_d[k] ) ;
            c2_d_yfield[k] = ( 2*dt ) / ( 2.*kappa_z_d[k] + dt*sigma_z_d[k] ) ;
            c3_d_xfield[k] = ( 2.*kappa_z_d[k] - dt*sigma_z_d[k] ) / ( 2.*kappa_z_d[k] + dt*sigma_z_d[k] ) ;
            c4_d_xfield[k] = ( 1. ) / ( 2.*kappa_z_d[k] + dt*sigma_z_d[k] ) ;
            c5_d_zfield[k] = ( 2.*kappa_z_d[k] + dt*sigma_z_d[k] ) ;
            c6_d_zfield[k] = ( 2.*kappa_z_d[k] - dt*sigma_z_d[k] ) ;
        }
    }
    else if (min_or_max==1){
        for ( int k=0 ; k<1 ; k++ ) {
            c1_p_yfield[k] = ( 2.*kappa_z_p[k] - dt*sigma_z_p[k] ) / ( 2.*kappa_z_p[k] + dt*sigma_z_p[k] ) ;
            c2_p_yfield[k] = ( 2*dt ) / ( 2.*kappa_z_p[k] + dt*sigma_z_p[k] ) ;
            c3_p_xfield[k] = ( 2.*kappa_z_p[k] - dt*sigma_z_p[k] ) / ( 2.*kappa_z_p[k] + dt*sigma_z_p[k] ) ;
            c4_p_xfield[k] = ( 1. ) / ( 2.*kappa_z_p[k] + dt*sigma_z_p[k] ) ;
            c5_p_zfield[k] = ( 2.*kappa_z_p[k] + dt*sigma_z_p[k] ) ;
            c6_p_zfield[k] = ( 2.*kappa_z_p[k] - dt*sigma_z_p[k] ) ;
        }

        for ( int k=0 ; k<1 ; k++ ) {
            c1_d_yfield[k] = ( 2.*kappa_z_d[k] - dt*sigma_z_d[k] ) / ( 2.*kappa_z_d[k] + dt*sigma_z_d[k] ) ;
            c2_d_yfield[k] = ( 2*dt ) / ( 2.*kappa_z_d[k] + dt*sigma_z_d[k] ) ;
            c3_d_xfield[k] = ( 2.*kappa_z_d[k] - dt*sigma_z_d[k] ) / ( 2.*kappa_z_d[k] + dt*sigma_z_d[k] ) ;
            c4_d_xfield[k] = ( 1. ) / ( 2.*kappa_z_d[k] + dt*sigma_z_d[k] ) ;
            c5_d_zfield[k] = ( 2.*kappa_z_d[k] + dt*sigma_z_d[k] ) ;
            c6_d_zfield[k] = ( 2.*kappa_z_d[k] - dt*sigma_z_d[k] ) ;
        }
    } // End Z
}

void PML_Solver2D_Yee::compute_E_from_D( ElectroMagn *fields, int iDim, int min_or_max, std::vector<unsigned int> dimPrim, unsigned int solvermin, unsigned int solvermax )
{
    const unsigned int nx_p = dimPrim[0];
    const unsigned int nx_d = dimPrim[0] + 1;
    const unsigned int ny_p = dimPrim[1];
    const unsigned int ny_d = dimPrim[1] + 1;
    
    ElectroMagnBC2D_PML* pml_fields = static_cast<ElectroMagnBC2D_PML*>( fields->emBoundCond[iDim*2+min_or_max] );
    Field2D* Ex_pml = NULL;
    Field2D* Ey_pml = NULL;
    Field2D* Ez_pml = NULL;
    Field2D* Hx_pml = NULL;
    Field2D* Hy_pml = NULL;
    Field2D* Hz_pml = NULL;
    Field2D* Dx_pml = NULL;
    Field2D* Dy_pml = NULL;
    Field2D* Dz_pml = NULL;

    Ex_pml = pml_fields->Ex_;
    Ey_pml = pml_fields->Ey_;
    Ez_pml = pml_fields->Ez_;
    Hx_pml = pml_fields->Hx_;
    Hy_pml = pml_fields->Hy_;
    Hz_pml = pml_fields->Hz_;
    Dx_pml = pml_fields->Dx_;
    Dy_pml = pml_fields->Dy_;
    Dz_pml = pml_fields->Dz_;

    if (min_or_max==0){
        isMin=true;
        isMax=false;
    }
    else if (min_or_max==1){
        isMin=false;
        isMax=true;
    }

    if (iDim == 0) {
        //Electric field Ex^(d,p,p) Remind that in PML, there no current
        for( unsigned int k=0 ; k<1 ; k++ ) {
            for( unsigned int i=solvermin ; i<(unsigned int)solvermax ; i++ ) {
                for( unsigned int j=0 ; j<ny_p ; j++ ) {
                    // Standard FDTD
                    // ( *Ex_pml )( i, j ) = + 1. * ( *Ex_pml )( i, j )
                    //                       + dt * ( ( *Hz_pml )( i, j+1 ) - ( *Hz_pml )( i, j ) )/dy;
                    // ( *Dx_pml )( i, j ) = 1*( *Ex_pml )( i, j );
                    // PML FDTD
                    Dx_pml_old = 1*( *Dx_pml )( i, j ) ;
                    ( *Dx_pml )( i, j ) = + c1_p_xfield[j] * ( *Dx_pml )( i, j )
                                          + c2_p_xfield[j] * ( ( *Hz_pml )( i, j+1 ) - ( *Hz_pml )( i, j ) )/dy;
                    ( *Ex_pml )( i, j ) = + c3_p_xfield[k] * ( *Ex_pml )( i, j )
                                          + c4_p_xfield[k] * (c5_d_xfield[i]*( *Dx_pml )( i, j ) - c6_d_xfield[i]*Dx_pml_old ) ;
                }
            }
        }
        //Electric field Ey^(p,d,p) Remind that in PML, there no current
        for( unsigned int k=0 ; k<1 ; k++ ) {
            for( unsigned int i=solvermin ; i<(unsigned int)solvermax ; i++ ) {
                for( unsigned int j=0 ; j<ny_d ; j++ ) {
                    // Standard FDTD
                    // ( *Ey_pml )( i, j ) = + 1. * ( *Ey_pml )( i, j )
                    //                       - dt * ( ( *Hz_pml )( i+1, j ) - ( *Hz_pml )( i, j ) )/dx;
                    // ( *Dy_pml )( i, j ) = 1*( *Ey_pml )( i, j );
                    // PML FDTD
                    Dy_pml_old = 1*( *Dy_pml )( i, j ) ;
                    ( *Dy_pml )( i, j ) = + c1_p_yfield[k] * ( *Dy_pml )( i, j )
                                          - c2_p_yfield[k] * ( ( *Hz_pml )( i+1, j ) - ( *Hz_pml )( i, j ) )/dx;
                    ( *Ey_pml )( i, j ) = + c3_p_yfield[i] * ( *Ey_pml )( i, j )
                                          + c4_p_yfield[i] * (c5_d_yfield[j]*( *Dy_pml )( i, j ) - c6_d_yfield[j]*Dy_pml_old ) ;
                }
            }
        }
        //Electric field Ez^(p,p,d) Remind that in PML, there no current
        for( unsigned int k=0 ; k<1 ; k++ ) {
            for( unsigned int i=solvermin ; i<(unsigned int)solvermax ; i++ ) {
                for( unsigned int j=0 ; j<ny_p ; j++ ) {
                    // Standard FDTD
                    // ( *Ez_pml )( i, j ) = + 1. * ( *Ez_pml )( i, j )
                    //                       - dt * ( ( ( *Hx_pml )( i, j+1 ) - ( *Hx_pml )( i, j ) )/dy - ( ( *Hy_pml )( i+1, j ) - ( *Hy_pml )( i, j ) )/dx );
                    // ( *Dz_pml )( i, j ) = 1*( *Ez_pml )( i, j );
                    // PML FDTD
                    Dz_pml_old = 1*( *Dz_pml )( i, j ) ;
                    ( *Dz_pml )( i, j ) = + c1_p_zfield[i] * ( *Dz_pml )( i, j )
                                          - c2_p_zfield[i] * ( ( ( *Hx_pml )( i, j+1 ) - ( *Hx_pml )( i, j ) )/dy - ( ( *Hy_pml )( i+1, j ) - ( *Hy_pml )( i, j ) )/dx );
                    ( *Ez_pml )( i, j ) = + c3_p_zfield[j] * ( *Ez_pml )( i, j )
                                          + c4_p_zfield[j] * (c5_d_zfield[k]*( *Dz_pml )( i, j ) - c6_d_zfield[k]*Dz_pml_old );
                }
            }
        }
    }
    else if (iDim == 1) {
        //Electric field Ex^(d,p,p) Remind that in PML, there no current
        for( unsigned int k=0 ; k<1 ; k++ ) {
            for( unsigned int i=0 ; i<nx_d ; i++ ) {
                for( unsigned int j=solvermin ; j<(unsigned int)solvermax ; j++ ) {
                    // Standard FDTD
                    // ( *Ex_pml )( i, j ) = + 1. * ( *Ex_pml )( i, j )
                    //                       + dt * ( ( *Hz_pml )( i, j+1 ) - ( *Hz_pml )( i, j ) )/dy;
                    // ( *Dx_pml )( i, j ) = 1*( *Ex_pml )( i, j );
                    // PML FDTD
                    Dx_pml_old = 1*( *Dx_pml )( i, j ) ;
                    ( *Dx_pml )( i, j ) = + c1_p_xfield[j] * ( *Dx_pml )( i, j )
                                          + c2_p_xfield[j] * ( ( *Hz_pml )( i, j+1 ) - ( *Hz_pml )( i, j ) )/dy;
                    ( *Ex_pml )( i, j ) = + c3_p_xfield[k] * ( *Ex_pml )( i, j )
                                          + c4_p_xfield[k] * (c5_d_xfield[i]*( *Dx_pml )( i, j ) - c6_d_xfield[i]*Dx_pml_old ) ;
                }
            }
        }
        //Electric field Ey^(p,d,p) Remind that in PML, there no current
        for( unsigned int k=0 ; k<1 ; k++ ) {
            for( unsigned int i=0 ; i<nx_p ; i++ ) {
                for( unsigned int j=solvermin ; j<(unsigned int)solvermax ; j++ ) {
                    // Standard FDTD
                    // ( *Ey_pml )( i, j ) = + 1. * ( *Ey_pml )( i, j )
                    //                       - dt * ( ( *Hz_pml )( i+1, j ) - ( *Hz_pml )( i, j ) )/dx;
                    // ( *Dy_pml )( i, j ) = 1*( *Ey_pml )( i, j );
                    // PML FDTD
                    Dy_pml_old = 1*( *Dy_pml )( i, j ) ;
                    ( *Dy_pml )( i, j ) = + c1_p_yfield[k] * ( *Dy_pml )( i, j )
                                          - c2_p_yfield[k] * ( ( *Hz_pml )( i+1, j ) - ( *Hz_pml )( i, j ) )/dx;
                    ( *Ey_pml )( i, j ) = + c3_p_yfield[i] * ( *Ey_pml )( i, j )
                                          + c4_p_yfield[i] * (c5_d_yfield[j]*( *Dy_pml )( i, j ) - c6_d_yfield[j]*Dy_pml_old ) ;
                }
            }
        }
        //Electric field Ez^(p,p,d) Remind that in PML, there no current
        for( unsigned int k=0 ; k<1 ; k++ ) {
            for( unsigned int i=0 ; i<nx_p ; i++ ) {
                for( unsigned int j=solvermin ; j<(unsigned int)solvermax ; j++ ) {
                    // Standard FDTD
                    // ( *Ez_pml )( i, j ) = + 1. * ( *Ez_pml )( i, j )
                    //                       - dt * ( ( ( *Hx_pml )( i, j+1 ) - ( *Hx_pml )( i, j ) )/dy - ( ( *Hy_pml )( i+1, j ) - ( *Hy_pml )( i, j ) )/dx );
                    // ( *Dz_pml )( i, j ) = 1*( *Ez_pml )( i, j );
                    // PML FDTD
                    Dz_pml_old = 1*( *Dz_pml )( i, j ) ;
                    ( *Dz_pml )( i, j ) = + c1_p_zfield[i] * ( *Dz_pml )( i, j )
                                          - c2_p_zfield[i] * ( ( ( *Hx_pml )( i, j+1 ) - ( *Hx_pml )( i, j ) )/dy - ( ( *Hy_pml )( i+1, j ) - ( *Hy_pml )( i, j ) )/dx );
                    ( *Ez_pml )( i, j ) = + c3_p_zfield[j] * ( *Ez_pml )( i, j )
                                          + c4_p_zfield[j] * (c5_d_zfield[k]*( *Dz_pml )( i, j ) - c6_d_zfield[k]*Dz_pml_old );
                }
            }
        }
    }
}

void PML_Solver2D_Yee::compute_H_from_B( ElectroMagn *fields, int iDim, int min_or_max, std::vector<unsigned int> dimPrim, unsigned int solvermin, unsigned int solvermax )
{
    const unsigned int nx_p = dimPrim[0];
    const unsigned int nx_d = dimPrim[0] + 1;
    const unsigned int ny_p = dimPrim[1];
    const unsigned int ny_d = dimPrim[1] + 1;
    
    ElectroMagnBC2D_PML* pml_fields = static_cast<ElectroMagnBC2D_PML*>( fields->emBoundCond[iDim*2+min_or_max] );
    Field2D* Ex_pml = NULL;
    Field2D* Ey_pml = NULL;
    Field2D* Ez_pml = NULL;
    Field2D* Hx_pml = NULL;
    Field2D* Hy_pml = NULL;
    Field2D* Hz_pml = NULL;
    Field2D* Bx_pml = NULL;
    Field2D* By_pml = NULL;
    Field2D* Bz_pml = NULL;

    Ex_pml = pml_fields->Ex_;
    Ey_pml = pml_fields->Ey_;
    Ez_pml = pml_fields->Ez_;
    Hx_pml = pml_fields->Hx_;
    Hy_pml = pml_fields->Hy_;
    Hz_pml = pml_fields->Hz_;
    Bx_pml = pml_fields->Bx_;
    By_pml = pml_fields->By_;
    Bz_pml = pml_fields->Bz_;

    if (min_or_max==0){
        isMin=true;
        isMax=false;
    }
    else if (min_or_max==1){
        isMin=false;
        isMax=true;
    }

    if (iDim==0){
        //Magnetic field Bx^(p,d,d) Remind that in PML, there no current
        for( unsigned int k=0 ; k<1 ; k++ ) {
            for( unsigned int i=solvermin ; i<(unsigned int)solvermax ; i++ ) {
                for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
                    // Standard FDTD
                    // ( *Bx_pml )( i, j ) = + 1 * ( *Bx_pml )( i, j )
                    //                       - dt * ( ( *Ez_pml )( i, j ) - ( *Ez_pml )( i, j-1 ) )/dy;
                    // ( *Hx_pml )( i, j ) = + 1 * ( *Bx_pml )( i, j );
                    // PML FDTD
                    Bx_pml_old = 1*( *Bx_pml )( i, j ) ;
                    ( *Bx_pml )( i, j ) = + c1_d_xfield[j] * ( *Bx_pml )( i, j )
                                          - c2_d_xfield[j] * ( ( *Ez_pml )( i, j ) - ( *Ez_pml )( i, j-1 ) )/dy;
                    ( *Hx_pml )( i, j ) = + c3_d_xfield[k] * ( *Hx_pml )( i, j )
                                          + c4_d_xfield[k] * (c5_p_xfield[i]*( *Bx_pml )( i, j ) - c6_p_xfield[i]*Bx_pml_old ) ;
                }
            }
        }
        //Magnetic field By^(d,p,d) Remind that in PML, there no current
        for( unsigned int k=0 ; k<1 ; k++ ) {
            for( unsigned int i=solvermin ; i<(unsigned int)solvermax ; i++ ) {
                for( unsigned int j=0 ; j<ny_p ; j++ ) {
                    // Standard FDTD
                    // ( *By_pml )( i, j ) = + 1 * ( *By_pml )( i, j )
                    //                       + dt * ( ( *Ez_pml )( i, j ) - ( *Ez_pml )( i-1, j ) )/dx;
                    // ( *Hy_pml )( i, j ) = + 1 * ( *By_pml )( i, j );
                    // PML FDTD
                    By_pml_old = 1*( *By_pml )( i, j ) ;
                    ( *By_pml )( i, j ) = + c1_d_yfield[k] * ( *By_pml )( i, j )
                                          + c2_d_yfield[k] * ( ( *Ez_pml )( i, j ) - ( *Ez_pml )( i-1, j ) )/dx;
                    ( *Hy_pml )( i, j ) = + c3_d_yfield[i] * ( *Hy_pml )( i, j )
                                          + c4_d_yfield[i] * (c5_p_yfield[j]*( *By_pml )( i, j ) - c6_p_yfield[j]*By_pml_old ) ;
                }
            }
        }
        //Magnetic field Bz^(d,d,p) Remind that in PML, there no current
        for( unsigned int k=0 ; k<1 ; k++ ) {
            for( unsigned int i=solvermin ; i<(unsigned int)solvermax ; i++ ) {
                for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
                    // Standard FDTD
                    // ( *Bz_pml )( i, j ) = + 1 * ( *Bz_pml )( i, j )
                    //                       + dt * ( ( ( *Ex_pml )( i, j ) - ( *Ex_pml )( i, j-1 ) )/dy - ( ( *Ey_pml )( i, j ) - ( *Ey_pml )( i-1, j ) )/dx );
                    // ( *Hz_pml )( i, j ) = + 1 * ( *Bz_pml )( i, j );
                    // PML FDTD
                    Bz_pml_old = 1*( *Bz_pml )( i, j ) ;
                    ( *Bz_pml )( i, j ) = + c1_d_zfield[i] * ( *Bz_pml )( i, j )
                                          + c2_d_zfield[i] * ( ( ( *Ex_pml )( i, j ) - ( *Ex_pml )( i, j-1 ) )/dy - ( ( *Ey_pml )( i, j ) - ( *Ey_pml )( i-1, j ) )/dx );
                    ( *Hz_pml )( i, j ) = + c3_d_zfield[j] * ( *Hz_pml )( i, j )
                                          + c4_d_zfield[j] * (c5_p_zfield[k]*( *Bz_pml )( i, j ) - c6_p_zfield[k]*Bz_pml_old );
                }
            }
        }
    }
    else if (iDim==1) {
        //Magnetic field Bx^(p,d,d) Remind that in PML, there no current
        for( unsigned int k=0 ; k<1 ; k++ ) {
            for( unsigned int i=0 ; i<nx_p ; i++ ) {
                for( unsigned int j=solvermin ; j<(unsigned int)solvermax ; j++ ) {
                    // Standard FDTD
                    // ( *Bx_pml )( i, j ) = + 1 * ( *Bx_pml )( i, j )
                    //                       - dt * ( ( *Ez_pml )( i, j ) - ( *Ez_pml )( i, j-1 ) )/dy;
                    // ( *Hx_pml )( i, j ) = + 1 * ( *Bx_pml )( i, j );
                    // PML FDTD
                    Bx_pml_old = 1*( *Bx_pml )( i, j ) ;
                    ( *Bx_pml )( i, j ) = + c1_d_xfield[j] * ( *Bx_pml )( i, j )
                                          - c2_d_xfield[j] * ( ( *Ez_pml )( i, j ) - ( *Ez_pml )( i, j-1 ) )/dy;
                    ( *Hx_pml )( i, j ) = + c3_d_xfield[k] * ( *Hx_pml )( i, j )
                                          + c4_d_xfield[k] * (c5_p_xfield[i]*( *Bx_pml )( i, j ) - c6_p_xfield[i]*Bx_pml_old ) ;
                }
            }
        }
        //Magnetic field By^(d,p,d) Remind that in PML, there no current
        for( unsigned int k=0 ; k<1 ; k++ ) {
            for( unsigned int i=1 ; i<nx_d-1 ; i++ ) {
                for( unsigned int j=solvermin ; j<(unsigned int)solvermax ; j++ ) {
                    // Standard FDTD
                    // ( *By_pml )( i, j ) = + 1 * ( *By_pml )( i, j )
                    //                       + dt * ( ( *Ez_pml )( i, j ) - ( *Ez_pml )( i-1, j ) )/dx;
                    // ( *Hy_pml )( i, j ) = + 1 * ( *By_pml )( i, j );
                    // PML FDTD
                    By_pml_old = 1*( *By_pml )( i, j ) ;
                    ( *By_pml )( i, j ) = + c1_d_yfield[k] * ( *By_pml )( i, j )
                                          + c2_d_yfield[k] * ( ( *Ez_pml )( i, j ) - ( *Ez_pml )( i-1, j ) )/dx;
                    ( *Hy_pml )( i, j ) = + c3_d_yfield[i] * ( *Hy_pml )( i, j )
                                          + c4_d_yfield[i] * (c5_p_yfield[j]*( *By_pml )( i, j ) - c6_p_yfield[j]*By_pml_old ) ;
                }
            }
        }
        //Magnetic field Bz^(d,d,p) Remind that in PML, there no current
        for( unsigned int k=0 ; k<1 ; k++ ) {
            for( unsigned int i=1 ; i<nx_d-1 ; i++ ) {
                for( unsigned int j=solvermin ; j<(unsigned int)solvermax ; j++ ) {
                    // Standard FDTD
                    // ( *Bz_pml )( i, j ) = + 1 * ( *Bz_pml )( i, j )
                    //                       + dt * ( ( ( *Ex_pml )( i, j ) - ( *Ex_pml )( i, j-1 ) )/dy - ( ( *Ey_pml )( i, j ) - ( *Ey_pml )( i-1, j ) )/dx );
                    // ( *Hz_pml )( i, j ) = + 1 * ( *Bz_pml )( i, j );
                    // PML FDTD
                    Bz_pml_old = 1*( *Bz_pml )( i, j ) ;
                    ( *Bz_pml )( i, j ) = + c1_d_zfield[i] * ( *Bz_pml )( i, j )
                                          + c2_d_zfield[i] * ( ( ( *Ex_pml )( i, j ) - ( *Ex_pml )( i, j-1 ) )/dy - ( ( *Ey_pml )( i, j ) - ( *Ey_pml )( i-1, j ) )/dx );
                    ( *Hz_pml )( i, j ) = + c3_d_zfield[j] * ( *Hz_pml )( i, j )
                                          + c4_d_zfield[j] * (c5_p_zfield[k]*( *Bz_pml )( i, j ) - c6_p_zfield[k]*Bz_pml_old );
                }
            }
        }
    }
}
