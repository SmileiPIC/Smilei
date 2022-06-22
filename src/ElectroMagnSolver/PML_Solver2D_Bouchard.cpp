#include "PML_Solver2D_Bouchard.h"
#include "ElectroMagn.h"
#include "ElectroMagnBC2D_PML.h"
#include "Field2D.h"
#include "Patch.h"

#include <algorithm>

PML_Solver2D_Bouchard::PML_Solver2D_Bouchard(Params &params)
    : Solver2D(params)
{
    //ERROR("Under development, not yet working");
    double dt = params.timestep;
    dx = params.cell_length[0];
    dy = params.cell_length[1];
    double dx_ov_dt  = dx/dt;
    //double dy_ov_dt  = dy/dt;
    //double dt_ov_dx  = dt/dx;
    //double dt_ov_dy  = dt/dy;
    if( dx!=dy ) {
        ERROR( "Bouchard solver requires the same cell-length in x and y directions" );
    }
    if( dx_ov_dt!=2 ) {
        WARNING( "Bouchard solver requires dx/dt = 2 (Magical Timestep)" );
    }

    // On the axes v_phi^max = 1.01c and is below c @ 0.54 kxdx/pi
    // So there could existe a numerical cherenkov emission at this point
    // On the diagonal v_phi^max = 1.01c and is below c @ 0.85 sqrt((kxdx)^2+(kydy)^2)
    double delta = 0.1222*(1-pow(2.,2))/4. ;
    double beta = -0.1727*(1-0.5*pow(2.,2)-4.*delta)/4. ;
    double alpha = 1-2.*beta-3.*delta;

    beta_xy = beta;
    beta_yx = beta;
    delta_y = delta;
    delta_x = delta;

    alpha_y =  alpha;
    alpha_x =  alpha;

    // Ax    = alpha_x*dt/dx;
    // Ay    = alpha_y*dt/dy;
    // Bx    = beta_xy*dt/dx;
    // By    = beta_yx*dt/dy;
    // Dx    = delta_x*dt/dy;
    // Dy    = delta_y*dt/dy;

    Ax    = alpha_x/dx;
    Ay    = alpha_y/dy;
    Bx    = beta_xy/dx;
    By    = beta_yx/dy;
    Dx    = delta_x/dy;
    Dy    = delta_y/dy;

    //Define here the value of coefficient kappa_x_max, power_kappa_x, sigma_x_max, power_sigma_x
    sigma_x_max = params.pml_parameters[0][0];
    kappa_x_max = params.pml_parameters[0][2];
    sigma_power_pml_x = params.pml_parameters[0][1];
    kappa_power_pml_x = params.pml_parameters[0][3];
    //Define here the value of coefficient kappa_y_max, power_kappa_y, sigma_y_max, power_sigma_y
    sigma_y_max = params.pml_parameters[1][0];
    kappa_y_max = params.pml_parameters[1][2];
    sigma_power_pml_y = params.pml_parameters[1][1];
    kappa_power_pml_y = params.pml_parameters[1][3];
}

PML_Solver2D_Bouchard::~PML_Solver2D_Bouchard()
{
}

void PML_Solver2D_Bouchard::operator() ( ElectroMagn* fields )
{
    ERROR( "This is not a solver for the main domain" );
}

void PML_Solver2D_Bouchard::setDomainSizeAndCoefficients( int iDim, int min_or_max, int ncells_pml_domain, int startpml, int* ncells_pml_min, int* ncells_pml_max, Patch* patch )
{
    if ( iDim == 0 ) {
        nx_p = ncells_pml_domain;
        nx_d = ncells_pml_domain+1;
    }
    else if ( iDim == 1 ) {
        ny_p = ncells_pml_domain;
        ny_d = ncells_pml_domain+1;
        nx_p += ncells_pml_min[0] + ncells_pml_max[0];
        nx_d += ncells_pml_min[0] + ncells_pml_max[0];
    }

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

    // Quote of the first primal grid-point where PML solver will be apply
    // The first dual grid-point is at getDomainLocalMax( 0 ) - 0.5*dx and getDomainLocalMax( 1 ) - 0.5*dr
    // xmax = patch->getDomainLocalMax( 0 );
    // ymax = patch->getDomainLocalMax( 1 );

    if ( iDim == 0 ) {
        // 3 cells (oversize) are vaccum so the PML media begin at r0 which is :
        // Eventually the size of PML media is :
        length_y_pml = 0 ;
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
        for ( int i=startpml; i<nx_p ; i++ ) {
            kappa_x_p[i] = 1. + (kappa_x_max - 1.) * pow( (i-startpml)*dx , kappa_power_pml_x ) / pow( length_x_pml , kappa_power_pml_x ) ;
            sigma_x_p[i] = sigma_x_max * pow( (i-startpml)*dx , sigma_power_pml_x ) / pow( length_x_pml , sigma_power_pml_x ) ;
        }
        // Y-direction
        for ( int j=0 ; j<ny_p ; j++ ) {
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
        for ( int i=startpml+1 ; i<nx_d ; i++ ) {
            kappa_x_d[i] = 1. + (kappa_x_max - 1.) * pow( (i-startpml-0.5)*dx , kappa_power_pml_x ) / pow( length_x_pml , kappa_power_pml_x ) ;
            sigma_x_d[i] = sigma_x_max * pow( (i-startpml-0.5)*dx , sigma_power_pml_x ) / pow( length_x_pml , sigma_power_pml_x ) ;
        }
        // Y-direction
        for ( int j=0 ; j<ny_d ; j++ ) {
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
        for ( int i=0 ; i<nx_p ; i++ ) {
            kappa_x_p[i] = 1. ;
            sigma_x_p[i] = 0. ;
        }
        if (ncells_pml_min[0] != 0 ){
            for ( int i=0 ; i<ncells_pml_min[0] ; i++ ) {
                kappa_x_p[i] = 1. + (kappa_x_max - 1.) * pow( ( ncells_pml_min[0] - 1 - i )*dx , kappa_power_pml_x ) / pow( length_x_pml_xmin , kappa_power_pml_x ) ;
                sigma_x_p[i] = sigma_x_max * pow( ( ncells_pml_min[0] - 1 - i )*dx , sigma_power_pml_x ) / pow( length_x_pml_xmin , sigma_power_pml_x ) ;
            }
        }
        if (ncells_pml_max[0] != 0 ){
            for ( int i=(nx_p-1)-(ncells_pml_max[0]-1) ; i<nx_p ; i++ ) {
                kappa_x_p[i] = 1. + (kappa_x_max - 1.) * pow( ( i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) )*dx , kappa_power_pml_x ) / pow( length_x_pml_xmax , kappa_power_pml_x ) ;
                sigma_x_p[i] = sigma_x_max * pow( (i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) )*dx , sigma_power_pml_x ) / pow( length_x_pml_xmax , sigma_power_pml_x ) ;
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
        for ( int j=startpml ; j<ny_p ; j++ ) {
            kappa_y_p[j] = 1. + (kappa_y_max - 1.) * pow( (j-startpml)*dy , kappa_power_pml_y ) / pow( length_y_pml , kappa_power_pml_y ) ;
            sigma_y_p[j] = sigma_y_max * pow( (j-startpml)*dy , sigma_power_pml_y ) / pow( length_y_pml , sigma_power_pml_y ) ;
        }
        // Z-direction
        for ( int k=0 ;k<1 ; k++ ) {
            kappa_z_p[k] = 1. ;
            sigma_z_p[k] = 0. ;
        }
        // Dual grid
        // X-direction
        for ( int i=0 ; i<nx_d ; i++ ) {
            kappa_x_d[i] = 1. ;
            sigma_x_d[i] = 0. ;
        }
        if (ncells_pml_min[0] != 0 ){
            for ( int i=0 ; i<ncells_pml_min[0] ; i++ ) {
                kappa_x_d[i] = 1. + (kappa_x_max - 1.) * pow( ( 0.5 + ncells_pml_min[0] - 1 - i )*dx , kappa_power_pml_x ) / pow( length_x_pml_xmin , kappa_power_pml_x ) ;
                sigma_x_d[i] = sigma_x_max * pow( ( 0.5 + ncells_pml_min[0] - 1 - i )*dx , sigma_power_pml_x ) / pow( length_x_pml_xmin , sigma_power_pml_x ) ;
            }
        }
        if (ncells_pml_max[0] != 0 ){
            for ( int i=(nx_p-1)-(ncells_pml_max[0]-1)+1 ; i<nx_d ; i++ ) {
                kappa_x_d[i] = 1. + (kappa_x_max - 1.) * pow( (i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) - 0.5 )*dx , kappa_power_pml_x ) / pow( length_x_pml_xmax , kappa_power_pml_x ) ;
                sigma_x_d[i] = sigma_x_max * pow( (i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) - 0.5 )*dx , sigma_power_pml_x ) / pow( length_x_pml_xmax , sigma_power_pml_x ) ;
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
        for ( int j=startpml+1 ; j<ny_d ; j++ ) {
            kappa_y_d[j] = 1. + (kappa_y_max - 1.) * pow( (j-startpml-0.5)*dy , kappa_power_pml_y ) / pow( length_y_pml , kappa_power_pml_y ) ;
            sigma_y_d[j] = sigma_y_max * pow( (j-startpml-0.5)*dy , sigma_power_pml_y ) / pow( length_y_pml , sigma_power_pml_y ) ;
        }
        // Z-direction
        for ( int k=0 ; k<1 ; k++ ) {
            kappa_z_d[k] = 1. ;
            sigma_z_d[k] = 0. ;
        }
    }

    if ((min_or_max==0)&&(iDim==0)){
        for ( int i=0 ; i<nx_p ; i++ ) {
            c1_p_zfield[i] = ( 2.*kappa_x_p[(nx_p-1)-i] - dt*sigma_x_p[(nx_p-1)-i] ) / ( 2.*kappa_x_p[(nx_p-1)-i] + dt*sigma_x_p[(nx_p-1)-i] ) ;
            c2_p_zfield[i] = ( 2*dt ) / ( 2.*kappa_x_p[(nx_p-1)-i] + dt*sigma_x_p[(nx_p-1)-i] ) ;
            c3_p_yfield[i] = ( 2.*kappa_x_p[(nx_p-1)-i] - dt*sigma_x_p[(nx_p-1)-i] ) / ( 2.*kappa_x_p[(nx_p-1)-i] + dt*sigma_x_p[(nx_p-1)-i] ) ;
            c4_p_yfield[i] = ( 1. ) / ( 2.*kappa_x_p[(nx_p-1)-i] + dt*sigma_x_p[(nx_p-1)-i] ) ;
            c5_p_xfield[i] = ( 2.*kappa_x_p[(nx_p-1)-i] + dt*sigma_x_p[(nx_p-1)-i] ) ;
            c6_p_xfield[i] = ( 2.*kappa_x_p[(nx_p-1)-i] - dt*sigma_x_p[(nx_p-1)-i] ) ;
        }

        for ( int i=0 ; i<nx_d ; i++ ) {
            c1_d_zfield[i] = ( 2.*kappa_x_d[(nx_d-1)-i] - dt*sigma_x_d[(nx_d-1)-i] ) / ( 2.*kappa_x_d[(nx_d-1)-i] + dt*sigma_x_d[(nx_d-1)-i] ) ;
            c2_d_zfield[i] = ( 2*dt ) / ( 2.*kappa_x_d[(nx_d-1)-i] + dt*sigma_x_d[(nx_d-1)-i] ) ;
            c3_d_yfield[i] = ( 2.*kappa_x_d[(nx_d-1)-i] - dt*sigma_x_d[(nx_d-1)-i] ) / ( 2.*kappa_x_d[(nx_d-1)-i] + dt*sigma_x_d[(nx_d-1)-i] ) ;
            c4_d_yfield[i] = ( 1. ) / ( 2.*kappa_x_d[(nx_d-1)-i] + dt*sigma_x_d[(nx_d-1)-i] ) ;
            c5_d_xfield[i] = ( 2.*kappa_x_d[(nx_d-1)-i] + dt*sigma_x_d[(nx_d-1)-i] ) ;
            c6_d_xfield[i] = ( 2.*kappa_x_d[(nx_d-1)-i] - dt*sigma_x_d[(nx_d-1)-i] ) ;
        }
    }
    else {
        for ( int i=0 ; i<nx_p ; i++ ) {
            c1_p_zfield[i] = ( 2.*kappa_x_p[i] - dt*sigma_x_p[i] ) / ( 2.*kappa_x_p[i] + dt*sigma_x_p[i] ) ;
            c2_p_zfield[i] = ( 2*dt ) / ( 2.*kappa_x_p[i] + dt*sigma_x_p[i] ) ;
            c3_p_yfield[i] = ( 2.*kappa_x_p[i] - dt*sigma_x_p[i] ) / ( 2.*kappa_x_p[i] + dt*sigma_x_p[i] ) ;
            c4_p_yfield[i] = ( 1. ) / ( 2.*kappa_x_p[i] + dt*sigma_x_p[i] ) ;
            c5_p_xfield[i] = ( 2.*kappa_x_p[i] + dt*sigma_x_p[i] ) ;
            c6_p_xfield[i] = ( 2.*kappa_x_p[i] - dt*sigma_x_p[i] ) ;
        }

        for ( int i=0 ; i<nx_d ; i++ ) {
            c1_d_zfield[i] = ( 2.*kappa_x_d[i] - dt*sigma_x_d[i] ) / ( 2.*kappa_x_d[i] + dt*sigma_x_d[i] ) ;
            c2_d_zfield[i] = ( 2*dt ) / ( 2.*kappa_x_d[i] + dt*sigma_x_d[i] ) ;
            c3_d_yfield[i] = ( 2.*kappa_x_d[i] - dt*sigma_x_d[i] ) / ( 2.*kappa_x_d[i] + dt*sigma_x_d[i] ) ;
            c4_d_yfield[i] = ( 1. ) / ( 2.*kappa_x_d[i] + dt*sigma_x_d[i] ) ;
            c5_d_xfield[i] = ( 2.*kappa_x_d[i] + dt*sigma_x_d[i] ) ;
            c6_d_xfield[i] = ( 2.*kappa_x_d[i] - dt*sigma_x_d[i] ) ;
        }
    } // End X

    if (min_or_max==0){
        for ( int j=0 ; j<ny_p ; j++ ) {
            c1_p_xfield[j] = ( 2.*kappa_y_p[(ny_p-1)-j] - dt*sigma_y_p[(ny_p-1)-j] ) / ( 2.*kappa_y_p[(ny_p-1)-j] + dt*sigma_y_p[(ny_p-1)-j] ) ;
            c2_p_xfield[j] = ( 2*dt ) / ( 2.*kappa_y_p[(ny_p-1)-j] + dt*sigma_y_p[(ny_p-1)-j] ) ;
            c3_p_zfield[j] = ( 2.*kappa_y_p[(ny_p-1)-j] - dt*sigma_y_p[(ny_p-1)-j] ) / ( 2.*kappa_y_p[(ny_p-1)-j] + dt*sigma_y_p[(ny_p-1)-j] ) ;
            c4_p_zfield[j] = ( 1. ) / ( 2.*kappa_y_p[(ny_p-1)-j] + dt*sigma_y_p[(ny_p-1)-j] ) ;
            c5_p_yfield[j] = ( 2.*kappa_y_p[(ny_p-1)-j] + dt*sigma_y_p[(ny_p-1)-j] ) ;
            c6_p_yfield[j] = ( 2.*kappa_y_p[(ny_p-1)-j] - dt*sigma_y_p[(ny_p-1)-j] ) ;
        }

        for ( int j=0 ; j<ny_d ; j++ ) {
            c1_d_xfield[j] = ( 2.*kappa_y_d[(ny_d-1)-j] - dt*sigma_y_d[(ny_d-1)-j] ) / ( 2.*kappa_y_d[(ny_d-1)-j] + dt*sigma_y_d[(ny_d-1)-j] ) ;
            c2_d_xfield[j] = ( 2*dt ) / ( 2.*kappa_y_d[(ny_d-1)-j] + dt*sigma_y_d[(ny_d-1)-j] ) ;
            c3_d_zfield[j] = ( 2.*kappa_y_d[(ny_d-1)-j] - dt*sigma_y_d[(ny_d-1)-j] ) / ( 2.*kappa_y_d[(ny_d-1)-j] + dt*sigma_y_d[(ny_d-1)-j] ) ;
            c4_d_zfield[j] = ( 1. ) / ( 2.*kappa_y_d[(ny_d-1)-j] + dt*sigma_y_d[(ny_d-1)-j] ) ;
            c5_d_yfield[j] = ( 2.*kappa_y_d[(ny_d-1)-j] + dt*sigma_y_d[(ny_d-1)-j] ) ;
            c6_d_yfield[j] = ( 2.*kappa_y_d[(ny_d-1)-j] - dt*sigma_y_d[(ny_d-1)-j] ) ;
        }
    }
    else if (min_or_max==1){
        for ( int j=0 ; j<ny_p ; j++ ) {
            c1_p_xfield[j] = ( 2.*kappa_y_p[j] - dt*sigma_y_p[j] ) / ( 2.*kappa_y_p[j] + dt*sigma_y_p[j] ) ;
            c2_p_xfield[j] = ( 2*dt ) / ( 2.*kappa_y_p[j] + dt*sigma_y_p[j] ) ;
            c3_p_zfield[j] = ( 2.*kappa_y_p[j] - dt*sigma_y_p[j] ) / ( 2.*kappa_y_p[j] + dt*sigma_y_p[j] ) ;
            c4_p_zfield[j] = ( 1. ) / ( 2.*kappa_y_p[j] + dt*sigma_y_p[j] ) ;
            c5_p_yfield[j] = ( 2.*kappa_y_p[j] + dt*sigma_y_p[j] ) ;
            c6_p_yfield[j] = ( 2.*kappa_y_p[j] - dt*sigma_y_p[j] ) ;
        }

        for ( int j=0 ; j<ny_d ; j++ ) {
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

void PML_Solver2D_Bouchard::compute_E_from_D( ElectroMagn *fields, int iDim, int min_or_max, int solvermin, int solvermax )
{
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
            for( unsigned int i=solvermin ; i<solvermax ; i++ ) {
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
            for( unsigned int i=solvermin ; i<solvermax ; i++ ) {
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
            for( unsigned int i=solvermin ; i<solvermax ; i++ ) {
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
                for( unsigned int j=solvermin ; j<solvermax ; j++ ) {
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
                for( unsigned int j=solvermin ; j<solvermax ; j++ ) {
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
                for( unsigned int j=solvermin ; j<solvermax ; j++ ) {
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

void PML_Solver2D_Bouchard::compute_H_from_B( ElectroMagn *fields, int iDim, int min_or_max, int solvermin, int solvermax )
{
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
            for( unsigned int i=solvermin ; i<solvermax ; i++ ) {
                for( unsigned int j=2 ; j<ny_d-2 ; j++ ) {
                    // Standard FDTD
                    // ( *Bx_pml )( i, j ) = + 1 * ( *Bx_pml )( i, j )
                    //                       - dt * ( ( *Ez_pml )( i, j ) - ( *Ez_pml )( i, j-1 ) )/dy;
                    // NS FDTD
                    // (*Bx2D)(i,j) += Ay * (( *Ez_pml )(i,j-1)-( *Ez_pml )(i,j))
                    //       + By * (( *Ez_pml )(i+1,j-1)-( *Ez_pml )(i+1,j) + ( *Ez_pml )(i-1,j-1)-( *Ez_pml )(i-1,j))
                    //       + Dy * (( *Ez_pml )(i,j-2)-( *Ez_pml )(i,j+1));
                    // ( *Hx_pml )( i, j ) = + 1 * ( *Bx_pml )( i, j );
                    // PML FDTD
                    Bx_pml_old = 1*( *Bx_pml )( i, j ) ;
                    ( *Bx_pml )( i, j ) = + c1_d_xfield[j] * ( *Bx_pml )( i, j )
                                          + c2_d_xfield[j] * (
                                              Ay * (( *Ez_pml )(i,j-1)-( *Ez_pml )(i,j))
                                            + By * (( *Ez_pml )(i+1,j-1)-( *Ez_pml )(i+1,j) + ( *Ez_pml )(i-1,j-1)-( *Ez_pml )(i-1,j))
                                            + Dy * (( *Ez_pml )(i,j-2)-( *Ez_pml )(i,j+1))
                                          ) ;
                    ( *Hx_pml )( i, j ) = + c3_d_xfield[k] * ( *Hx_pml )( i, j )
                                          + c4_d_xfield[k] * (c5_p_xfield[i]*( *Bx_pml )( i, j ) - c6_p_xfield[i]*Bx_pml_old ) ;
                }
            }
        }
        //Magnetic field By^(d,p,d) Remind that in PML, there no current
        for( unsigned int k=0 ; k<1 ; k++ ) {
            for( unsigned int i=solvermin ; i<solvermax ; i++ ) {
                for( unsigned int j=1 ; j<ny_p-1 ; j++ ) {
                    // Standard FDTD
                    // ( *By_pml )( i, j ) = + 1 * ( *By_pml )( i, j )
                    //                       + dt * ( ( *Ez_pml )( i, j ) - ( *Ez_pml )( i-1, j ) )/dx;
                    // NS FDTD
                    // (*By2D)(i,j) += Ax * (( *Ez_pml )(i,j) - ( *Ez_pml )(i-1,j))
                    //       + Bx * (( *Ez_pml )(i,j+1)-( *Ez_pml )(i-1,j+1) + ( *Ez_pml )(i,j-1)-( *Ez_pml )(i-1,j-1))
                    //       + Dx * (( *Ez_pml )(i+1,j) - ( *Ez_pml )(i-2,j));
                    // ( *Hy_pml )( i, j ) = + 1 * ( *By_pml )( i, j );
                    // PML FDTD
                    By_pml_old = 1*( *By_pml )( i, j ) ;
                    ( *By_pml )( i, j ) = + c1_d_yfield[k] * ( *By_pml )( i, j )
                                          + c2_d_yfield[k] * (
                                              Ax * (( *Ez_pml )(i,j) - ( *Ez_pml )(i-1,j))
                                            + Bx * (( *Ez_pml )(i,j+1)-( *Ez_pml )(i-1,j+1) + ( *Ez_pml )(i,j-1)-( *Ez_pml )(i-1,j-1))
                                            + Dx * (( *Ez_pml )(i+1,j) - ( *Ez_pml )(i-2,j))
                                          );
                    ( *Hy_pml )( i, j ) = + c3_d_yfield[i] * ( *Hy_pml )( i, j )
                                          + c4_d_yfield[i] * (c5_p_yfield[j]*( *By_pml )( i, j ) - c6_p_yfield[j]*By_pml_old ) ;
                }
            }
        }
        //Magnetic field Bz^(d,d,p) Remind that in PML, there no current
        for( unsigned int k=0 ; k<1 ; k++ ) {
            for( unsigned int i=solvermin ; i<solvermax ; i++ ) {
                for( unsigned int j=2 ; j<ny_d-2 ; j++ ) {
                    // Standard FDTD
                    // ( *Bz_pml )( i, j ) = + 1 * ( *Bz_pml )( i, j )
                    //                       + dt * ( ( ( *Ex_pml )( i, j ) - ( *Ex_pml )( i, j-1 ) )/dy - ( ( *Ey_pml )( i, j ) - ( *Ey_pml )( i-1, j ) )/dx );
                    // NS FDTD
                    // (*Bz2D)(i,j) += Ay * (( *Ex_pml )(i,j)-( *Ex_pml )(i,j-1))
                    //             + By * (( *Ex_pml )(i+1,j)-( *Ex_pml )(i+1,j-1) + ( *Ex_pml )(i-1,j)-( *Ex_pml )(i-1,j-1))
                    //             + Dy * (( *Ex_pml )(i,j+1)-( *Ex_pml )(i,j-2))
                    //             + Ax * (( *Ey_pml )(i-1,j)-( *Ey_pml )(i,j))
                    //             + Bx * (( *Ey_pml )(i-1,j+1)-( *Ey_pml )(i,j+1) + ( *Ey_pml )(i-1,j-1)-( *Ey_pml )(i,j-1))
                    //             + Dx * (( *Ey_pml )(i-2,j)-( *Ey_pml )(i+1,j));
                    // ( *Hz_pml )( i, j ) = + 1 * ( *Bz_pml )( i, j );
                    // PML FDTD
                    Bz_pml_old = 1*( *Bz_pml )( i, j ) ;
                    ( *Bz_pml )( i, j ) = + c1_d_zfield[i] * ( *Bz_pml )( i, j )
                                          + c2_d_zfield[i] * (
                                              Ay * (( *Ex_pml )(i,j)-( *Ex_pml )(i,j-1))
                                            + By * (( *Ex_pml )(i+1,j)-( *Ex_pml )(i+1,j-1) + ( *Ex_pml )(i-1,j)-( *Ex_pml )(i-1,j-1))
                                            + Dy * (( *Ex_pml )(i,j+1)-( *Ex_pml )(i,j-2))
                                            + Ax * (( *Ey_pml )(i-1,j)-( *Ey_pml )(i,j))
                                            + Bx * (( *Ey_pml )(i-1,j+1)-( *Ey_pml )(i,j+1) + ( *Ey_pml )(i-1,j-1)-( *Ey_pml )(i,j-1))
                                            + Dx * (( *Ey_pml )(i-2,j)-( *Ey_pml )(i+1,j))
                                          );
                    ( *Hz_pml )( i, j ) = + c3_d_zfield[j] * ( *Hz_pml )( i, j )
                                          + c4_d_zfield[j] * (c5_p_zfield[k]*( *Bz_pml )( i, j ) - c6_p_zfield[k]*Bz_pml_old );
                }
            }
        }
    }
    else if (iDim==1) {
        //Magnetic field Bx^(p,d,d) Remind that in PML, there no current
        for( unsigned int k=0 ; k<1 ; k++ ) {
            for( unsigned int i=1 ; i<nx_p-1 ; i++ ) {
                for( unsigned int j=solvermin ; j<solvermax ; j++ ) {
                    // Standard FDTD
                    // ( *Bx_pml )( i, j ) = + 1 * ( *Bx_pml )( i, j )
                    //                       - dt * ( ( *Ez_pml )( i, j ) - ( *Ez_pml )( i, j-1 ) )/dy;
                    // ( *Hx_pml )( i, j ) = + 1 * ( *Bx_pml )( i, j );
                    // PML FDTD
                    Bx_pml_old = 1*( *Bx_pml )( i, j ) ;
                    ( *Bx_pml )( i, j ) = + c1_d_xfield[j] * ( *Bx_pml )( i, j )
                                          + c2_d_xfield[j] * (
                                              Ay * (( *Ez_pml )(i,j-1)-( *Ez_pml )(i,j))
                                            + By * (( *Ez_pml )(i+1,j-1)-( *Ez_pml )(i+1,j) + ( *Ez_pml )(i-1,j-1)-( *Ez_pml )(i-1,j))
                                            + Dy * (( *Ez_pml )(i,j-2)-( *Ez_pml )(i,j+1))
                                          ) ;
                    ( *Hx_pml )( i, j ) = + c3_d_xfield[k] * ( *Hx_pml )( i, j )
                                          + c4_d_xfield[k] * (c5_p_xfield[i]*( *Bx_pml )( i, j ) - c6_p_xfield[i]*Bx_pml_old ) ;
                }
            }
        }
        //Magnetic field By^(d,p,d) Remind that in PML, there no current
        for( unsigned int k=0 ; k<1 ; k++ ) {
            for( unsigned int i=2 ; i<nx_d-2 ; i++ ) {
                for( unsigned int j=solvermin ; j<solvermax ; j++ ) {
                    // Standard FDTD
                    // ( *By_pml )( i, j ) = + 1 * ( *By_pml )( i, j )
                    //                       + dt * ( ( *Ez_pml )( i, j ) - ( *Ez_pml )( i-1, j ) )/dx;
                    // ( *Hy_pml )( i, j ) = + 1 * ( *By_pml )( i, j );
                    // PML FDTD
                    By_pml_old = 1*( *By_pml )( i, j ) ;
                    ( *By_pml )( i, j ) = + c1_d_yfield[k] * ( *By_pml )( i, j )
                                          + c2_d_yfield[k] * (
                                              Ax * (( *Ez_pml )(i,j) - ( *Ez_pml )(i-1,j))
                                            + Bx * (( *Ez_pml )(i,j+1)-( *Ez_pml )(i-1,j+1) + ( *Ez_pml )(i,j-1)-( *Ez_pml )(i-1,j-1))
                                            + Dx * (( *Ez_pml )(i+1,j) - ( *Ez_pml )(i-2,j))
                                          );
                    ( *Hy_pml )( i, j ) = + c3_d_yfield[i] * ( *Hy_pml )( i, j )
                                          + c4_d_yfield[i] * (c5_p_yfield[j]*( *By_pml )( i, j ) - c6_p_yfield[j]*By_pml_old ) ;
                }
            }
        }
        //Magnetic field Bz^(d,d,p) Remind that in PML, there no current
        for( unsigned int k=0 ; k<1 ; k++ ) {
            for( unsigned int i=2 ; i<nx_d-2 ; i++ ) {
                for( unsigned int j=solvermin ; j<solvermax ; j++ ) {
                    // Standard FDTD
                    // ( *Bz_pml )( i, j ) = + 1 * ( *Bz_pml )( i, j )
                    //                       + dt * ( ( ( *Ex_pml )( i, j ) - ( *Ex_pml )( i, j-1 ) )/dy - ( ( *Ey_pml )( i, j ) - ( *Ey_pml )( i-1, j ) )/dx );
                    // ( *Hz_pml )( i, j ) = + 1 * ( *Bz_pml )( i, j );
                    // PML FDTD
                    Bz_pml_old = 1*( *Bz_pml )( i, j ) ;
                    ( *Bz_pml )( i, j ) = + c1_d_zfield[i] * ( *Bz_pml )( i, j )
                                          + c2_d_zfield[i] * (
                                              Ay * (( *Ex_pml )(i,j)-( *Ex_pml )(i,j-1))
                                            + By * (( *Ex_pml )(i+1,j)-( *Ex_pml )(i+1,j-1) + ( *Ex_pml )(i-1,j)-( *Ex_pml )(i-1,j-1))
                                            + Dy * (( *Ex_pml )(i,j+1)-( *Ex_pml )(i,j-2))
                                            + Ax * (( *Ey_pml )(i-1,j)-( *Ey_pml )(i,j))
                                            + Bx * (( *Ey_pml )(i-1,j+1)-( *Ey_pml )(i,j+1) + ( *Ey_pml )(i-1,j-1)-( *Ey_pml )(i,j-1))
                                            + Dx * (( *Ey_pml )(i-2,j)-( *Ey_pml )(i+1,j))
                                          );
                    ( *Hz_pml )( i, j ) = + c3_d_zfield[j] * ( *Hz_pml )( i, j )
                                          + c4_d_zfield[j] * (c5_p_zfield[k]*( *Bz_pml )( i, j ) - c6_p_zfield[k]*Bz_pml_old );
                }
            }
        }
    }
}
