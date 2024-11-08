#include "PML_Solver3D_Bouchard.h"
#include "ElectroMagn.h"
#include "ElectroMagnBC3D_PML.h"
#include "Field3D.h"
#include "Patch.h"

PML_Solver3D_Bouchard::PML_Solver3D_Bouchard( Params &params ):
    Solver3D( params ),
    pml_sigma_( 3, NULL ),
    pml_kappa_( 3, NULL )
{
    //ERROR("Under development, not yet working");
    double dt = params.timestep;
    dx = params.cell_length[0];
    dy = params.cell_length[1];
    dz = params.cell_length[2];
    double dx_ov_dt  = dx/dt;
    //double dy_ov_dt  = dy/dt;
    //double dz_ov_dt  = dz/dt;
    //double dt_ov_dx  = dt/dx;
    //double dt_ov_dy  = dt/dy;
    //double dt_ov_dz  = dt/dz;
    //Not necessary to have dx=dy=dz, but dispersion law are modify
    //In particular if dz >> dx,dy then solver become like the 2d solver
    //if( (dx!=dy)||(dx!=dz)||(dy!=dz) ) {
    //    ERROR( "Bouchard solver requires the same cell-length in x, y and z directions" );
    //}
    if( dx_ov_dt!=2 ) {
        WARNING( "Bouchard solver requires dx/dt = 2 (Magical Timestep)" );
    }

    double delta = -0.0916500000000000;//0.1222*(1-pow (2.,2))/4. ;
    double beta  =  0.027347045;//-0.1727*(1-0.5*pow (2.,2)-4.*delta)/4. ;
    double alpha =  1.16556182;//1-4.*beta-3.*delta ;

    delta_x = delta ;
    delta_y = delta ;
    delta_z = delta ;
    beta_xy = beta ;
    beta_yx = beta ;
    beta_xz = beta ;
    beta_zx = beta ;
    beta_yz = beta ;
    beta_zy = beta ;
    alpha_x = alpha ;
    alpha_y = alpha ;
    alpha_z = alpha ;

    Ax  = alpha_x/dx;
    Ay  = alpha_y/dy;
    Az  = alpha_z/dz;
    Bxy = beta_xy/dx;
    Byx = beta_yx/dy;
    Bxz = beta_xz/dx;
    Bzx = beta_zx/dz;
    Byz = beta_yz/dy;
    Bzy = beta_zy/dz;
    Dx  = delta_x/dx;
    Dy  = delta_y/dy;
    Dz  = delta_z/dz;

    ////Define here the value of coefficient kappa_x_max, power_kappa_x, sigma_x_max, power_sigma_x
    //sigma_x_max = params.pml_sigma_parameters[0][0];
    //kappa_x_max = params.pml_kappa_parameters[0][0];
    //sigma_power_pml_x = params.pml_sigma_parameters[0][1];
    //kappa_power_pml_x = params.pml_kappa_parameters[0][1];
    ////Define here the value of coefficient kappa_y_max, power_kappa_y, sigma_y_max, power_sigma_y
    //sigma_y_max = params.pml_sigma_parameters[1][0];
    //kappa_y_max = params.pml_kappa_parameters[1][0];
    //sigma_power_pml_y = params.pml_sigma_parameters[1][1];
    //kappa_power_pml_y = params.pml_kappa_parameters[1][1];
    ////Define here the value of coefficient kappa_z_max, power_kappa_z, sigma_z_max, power_sigma_z
    //sigma_z_max = params.pml_sigma_parameters[2][0];
    //kappa_z_max = params.pml_kappa_parameters[2][0];
    //sigma_power_pml_z = params.pml_sigma_parameters[2][1];
    //kappa_power_pml_z = params.pml_kappa_parameters[2][1];  

    std::vector<PyObject *> prof;
    if( PyTools::extract_pyProfiles( "pml_sigma", "Main", 0, prof )){
        if( prof.size() == 0 or prof.size() == 2 ){
            ERROR(" in pml_sigma, expecting a list of 1 or 3 profiles.");
        }
    // extracted profile // number of variables of the function // name of the profile extracted // params // try numpy ?? // try file ?? // time variable ??
        pml_sigma_[0] = new Profile( prof[0], 1, "pml_sigma_x_profile", params, true, false, false );
        if( prof.size() == 1){ 
            pml_sigma_[1] = new Profile( prof[0], 1, "pml_sigma_y_profile", params, true, false, false );
            pml_sigma_[2] = new Profile( prof[0], 1, "pml_sigma_z_profile", params, true, false, false );
        } else {
            pml_sigma_[1] = new Profile( prof[1], 1, "pml_sigma_y_profile", params, true, false, false );
            pml_sigma_[2] = new Profile( prof[2], 1, "pml_sigma_z_profile", params, true, false, false );
        }
    }
    if( PyTools::extract_pyProfiles( "pml_kappa", "Main", 0, prof )){
        if( prof.size() == 0 or prof.size() == 2 ){
            ERROR(" in pml_kappa, expecting a list of 1 or 3 profiles.");
        }
        pml_kappa_[0] = new Profile( prof[0], 1, "pml_kappa_x_profile", params, true, false, false );
        if( prof.size() == 1){ 
            pml_kappa_[1] = new Profile( prof[0], 1, "pml_kappa_y_profile", params, true, false, false );
            pml_kappa_[2] = new Profile( prof[0], 1, "pml_kappa_z_profile", params, true, false, false );
        } else {
            pml_kappa_[1] = new Profile( prof[1], 1, "pml_kappa_y_profile", params, true, false, false );
            pml_kappa_[2] = new Profile( prof[2], 1, "pml_kappa_z_profile", params, true, false, false );
        }
    }

}

PML_Solver3D_Bouchard::~PML_Solver3D_Bouchard()
{
    for( unsigned int i=0; i<pml_sigma_.size(); i++ ) {
        delete pml_sigma_[i];
    }
    for( unsigned int i=0; i<pml_kappa_.size(); i++ ) {
        delete pml_kappa_[i];
    }
}

void PML_Solver3D_Bouchard::operator()( ElectroMagn * )
{
    ERROR( "This is not a solver for the main domain" );
}

void PML_Solver3D_Bouchard::setDomainSizeAndCoefficients( int iDim, int min_or_max, std::vector<unsigned int> dimPrim, int ncells_pml_domain, int startpml, int* ncells_pml_min, int* ncells_pml_max, Patch* )
{
    const unsigned int nx_p = dimPrim[0];
    const unsigned int nx_d = dimPrim[0] + 1;
    const unsigned int ny_p = dimPrim[1];
    const unsigned int ny_d = dimPrim[1] + 1;
    const unsigned int nz_p = dimPrim[2];
    const unsigned int nz_d = dimPrim[2] + 1;

    //PML Coeffs Kappa,Sigma ...
    //Primal
    kappa_x_p.resize( nx_p );
    sigma_x_p.resize( nx_p );
    kappa_y_p.resize( ny_p );
    sigma_y_p.resize( ny_p );
    kappa_z_p.resize( nz_p );
    sigma_z_p.resize( nz_p );
    //Dual
    kappa_x_d.resize( nx_d );
    sigma_x_d.resize( nx_d );
    kappa_y_d.resize( ny_d );
    sigma_y_d.resize( ny_d );
    kappa_z_d.resize( nz_d );
    sigma_z_d.resize( nz_d );

    // While min and max of each dim have got the same ncells_pml, coefficient are the same
    // for min and max PML

    //Coeffs for Ex,Dx,Hx,Bx
    //Primal
    c1_p_xfield.resize( ny_p ); // j-dependent
    c2_p_xfield.resize( ny_p ); // j-dependent
    c3_p_xfield.resize( nz_p ); // k-dependent
    c4_p_xfield.resize( nz_p ); // k-dependent
    c5_p_xfield.resize( nx_p ); // i-dependent
    c6_p_xfield.resize( nx_p ); // i-dependent
    //Dual
    c1_d_xfield.resize( ny_d ); // j-dependent
    c2_d_xfield.resize( ny_d ); // j-dependent
    c3_d_xfield.resize( nz_d ); // j-dependent
    c4_d_xfield.resize( nz_d ); // j-dependent
    c5_d_xfield.resize( nx_d ); // i-dependent
    c6_d_xfield.resize( nx_d ); // i-dependent
    //Coefs for Ey,Dy,Hy,By
    //Primal
    c1_p_yfield.resize( nz_p ); // k-dependent
    c2_p_yfield.resize( nz_p ); // k-dependent
    c3_p_yfield.resize( nx_p ); // i-dependent
    c4_p_yfield.resize( nx_p ); // i-dependent
    c5_p_yfield.resize( ny_p ); // j-dependent
    c6_p_yfield.resize( ny_p ); // j-dependent
    //Dual
    c1_d_yfield.resize( nz_d ); // k-dependent
    c2_d_yfield.resize( nz_d ); // k-dependent
    c3_d_yfield.resize( nx_d ); // i-dependent
    c4_d_yfield.resize( nx_d ); // i-dependent
    c5_d_yfield.resize( ny_d ); // k-dependent
    c6_d_yfield.resize( ny_d ); // k-dependent
    //Coefs for Ez,Dz,Hz,Bz
    //Primal
    c1_p_zfield.resize( nx_p ); // j-dependent
    c2_p_zfield.resize( nx_p ); // j-dependent
    c3_p_zfield.resize( ny_p ); // i-dependent
    c4_p_zfield.resize( ny_p ); // i-dependent
    c5_p_zfield.resize( nz_p ); // j-dependent
    c6_p_zfield.resize( nz_p ); // j-dependent
    //Dual
    c1_d_zfield.resize( nx_d ); // i-dependent
    c2_d_zfield.resize( nx_d ); // i-dependent
    c3_d_zfield.resize( ny_d ); // j-dependent
    c4_d_zfield.resize( ny_d ); // j-dependent
    c5_d_zfield.resize( nz_d ); // k-dependent
    c6_d_zfield.resize( nz_d ); // k-dependent
 
    if ( iDim == 0 ) {
        // 3 cells (oversize) are vaccum so the PML media begin at y0 which is :
        // Eventually the size of PML media is :
        length_x_pml = (ncells_pml_domain-startpml+0.5)*dx ;
        length_y_pml = 0. ;
        length_z_pml = 0. ;
        // Primal grid
        // X-direction
        // Params for first cell of PML-patch (vacuum) i = 0,1,2
        for ( int i=0 ; i<startpml ; i++ ) {
            // Coeffs for the first cell
            kappa_x_p[i] = 1. ;
            sigma_x_p[i] = 0. ;
        }
        // Params for other cells (PML Media) when i>=3
        for( int i = startpml ; i< (int) nx_p ; i++ ) {
            //kappa_x_p[i] = 1. + (kappa_x_max - 1.) * pow( (i-startpml)*dx , kappa_power_pml_x ) / pow( length_x_pml , kappa_power_pml_x ) ;
            //sigma_x_p[i] = sigma_x_max * pow( (i-startpml)*dx , sigma_power_pml_x ) / pow( length_x_pml , sigma_power_pml_x ) ;
            kappa_x_p[i] = pml_kappa_[0]->valueAt((i-startpml)*dx/length_x_pml);
            sigma_x_p[i] = pml_sigma_[0]->valueAt((i-startpml)*dx/length_x_pml);
        }
        // Y-direction
        for( unsigned int j = 0 ; j<ny_p ; j++ ) {
            kappa_y_p[j] = 1. ;
            sigma_y_p[j] = 0. ;
        }
        // Z-direction
        for( unsigned int k = 0 ; k<nz_p ; k++ ) {
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
            //kappa_x_d[i] = 1. + (kappa_x_max - 1.) * pow( (i-startpml-0.5)*dx , kappa_power_pml_x ) / pow( length_x_pml , kappa_power_pml_x ) ;
            //sigma_x_d[i] = sigma_x_max * pow( (i-startpml-0.5)*dx , sigma_power_pml_x ) / pow( length_x_pml , sigma_power_pml_x ) ;
            kappa_x_d[i] = pml_kappa_[0]->valueAt((i-startpml-0.5)*dx/length_x_pml);
            sigma_x_d[i] = pml_sigma_[0]->valueAt((i-startpml-0.5)*dx/length_x_pml);
        }
        // Y-direction
        for( unsigned int j = 0 ; j<ny_d ; j++ ) {
            kappa_y_d[j] = 1. ;
            sigma_y_d[j] = 0. ;
        }
        // Z-direction
        for( unsigned int k = 0 ; k<nz_d ; k++ ) {
            kappa_z_d[k] = 1. ;
            sigma_z_d[k] = 0. ;
        }
    }

    if ( iDim == 1 ) {
        // 3 cells are vaccum so the PML media begin at y0 which is :
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
                //kappa_x_p[i] = 1. + (kappa_x_max - 1.) * pow( ( ncells_pml_min[0] - 1 - i )*dx , kappa_power_pml_x ) / pow( length_x_pml_xmin , kappa_power_pml_x ) ;
                //sigma_x_p[i] = sigma_x_max * pow( ( ncells_pml_min[0] - 1 - i )*dx , sigma_power_pml_x ) / pow( length_x_pml_xmin , sigma_power_pml_x ) ;
                kappa_x_p[i] = pml_kappa_[0]->valueAt((ncells_pml_min[0] - 1 - i)*dx/length_x_pml_xmin);
                sigma_x_p[i] = pml_sigma_[0]->valueAt((ncells_pml_min[0] - 1 - i)*dx/length_x_pml_xmin);
            }
        }
        if (ncells_pml_max[0] != 0 ){
            for( int i = (nx_p-1)-(ncells_pml_max[0]-1) ; i< (int) nx_p ; i++ ) {
                //kappa_x_p[i] = 1. + (kappa_x_max - 1.) * pow( ( i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) )*dx , kappa_power_pml_x ) / pow( length_x_pml_xmax , kappa_power_pml_x ) ;
                //sigma_x_p[i] = sigma_x_max * pow( (i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) )*dx , sigma_power_pml_x ) / pow( length_x_pml_xmax , sigma_power_pml_x ) ;
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
            //kappa_y_p[j] = 1. + (kappa_y_max - 1.) * pow( (j-startpml)*dy , kappa_power_pml_y ) / pow( length_y_pml , kappa_power_pml_y ) ;
            //sigma_y_p[j] = sigma_y_max * pow( (j-startpml)*dy , sigma_power_pml_y ) / pow( length_y_pml , sigma_power_pml_y ) ;
            kappa_y_p[j] = pml_kappa_[1]->valueAt((j-startpml)*dy/length_y_pml);
            sigma_y_p[j] = pml_sigma_[1]->valueAt((j-startpml)*dy/length_y_pml);
        }
        // Z-direction
        for( unsigned int k = 0 ;k<nz_p ; k++ ) {
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
                //kappa_x_d[i] = 1. + (kappa_x_max - 1.) * pow( ( 0.5 + ncells_pml_min[0] - 1 - i )*dx , kappa_power_pml_x ) / pow( length_x_pml_xmin , kappa_power_pml_x ) ;
                //sigma_x_d[i] = sigma_x_max * pow( ( 0.5 + ncells_pml_min[0] - 1 - i )*dx , sigma_power_pml_x ) / pow( length_x_pml_xmin , sigma_power_pml_x ) ;
                kappa_x_d[i] = pml_kappa_[0]->valueAt(( ncells_pml_min[0] - 0.5 - i )*dx/length_x_pml_xmin);
                sigma_x_d[i] = pml_sigma_[0]->valueAt(( ncells_pml_min[0] - 0.5 - i )*dx/length_x_pml_xmin);
            }
        }
        if (ncells_pml_max[0] != 0 ){
            for( int i = (nx_p-1)-(ncells_pml_max[0]-1)+1 ; i< (int) nx_d ; i++ ) {
                //kappa_x_d[i] = 1. + (kappa_x_max - 1.) * pow( (i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) - 0.5 )*dx , kappa_power_pml_x ) / pow( length_x_pml_xmax , kappa_power_pml_x ) ;
                //sigma_x_d[i] = sigma_x_max * pow( (i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) - 0.5 )*dx , sigma_power_pml_x ) / pow( length_x_pml_xmax , sigma_power_pml_x ) ;
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
            //kappa_y_d[j] = 1. + (kappa_y_max - 1.) * pow( (j-startpml-0.5)*dy , kappa_power_pml_y ) / pow( length_y_pml , kappa_power_pml_y ) ;
            //sigma_y_d[j] = sigma_y_max * pow( (j-startpml-0.5)*dy , sigma_power_pml_y ) / pow( length_y_pml , sigma_power_pml_y ) ;
            kappa_y_d[j] = pml_kappa_[1]->valueAt((j-startpml-0.5)*dy/length_y_pml);
            sigma_y_d[j] = pml_sigma_[1]->valueAt((j-startpml-0.5)*dy/length_y_pml);
        }
        // Z-direction
        for( unsigned int k = 0 ; k<nz_d ; k++ ) {
            kappa_z_d[k] = 1. ;
            sigma_z_d[k] = 0. ;
        }
    }

    if ( iDim == 2 ) {
        // 3 cells are vaccum so the PML media begin at z0 which is :
        // Eventually the size of PML media is :
        length_x_pml_xmax = (ncells_pml_max[0]+0.5)*dx ;
        length_x_pml_xmin = (ncells_pml_min[0]+0.5)*dx ;
        length_y_pml_ymax = (ncells_pml_max[1]+0.5)*dy ;
        length_y_pml_ymin = (ncells_pml_min[1]+0.5)*dy ;
        length_z_pml = (ncells_pml_domain-startpml+0.5)*dz ;
        // Primal grid
        // X-direction
        for( unsigned int i = 0 ; i<nx_p ; i++ ) {
            kappa_x_p[i] = 1. ;
            sigma_x_p[i] = 0. ;
        }
        if (ncells_pml_min[0] != 0 ){
            for ( int i=0 ; i<ncells_pml_min[0] ; i++ ) {
                //kappa_x_p[i] = 1. + (kappa_x_max - 1.) * pow( ( ncells_pml_min[0] - 1 - i )*dx , kappa_power_pml_x ) / pow( length_x_pml_xmin , kappa_power_pml_x ) ;
                //sigma_x_p[i] = sigma_x_max * pow( ( ncells_pml_min[0] - 1 - i )*dx , sigma_power_pml_x ) / pow( length_x_pml_xmin , sigma_power_pml_x ) ;
                kappa_x_p[i] = pml_kappa_[0]->valueAt(( ncells_pml_min[0] - 1 - i )*dx/length_x_pml_xmin);
                sigma_x_p[i] = pml_sigma_[0]->valueAt(( ncells_pml_min[0] - 1 - i )*dx/length_x_pml_xmin);
            }
        }
        if (ncells_pml_max[0] != 0 ){
            for( int i = (nx_p-1)-(ncells_pml_max[0]-1) ; i< (int) nx_p ; i++ ) {
                //kappa_x_p[i] = 1. + (kappa_x_max - 1.) * pow( ( i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) )*dx , kappa_power_pml_x ) / pow( length_x_pml_xmax , kappa_power_pml_x ) ;
                //sigma_x_p[i] = sigma_x_max * pow( (i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) )*dx , sigma_power_pml_x ) / pow( length_x_pml_xmax , sigma_power_pml_x ) ;
                kappa_x_p[i] = pml_kappa_[0]->valueAt(( i - nx_p+ncells_pml_max[0] )*dx/length_x_pml_xmax);
                sigma_x_p[i] = pml_sigma_[0]->valueAt(( i - nx_p+ncells_pml_max[0] )*dx/length_x_pml_xmax);
            }
        }
        // Y-direction
        for( unsigned int j = 0 ; j<ny_p ; j++ ) {
            kappa_y_p[j] = 1. ;
            sigma_y_p[j] = 0. ;
        }
        if (ncells_pml_min[1] != 0 ){
            for ( int j=0 ; j<ncells_pml_min[1] ; j++ ) {
                //kappa_y_p[j] = 1. + (kappa_y_max - 1.) * pow( ( ncells_pml_min[1] - 1 - j )*dy , kappa_power_pml_y ) / pow( length_y_pml_ymin , kappa_power_pml_y ) ;
                //sigma_y_p[j] = sigma_y_max * pow( ( ncells_pml_min[1] - 1 - j )*dy , sigma_power_pml_y ) / pow( length_y_pml_ymin , sigma_power_pml_y ) ;
                kappa_y_p[j] = pml_kappa_[1]->valueAt(( ncells_pml_min[1] - 1 - j )*dy/length_y_pml_ymin);
                sigma_y_p[j] = pml_sigma_[1]->valueAt(( ncells_pml_min[1] - 1 - j )*dy/length_y_pml_ymin);
            }
        }
        if (ncells_pml_max[1] != 0 ){
            for( int j = (ny_p-1)-(ncells_pml_max[1]-1) ; j< (int) ny_p ; j++ ) {
                //kappa_y_p[j] = 1. + (kappa_y_max - 1.) * pow( ( j - ( (ny_p-1)-(ncells_pml_max[1]-1) ) )*dy , kappa_power_pml_y ) / pow( length_y_pml_ymax , kappa_power_pml_y ) ;
                //sigma_y_p[j] = sigma_y_max * pow( (j - ( (ny_p-1)-(ncells_pml_max[1]-1) ) )*dy , sigma_power_pml_y ) / pow( length_y_pml_ymax , sigma_power_pml_y ) ;
                kappa_y_p[j] = pml_kappa_[1]->valueAt(( j - ny_p+ncells_pml_max[1] )*dy/length_y_pml_ymax);
                sigma_y_p[j] = pml_sigma_[1]->valueAt(( j - ny_p+ncells_pml_max[1] )*dy/length_y_pml_ymax);
            }
        }
        // Z-direction
        // Params for first cell of PML-patch (vacuum) i = 0,1,2
        for ( int k=0 ; k<startpml ; k++ ) {
            // Coeffs for the first cell
            kappa_z_p[k] = 1. ;
            sigma_z_p[k] = 0. ;
        }
        // Params for other cells (PML Media) when j>=3
        for( int k = startpml ; k< (int) nz_p ; k++ ) {
            //kappa_z_p[k] = 1. + (kappa_z_max - 1.) * pow( (k-startpml)*dz , kappa_power_pml_z ) / pow( length_z_pml , kappa_power_pml_z ) ;
            //sigma_z_p[k] = sigma_z_max * pow( (k-startpml)*dz , sigma_power_pml_z ) / pow( length_z_pml , sigma_power_pml_z ) ;
            kappa_z_p[k] = pml_kappa_[2]->valueAt((k-startpml)*dz/length_z_pml);
            sigma_z_p[k] = pml_sigma_[2]->valueAt((k-startpml)*dz/length_z_pml);
        }
        // Dual grid
        // X-direction
        for( unsigned int i = 0 ; i<nx_d ; i++ ) {
            kappa_x_d[i] = 1. ;
            sigma_x_d[i] = 0. ;
        }
        if (ncells_pml_min[0] != 0 ){
            for ( int i=0 ; i<ncells_pml_min[0] ; i++ ) {
                //kappa_x_d[i] = 1. + (kappa_x_max - 1.) * pow( ( 0.5 + ncells_pml_min[0] - 1 - i )*dx , kappa_power_pml_x ) / pow( length_x_pml_xmin , kappa_power_pml_x ) ;
                //sigma_x_d[i] = sigma_x_max * pow( ( 0.5 + ncells_pml_min[0] - 1 - i )*dx , sigma_power_pml_x ) / pow( length_x_pml_xmin , sigma_power_pml_x ) ;
                kappa_x_d[i] = pml_kappa_[0]->valueAt(( ncells_pml_min[0] - 0.5 - i )*dx/length_x_pml_xmin);
                sigma_x_d[i] = pml_sigma_[0]->valueAt(( ncells_pml_min[0] - 0.5 - i )*dx/length_x_pml_xmin);
            }
        }
        if (ncells_pml_max[0] != 0 ){
            for( int i = (nx_p-1)-(ncells_pml_max[0]-1)+1 ; i< (int) nx_d ; i++ ) {
                //kappa_x_d[i] = 1. + (kappa_x_max - 1.) * pow( (i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) - 0.5 )*dx , kappa_power_pml_x ) / pow( length_x_pml_xmax , kappa_power_pml_x ) ;
                //sigma_x_d[i] = sigma_x_max * pow( (i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) - 0.5 )*dx , sigma_power_pml_x ) / pow( length_x_pml_xmax , sigma_power_pml_x ) ;
                kappa_x_d[i] = pml_kappa_[0]->valueAt((i - nx_p + ncells_pml_max[0] - 0.5 )*dx/length_x_pml_xmax);
                sigma_x_d[i] = pml_sigma_[0]->valueAt((i - nx_p  + ncells_pml_max[0] - 0.5)*dx/length_x_pml_xmax);
            }
        }
        // Y-direction
        for( unsigned int j = 0 ; j<ny_d ; j++ ) {
            kappa_y_d[j] = 1. ;
            sigma_y_d[j] = 0. ;
        }
        if (ncells_pml_min[1] != 0 ){
            for ( int j=0 ; j<ncells_pml_min[1] ; j++ ) {
                //kappa_y_d[j] = 1. + (kappa_y_max - 1.) * pow( ( 0.5 + ncells_pml_min[1] - 1 - j )*dy , kappa_power_pml_y ) / pow( length_y_pml_ymin , kappa_power_pml_y ) ;
                //sigma_y_d[j] = sigma_y_max * pow( ( 0.5 + ncells_pml_min[1] - 1 - j )*dy , sigma_power_pml_y ) / pow( length_y_pml_ymin , sigma_power_pml_y ) ;
                kappa_y_d[j] = pml_kappa_[1]->valueAt(( ncells_pml_min[1] - 0.5 - j )*dy/length_y_pml_ymin);
                sigma_y_d[j] = pml_sigma_[1]->valueAt(( ncells_pml_min[1] - 0.5 - j )*dy/length_y_pml_ymin);
            }
        }
        if (ncells_pml_max[1] != 0 ){
            for( int j = (ny_p-1)-(ncells_pml_max[1]-1)+1 ; j< (int) ny_d ; j++ ) {
                //kappa_y_d[j] = 1. + (kappa_y_max - 1.) * pow( (j - ( (ny_p-1)-(ncells_pml_max[1]-1) ) - 0.5 )*dy , kappa_power_pml_y ) / pow( length_y_pml_ymax , kappa_power_pml_y ) ;
                //sigma_y_d[j] = sigma_y_max * pow( (j - ( (ny_p-1)-(ncells_pml_max[1]-1) ) - 0.5 )*dy , sigma_power_pml_y ) / pow( length_y_pml_ymax , sigma_power_pml_y ) ;
                kappa_y_d[j] = pml_kappa_[1]->valueAt((j - ny_p + ncells_pml_max[1] - 0.5)*dy/length_y_pml_ymax);
                sigma_y_d[j] = pml_sigma_[1]->valueAt((j - ny_p + ncells_pml_max[1] - 0.5)*dy/length_y_pml_ymax);
            }
        }
        // Z-direction
        // Params for first cell of PML-patch (vacuum) i = 0,1,2,3
        for ( int k=0 ; k<startpml+1 ; k++ ) {
            // Coeffs for the first cell
            kappa_z_d[k] = 1. ;
            sigma_z_d[k] = 0. ;
        }
        // Params for other cells (PML Media) when j>=4
        for( int k = startpml+1 ; k< (int) nz_d ; k++ ) {
            //kappa_z_d[k] = 1. + (kappa_z_max - 1.) * pow( (k-startpml-0.5)*dz , kappa_power_pml_z ) / pow( length_z_pml , kappa_power_pml_z ) ;
            //sigma_z_d[k] = sigma_z_max * pow( (k-startpml-0.5)*dz , sigma_power_pml_z ) / pow( length_z_pml , sigma_power_pml_z ) ;
            kappa_z_d[k] = pml_kappa_[2]->valueAt((k-startpml-0.5)*dz/length_z_pml);
            sigma_z_d[k] = pml_sigma_[2]->valueAt((k-startpml-0.5)*dz/length_z_pml);
        }
    }

    //Warning ncells_pml_max[0] and ncells_pml_max[1] could be different

    //Coefficients for PML in Xmin and Xmax
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

    //Coefficients for PML in Ymin and Ymax
    if ((min_or_max==0)&&(iDim==1)){
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
    else {
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

    //Coefficients for PML in Zmin and Zmax
    if (min_or_max==0){
        // for( int k = 0 ; k< (int) nz_p ; k++ ) {
        //     c1_p_yfield[k] = ( 2.*kappa_z_p[k] - dt*sigma_z_p[k] ) / ( 2.*kappa_z_p[k] + dt*sigma_z_p[k] ) ;
        //     c2_p_yfield[k] = ( 2*dt ) / ( 2.*kappa_z_p[k] + dt*sigma_z_p[k] ) ;
        //     c3_p_xfield[k] = ( 2.*kappa_z_p[k] - dt*sigma_z_p[k] ) / ( 2.*kappa_z_p[k] + dt*sigma_z_p[k] ) ;
        //     c4_p_xfield[k] = ( 1. ) / ( 2.*kappa_z_p[k] + dt*sigma_z_p[k] ) ;
        //     c5_p_zfield[k] = ( 2.*kappa_z_p[k] + dt*sigma_z_p[k] ) ;
        //     c6_p_zfield[k] = ( 2.*kappa_z_p[k] - dt*sigma_z_p[k] ) ;
        // }

        // for( int k = 0 ; k< (int) nz_d ; k++ ) {
        //     c1_d_yfield[k] = ( 2.*kappa_z_d[k] - dt*sigma_z_d[k] ) / ( 2.*kappa_z_d[k] + dt*sigma_z_d[k] ) ;
        //     c2_d_yfield[k] = ( 2*dt ) / ( 2.*kappa_z_d[k] + dt*sigma_z_d[k] ) ;
        //     c3_d_xfield[k] = ( 2.*kappa_z_d[k] - dt*sigma_z_d[k] ) / ( 2.*kappa_z_d[k] + dt*sigma_z_d[k] ) ;
        //     c4_d_xfield[k] = ( 1. ) / ( 2.*kappa_z_d[k] + dt*sigma_z_d[k] ) ;
        //     c5_d_zfield[k] = ( 2.*kappa_z_d[k] + dt*sigma_z_d[k] ) ;
        //     c6_d_zfield[k] = ( 2.*kappa_z_d[k] - dt*sigma_z_d[k] ) ;
        // }
        for( int k = 0 ; k< (int) nz_p ; k++ ) {
            c1_p_yfield[k] = ( 2.*kappa_z_p[(nz_p-1)-k] - dt*sigma_z_p[(nz_p-1)-k] ) / ( 2.*kappa_z_p[(nz_p-1)-k] + dt*sigma_z_p[(nz_p-1)-k] ) ;
            c2_p_yfield[k] = ( 2*dt ) / ( 2.*kappa_z_p[(nz_p-1)-k] + dt*sigma_z_p[(nz_p-1)-k] ) ;
            c3_p_xfield[k] = ( 2.*kappa_z_p[(nz_p-1)-k] - dt*sigma_z_p[(nz_p-1)-k] ) / ( 2.*kappa_z_p[(nz_p-1)-k] + dt*sigma_z_p[(nz_p-1)-k] ) ;
            c4_p_xfield[k] = ( 1. ) / ( 2.*kappa_z_p[(nz_p-1)-k] + dt*sigma_z_p[(nz_p-1)-k] ) ;
            c5_p_zfield[k] = ( 2.*kappa_z_p[(nz_p-1)-k] + dt*sigma_z_p[(nz_p-1)-k] ) ;
            c6_p_zfield[k] = ( 2.*kappa_z_p[(nz_p-1)-k] - dt*sigma_z_p[(nz_p-1)-k] ) ;
        }

        for( int k = 0 ; k< (int) nz_d ; k++ ) {
            c1_d_yfield[k] = ( 2.*kappa_z_d[(nz_d-1)-k] - dt*sigma_z_d[(nz_d-1)-k] ) / ( 2.*kappa_z_d[(nz_d-1)-k] + dt*sigma_z_d[(nz_d-1)-k] ) ;
            c2_d_yfield[k] = ( 2*dt ) / ( 2.*kappa_z_d[(nz_d-1)-k] + dt*sigma_z_d[(nz_d-1)-k] ) ;
            c3_d_xfield[k] = ( 2.*kappa_z_d[(nz_d-1)-k] - dt*sigma_z_d[(nz_d-1)-k] ) / ( 2.*kappa_z_d[(nz_d-1)-k] + dt*sigma_z_d[(nz_d-1)-k] ) ;
            c4_d_xfield[k] = ( 1. ) / ( 2.*kappa_z_d[(nz_d-1)-k] + dt*sigma_z_d[(nz_d-1)-k] ) ;
            c5_d_zfield[k] = ( 2.*kappa_z_d[(nz_d-1)-k] + dt*sigma_z_d[(nz_d-1)-k] ) ;
            c6_d_zfield[k] = ( 2.*kappa_z_d[(nz_d-1)-k] - dt*sigma_z_d[(nz_d-1)-k] ) ;
        }
    }
    else if (min_or_max==1){
        for( int k = 0 ; k< (int) nz_p ; k++ ) {
            c1_p_yfield[k] = ( 2.*kappa_z_p[k] - dt*sigma_z_p[k] ) / ( 2.*kappa_z_p[k] + dt*sigma_z_p[k] ) ;
            c2_p_yfield[k] = ( 2*dt ) / ( 2.*kappa_z_p[k] + dt*sigma_z_p[k] ) ;
            c3_p_xfield[k] = ( 2.*kappa_z_p[k] - dt*sigma_z_p[k] ) / ( 2.*kappa_z_p[k] + dt*sigma_z_p[k] ) ;
            c4_p_xfield[k] = ( 1. ) / ( 2.*kappa_z_p[k] + dt*sigma_z_p[k] ) ;
            c5_p_zfield[k] = ( 2.*kappa_z_p[k] + dt*sigma_z_p[k] ) ;
            c6_p_zfield[k] = ( 2.*kappa_z_p[k] - dt*sigma_z_p[k] ) ;
        }

        for( int k = 0 ; k< (int) nz_d ; k++ ) {
            c1_d_yfield[k] = ( 2.*kappa_z_d[k] - dt*sigma_z_d[k] ) / ( 2.*kappa_z_d[k] + dt*sigma_z_d[k] ) ;
            c2_d_yfield[k] = ( 2*dt ) / ( 2.*kappa_z_d[k] + dt*sigma_z_d[k] ) ;
            c3_d_xfield[k] = ( 2.*kappa_z_d[k] - dt*sigma_z_d[k] ) / ( 2.*kappa_z_d[k] + dt*sigma_z_d[k] ) ;
            c4_d_xfield[k] = ( 1. ) / ( 2.*kappa_z_d[k] + dt*sigma_z_d[k] ) ;
            c5_d_zfield[k] = ( 2.*kappa_z_d[k] + dt*sigma_z_d[k] ) ;
            c6_d_zfield[k] = ( 2.*kappa_z_d[k] - dt*sigma_z_d[k] ) ;
        }
    } // End Z
}

void PML_Solver3D_Bouchard::compute_E_from_D( ElectroMagn *fields, int iDim, int min_or_max, std::vector<unsigned int> dimPrim, unsigned int solvermin, unsigned int solvermax )
{
    const unsigned int nx_p = dimPrim[0];
    const unsigned int nx_d = dimPrim[0] + 1;
    const unsigned int ny_p = dimPrim[1];
    const unsigned int ny_d = dimPrim[1] + 1;
    const unsigned int nz_p = dimPrim[2];
    const unsigned int nz_d = dimPrim[2] + 1;
    
    ElectroMagnBC3D_PML* pml_fields = static_cast<ElectroMagnBC3D_PML*>( fields->emBoundCond[iDim*2+min_or_max] );
    Field3D* Ex_pml = NULL;
    Field3D* Ey_pml = NULL;
    Field3D* Ez_pml = NULL;
    Field3D* Hx_pml = NULL;
    Field3D* Hy_pml = NULL;
    Field3D* Hz_pml = NULL;
    Field3D* Dx_pml = NULL;
    Field3D* Dy_pml = NULL;
    Field3D* Dz_pml = NULL;

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
        for( unsigned int i=solvermin ; i<(unsigned int)solvermax ; i++ ) {
            for( unsigned int j=0 ; j<ny_p ; j++ ) {
                for( unsigned int k=0 ; k<nz_p ; k++ ) {
                    // Standard FDTD
                    // ( *Ex_pml )( i, j, k ) = + 1. * ( *Ex_pml )( i, j, k )
                    //                          + dt/dy * ( ( *Hz_pml )( i, j+1, k ) - ( *Hz_pml )( i, j, k ) )
                    //                          - dt/dz * ( ( *Hy_pml )( i, j, k+1 ) - ( *Hy_pml )( i, j, k ) );
                    // ( *Dx_pml )( i, j, k ) = 1*( *Ex_pml )( i, j, k );
                    // PML FDTD
                    Dx_pml_old = 1*( *Dx_pml )( i, j, k ) ;
                    ( *Dx_pml )( i, j, k ) = + c1_p_xfield[j] * ( *Dx_pml )( i, j, k )
                                             + c2_p_xfield[j]/dy * ( ( *Hz_pml )( i, j+1, k ) - ( *Hz_pml )( i, j, k ) )
                                             - c2_p_xfield[j]/dz * ( ( *Hy_pml )( i, j, k+1 ) - ( *Hy_pml )( i, j, k ) ) ;
                    ( *Ex_pml )( i, j, k ) = + c3_p_xfield[k] * ( *Ex_pml )( i, j, k )
                                             + c4_p_xfield[k] * ( c5_d_xfield[i]*( *Dx_pml )( i, j, k ) - c6_d_xfield[i]*Dx_pml_old ) ;
                }
            }
        }
        //Electric field Ey^(p,d,p) Remind that in PML, there no current
        for( unsigned int i=solvermin ; i<(unsigned int)solvermax ; i++ ) {
            for( unsigned int j=0 ; j<ny_d ; j++ ) {
                for( unsigned int k=0 ; k<nz_p ; k++ ) {
                    // Standard FDTD
                    // ( *Ey_pml )( i, j , k) = + 1. * ( *Ey_pml )( i, j, k )
                    //                       - dt/dx * ( ( *Hz_pml )( i+1, j, k ) - ( *Hz_pml )( i, j, k ) )
                    //                       + dt/dz * ( ( *Hx_pml )( i, j, k+1 ) - ( *Hx_pml )( i, j, k ) );
                    // ( *Dy_pml )( i, j, k ) = 1*( *Ey_pml )( i, j, k );
                    // PML FDTD
                    Dy_pml_old = 1*( *Dy_pml )( i, j, k ) ;
                    ( *Dy_pml )( i, j, k ) = + c1_p_yfield[k] * ( *Dy_pml )( i, j, k )
                                             - c2_p_yfield[k]/dx * ( ( *Hz_pml )( i+1, j, k ) - ( *Hz_pml )( i, j, k ) )
                                             + c2_p_yfield[k]/dz * ( ( *Hx_pml )( i, j, k+1 ) - ( *Hx_pml )( i, j, k ) ) ;
                    ( *Ey_pml )( i, j, k ) = + c3_p_yfield[i] * ( *Ey_pml )( i, j, k )
                                             + c4_p_yfield[i] * (c5_d_yfield[j]*( *Dy_pml )( i, j,k ) - c6_d_yfield[j]*Dy_pml_old ) ;
                }
            }
        }
        //Electric field Ez^(p,p,d) Remind that in PML, there no current
        for( unsigned int i=solvermin ; i<(unsigned int)solvermax ; i++ ) {
            for( unsigned int j=0 ; j<ny_p ; j++ ) {
                for( unsigned int k=0 ; k<nz_d ; k++ ) {
                    // Standard FDTD
                    // ( *Ez_pml )( i, j, k ) = + 1. * ( *Ez_pml )( i, j, k )
                    //                       - dt/dy * ( ( *Hx_pml )( i, j+1, k ) - ( *Hx_pml )( i, j, k ) )
                    //                       + dt/dx * ( ( *Hy_pml )( i+1, j, k ) - ( *Hy_pml )( i, j, k ) );
                    // ( *Dz_pml )( i, j, k ) = 1*( *Ez_pml )( i, j, k );
                    // PML FDTD
                    Dz_pml_old = 1*( *Dz_pml )( i, j, k ) ;
                    ( *Dz_pml )( i, j, k ) = + c1_p_zfield[i] * ( *Dz_pml )( i, j, k )
                                             - c2_p_zfield[i]/dy * ( ( *Hx_pml )( i, j+1, k ) - ( *Hx_pml )( i, j, k ) )
                                             + c2_p_zfield[i]/dx * ( ( *Hy_pml )( i+1, j, k ) - ( *Hy_pml )( i, j, k ) );
                    ( *Ez_pml )( i, j, k ) = + c3_p_zfield[j] * ( *Ez_pml )( i, j, k )
                                             + c4_p_zfield[j] * (c5_d_zfield[k]*( *Dz_pml )( i, j, k ) - c6_d_zfield[k]*Dz_pml_old );
                }
            }
        }
    }
    else if (iDim == 1) {
        //Electric field Ex^(d,p,p) Remind that in PML, there no current
        for( unsigned int i=0 ; i<nx_d ; i++ ) {
            for( unsigned int j=solvermin ; j<(unsigned int)solvermax ; j++ ) {
                for( unsigned int k=0 ; k<nz_p ; k++ ) {
                    // Standard FDTD
                    // ( *Ex_pml )( i, j, k ) = + 1. * ( *Ex_pml )( i, j, k )
                    //                          + dt/dy * ( ( *Hz_pml )( i, j+1, k ) - ( *Hz_pml )( i, j, k ) )
                    //                          - dt/dz * ( ( *Hy_pml )( i, j, k+1 ) - ( *Hy_pml )( i, j, k ) );
                    // ( *Dx_pml )( i, j, k ) = 1*( *Ex_pml )( i, j, k );
                    // PML FDTD
                    Dx_pml_old = 1*( *Dx_pml )( i, j, k ) ;
                    ( *Dx_pml )( i, j, k ) = + c1_p_xfield[j] * ( *Dx_pml )( i, j, k )
                                             + c2_p_xfield[j]/dy * ( ( *Hz_pml )( i, j+1, k ) - ( *Hz_pml )( i, j, k ) )
                                             - c2_p_xfield[j]/dz * ( ( *Hy_pml )( i, j, k+1 ) - ( *Hy_pml )( i, j, k ) ) ;
                    ( *Ex_pml )( i, j, k ) = + c3_p_xfield[k] * ( *Ex_pml )( i, j, k )
                                             + c4_p_xfield[k] * ( c5_d_xfield[i]*( *Dx_pml )( i, j, k ) - c6_d_xfield[i]*Dx_pml_old ) ;
                }
            }
        }
        //Electric field Ey^(p,d,p) Remind that in PML, there no current
        for( unsigned int i=0 ; i<nx_p ; i++ ) {
            for( unsigned int j=solvermin ; j<(unsigned int)solvermax ; j++ ) {
                for( unsigned int k=0 ; k<nz_p ; k++ ) {
                    // Standard FDTD
                    // ( *Ey_pml )( i, j , k) = + 1. * ( *Ey_pml )( i, j, k )
                    //                       - dt/dx * ( ( *Hz_pml )( i+1, j, k ) - ( *Hz_pml )( i, j, k ) )
                    //                       + dt/dz * ( ( *Hx_pml )( i, j, k+1 ) - ( *Hx_pml )( i, j, k ) );
                    // ( *Dy_pml )( i, j, k ) = 1*( *Ey_pml )( i, j, k );
                    // PML FDTD
                    Dy_pml_old = 1*( *Dy_pml )( i, j, k ) ;
                    ( *Dy_pml )( i, j, k ) = + c1_p_yfield[k] * ( *Dy_pml )( i, j, k )
                                             - c2_p_yfield[k]/dx * ( ( *Hz_pml )( i+1, j, k ) - ( *Hz_pml )( i, j, k ) )
                                             + c2_p_yfield[k]/dz * ( ( *Hx_pml )( i, j, k+1 ) - ( *Hx_pml )( i, j, k ) ) ;
                    ( *Ey_pml )( i, j, k ) = + c3_p_yfield[i] * ( *Ey_pml )( i, j, k )
                                             + c4_p_yfield[i] * (c5_d_yfield[j]*( *Dy_pml )( i, j,k ) - c6_d_yfield[j]*Dy_pml_old ) ;
                }
            }
        }
        //Electric field Ez^(p,p,d) Remind that in PML, there no current
        for( unsigned int i=0 ; i<nx_p ; i++ ) {
            for( unsigned int j=solvermin ; j<(unsigned int)solvermax ; j++ ) {
                for( unsigned int k=0 ; k<nz_d ; k++ ) {
                    // Standard FDTD
                    // ( *Ez_pml )( i, j, k ) = + 1. * ( *Ez_pml )( i, j, k )
                    //                       - dt/dy * ( ( *Hx_pml )( i, j+1, k ) - ( *Hx_pml )( i, j, k ) )
                    //                       + dt/dx * ( ( *Hy_pml )( i+1, j, k ) - ( *Hy_pml )( i, j, k ) );
                    // ( *Dz_pml )( i, j, k ) = 1*( *Ez_pml )( i, j, k );
                    // PML FDTD
                    Dz_pml_old = 1*( *Dz_pml )( i, j, k ) ;
                    ( *Dz_pml )( i, j, k ) = + c1_p_zfield[i] * ( *Dz_pml )( i, j, k )
                                             - c2_p_zfield[i]/dy * ( ( *Hx_pml )( i, j+1, k ) - ( *Hx_pml )( i, j, k ) )
                                             + c2_p_zfield[i]/dx * ( ( *Hy_pml )( i+1, j, k ) - ( *Hy_pml )( i, j, k ) );
                    ( *Ez_pml )( i, j, k ) = + c3_p_zfield[j] * ( *Ez_pml )( i, j, k )
                                             + c4_p_zfield[j] * (c5_d_zfield[k]*( *Dz_pml )( i, j, k ) - c6_d_zfield[k]*Dz_pml_old );
                }
            }
        }
    }
    else if (iDim == 2) {
        //Electric field Ex^(d,p,p) Remind that in PML, there no current
        for( unsigned int i=0 ; i<nx_d ; i++ ) {
            for( unsigned int j=0 ; j<ny_p ; j++ ) {
                for( unsigned int k=solvermin ; k<(unsigned int)solvermax ; k++ ) {
                    // Standard FDTD
                    // ( *Ex_pml )( i, j, k ) = + 1. * ( *Ex_pml )( i, j, k )
                    //                          + dt/dy * ( ( *Hz_pml )( i, j+1, k ) - ( *Hz_pml )( i, j, k ) )
                    //                          - dt/dz * ( ( *Hy_pml )( i, j, k+1 ) - ( *Hy_pml )( i, j, k ) );
                    // ( *Dx_pml )( i, j, k ) = 1*( *Ex_pml )( i, j, k );
                    // PML FDTD
                    Dx_pml_old = 1*( *Dx_pml )( i, j, k ) ;
                    ( *Dx_pml )( i, j, k ) = + c1_p_xfield[j] * ( *Dx_pml )( i, j, k )
                                             + c2_p_xfield[j]/dy * ( ( *Hz_pml )( i, j+1, k ) - ( *Hz_pml )( i, j, k ) )
                                             - c2_p_xfield[j]/dz * ( ( *Hy_pml )( i, j, k+1 ) - ( *Hy_pml )( i, j, k ) ) ;
                    ( *Ex_pml )( i, j, k ) = + c3_p_xfield[k] * ( *Ex_pml )( i, j, k )
                                             + c4_p_xfield[k] * ( c5_d_xfield[i]*( *Dx_pml )( i, j, k ) - c6_d_xfield[i]*Dx_pml_old ) ;
                }
            }
        }
        //Electric field Ey^(p,d,p) Remind that in PML, there no current
        for( unsigned int i=0 ; i<nx_p ; i++ ) {
            for( unsigned int j=0 ; j<ny_d ; j++ ) {
                for( unsigned int k=solvermin ; k<(unsigned int)solvermax ; k++ ) {
                    // Standard FDTD
                    // ( *Ey_pml )( i, j , k) = + 1. * ( *Ey_pml )( i, j, k )
                    //                       - dt/dx * ( ( *Hz_pml )( i+1, j, k ) - ( *Hz_pml )( i, j, k ) )
                    //                       + dt/dz * ( ( *Hx_pml )( i, j, k+1 ) - ( *Hx_pml )( i, j, k ) );
                    // ( *Dy_pml )( i, j, k ) = 1*( *Ey_pml )( i, j, k );
                    // PML FDTD
                    Dy_pml_old = 1*( *Dy_pml )( i, j, k ) ;
                    ( *Dy_pml )( i, j, k ) = + c1_p_yfield[k] * ( *Dy_pml )( i, j, k )
                                             - c2_p_yfield[k]/dx * ( ( *Hz_pml )( i+1, j, k ) - ( *Hz_pml )( i, j, k ) )
                                             + c2_p_yfield[k]/dz * ( ( *Hx_pml )( i, j, k+1 ) - ( *Hx_pml )( i, j, k ) ) ;
                    ( *Ey_pml )( i, j, k ) = + c3_p_yfield[i] * ( *Ey_pml )( i, j, k )
                                             + c4_p_yfield[i] * (c5_d_yfield[j]*( *Dy_pml )( i, j,k ) - c6_d_yfield[j]*Dy_pml_old ) ;
                }
            }
        }
        //Electric field Ez^(p,p,d) Remind that in PML, there no current
        for( unsigned int i=0 ; i<nx_p ; i++ ) {
            for( unsigned int j=0 ; j<ny_p ; j++ ) {
                for( unsigned int k=solvermin ; k<(unsigned int)solvermax ; k++ ) {
                    // Standard FDTD
                    // ( *Ez_pml )( i, j, k ) = + 1. * ( *Ez_pml )( i, j, k )
                    //                       - dt/dy * ( ( *Hx_pml )( i, j+1, k ) - ( *Hx_pml )( i, j, k ) )
                    //                       + dt/dx * ( ( *Hy_pml )( i+1, j, k ) - ( *Hy_pml )( i, j, k ) );
                    // ( *Dz_pml )( i, j, k ) = 1*( *Ez_pml )( i, j, k );
                    // PML FDTD
                    Dz_pml_old = 1*( *Dz_pml )( i, j, k ) ;
                    ( *Dz_pml )( i, j, k ) = + c1_p_zfield[i] * ( *Dz_pml )( i, j, k )
                                             - c2_p_zfield[i]/dy * ( ( *Hx_pml )( i, j+1, k ) - ( *Hx_pml )( i, j, k ) )
                                             + c2_p_zfield[i]/dx * ( ( *Hy_pml )( i+1, j, k ) - ( *Hy_pml )( i, j, k ) );
                    ( *Ez_pml )( i, j, k ) = + c3_p_zfield[j] * ( *Ez_pml )( i, j, k )
                                             + c4_p_zfield[j] * (c5_d_zfield[k]*( *Dz_pml )( i, j, k ) - c6_d_zfield[k]*Dz_pml_old );
                }
            }
        }
    }
}

void PML_Solver3D_Bouchard::compute_H_from_B( ElectroMagn *fields, int iDim, int min_or_max, std::vector<unsigned int> dimPrim, unsigned int solvermin, unsigned int solvermax )
{
    const unsigned int nx_p = dimPrim[0];
    const unsigned int nx_d = dimPrim[0] + 1;
    const unsigned int ny_p = dimPrim[1];
    const unsigned int ny_d = dimPrim[1] + 1;
    const unsigned int nz_p = dimPrim[2];
    const unsigned int nz_d = dimPrim[2] + 1;
    
    ElectroMagnBC3D_PML* pml_fields = static_cast<ElectroMagnBC3D_PML*>( fields->emBoundCond[iDim*2+min_or_max] );
    Field3D* Ex_pml = NULL;
    Field3D* Ey_pml = NULL;
    Field3D* Ez_pml = NULL;
    Field3D* Hx_pml = NULL;
    Field3D* Hy_pml = NULL;
    Field3D* Hz_pml = NULL;
    Field3D* Bx_pml = NULL;
    Field3D* By_pml = NULL;
    Field3D* Bz_pml = NULL;

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
        for( unsigned int i=solvermin ; i<(unsigned int)solvermax ; i++ ) {
            for( unsigned int j=2 ; j<ny_d-2 ; j++ ) {
                for( unsigned int k=2 ; k<nz_d-2 ; k++ ) {
                    // Standard FDTD
                    // ( *Bx_pml )( i, j, k ) = + 1 * ( *Bx_pml )( i, j, k )
                    //                          - dt/dy * ( ( *Ez_pml )( i, j, k ) - ( *Ez_pml )( i, j-1, k ) )
                    //                          + dt/dz * ( ( *Ey_pml )( i, j, k ) - ( *Ey_pml )( i, j, k-1 ) ) ;
                    // NS FDTD
                    // ( *Bx_pml )( i, j, k ) = + 1 * ( *Bx_pml )( i, j, k ) + dt * (
                    //     + Az *  ( ( *Ey_pml )( i, j, k )    - ( *Ey_pml )( i, j, k-1 ) )
                    //     + Bzx * ( ( *Ey_pml )( i+1, j, k )  - ( *Ey_pml )( i+1, j, k-1 ) + ( *Ey_pml )( i-1, j, k )-( *Ey_pml )( i-1, j, k-1 ) )
                    //     + Bzy * ( ( *Ey_pml )( i, j+1, k )  - ( *Ey_pml )( i, j+1, k-1 ) + ( *Ey_pml )( i, j-1, k )-( *Ey_pml )( i, j-1, k-1 ) )
                    //     +  Dz * ( ( *Ey_pml )( i, j, k+1 )  - ( *Ey_pml )( i, j, k-2) )
                    //     -  Ay * ( ( *Ez_pml )( i,  j, k )   - ( *Ez_pml )( i,  j-1, k ) )
                    //     - Byx * ( ( *Ez_pml )( i+1, j, k )  - ( *Ez_pml )( i+1, j-1, k ) + ( *Ez_pml )( i-1, j, k )-( *Ez_pml )( i-1, j-1, k ) )
                    //     - Byz * ( ( *Ez_pml )( i, j, k+1 )  - ( *Ez_pml )( i, j-1, k+1 ) + ( *Ez_pml )( i, j, k-1 )-( *Ez_pml )( i, j-1, k-1 ) )
                    //     - Dy *  ( ( *Ez_pml )( i, j+1, k )  - ( *Ez_pml )( i, j-2, k ) )
                    // );
                    // ( *Hx_pml )( i, j, k ) = + 1 * ( *Bx_pml )( i, j, k );
                    // PML FDTD
                    Bx_pml_old = 1*( *Bx_pml )( i, j, k ) ;
                    ( *Bx_pml )( i, j, k ) = + c1_d_xfield[j] * ( *Bx_pml )( i, j, k )
                                             + c2_d_xfield[j] * (
                                                + Az *  ( ( *Ey_pml )( i, j, k )    - ( *Ey_pml )( i, j, k-1 ) )
                                                + Bzx * ( ( *Ey_pml )( i+1, j, k )  - ( *Ey_pml )( i+1, j, k-1 ) + ( *Ey_pml )( i-1, j, k )-( *Ey_pml )( i-1, j, k-1 ) )
                                                + Bzy * ( ( *Ey_pml )( i, j+1, k )  - ( *Ey_pml )( i, j+1, k-1 ) + ( *Ey_pml )( i, j-1, k )-( *Ey_pml )( i, j-1, k-1 ) )
                                                +  Dz * ( ( *Ey_pml )( i, j, k+1 )  - ( *Ey_pml )( i, j, k-2) )
                                                -  Ay * ( ( *Ez_pml )( i,  j, k )   - ( *Ez_pml )( i,  j-1, k ) )
                                                - Byx * ( ( *Ez_pml )( i+1, j, k )  - ( *Ez_pml )( i+1, j-1, k ) + ( *Ez_pml )( i-1, j, k )-( *Ez_pml )( i-1, j-1, k ) )
                                                - Byz * ( ( *Ez_pml )( i, j, k+1 )  - ( *Ez_pml )( i, j-1, k+1 ) + ( *Ez_pml )( i, j, k-1 )-( *Ez_pml )( i, j-1, k-1 ) )
                                                - Dy *  ( ( *Ez_pml )( i, j+1, k )  - ( *Ez_pml )( i, j-2, k ) )
                                             ) ;
                    ( *Hx_pml )( i, j, k ) = + c3_d_xfield[k] * ( *Hx_pml )( i, j, k )
                                             + c4_d_xfield[k] * (c5_p_xfield[i]*( *Bx_pml )( i, j, k ) - c6_p_xfield[i]*Bx_pml_old ) ;
                }
            }
        }
        //Magnetic field By^(d,p,d) Remind that in PML, there no current
        for( unsigned int i=solvermin ; i<(unsigned int)solvermax ; i++ ) {
            for( unsigned int j=1 ; j<ny_p-1 ; j++ ) {
                for( unsigned int k=2 ; k<nz_d-2 ; k++ ) {
                    // Standard FDTD
                    // ( *By_pml )( i, j, k ) = + 1 * ( *By_pml )( i, j, k )
                    //                          + dt/dx * ( ( *Ez_pml )( i, j, k ) - ( *Ez_pml )( i-1, j, k ) )
                    //                          - dt/dz * ( ( *Ex_pml )( i, j, k ) - ( *Ex_pml )( i, j, k-1 ) );
                    // NS FDTD
                    // ( *By_pml )( i, j, k ) = + 1 * ( *By_pml )( i, j, k ) + dt * (
                    //     + Ax *  ( ( *Ez_pml )( i,  j, k )   - ( *Ez_pml )( i-1, j, k ) )
                    //     + Bxy * ( ( *Ez_pml )( i,  j+1, k ) - ( *Ez_pml )( i-1, j+1, k ) + ( *Ez_pml )( i, j-1, k )-( *Ez_pml )( i-1, j-1, k ) )
                    //     + Bxz * ( ( *Ez_pml )( i, j, k+1 )  - ( *Ez_pml )( i-1, j, k+1 ) + ( *Ez_pml )( i, j, k-1 )-( *Ez_pml )( i-1, j, k-1 ) )
                    //     +  Dx * ( ( *Ez_pml )( i+1, j, k )  - ( *Ez_pml )( i-2, j, k ) )
                    //     -  Az * ( ( *Ex_pml )( i, j, k )    - ( *Ex_pml )( i, j, k-1 ) )
                    //     - Bzx * ( ( *Ex_pml )( i+1, j, k )  - ( *Ex_pml )( i+1, j, k-1 ) + ( *Ex_pml )( i-1, j, k )-( *Ex_pml )( i-1, j, k-1 ) )
                    //     - Bzy * ( ( *Ex_pml )( i, j+1, k )  - ( *Ex_pml )( i, j+1, k-1 ) + ( *Ex_pml )( i, j-1, k )-( *Ex_pml )( i, j-1, k-1 ) )
                    //     -  Dz * ( ( *Ex_pml )( i, j, k+1 )  - ( *Ex_pml )( i, j, k-2) )
                    // );
                    // ( *Hy_pml )( i, j, k ) = + 1 * ( *By_pml )( i, j, k );
                    // PML FDTD
                    By_pml_old = 1*( *By_pml )( i, j, k ) ;
                    ( *By_pml )( i, j, k ) = + c1_d_yfield[k] * ( *By_pml )( i, j, k )
                                             + c2_d_yfield[k] * (
                                                + Ax *  ( ( *Ez_pml )( i,  j, k )   - ( *Ez_pml )( i-1, j, k ) )
                                                + Bxy * ( ( *Ez_pml )( i,  j+1, k ) - ( *Ez_pml )( i-1, j+1, k ) + ( *Ez_pml )( i, j-1, k )-( *Ez_pml )( i-1, j-1, k ) )
                                                + Bxz * ( ( *Ez_pml )( i, j, k+1 )  - ( *Ez_pml )( i-1, j, k+1 ) + ( *Ez_pml )( i, j, k-1 )-( *Ez_pml )( i-1, j, k-1 ) )
                                                +  Dx * ( ( *Ez_pml )( i+1, j, k )  - ( *Ez_pml )( i-2, j, k ) )
                                                -  Az * ( ( *Ex_pml )( i, j, k )    - ( *Ex_pml )( i, j, k-1 ) )
                                                - Bzx * ( ( *Ex_pml )( i+1, j, k )  - ( *Ex_pml )( i+1, j, k-1 ) + ( *Ex_pml )( i-1, j, k )-( *Ex_pml )( i-1, j, k-1 ) )
                                                - Bzy * ( ( *Ex_pml )( i, j+1, k )  - ( *Ex_pml )( i, j+1, k-1 ) + ( *Ex_pml )( i, j-1, k )-( *Ex_pml )( i, j-1, k-1 ) )
                                                -  Dz * ( ( *Ex_pml )( i, j, k+1 )  - ( *Ex_pml )( i, j, k-2) )
                                             ) ;
                    ( *Hy_pml )( i, j, k ) = + c3_d_yfield[i] * ( *Hy_pml )( i, j, k )
                                             + c4_d_yfield[i] * (c5_p_yfield[j]*( *By_pml )( i, j, k ) - c6_p_yfield[j]*By_pml_old ) ;
                }
            }
        }
        //Magnetic field Bz^(d,d,p) Remind that in PML, there no current
        for( unsigned int i=solvermin ; i<(unsigned int)solvermax ; i++ ) {
            for( unsigned int j=2 ; j<ny_d-2 ; j++ ) {
                for( unsigned int k=1 ; k<nz_p-1 ; k++ ) {
                    // Standard FDTD
                    // ( *Bz_pml )( i, j, k ) = + 1 * ( *Bz_pml )( i, j, k )
                    //                          + dt/dy * ( ( *Ex_pml )( i, j, k ) - ( *Ex_pml )( i, j-1, k ) )
                    //                          - dt/dx * ( ( *Ey_pml )( i, j, k ) - ( *Ey_pml )( i-1, j, k ) );
                    // NS FDTD
                    // ( *Bz_pml )( i, j, k ) = + 1 * ( *Bz_pml )( i, j, k ) + dt * (
                    //     + Ay *  ( ( *Ex_pml )( i, j, k )  -( *Ex_pml )( i, j-1, k ) )
                    //     + Byx * ( ( *Ex_pml )( i+1, j, k )-( *Ex_pml )( i+1, j-1, k ) + ( *Ex_pml )( i-1, j, k )-( *Ex_pml )( i-1, j-1, k ))
                    //     + Byz * ( ( *Ex_pml )( i, j, k+1 )-( *Ex_pml )( i, j-1, k+1 ) + ( *Ex_pml )( i, j, k-1 )-( *Ex_pml )( i, j-1, k-1 ))
                    //     + Dy  * ( ( *Ex_pml )( i, j+1, k )-( *Ex_pml )( i, j-2, k ) )
                    //     -  Ax * ( ( *Ey_pml )( i, j, k )  -( *Ey_pml )( i-1, j, k ) )
                    //     - Bxy * ( ( *Ey_pml )( i, j+1, k )-( *Ey_pml )( i-1, j+1, k ) + ( *Ey_pml )( i, j-1, k )-( *Ey_pml )( i-1, j-1, k ))
                    //     - Bxz * ( ( *Ey_pml )( i, j, k+1 )-( *Ey_pml )( i-1, j, k+1 ) + ( *Ey_pml )( i, j, k-1 )-( *Ey_pml )( i-1, j, k-1 ))
                    //     - Dx  * ( ( *Ey_pml )( i+1, j, k )-( *Ey_pml )( i-2, j, k ) )
                    // );
                    // ( *Hz_pml )( i, j, k ) = + 1 * ( *Bz_pml )( i, j, k );
                    // PML FDTD
                    Bz_pml_old = 1*( *Bz_pml )( i, j, k ) ;
                    ( *Bz_pml )( i, j, k ) = + c1_d_zfield[i] * ( *Bz_pml )( i, j, k )
                                             + c2_d_zfield[i] * (
                                                + Ay *  ( ( *Ex_pml )( i, j, k )  -( *Ex_pml )( i, j-1, k ) )
                                                + Byx * ( ( *Ex_pml )( i+1, j, k )-( *Ex_pml )( i+1, j-1, k ) + ( *Ex_pml )( i-1, j, k )-( *Ex_pml )( i-1, j-1, k ))
                                                + Byz * ( ( *Ex_pml )( i, j, k+1 )-( *Ex_pml )( i, j-1, k+1 ) + ( *Ex_pml )( i, j, k-1 )-( *Ex_pml )( i, j-1, k-1 ))
                                                + Dy  * ( ( *Ex_pml )( i, j+1, k )-( *Ex_pml )( i, j-2, k ) )
                                                -  Ax * ( ( *Ey_pml )( i, j, k )  -( *Ey_pml )( i-1, j, k ) )
                                                - Bxy * ( ( *Ey_pml )( i, j+1, k )-( *Ey_pml )( i-1, j+1, k ) + ( *Ey_pml )( i, j-1, k )-( *Ey_pml )( i-1, j-1, k ))
                                                - Bxz * ( ( *Ey_pml )( i, j, k+1 )-( *Ey_pml )( i-1, j, k+1 ) + ( *Ey_pml )( i, j, k-1 )-( *Ey_pml )( i-1, j, k-1 ))
                                                - Dx  * ( ( *Ey_pml )( i+1, j, k )-( *Ey_pml )( i-2, j, k ) )
                                             ) ;
                    ( *Hz_pml )( i, j, k ) = + c3_d_zfield[j] * ( *Hz_pml )( i, j, k )
                                             + c4_d_zfield[j] * (c5_p_zfield[k]*( *Bz_pml )( i, j, k ) - c6_p_zfield[k]*Bz_pml_old );
                }
            }
        }
    }
    else if (iDim==1) {
        //Magnetic field Bx^(p,d,d) Remind that in PML, there no current
        for( unsigned int i=1 ; i<nx_p-1 ; i++ ) {
            for( unsigned int j=solvermin ; j<(unsigned int)solvermax ; j++ ) {
                for( unsigned int k=2 ; k<nz_d-2 ; k++ ) {
                    // Standard FDTD
                    // ( *Bx_pml )( i, j, k ) = + 1 * ( *Bx_pml )( i, j, k )
                    //                          - dt/dy * ( ( *Ez_pml )( i, j, k ) - ( *Ez_pml )( i, j-1, k ) )
                    //                          + dt/dz * ( ( *Ey_pml )( i, j, k ) - ( *Ey_pml )( i, j, k-1 ) ) ;
                    // ( *Hx_pml )( i, j, k ) = + 1 * ( *Bx_pml )( i, j, k );
                    // PML FDTD
                    Bx_pml_old = 1*( *Bx_pml )( i, j, k ) ;
                    ( *Bx_pml )( i, j, k ) = + c1_d_xfield[j] * ( *Bx_pml )( i, j, k )
                                             + c2_d_xfield[j] * (
                                                + Az *  ( ( *Ey_pml )( i, j, k )    - ( *Ey_pml )( i, j, k-1 ) )
                                                + Bzx * ( ( *Ey_pml )( i+1, j, k )  - ( *Ey_pml )( i+1, j, k-1 ) + ( *Ey_pml )( i-1, j, k )-( *Ey_pml )( i-1, j, k-1 ) )
                                                + Bzy * ( ( *Ey_pml )( i, j+1, k )  - ( *Ey_pml )( i, j+1, k-1 ) + ( *Ey_pml )( i, j-1, k )-( *Ey_pml )( i, j-1, k-1 ) )
                                                +  Dz * ( ( *Ey_pml )( i, j, k+1 )  - ( *Ey_pml )( i, j, k-2) )
                                                -  Ay * ( ( *Ez_pml )( i,  j, k )   - ( *Ez_pml )( i,  j-1, k ) )
                                                - Byx * ( ( *Ez_pml )( i+1, j, k )  - ( *Ez_pml )( i+1, j-1, k ) + ( *Ez_pml )( i-1, j, k )-( *Ez_pml )( i-1, j-1, k ) )
                                                - Byz * ( ( *Ez_pml )( i, j, k+1 )  - ( *Ez_pml )( i, j-1, k+1 ) + ( *Ez_pml )( i, j, k-1 )-( *Ez_pml )( i, j-1, k-1 ) )
                                                - Dy *  ( ( *Ez_pml )( i, j+1, k )  - ( *Ez_pml )( i, j-2, k ) )
                                             ) ;
                    ( *Hx_pml )( i, j, k ) = + c3_d_xfield[k] * ( *Hx_pml )( i, j, k )
                                             + c4_d_xfield[k] * (c5_p_xfield[i]*( *Bx_pml )( i, j, k ) - c6_p_xfield[i]*Bx_pml_old ) ;
                }
            }
        }
        //Magnetic field By^(d,p,d) Remind that in PML, there no current
        for( unsigned int i=2 ; i<nx_d-2 ; i++ ) {
            for( unsigned int j=solvermin ; j<(unsigned int)solvermax ; j++ ) {
                for( unsigned int k=2 ; k<nz_d-2 ; k++ ) {
                    // Standard FDTD
                    // ( *By_pml )( i, j, k ) = + 1 * ( *By_pml )( i, j, k )
                    //                          + dt/dx * ( ( *Ez_pml )( i, j, k ) - ( *Ez_pml )( i-1, j, k ) )
                    //                          - dt/dz * ( ( *Ex_pml )( i, j, k ) - ( *Ex_pml )( i, j, k-1 ) );
                    // ( *Hy_pml )( i, j, k ) = + 1 * ( *By_pml )( i, j, k );
                    // PML FDTD
                    By_pml_old = 1*( *By_pml )( i, j, k ) ;
                    ( *By_pml )( i, j, k ) = + c1_d_yfield[k] * ( *By_pml )( i, j, k )
                                             + c2_d_yfield[k] * (
                                                + Ax *  ( ( *Ez_pml )( i,  j, k )   - ( *Ez_pml )( i-1, j, k ) )
                                                + Bxy * ( ( *Ez_pml )( i,  j+1, k ) - ( *Ez_pml )( i-1, j+1, k ) + ( *Ez_pml )( i, j-1, k )-( *Ez_pml )( i-1, j-1, k ) )
                                                + Bxz * ( ( *Ez_pml )( i, j, k+1 )  - ( *Ez_pml )( i-1, j, k+1 ) + ( *Ez_pml )( i, j, k-1 )-( *Ez_pml )( i-1, j, k-1 ) )
                                                +  Dx * ( ( *Ez_pml )( i+1, j, k )  - ( *Ez_pml )( i-2, j, k ) )
                                                -  Az * ( ( *Ex_pml )( i, j, k )    - ( *Ex_pml )( i, j, k-1 ) )
                                                - Bzx * ( ( *Ex_pml )( i+1, j, k )  - ( *Ex_pml )( i+1, j, k-1 ) + ( *Ex_pml )( i-1, j, k )-( *Ex_pml )( i-1, j, k-1 ) )
                                                - Bzy * ( ( *Ex_pml )( i, j+1, k )  - ( *Ex_pml )( i, j+1, k-1 ) + ( *Ex_pml )( i, j-1, k )-( *Ex_pml )( i, j-1, k-1 ) )
                                                -  Dz * ( ( *Ex_pml )( i, j, k+1 )  - ( *Ex_pml )( i, j, k-2) )
                                             ) ;
                    ( *Hy_pml )( i, j, k ) = + c3_d_yfield[i] * ( *Hy_pml )( i, j, k )
                                             + c4_d_yfield[i] * (c5_p_yfield[j]*( *By_pml )( i, j, k ) - c6_p_yfield[j]*By_pml_old ) ;
                }
            }
        }
        //Magnetic field Bz^(d,d,p) Remind that in PML, there no current
        for( unsigned int i=2 ; i<nx_d-2 ; i++ ) {
            for( unsigned int j=solvermin ; j<(unsigned int)solvermax ; j++ ) {
                for( unsigned int k=1 ; k<nz_p-1 ; k++ ) {
                    // Standard FDTD
                    // ( *Bz_pml )( i, j, k ) = + 1 * ( *Bz_pml )( i, j, k )
                    //                          + dt/dy * ( ( *Ex_pml )( i, j, k ) - ( *Ex_pml )( i, j-1, k ) )
                    //                          - dt/dx * ( ( *Ey_pml )( i, j, k ) - ( *Ey_pml )( i-1, j, k ) );
                    // ( *Hz_pml )( i, j, k ) = + 1 * ( *Bz_pml )( i, j, k );
                    // PML FDTD
                    Bz_pml_old = 1*( *Bz_pml )( i, j, k ) ;
                    ( *Bz_pml )( i, j, k ) = + c1_d_zfield[i] * ( *Bz_pml )( i, j, k )
                                             + c2_d_zfield[i] * (
                                                + Ay *  ( ( *Ex_pml )( i, j, k )  -( *Ex_pml )( i, j-1, k ) )
                                                + Byx * ( ( *Ex_pml )( i+1, j, k )-( *Ex_pml )( i+1, j-1, k ) + ( *Ex_pml )( i-1, j, k )-( *Ex_pml )( i-1, j-1, k ))
                                                + Byz * ( ( *Ex_pml )( i, j, k+1 )-( *Ex_pml )( i, j-1, k+1 ) + ( *Ex_pml )( i, j, k-1 )-( *Ex_pml )( i, j-1, k-1 ))
                                                + Dy  * ( ( *Ex_pml )( i, j+1, k )-( *Ex_pml )( i, j-2, k ) )
                                                -  Ax * ( ( *Ey_pml )( i, j, k )  -( *Ey_pml )( i-1, j, k ) )
                                                - Bxy * ( ( *Ey_pml )( i, j+1, k )-( *Ey_pml )( i-1, j+1, k ) + ( *Ey_pml )( i, j-1, k )-( *Ey_pml )( i-1, j-1, k ))
                                                - Bxz * ( ( *Ey_pml )( i, j, k+1 )-( *Ey_pml )( i-1, j, k+1 ) + ( *Ey_pml )( i, j, k-1 )-( *Ey_pml )( i-1, j, k-1 ))
                                                - Dx  * ( ( *Ey_pml )( i+1, j, k )-( *Ey_pml )( i-2, j, k ) )
                                             ) ;
                    ( *Hz_pml )( i, j, k ) = + c3_d_zfield[j] * ( *Hz_pml )( i, j, k )
                                             + c4_d_zfield[j] * (c5_p_zfield[k]*( *Bz_pml )( i, j, k ) - c6_p_zfield[k]*Bz_pml_old );
                }
            }
        }
    }
    else if (iDim==2) {
        //Magnetic field Bx^(p,d,d) Remind that in PML, there no current
        for( unsigned int i=1 ; i<nx_p-1 ; i++ ) {
            for( unsigned int j=2 ; j<ny_d-2 ; j++ ) {
                for( unsigned int k=solvermin ; k<(unsigned int)solvermax ; k++ ) {
                    // Standard FDTD
                    // ( *Bx_pml )( i, j, k ) = + 1 * ( *Bx_pml )( i, j, k )
                    //                          - dt/dy * ( ( *Ez_pml )( i, j, k ) - ( *Ez_pml )( i, j-1, k ) )
                    //                          + dt/dz * ( ( *Ey_pml )( i, j, k ) - ( *Ey_pml )( i, j, k-1 ) ) ;
                    // ( *Hx_pml )( i, j, k ) = + 1 * ( *Bx_pml )( i, j, k );
                    // PML FDTD
                    Bx_pml_old = 1*( *Bx_pml )( i, j, k ) ;
                    ( *Bx_pml )( i, j, k ) = + c1_d_xfield[j] * ( *Bx_pml )( i, j, k )
                                             + c2_d_xfield[j] * (
                                                + Az *  ( ( *Ey_pml )( i, j, k )    - ( *Ey_pml )( i, j, k-1 ) )
                                                + Bzx * ( ( *Ey_pml )( i+1, j, k )  - ( *Ey_pml )( i+1, j, k-1 ) + ( *Ey_pml )( i-1, j, k )-( *Ey_pml )( i-1, j, k-1 ) )
                                                + Bzy * ( ( *Ey_pml )( i, j+1, k )  - ( *Ey_pml )( i, j+1, k-1 ) + ( *Ey_pml )( i, j-1, k )-( *Ey_pml )( i, j-1, k-1 ) )
                                                +  Dz * ( ( *Ey_pml )( i, j, k+1 )  - ( *Ey_pml )( i, j, k-2) )
                                                -  Ay * ( ( *Ez_pml )( i,  j, k )   - ( *Ez_pml )( i,  j-1, k ) )
                                                - Byx * ( ( *Ez_pml )( i+1, j, k )  - ( *Ez_pml )( i+1, j-1, k ) + ( *Ez_pml )( i-1, j, k )-( *Ez_pml )( i-1, j-1, k ) )
                                                - Byz * ( ( *Ez_pml )( i, j, k+1 )  - ( *Ez_pml )( i, j-1, k+1 ) + ( *Ez_pml )( i, j, k-1 )-( *Ez_pml )( i, j-1, k-1 ) )
                                                - Dy *  ( ( *Ez_pml )( i, j+1, k )  - ( *Ez_pml )( i, j-2, k ) )
                                             ) ;
                    ( *Hx_pml )( i, j, k ) = + c3_d_xfield[k] * ( *Hx_pml )( i, j, k )
                                             + c4_d_xfield[k] * (c5_p_xfield[i]*( *Bx_pml )( i, j, k ) - c6_p_xfield[i]*Bx_pml_old ) ;
                }
            }
        }
        //Magnetic field By^(d,p,d) Remind that in PML, there no current
        for( unsigned int i=2 ; i<nx_d-2 ; i++ ) {
            for( unsigned int j=1 ; j<ny_p-1 ; j++ ) {
                for( unsigned int k=solvermin ; k<(unsigned int)solvermax ; k++ ) {
                    // Standard FDTD
                    // ( *By_pml )( i, j, k ) = + 1 * ( *By_pml )( i, j, k )
                    //                          + dt/dx * ( ( *Ez_pml )( i, j, k ) - ( *Ez_pml )( i-1, j, k ) )
                    //                          - dt/dz * ( ( *Ex_pml )( i, j, k ) - ( *Ex_pml )( i, j, k-1 ) );
                    // ( *Hy_pml )( i, j, k ) = + 1 * ( *By_pml )( i, j, k );
                    // PML FDTD
                    By_pml_old = 1*( *By_pml )( i, j, k ) ;
                    ( *By_pml )( i, j, k ) = + c1_d_yfield[k] * ( *By_pml )( i, j, k )
                                             + c2_d_yfield[k] * (
                                                + Ax *  ( ( *Ez_pml )( i,  j, k )   - ( *Ez_pml )( i-1, j, k ) )
                                                + Bxy * ( ( *Ez_pml )( i,  j+1, k ) - ( *Ez_pml )( i-1, j+1, k ) + ( *Ez_pml )( i, j-1, k )-( *Ez_pml )( i-1, j-1, k ) )
                                                + Bxz * ( ( *Ez_pml )( i, j, k+1 )  - ( *Ez_pml )( i-1, j, k+1 ) + ( *Ez_pml )( i, j, k-1 )-( *Ez_pml )( i-1, j, k-1 ) )
                                                +  Dx * ( ( *Ez_pml )( i+1, j, k )  - ( *Ez_pml )( i-2, j, k ) )
                                                -  Az * ( ( *Ex_pml )( i, j, k )    - ( *Ex_pml )( i, j, k-1 ) )
                                                - Bzx * ( ( *Ex_pml )( i+1, j, k )  - ( *Ex_pml )( i+1, j, k-1 ) + ( *Ex_pml )( i-1, j, k )-( *Ex_pml )( i-1, j, k-1 ) )
                                                - Bzy * ( ( *Ex_pml )( i, j+1, k )  - ( *Ex_pml )( i, j+1, k-1 ) + ( *Ex_pml )( i, j-1, k )-( *Ex_pml )( i, j-1, k-1 ) )
                                                -  Dz * ( ( *Ex_pml )( i, j, k+1 )  - ( *Ex_pml )( i, j, k-2) )
                                             ) ;
                    ( *Hy_pml )( i, j, k ) = + c3_d_yfield[i] * ( *Hy_pml )( i, j, k )
                                             + c4_d_yfield[i] * (c5_p_yfield[j]*( *By_pml )( i, j, k ) - c6_p_yfield[j]*By_pml_old ) ;
                }
            }
        }
        //Magnetic field Bz^(d,d,p) Remind that in PML, there no current
        for( unsigned int i=2 ; i<nx_d-2 ; i++ ) {
            for( unsigned int j=2 ; j<ny_d-2 ; j++ ) {
                for( unsigned int k=solvermin ; k<(unsigned int)solvermax ; k++ ) {
                    // Standard FDTD
                    // ( *Bz_pml )( i, j, k ) = + 1 * ( *Bz_pml )( i, j, k )
                    //                          + dt/dy * ( ( *Ex_pml )( i, j, k ) - ( *Ex_pml )( i, j-1, k ) )
                    //                          - dt/dx * ( ( *Ey_pml )( i, j, k ) - ( *Ey_pml )( i-1, j, k ) );
                    // ( *Hz_pml )( i, j, k ) = + 1 * ( *Bz_pml )( i, j, k );
                    // PML FDTD
                    Bz_pml_old = 1*( *Bz_pml )( i, j, k ) ;
                    ( *Bz_pml )( i, j, k ) = + c1_d_zfield[i] * ( *Bz_pml )( i, j, k )
                                             + c2_d_zfield[i] * (
                                                + Ay *  ( ( *Ex_pml )( i, j, k )  -( *Ex_pml )( i, j-1, k ) )
                                                + Byx * ( ( *Ex_pml )( i+1, j, k )-( *Ex_pml )( i+1, j-1, k ) + ( *Ex_pml )( i-1, j, k )-( *Ex_pml )( i-1, j-1, k ))
                                                + Byz * ( ( *Ex_pml )( i, j, k+1 )-( *Ex_pml )( i, j-1, k+1 ) + ( *Ex_pml )( i, j, k-1 )-( *Ex_pml )( i, j-1, k-1 ))
                                                + Dy  * ( ( *Ex_pml )( i, j+1, k )-( *Ex_pml )( i, j-2, k ) )
                                                -  Ax * ( ( *Ey_pml )( i, j, k )  -( *Ey_pml )( i-1, j, k ) )
                                                - Bxy * ( ( *Ey_pml )( i, j+1, k )-( *Ey_pml )( i-1, j+1, k ) + ( *Ey_pml )( i, j-1, k )-( *Ey_pml )( i-1, j-1, k ))
                                                - Bxz * ( ( *Ey_pml )( i, j, k+1 )-( *Ey_pml )( i-1, j, k+1 ) + ( *Ey_pml )( i, j, k-1 )-( *Ey_pml )( i-1, j, k-1 ))
                                                - Dx  * ( ( *Ey_pml )( i+1, j, k )-( *Ey_pml )( i-2, j, k ) )
                                             ) ;
                    ( *Hz_pml )( i, j, k ) = + c3_d_zfield[j] * ( *Hz_pml )( i, j, k )
                                             + c4_d_zfield[j] * (c5_p_zfield[k]*( *Bz_pml )( i, j, k ) - c6_p_zfield[k]*Bz_pml_old );
                }
            }
        }
    }
}
