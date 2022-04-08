#include "PML_Solver3D_Envelope.h"
#include "ElectroMagn.h"
#include "LaserEnvelope.h"
#include "EnvelopeBC3D_PML.h"
#include "Field3D.h"
#include "cField3D.h"
#include "Patch.h"
#include <complex>
#include <algorithm>

PML_Solver3D_Envelope::PML_Solver3D_Envelope( Params &params )
    : Solver3D( params )
{
    // // X-PML
    // kappa_x_max = 1.0 ;
    // sigma_x_max = 0.0 ; // 1.16 for 20 cells ; // 1.36 for 10 cells ;
    // alpha_x_max = 0.0 ;
    // power_pml_kappa_x = 3.;
    // power_pml_sigma_x = 2.;
    // power_pml_alpha_x = 1.;
    // alpha_cx = 0. ; // Try to use a more practical timestep !
    // // Y-PML
    // kappa_y_max = 1.0 ;
    // sigma_y_max = 0.0 ; // 2.32 ; // 2.32 for 20 cells ; // 2.72 for 10 cells ;
    // alpha_y_max = 0.0 ;
    // power_pml_kappa_y = 3.;
    // power_pml_sigma_y = 2.;
    // power_pml_alpha_y = 1.;
    // alpha_cy = 0. ; // 0.25 ; // 0.8 for transverse is ok
    // // Z-PML
    // kappa_z_max = 1.0 ;
    // sigma_z_max = 0.0 ; // 2.32 ; // 2.32 for 20 cells ; // 2.72 for 10 cells ;
    // alpha_z_max = 0.0 ;
    // power_pml_kappa_z = 3.;
    // power_pml_sigma_z = 2.;
    // power_pml_alpha_z = 1.;
    // alpha_cz = 0. ; // 0.25 ; // 0.8 for transverse is ok

    // X-PML
    kappa_x_max = 1.0 ;
    sigma_x_max = 1.36 ; // 1.16 for 20 cells ; // 1.36 for 10 cells ;
    alpha_x_max = 0.0 ;
    power_pml_kappa_x = 3.;
    power_pml_sigma_x = 2.;
    power_pml_alpha_x = 1.;
    alpha_cx = 1.01 ; // Try to use a more practical timestep !
    // Y-PML
    kappa_y_max = 1. ;
    sigma_y_max = 1.8 ; // 2.32 ; // 2.32 for 20 cells ; // 2.72 for 10 cells ;
    alpha_y_max = 0.0 ;
    power_pml_kappa_y = 3.;
    power_pml_sigma_y = 2.;
    power_pml_alpha_y = 1.;
    alpha_cy = 0.10 ; // 0.25 ; // 0.8 for transverse is ok
    // Z-PML
    kappa_z_max = 1. ;
    sigma_z_max = 1.8 ; // 2.32 ; // 2.32 for 20 cells ; // 2.72 for 10 cells ;
    alpha_z_max = 0.0 ;
    power_pml_kappa_z = 3.;
    power_pml_sigma_z = 2.;
    power_pml_alpha_z = 1.;
    alpha_cz = 0.10 ; // 0.25 ; // 0.8 for transverse is ok
}

PML_Solver3D_Envelope::~PML_Solver3D_Envelope()
{
}

void PML_Solver3D_Envelope::operator()( ElectroMagn *fields )
{
    ERROR( "This is not a solver for the main domain" );

    //cField2D *A_n        = static_cast<cField2D *>( envelope->A_ );   // the envelope at timestep n
    //cField2D *A_nm1      = static_cast<cField2D *>( envelope->A0_ );  // the envelope at timestep n-1
}

void PML_Solver3D_Envelope::setDomainSizeAndCoefficients( int iDim, int min_or_max, int ncells_pml_domain, int startpml, int* ncells_pml_min, int* ncells_pml_max, Patch* patch )
{
    if ( iDim == 0 ) {
        nx_p = ncells_pml_domain;
    }
    else if ( iDim == 1 ) {
        ny_p = ncells_pml_domain;
        nx_p += ncells_pml_min[0]-1*(patch->isXmin()) + ncells_pml_max[0]-1*(patch->isXmax());
    }
    else if ( iDim == 2 ) {
        nz_p = ncells_pml_domain;
        nx_p += ncells_pml_min[0]-1*(patch->isXmin()) + ncells_pml_max[0]-1*(patch->isXmax());
        ny_p += ncells_pml_min[1]-1*(patch->isYmin()) + ncells_pml_max[1]-1*(patch->isYmax());
    }

    //PML Coeffs Kappa,Sigma ...
    //Primal
    kappa_x_p.resize( nx_p );
    sigma_x_p.resize( nx_p );
    alpha_x_p.resize( nx_p );
    kappa_prime_x_p.resize( nx_p );
    sigma_prime_x_p.resize( nx_p );
    alpha_prime_x_p.resize( nx_p );
    kappa_y_p.resize( ny_p );
    sigma_y_p.resize( ny_p );
    alpha_y_p.resize( ny_p );
    kappa_prime_y_p.resize( ny_p );
    sigma_prime_y_p.resize( ny_p );
    alpha_prime_y_p.resize( ny_p );
    kappa_z_p.resize( nz_p );
    sigma_z_p.resize( nz_p );
    alpha_z_p.resize( nz_p );
    kappa_prime_z_p.resize( nz_p );
    sigma_prime_z_p.resize( nz_p );
    alpha_prime_z_p.resize( nz_p );

    if ( iDim == 0 ) {
        // 3 cells (oversize) are vaccum so the PML media begin at r0 which is :
        // Eventually the size of PML media is :
        length_z_pml = 0. ;
        length_y_pml = 0. ;
        length_x_pml = (ncells_pml_domain-startpml+0.5)*dx ;
        // Primal grid
        // X-direction
        // Params for first cell of PML-patch (vacuum) i = 0,1,2
        for ( int i=0 ; i<startpml ; i++ ) {
            // Coeffs for the first cell
            kappa_x_p[i] = 1. ;
            sigma_x_p[i] = 0. ;
            alpha_x_p[i] = 0. ;
            kappa_prime_x_p[i] = 0. ;
            sigma_prime_x_p[i] = 0. ;
            alpha_prime_x_p[i] = 0. ;
        }
        for ( int i=startpml; i<nx_p ; i++ ) {
            // Parameters
            kappa_x_p[i] = 1. - (kappa_x_max - 1.) * pow( (i-startpml)*dx , power_pml_kappa_x ) / pow( length_x_pml , power_pml_kappa_x ) ;
            sigma_x_p[i] = sigma_x_max * pow( (i-startpml)*dx , power_pml_sigma_x ) / pow( length_x_pml , power_pml_sigma_x ) ;
            alpha_x_p[i] = alpha_cx + alpha_x_max * (1. - pow( (i-startpml)*dx , power_pml_alpha_x ) / pow( length_x_pml , power_pml_alpha_x ) );
            // Derivatives
            kappa_prime_x_p[i] = -(kappa_x_max - 1.) * power_pml_kappa_x * pow( (i-startpml)*dx , power_pml_kappa_x-1 ) / pow( length_x_pml , power_pml_kappa_x ) ;
            sigma_prime_x_p[i] = sigma_x_max * power_pml_sigma_x * pow( (i-startpml)*dx , power_pml_sigma_x-1 ) / pow( length_x_pml , power_pml_sigma_x ) ;
            alpha_prime_x_p[i] = -alpha_x_max * power_pml_alpha_x * pow( (i-startpml)*dx , power_pml_alpha_x-1 ) / pow( length_x_pml , power_pml_alpha_x ) ;
        }
        if (min_or_max==0) {
            std::reverse(kappa_x_p.begin(), kappa_x_p.end());
            std::reverse(sigma_x_p.begin(), sigma_x_p.end());
            std::reverse(alpha_x_p.begin(), alpha_x_p.end());
            std::reverse(kappa_prime_x_p.begin(), kappa_prime_x_p.end());
            std::reverse(sigma_prime_x_p.begin(), sigma_prime_x_p.end());
            std::reverse(alpha_prime_x_p.begin(), alpha_prime_x_p.end());
            for (int i=0 ; i<nx_p ; i++){
                // Due to SMILEI convention for propagating wave
                kappa_x_p[i] *= +1;
                sigma_x_p[i] *= -1;
                alpha_x_p[i] *= -1;
                // Due to SMILEI convention for propagating wave
                // For the envelope it's not a good solution for the min value !
                // kappa_prime_x_p[i] *= -1.;
                // sigma_prime_x_p[i] *= +1.;
                // alpha_prime_x_p[i] *= +1.;
                // Upper solution make xmin diverge quickly
                kappa_prime_x_p[i] *= +1;
                sigma_prime_x_p[i] *= -1;
                alpha_prime_x_p[i] *= -1;
            }
        }
        if (min_or_max==1) {
            for (int i=0 ; i<nx_p ; i++){
                // Due to SMILEI convention for propagating wave
                kappa_x_p[i] *= +1;
                sigma_x_p[i] *= -1;
                alpha_x_p[i] *= -1;
                // Due to SMILEI convention for propagating wave
                kappa_prime_x_p[i] *= +1;
                sigma_prime_x_p[i] *= -1;
                alpha_prime_x_p[i] *= -1;
            }
        }
        // Y-direction
        for ( int j=0 ; j<ny_p ; j++ ) {
            kappa_y_p[j] = 1. ;
            sigma_y_p[j] = 0. ;
        }
        // Z-direction
        for ( int k=0 ; k<nz_p ; k++ ) {
            kappa_z_p[k] = 1. ;
            sigma_z_p[k] = 0. ;
        }
    }
    if ( iDim == 1 ) {
        // 3 cells are vaccum so the PML media begin at r0 which is :
        // Eventually the size of PML media is :
        length_y_pml = (ncells_pml_domain-startpml+0.5)*dy ;
        length_x_pml_xmax = (ncells_pml_max[0]+0.5)*dx ;
        length_x_pml_xmin = (ncells_pml_min[0]+0.5)*dx ;
        for ( int i=0 ; i<nx_p ; i++ ) {
            kappa_x_p[i] = 1. ;
            sigma_x_p[i] = 0. ;
            alpha_x_p[i] = 0. ;
            kappa_prime_x_p[i] = 0. ;
            sigma_prime_x_p[i] = 0. ;
            alpha_prime_x_p[i] = 0. ;
        }
        if (ncells_pml_min[0] != 0 ){
            for ( int i=0 ; i<ncells_pml_min[0] ; i++ ) {
                // Parameters
                kappa_x_p[i] = 1. - (kappa_x_max - 1.) * pow( ( ncells_pml_min[0] - 1 - i )*dx , power_pml_kappa_x ) / pow( length_x_pml_xmin , power_pml_kappa_x ) ;
                sigma_x_p[i] = sigma_x_max * pow( ( ncells_pml_min[0] - 1 - i )*dx , power_pml_sigma_x ) / pow( length_x_pml_xmin , power_pml_sigma_x ) ;
                alpha_x_p[i] = alpha_cx + alpha_x_max * (1. - pow( ( ncells_pml_min[0] - 1 - i )*dx , power_pml_alpha_x ) / pow( length_x_pml_xmin , power_pml_alpha_x ) );
                // Derivatives
                kappa_prime_x_p[i] = - (kappa_x_max - 1.) * power_pml_kappa_x * pow( ( ncells_pml_min[0] - 1 - i )*dx , power_pml_kappa_x-1 ) / pow( length_x_pml_xmin , power_pml_kappa_x ) ;
                sigma_prime_x_p[i] = sigma_x_max * power_pml_sigma_x * pow( ( ncells_pml_min[0] - 1 - i )*dx , power_pml_sigma_x-1 ) / pow( length_x_pml_xmin , power_pml_sigma_x ) ;
                alpha_prime_x_p[i] = -alpha_x_max * power_pml_alpha_x * pow( ( ncells_pml_min[0] - 1 - i )*dx , power_pml_alpha_x-1 ) / pow( length_x_pml_xmin , power_pml_alpha_x ) ;
                // Convention Envelop Smilei
                kappa_x_p[i] *= +1 ;
                sigma_x_p[i] *= -1 ;
                alpha_x_p[i] *= -1 ;
                // kappa_prime_x_p[i] *= -1 ;
                // sigma_prime_x_p[i] *= +1 ;
                // alpha_prime_x_p[i] *= +1 ;
                kappa_prime_x_p[i] *= +1 ;
                sigma_prime_x_p[i] *= -1 ;
                alpha_prime_x_p[i] *= -1 ;
            }
        }
        if (ncells_pml_max[0] != 0 ){
            for ( int i=(nx_p-1)-(ncells_pml_max[0]-1) ; i<nx_p ; i++ ) { // La aussi, il y a 2 cellules de trop pour les pml xmax avec 1 seul patch
                // Parameters
                kappa_x_p[i] = 1. - (kappa_x_max - 1.) * pow( ( i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) )*dx , power_pml_kappa_x ) / pow( length_x_pml_xmax , power_pml_kappa_x ) ;
                sigma_x_p[i] = sigma_x_max * pow( (i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) )*dx , power_pml_sigma_x ) / pow( length_x_pml_xmax, power_pml_sigma_x ) ;
                alpha_x_p[i] = alpha_cx + alpha_x_max * (1. - pow( ( i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) )*dx , power_pml_alpha_x ) / pow( length_x_pml_xmin , power_pml_alpha_x ) );
                // Derivatives
                kappa_prime_x_p[i] = - (kappa_x_max - 1.) * power_pml_kappa_x * pow( ( i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) )*dx , power_pml_kappa_x-1 ) / pow( length_x_pml_xmax , power_pml_kappa_x ) ;
                sigma_prime_x_p[i] = sigma_x_max * power_pml_sigma_x * pow( ( i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) )*dx , power_pml_sigma_x-1 ) / pow( length_x_pml_xmax , power_pml_sigma_x ) ;
                alpha_prime_x_p[i] = -alpha_x_max * power_pml_alpha_x * pow( ( i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) )*dx , power_pml_alpha_x-1 ) / pow( length_x_pml_xmax , power_pml_alpha_x ) ;
                // Convention Envelop Smilei
                kappa_x_p[i] *= +1 ;
                sigma_x_p[i] *= -1 ;
                alpha_x_p[i] *= -1 ;
                kappa_prime_x_p[i] *= +1 ;
                sigma_prime_x_p[i] *= -1 ;
                alpha_prime_x_p[i] *= -1 ;
            }
        }
        // Y-direction
        for ( int j=0 ; j<startpml ; j++ ) {
            // Coeffs for the first cell
            kappa_y_p[j] = 1. ;
            sigma_y_p[j] = 0. ;
            alpha_y_p[j] = 0. ;
            kappa_prime_y_p[j] = 0. ;
            sigma_prime_y_p[j] = 0. ;
            alpha_prime_y_p[j] = 0. ;
        }
        // Params for other cells (PML Media) when i>=3
        for ( int j=startpml; j<ny_p ; j++ ) {
            // Parameters
            kappa_y_p[j] = 1. + (kappa_y_max - 1.) * pow( (j-startpml)*dy , power_pml_kappa_y ) / pow( length_y_pml , power_pml_kappa_y ) ;
            sigma_y_p[j] = sigma_y_max * pow( (j-startpml)*dy , power_pml_sigma_y ) / pow( length_y_pml , power_pml_sigma_y ) ;
            alpha_y_p[j] = alpha_cy + alpha_y_max * (1. - pow( (j-startpml)*dy , power_pml_alpha_y ) / pow( length_y_pml , power_pml_alpha_y ) );
            // Derivatives
            kappa_prime_y_p[j] = (kappa_y_max - 1.) * power_pml_kappa_y * pow( (j-startpml)*dy , power_pml_kappa_y-1 ) / pow( length_y_pml , power_pml_kappa_y ) ;
            sigma_prime_y_p[j] = sigma_y_max * power_pml_sigma_y * pow( (j-startpml)*dy , power_pml_sigma_y-1 ) / pow( length_y_pml , power_pml_sigma_y ) ;
            alpha_prime_y_p[j] = -alpha_y_max * power_pml_alpha_y * pow( (j-startpml)*dy , power_pml_alpha_y-1 ) / pow( length_y_pml , power_pml_alpha_y ) ;
        }
        if (min_or_max==0) {
            std::reverse(kappa_y_p.begin(), kappa_y_p.end());
            std::reverse(sigma_y_p.begin(), sigma_y_p.end());
            std::reverse(alpha_y_p.begin(), alpha_y_p.end());
            std::reverse(kappa_prime_y_p.begin(), kappa_prime_y_p.end());
            std::reverse(sigma_prime_y_p.begin(), sigma_prime_y_p.end());
            std::reverse(alpha_prime_y_p.begin(), alpha_prime_y_p.end());
            for (int j=0 ; j<ny_p ; j++){
                // Due to SMILEI convention for propagating wave
                kappa_y_p[j] *= +1;
                sigma_y_p[j] *= -1;
                alpha_y_p[j] *= -1;
                // Due to SMILEI convention for propagating wave
                kappa_prime_y_p[j] *= -1.;
                sigma_prime_y_p[j] *= +1.;
                alpha_prime_y_p[j] *= +1.;
                // kappa_prime_y_p[j] *= +1.;
                // sigma_prime_y_p[j] *= -1.;
                // alpha_prime_y_p[j] *= -1.;
            }
        }
        if (min_or_max==1) {
            for (int j=0 ; j<ny_p ; j++){
                // Due to SMILEI convention for propagating wave
                kappa_y_p[j] *= +1;
                sigma_y_p[j] *= -1;
                alpha_y_p[j] *= -1;
                // Due to SMILEI convention for propagating wave
                kappa_prime_y_p[j] *= +1;
                sigma_prime_y_p[j] *= -1;
                alpha_prime_y_p[j] *= -1;
            }
        }
        // Z-direction
        for ( int k=0 ; k<nz_p ; k++ ) {
            kappa_z_p[k] = 1. ;
            sigma_z_p[k] = 0. ;
        }
    }
    if ( iDim == 2 ) {
        // 3 cells are vaccum so the PML media begin at r0 which is :
        // Eventually the size of PML media is :
        length_z_pml = (ncells_pml_domain-startpml+0.5)*dz ;
        length_x_pml_xmax = (ncells_pml_max[0]+0.5)*dx ;
        length_x_pml_xmin = (ncells_pml_min[0]+0.5)*dx ;
        length_y_pml_ymax = (ncells_pml_max[1]+0.5)*dy ;
        length_y_pml_ymin = (ncells_pml_min[1]+0.5)*dy ;
        for ( int i=0 ; i<nx_p ; i++ ) {
            kappa_x_p[i] = 1. ;
            sigma_x_p[i] = 0. ;
            alpha_x_p[i] = 0. ;
            kappa_prime_x_p[i] = 0. ;
            sigma_prime_x_p[i] = 0. ;
            alpha_prime_x_p[i] = 0. ;
        }
        if (ncells_pml_min[0] != 0 ){
            for ( int i=0 ; i<ncells_pml_min[0] ; i++ ) {
                // Parameters
                kappa_x_p[i] = 1. - (kappa_x_max - 1.) * pow( ( ncells_pml_min[0] - 1 - i )*dx , power_pml_kappa_x ) / pow( length_x_pml_xmin , power_pml_kappa_x ) ;
                sigma_x_p[i] = sigma_x_max * pow( ( ncells_pml_min[0] - 1 - i )*dx , power_pml_sigma_x ) / pow( length_x_pml_xmin , power_pml_sigma_x ) ;
                alpha_x_p[i] = alpha_cx + alpha_x_max * (1. - pow( ( ncells_pml_min[0] - 1 - i )*dx , power_pml_alpha_x ) / pow( length_x_pml_xmin , power_pml_alpha_x ) );
                // Derivatives
                kappa_prime_x_p[i] = - (kappa_x_max - 1.) * power_pml_kappa_x * pow( ( ncells_pml_min[0] - 1 - i )*dx , power_pml_kappa_x-1 ) / pow( length_x_pml_xmin , power_pml_kappa_x ) ;
                sigma_prime_x_p[i] = sigma_x_max * power_pml_sigma_x * pow( ( ncells_pml_min[0] - 1 - i )*dx , power_pml_sigma_x-1 ) / pow( length_x_pml_xmin , power_pml_sigma_x ) ;
                alpha_prime_x_p[i] = -alpha_x_max * power_pml_alpha_x * pow( ( ncells_pml_min[0] - 1 - i )*dx , power_pml_alpha_x-1 ) / pow( length_x_pml_xmin , power_pml_alpha_x ) ;
                // Convention Envelop Smilei
                kappa_x_p[i] *= +1 ;
                sigma_x_p[i] *= -1 ;
                alpha_x_p[i] *= -1 ;
                // kappa_prime_x_p[i] *= -1 ;
                // sigma_prime_x_p[i] *= +1 ;
                // alpha_prime_x_p[i] *= +1 ;
                kappa_prime_x_p[i] *= +1 ;
                sigma_prime_x_p[i] *= -1 ;
                alpha_prime_x_p[i] *= -1 ;
            }
        }
        if (ncells_pml_max[0] != 0 ){
            for ( int i=(nx_p-1)-(ncells_pml_max[0]-1) ; i<nx_p ; i++ ) { // La aussi, il y a 2 cellules de trop pour les pml xmax avec 1 seul patch
                // Parameters
                kappa_x_p[i] = 1. - (kappa_x_max - 1.) * pow( ( i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) )*dx , power_pml_kappa_x ) / pow( length_x_pml_xmax , power_pml_kappa_x ) ;
                sigma_x_p[i] = sigma_x_max * pow( (i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) )*dx , power_pml_sigma_x ) / pow( length_x_pml_xmax, power_pml_sigma_x ) ;
                alpha_x_p[i] = alpha_cx + alpha_x_max * (1. - pow( ( i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) )*dx , power_pml_alpha_x ) / pow( length_x_pml_xmin , power_pml_alpha_x ) );
                // Derivatives
                kappa_prime_x_p[i] = - (kappa_x_max - 1.) * power_pml_kappa_x * pow( ( i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) )*dx , power_pml_kappa_x-1 ) / pow( length_x_pml_xmax , power_pml_kappa_x ) ;
                sigma_prime_x_p[i] = sigma_x_max * power_pml_sigma_x * pow( ( i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) )*dx , power_pml_sigma_x-1 ) / pow( length_x_pml_xmax , power_pml_sigma_x ) ;
                alpha_prime_x_p[i] = -alpha_x_max * power_pml_alpha_x * pow( ( i - ( (nx_p-1)-(ncells_pml_max[0]-1) ) )*dx , power_pml_alpha_x-1 ) / pow( length_x_pml_xmax , power_pml_alpha_x ) ;
                // Convention Envelop Smilei
                kappa_x_p[i] *= +1 ;
                sigma_x_p[i] *= -1 ;
                alpha_x_p[i] *= -1 ;
                kappa_prime_x_p[i] *= +1 ;
                sigma_prime_x_p[i] *= -1 ;
                alpha_prime_x_p[i] *= -1 ;
            }
        }
        // Y-direction
        for ( int j=0 ; j<ny_p ; j++ ) {
            // Coeffs for the first cell
            kappa_y_p[j] = 1. ;
            sigma_y_p[j] = 0. ;
            alpha_y_p[j] = 0. ;
            kappa_prime_y_p[j] = 0. ;
            sigma_prime_y_p[j] = 0. ;
            alpha_prime_y_p[j] = 0. ;
        }
        //Params for other cells (PML Media) when i>=3
        if (ncells_pml_min[1] != 0 ){
            for ( int j=0 ; j<ncells_pml_min[1] ; j++ ) {
                // Parameters
                kappa_y_p[j] = 1. - (kappa_y_max - 1.) * pow( ( ncells_pml_min[1] - 1 - j )*dy , power_pml_kappa_y ) / pow( length_y_pml_ymin , power_pml_kappa_y ) ;
                sigma_y_p[j] = sigma_y_max * pow( ( ncells_pml_min[1] - 1 - j )*dy , power_pml_sigma_y ) / pow( length_y_pml_ymin , power_pml_sigma_y ) ;
                alpha_y_p[j] = alpha_cy + alpha_y_max * (1. - pow( ( ncells_pml_min[1] - 1 - j )*dy , power_pml_alpha_y ) / pow( length_y_pml_ymin , power_pml_alpha_y ) );
                // Derivatives
                kappa_prime_y_p[j] = - (kappa_y_max - 1.) * power_pml_kappa_y * pow( ( ncells_pml_min[1] - 1 - j )*dy , power_pml_kappa_y-1 ) / pow( length_y_pml_ymin , power_pml_kappa_y ) ;
                sigma_prime_y_p[j] = sigma_y_max * power_pml_sigma_y * pow( ( ncells_pml_min[1] - 1 - j )*dy , power_pml_sigma_y-1 ) / pow( length_y_pml_ymin , power_pml_sigma_y ) ;
                alpha_prime_y_p[j] = -alpha_y_max * power_pml_alpha_y * pow( ( ncells_pml_min[1] - 1 - j )*dy , power_pml_alpha_y-1 ) / pow( length_y_pml_ymin , power_pml_alpha_y ) ;
                // Convention Envelop Smilei
                kappa_y_p[j] *= +1 ;
                sigma_y_p[j] *= -1 ;
                alpha_y_p[j] *= -1 ;
                // kappa_prime_y_p[j] *= -1 ;
                // sigma_prime_y_p[j] *= +1 ;
                // alpha_prime_y_p[j] *= +1 ;
                kappa_prime_y_p[j] *= +1 ;
                sigma_prime_y_p[j] *= -1 ;
                alpha_prime_y_p[j] *= -1 ;
            }
        }
        if (ncells_pml_max[1] != 0 ){
            for ( int j=(ny_p-1)-(ncells_pml_max[1]-1) ; j<ny_p ; j++ ) { // La aussi, il y a 2 cellules de trop pour les pml xmax avec 1 seul patch
                // Parameters
                kappa_y_p[j] = 1. - (kappa_y_max - 1.) * pow( ( j - ( (ny_p-1)-(ncells_pml_max[1]-1) ) )*dy , power_pml_kappa_y ) / pow( length_y_pml_ymax , power_pml_kappa_y ) ;
                sigma_y_p[j] = sigma_y_max * pow( (j - ( (ny_p-1)-(ncells_pml_max[1]-1) ) )*dy , power_pml_sigma_y ) / pow( length_y_pml_ymax, power_pml_sigma_y ) ;
                alpha_y_p[j] = alpha_cy + alpha_y_max * (1. - pow( ( j - ( (ny_p-1)-(ncells_pml_max[1]-1) ) )*dy , power_pml_alpha_y ) / pow( length_y_pml_ymin , power_pml_alpha_y ) );
                // Derivatives
                kappa_prime_y_p[j] = - (kappa_y_max - 1.) * power_pml_kappa_y * pow( ( j - ( (ny_p-1)-(ncells_pml_max[1]-1) ) )*dy , power_pml_kappa_y-1 ) / pow( length_y_pml_ymax , power_pml_kappa_y ) ;
                sigma_prime_y_p[j] = sigma_y_max * power_pml_sigma_y * pow( ( j - ( (ny_p-1)-(ncells_pml_max[1]-1) ) )*dy , power_pml_sigma_y-1 ) / pow( length_y_pml_ymax , power_pml_sigma_y ) ;
                alpha_prime_y_p[j] = -alpha_y_max * power_pml_alpha_y * pow( ( j - ( (ny_p-1)-(ncells_pml_max[1]-1) ) )*dy , power_pml_alpha_y-1 ) / pow( length_y_pml_ymax , power_pml_alpha_y ) ;
                // Convention Envelop Smilei
                kappa_y_p[j] *= +1 ;
                sigma_y_p[j] *= -1 ;
                alpha_y_p[j] *= -1 ;
                kappa_prime_y_p[j] *= +1 ;
                sigma_prime_y_p[j] *= -1 ;
                alpha_prime_y_p[j] *= -1 ;
            }
        }
        // Z-direction
        for ( int k=0 ; k<startpml ; k++ ) {
            // Coeffs for the first cell
            kappa_z_p[k] = 1. ;
            sigma_z_p[k] = 0. ;
            alpha_z_p[k] = 0. ;
            kappa_prime_z_p[k] = 0. ;
            sigma_prime_z_p[k] = 0. ;
            alpha_prime_z_p[k] = 0. ;
        }
        // Params for other cells (PML Media) when i>=3
        for ( int k=startpml; k<nz_p ; k++ ) {
            // Parameters
            kappa_z_p[k] = 1. + (kappa_z_max - 1.) * pow( (k-startpml)*dz , power_pml_kappa_z ) / pow( length_z_pml , power_pml_kappa_z ) ;
            sigma_z_p[k] = sigma_z_max * pow( (k-startpml)*dz , power_pml_sigma_z ) / pow( length_z_pml , power_pml_sigma_z ) ;
            alpha_z_p[k] = alpha_cz + alpha_z_max * (1. - pow( (k-startpml)*dz , power_pml_alpha_z ) / pow( length_z_pml , power_pml_alpha_z ) );
            // Derivatives
            kappa_prime_z_p[k] = (kappa_z_max - 1.) * power_pml_kappa_z * pow( (k-startpml)*dz , power_pml_kappa_z-1 ) / pow( length_z_pml , power_pml_kappa_z ) ;
            sigma_prime_z_p[k] = sigma_z_max * power_pml_sigma_z * pow( (k-startpml)*dz , power_pml_sigma_z-1 ) / pow( length_z_pml , power_pml_sigma_z ) ;
            alpha_prime_z_p[k] = -alpha_z_max * power_pml_alpha_z * pow( (k-startpml)*dz , power_pml_alpha_z-1 ) / pow( length_z_pml , power_pml_alpha_z ) ;
        }
        if (min_or_max==0) {
            std::reverse(kappa_z_p.begin(), kappa_z_p.end());
            std::reverse(sigma_z_p.begin(), sigma_z_p.end());
            std::reverse(alpha_z_p.begin(), alpha_z_p.end());
            std::reverse(kappa_prime_z_p.begin(), kappa_prime_z_p.end());
            std::reverse(sigma_prime_z_p.begin(), sigma_prime_z_p.end());
            std::reverse(alpha_prime_z_p.begin(), alpha_prime_z_p.end());
            for (int k=0 ; k<nz_p ; k++){
                // Due to SMILEI convention for propagating wave
                kappa_z_p[k] *= +1;
                sigma_z_p[k] *= -1;
                alpha_z_p[k] *= -1;
                // Due to SMILEI convention for propagating wave
                kappa_prime_z_p[k] *= -1.;
                sigma_prime_z_p[k] *= +1.;
                alpha_prime_z_p[k] *= +1.;
                // kappa_prime_z_p[k] *= +1.;
                // sigma_prime_z_p[k] *= -1.;
                // alpha_prime_z_p[k] *= -1.;
            }
        }
        if (min_or_max==1) {
            for (int k=0 ; k<nz_p ; k++){
                // Due to SMILEI convention for propagating wave
                kappa_z_p[k] *= +1;
                sigma_z_p[k] *= -1;
                alpha_z_p[k] *= -1;
                // Due to SMILEI convention for propagating wave
                kappa_prime_z_p[k] *= +1;
                sigma_prime_z_p[k] *= -1;
                alpha_prime_z_p[k] *= -1;
            }
        }
    }
}

void PML_Solver3D_Envelope::compute_A_from_G( LaserEnvelope *envelope, int iDim, int min_or_max, int solvermin, int solvermax )
{
    EnvelopeBC3D_PML* pml_fields = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[iDim*2+min_or_max] );

    cField3D* A_nm1_pml = NULL;
    cField3D* u1_nm1_x_pml = NULL;
    cField3D* u2_nm1_x_pml = NULL;
    cField3D* u3_nm1_x_pml = NULL;
    cField3D* u1_nm1_y_pml = NULL;
    cField3D* u2_nm1_y_pml = NULL;
    cField3D* u3_nm1_y_pml = NULL;
    cField3D* u1_nm1_z_pml = NULL;
    cField3D* u2_nm1_z_pml = NULL;
    cField3D* u3_nm1_z_pml = NULL;

    cField3D* A_n_pml = NULL;

    cField3D* A_np1_pml = NULL;
    cField3D* u1_np1_x_pml = NULL;
    cField3D* u2_np1_x_pml = NULL;
    cField3D* u3_np1_x_pml = NULL;
    cField3D* u1_np1_y_pml = NULL;
    cField3D* u2_np1_y_pml = NULL;
    cField3D* u3_np1_y_pml = NULL;
    cField3D* u1_np1_z_pml = NULL;
    cField3D* u2_np1_z_pml = NULL;
    cField3D* u3_np1_z_pml = NULL;


    A_nm1_pml = pml_fields->A_nm1_;
    u1_nm1_x_pml = pml_fields->u1_nm1_x_;
    u2_nm1_x_pml = pml_fields->u2_nm1_x_;
    u3_nm1_x_pml = pml_fields->u3_nm1_x_;
    u1_nm1_y_pml = pml_fields->u1_nm1_y_;
    u2_nm1_y_pml = pml_fields->u2_nm1_y_;
    u3_nm1_y_pml = pml_fields->u3_nm1_y_;
    u1_nm1_z_pml = pml_fields->u1_nm1_z_;
    u2_nm1_z_pml = pml_fields->u2_nm1_z_;
    u3_nm1_z_pml = pml_fields->u3_nm1_z_;

    A_n_pml = pml_fields->A_n_;

    A_np1_pml = pml_fields->A_np1_;
    u1_np1_x_pml = pml_fields->u1_np1_x_;
    u2_np1_x_pml = pml_fields->u2_np1_x_;
    u3_np1_x_pml = pml_fields->u3_np1_x_;
    u1_np1_y_pml = pml_fields->u1_np1_y_;
    u2_np1_y_pml = pml_fields->u2_np1_y_;
    u3_np1_y_pml = pml_fields->u3_np1_y_;
    u1_np1_z_pml = pml_fields->u1_np1_z_;
    u2_np1_z_pml = pml_fields->u2_np1_z_;
    u3_np1_z_pml = pml_fields->u3_np1_z_;

    // Auxiliary Quantities
    std::complex<double> i1 = std::complex<double>( 0., 1. ); // imaginary unit
    double k0 = 1.; // laser wavenumber
    std::complex<double> source_term_x ;
    std::complex<double> source_term_y ;
    std::complex<double> source_term_z ;

    if (iDim == 0) {
        // A (p,p,p) Remind that in PML, there no current
        for( unsigned int i=solvermin ; i<solvermax; i++ ) { // x loop
            for( unsigned int j=1 ; j < ny_p-1 ; j++ ) { // y loop
                for( unsigned int k=1 ; k < nz_p-1 ; k++ ) {
                    // ====
                    // STD Solver for propagation in vacuum
                    // ====
                    // ( *A_np1_pml )( i, j, k ) = - (1+i1*k0*dt) * ( *A_nm1_pml )( i, j, k ) ;
                    // ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + 2. * ( *A_n_pml )( i, j, k ) ;
                    // ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + 2.*i1*k0*dt*dt*( ( *A_n_pml )( i+1, j, k )-( *A_n_pml )( i-1, j, k ) )/(2.*dx) ;
                    // ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + dt*dt*( ( *A_n_pml )( i-1, j, k )-2.*( *A_n_pml )( i, j, k )+( *A_n_pml )( i+1, j, k ) )/(dx*dx) ;
                    // ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + dt*dt*( ( *A_n_pml )( i, j-1, k )-2.*( *A_n_pml )( i, j, k )+( *A_n_pml )( i, j+1, k ) )/(dy*dy) ;
                    // ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + dt*dt*( ( *A_n_pml )( i, j, k-1 )-2.*( *A_n_pml )( i, j, k )+( *A_n_pml )( i, j, k+1 ) )/(dz*dz) ;
                    // ( *A_np1_pml )( i, j, k ) = ( ( 1+i1*k0*dt) / (1+k0*k0*dt*dt) )*( *A_np1_pml )( i, j, k );
                    // ====
                    // ADE-CFS-PML
                    // ====
                    // 1. update u3
                    // 2. update u2
                    // 3. update u1
                    // 4. use u1 like a source terme in the standard enveloppe scheme
                    // ====
                    // Envelope hypothesis with the carrier wave which propagate along x :
                    // ----
                    // dA/dx = dA/dx + ik0 A
                    std::complex<double> dA_over_dx_fdtd = ( ( *A_n_pml )( i+1, j, k )-( *A_n_pml )( i-1, j, k ) )/(2.*dx) ;
                    std::complex<double> dA_over_dx = dA_over_dx_fdtd
                                                      + i1*k0*( *A_n_pml )( i, j, k ) ;
                    // d2A/dx^2 = d2A/dx^2 + 2ik0 dA/dx - k0^2 A
                    std::complex<double> d2A_over_dx2_fdtd = ( ( *A_n_pml )( i-1, j, k )-2.*( *A_n_pml )( i, j, k )+( *A_n_pml )( i+1, j, k ) )/(dx*dx) ;
                    std::complex<double> d2A_over_dx2 = d2A_over_dx2_fdtd
                                                        + 2.*i1*k0*( ( *A_n_pml )( i+1, j, k )-( *A_n_pml )( i-1, j, k ) )/(2.*dx)
                                                        - k0*k0*( *A_n_pml )( i, j, k ) ;
                    // d2A/dy^2 = d2A/dy^2
                    std::complex<double> d2A_over_dy2 = ( ( *A_n_pml )( i, j-1, k )-2.*( *A_n_pml )( i, j, k )+( *A_n_pml )( i, j+1, k ) )/(dy*dy) ;
                    // d2A/dz^2 = d2A/dz^2
                    std::complex<double> d2A_over_dz2 = ( ( *A_n_pml )( i, j, k-1 )-2.*( *A_n_pml )( i, j, k )+( *A_n_pml )( i, j, k+1 ) )/(dz*dz) ;
                    // ====
                    // ADE update
                    // ----
                    // CFS-PML Block
                    // 1. update u3
                    ( *u3_np1_x_pml )( i, j, k ) = -kappa_prime_x_p[i]*sigma_x_p[i] ;
                    ( *u3_np1_x_pml )( i, j, k ) = ( *u3_np1_x_pml )( i, j, k ) + sigma_prime_x_p[i]*kappa_x_p[i] ;
                    ( *u3_np1_x_pml )( i, j, k ) = ( *u3_np1_x_pml )( i, j, k ) + alpha_prime_x_p[i]*pow(kappa_x_p[i],2) ;
                    ( *u3_np1_x_pml )( i, j, k ) = ( *u3_np1_x_pml )( i, j, k ) * pow(sigma_x_p[i],2) * dA_over_dx / pow(kappa_x_p[i],4) ;
                    // time operation on u3 : Be carefull, u3 has to be considered like an envelop * a carrier wave
                    ( *u3_np1_x_pml )( i, j, k ) = ( ( *u3_np1_x_pml )( i, j, k ) - ( *u3_nm1_x_pml )( i, j, k )*( 1. + 0.5*dt*( i1*k0 + alpha_x_p[i] + sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) * dt / ( 0.5*dt*(i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    //( *u3_np1_x_pml )( i, j, k ) = ( ( *u3_np1_x_pml )( i, j, k ) - ( *u3_nm1_x_pml )( i, j, k )*( 1. + 1.0*dt*( i1*k0 + alpha_x_p[i] + sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) * dt / ( 1.0*dt*(i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    // 2. update u2
                    ( *u2_np1_x_pml )( i, j, k ) = (2.*sigma_prime_x_p[i]*kappa_x_p[i]+pow(kappa_x_p[i],2)*alpha_prime_x_p[i]-3.*kappa_prime_x_p[i]*sigma_x_p[i])*dA_over_dx ;
                    ( *u2_np1_x_pml )( i, j, k ) = ( *u2_np1_x_pml )( i, j, k ) + sigma_x_p[i]*kappa_x_p[i]*d2A_over_dx2 ;
                    ( *u2_np1_x_pml )( i, j, k ) = ( *u2_np1_x_pml )( i, j, k ) * sigma_x_p[i] ;
                    ( *u2_np1_x_pml )( i, j, k ) = ( *u2_np1_x_pml )( i, j, k ) - pow(kappa_x_p[i],3)*0.5*( ( *u3_np1_x_pml )( i, j, k ) + ( *u3_nm1_x_pml )( i, j, k ) ) ;
                    //( *u2_np1_x_pml )( i, j, k ) = ( *u2_np1_x_pml )( i, j, k ) - pow(kappa_x_p[i],3)*( *u3_np1_x_pml )( i, j, k ) ;
                    ( *u2_np1_x_pml )( i, j, k ) = ( *u2_np1_x_pml )( i, j, k ) / pow(kappa_x_p[i],4) ;
                    // time operation on u2 : Be carefull, u2 has to be considered like an envelop * a carrier wave
                    ( *u2_np1_x_pml )( i, j, k ) = ( ( *u2_np1_x_pml )( i, j, k ) - ( *u2_nm1_x_pml )( i, j, k )*( 1. + 0.5*dt*( i1*k0 + alpha_x_p[i] + sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) * dt / ( 0.5*dt*(i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    //( *u2_np1_x_pml )( i, j, k ) = ( ( *u2_np1_x_pml )( i, j, k ) - ( *u2_nm1_x_pml )( i, j, k )*( 1. + 1.0*dt*( i1*k0 + alpha_x_p[i] + sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) * dt / ( 1.0*dt*(i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    // 3. update u1
                    ( *u1_np1_x_pml )( i, j, k ) = ( sigma_prime_x_p[i]*kappa_x_p[i] - 3*kappa_prime_x_p[i]*sigma_x_p[i] ) * dA_over_dx ;
                    ( *u1_np1_x_pml )( i, j, k ) = ( *u1_np1_x_pml )( i, j, k ) + 2.*sigma_x_p[i]*kappa_x_p[i]*d2A_over_dx2 ;
                    ( *u1_np1_x_pml )( i, j, k ) = ( *u1_np1_x_pml )( i, j, k ) - pow(kappa_x_p[i],3)*0.5*( ( *u2_np1_x_pml )( i, j, k ) + ( *u2_nm1_x_pml )( i, j, k ) ) ;
                    //( *u1_np1_x_pml )( i, j, k ) = ( *u1_np1_x_pml )( i, j, k ) - pow(kappa_x_p[i],3)*( *u2_np1_x_pml )( i, j, k ) ;
                    ( *u1_np1_x_pml )( i, j, k ) = ( *u1_np1_x_pml )( i, j, k ) / pow(kappa_x_p[i],4) ;
                    // time operation on u1 : Be carefull, u1 has to be considered like an envelop * a carrier wave
                    ( *u1_np1_x_pml )( i, j, k ) = ( ( *u1_np1_x_pml )( i, j, k ) - ( *u1_nm1_x_pml )( i, j, k )*( 1. + 0.5*dt*( i1*k0 + alpha_x_p[i] + sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) * dt / ( 0.5*dt*(i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    //( *u1_np1_x_pml )( i, j, k ) = ( ( *u1_np1_x_pml )( i, j, k ) - ( *u1_nm1_x_pml )( i, j, k )*( 1. + 1.0*dt*( i1*k0 + alpha_x_p[i] + sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) * dt / ( 1.0*dt*(i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    // ----
                    // Envelop udpate with correction/source terms
                    // ----
                    // 4.a update A : Correction/source terms
                    source_term_x = ( kappa_x_p[i] - pow(kappa_x_p[i],3) )*d2A_over_dx2 ;
                    source_term_x = source_term_x - kappa_prime_x_p[i]*dA_over_dx ;
                    source_term_x = source_term_x - pow(kappa_x_p[i],3)*0.5*( ( *u1_np1_x_pml )( i, j, k ) + ( *u1_nm1_x_pml )( i, j, k ) ) ;
                    source_term_x = dt*dt*source_term_x / pow(kappa_x_p[i],3) ;
                    // ----
                    ( *A_np1_pml )( i, j, k ) = 1.*source_term_x ;
                    //( *A_np1_pml )( i, j, k ) = 0;
                    // 4.b standard envelope FDTD
                    ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + dt*dt*d2A_over_dz2 ;
                    ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + dt*dt*d2A_over_dy2 ;
                    ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + dt*dt*d2A_over_dx2 ;
                    ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + dt*dt*k0*k0*( *A_n_pml )( i, j, k ) ;
                    ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) - (1.+i1*k0*dt) * ( *A_nm1_pml )( i, j, k ) ;
                    ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + 2.*( *A_n_pml )( i, j, k ) ;
                    ( *A_np1_pml )( i, j, k ) = ( ( 1.+i1*k0*dt) / (1.+k0*k0*dt*dt) )*( *A_np1_pml )( i, j, k );
                } // end y loop
            } // end x loop
        } // end z loop

        for( unsigned int i=0 ; i<nx_p ; i++ ) { // x loop
            for( unsigned int j=0 ; j < ny_p ; j++ ) { // y loop
                for( unsigned int k=0 ; k < nz_p ; k++ ) { // z loop
                    // final back-substitution
                    // Auxillary Variable
                    // X-PML-ADE
                    ( *u3_nm1_x_pml )( i, j, k )        = 1.*( *u3_np1_x_pml )( i, j, k );
                    ( *u2_nm1_x_pml )( i, j, k )        = 1.*( *u2_np1_x_pml )( i, j, k );
                    ( *u1_nm1_x_pml )( i, j, k )        = 1.*( *u1_np1_x_pml )( i, j, k );
                    // A-field
                    ( *A_nm1_pml )( i, j, k )       = 1.*( *A_n_pml )( i, j, k );
                    ( *A_n_pml )( i, j, k )         = 1.*( *A_np1_pml )( i, j, k );
                } // end z loop
            } // end y loop
        } // end x loop
    }

    else if (iDim == 1) {
        // A (p,p,p) Remind that in PML, there no current
        for( unsigned int i=1 ; i<nx_p-1; i++ ) { // x loop
            for( unsigned int j=solvermin ; j < solvermax ; j++ ) { // y loop
                for( unsigned int k=1 ; k < nz_p-1 ; k++ ) { // z loop
                    // ====
                    // STD Solver for propagation in vacuum
                    // ====
                    // ( *A_np1_pml )( i, j, k ) = - (1+i1*k0*dt) * ( *A_nm1_pml )( i, j, k ) ;
                    // ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + 2. * ( *A_n_pml )( i, j, k ) ;
                    // ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + 2.*i1*k0*dt*dt*( ( *A_n_pml )( i+1, j, k )-( *A_n_pml )( i-1, j, k ) )/(2.*dx) ;
                    // ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + dt*dt*( ( *A_n_pml )( i-1, j, k )-2.*( *A_n_pml )( i, j, k )+( *A_n_pml )( i+1, j, k ) )/(dx*dx) ;
                    // ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + dt*dt*( ( *A_n_pml )( i, j-1, k )-2.*( *A_n_pml )( i, j, k )+( *A_n_pml )( i, j+1, k ) )/(dy*dy) ;
                    // ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + dt*dt*( ( *A_n_pml )( i, j, k-1 )-2.*( *A_n_pml )( i, j, k )+( *A_n_pml )( i, j, k+1 ) )/(dz*dz) ;
                    // ( *A_np1_pml )( i, j, k ) = ( ( 1+i1*k0*dt) / (1+k0*k0*dt*dt) )*( *A_np1_pml )( i, j, k );
                    // ====
                    // ADE-CFS-PML
                    // ====
                    // 1. update u3
                    // 2. update u2
                    // 3. update u1
                    // 4. use u1 like a source terme in the standard enveloppe scheme
                    // ====
                    // Envelope hypothesis with the carrier wave which propagate along x :
                    // ----
                    // dA/dx = dA/dx + ik0 A
                    std::complex<double> dA_over_dx_fdtd = ( ( *A_n_pml )( i+1, j, k )-( *A_n_pml )( i-1, j, k ) )/(2.*dx) ;
                    std::complex<double> dA_over_dx = dA_over_dx_fdtd
                                                      + i1*k0*( *A_n_pml )( i, j, k ) ;
                    // d2A/dx^2 = d2A/dx^2 + 2ik0 dA/dx - k0^2 A
                    std::complex<double> d2A_over_dx2_fdtd = ( ( *A_n_pml )( i-1, j, k )-2.*( *A_n_pml )( i, j, k )+( *A_n_pml )( i+1, j, k ) )/(dx*dx) ;
                    std::complex<double> d2A_over_dx2 = d2A_over_dx2_fdtd
                                                        + 2.*i1*k0*( ( *A_n_pml )( i+1, j, k )-( *A_n_pml )( i-1, j, k ) )/(2.*dx)
                                                        - k0*k0*( *A_n_pml )( i, j, k ) ;
                    // dA/dy = dA/dy
                    std::complex<double> dA_over_dy = ( ( *A_n_pml )( i, j+1, k )-( *A_n_pml )( i, j-1, k ) )/(2.*dy) ;
                    // d2A/dy^2 = d2A/dy^2
                    std::complex<double> d2A_over_dy2 = ( ( *A_n_pml )( i, j-1, k )-2.*( *A_n_pml )( i, j, k )+( *A_n_pml )( i, j+1, k ) )/(dy*dy) ;
                    // d2A/dz^2 = d2A/dz^2
                    std::complex<double> d2A_over_dz2 = ( ( *A_n_pml )( i, j, k-1 )-2.*( *A_n_pml )( i, j, k )+( *A_n_pml )( i, j, k+1 ) )/(dz*dz) ;
                    // ====
                    // ADE update
                    // ----
                    // La modification des opérateurs différentiels est independant les uns des autres
                    // Il ammène a résoudres plusieurs ADE et simplement d'additionner le resultat dans les regions ou les coefficients ne sont pas nul
                    // Ainsi il faut resoudre u3_x -> u2_x -> u1_x -> source_terme_x
                    // Puis                   u3_y -> u2_y -> u1_y -> source_terme_y
                    // Enfin                  Equations diff avec les termes sources !
                    // ----
                    // 1. update u3
                    ( *u3_np1_x_pml )( i, j, k ) = -kappa_prime_x_p[i]*sigma_x_p[i] ;
                    ( *u3_np1_x_pml )( i, j, k ) = ( *u3_np1_x_pml )( i, j, k ) + sigma_prime_x_p[i]*kappa_x_p[i] ;
                    ( *u3_np1_x_pml )( i, j, k ) = ( *u3_np1_x_pml )( i, j, k ) + alpha_prime_x_p[i]*pow(kappa_x_p[i],2) ;
                    ( *u3_np1_x_pml )( i, j, k ) = ( *u3_np1_x_pml )( i, j, k ) * pow(sigma_x_p[i],2) * dA_over_dx / pow(kappa_x_p[i],4) ;
                    // time operation on u3 : Be carefull, u3 has to be considered like an envelop * a carrier wave
                    ( *u3_np1_x_pml )( i, j, k ) = ( ( *u3_np1_x_pml )( i, j, k ) - ( *u3_nm1_x_pml )( i, j, k )*( 1. + 0.5*dt*( i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) * dt / ( 0.5*dt*(i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    // 2. update u2
                    ( *u2_np1_x_pml )( i, j, k ) = (2.*sigma_prime_x_p[i]*kappa_x_p[i]+pow(kappa_x_p[i],2)*alpha_prime_x_p[i]-3.*kappa_prime_x_p[i]*sigma_x_p[i])*dA_over_dx ;
                    ( *u2_np1_x_pml )( i, j, k ) = ( *u2_np1_x_pml )( i, j, k ) + sigma_x_p[i]*kappa_x_p[i]*d2A_over_dx2 ;
                    ( *u2_np1_x_pml )( i, j, k ) = ( *u2_np1_x_pml )( i, j, k ) * sigma_x_p[i] ;
                    ( *u2_np1_x_pml )( i, j, k ) = ( *u2_np1_x_pml )( i, j, k ) - pow(kappa_x_p[i],3)*0.5*( ( *u3_np1_x_pml )( i, j, k ) + ( *u3_nm1_x_pml )( i, j, k ) ) ;
                    ( *u2_np1_x_pml )( i, j, k ) = ( *u2_np1_x_pml )( i, j, k ) / pow(kappa_x_p[i],4) ;
                    // time operation on u2 : Be carefull, u2 has to be considered like an envelop * a carrier wave
                    ( *u2_np1_x_pml )( i, j, k ) = ( ( *u2_np1_x_pml )( i, j, k ) - ( *u2_nm1_x_pml )( i, j, k )*( 1. + 0.5*dt*( i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) * dt / ( 0.5*dt*(i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    // 3. update u1
                    ( *u1_np1_x_pml )( i, j, k ) = ( sigma_prime_x_p[i]*kappa_x_p[i] - 3*kappa_prime_x_p[i]*sigma_x_p[i] ) * dA_over_dx ;
                    ( *u1_np1_x_pml )( i, j, k ) = ( *u1_np1_x_pml )( i, j, k ) + 2.*sigma_x_p[i]*kappa_x_p[i]*d2A_over_dx2 ;
                    ( *u1_np1_x_pml )( i, j, k ) = ( *u1_np1_x_pml )( i, j, k ) - pow(kappa_x_p[i],3)*0.5*( ( *u2_np1_x_pml )( i, j, k ) + ( *u2_nm1_x_pml )( i, j, k ) ) ;
                    ( *u1_np1_x_pml )( i, j, k ) = ( *u1_np1_x_pml )( i, j, k ) / pow(kappa_x_p[i],4) ;
                    // time operation on u1 : Be carefull, u1 has to be considered like an envelop * a carrier wave
                    ( *u1_np1_x_pml )( i, j, k ) = ( ( *u1_np1_x_pml )( i, j, k ) - ( *u1_nm1_x_pml )( i, j, k )*( 1. + 0.5*dt*( i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) *dt / ( 0.5*dt*(i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    // 1. update u3
                    ( *u3_np1_y_pml )( i, j, k ) = -kappa_prime_y_p[j]*sigma_y_p[j] ;
                    ( *u3_np1_y_pml )( i, j, k ) = ( *u3_np1_y_pml )( i, j, k ) + sigma_prime_y_p[j]*kappa_y_p[j] ;
                    ( *u3_np1_y_pml )( i, j, k ) = ( *u3_np1_y_pml )( i, j, k ) + alpha_prime_y_p[j]*pow(kappa_y_p[j],2) ;
                    ( *u3_np1_y_pml )( i, j, k ) = ( *u3_np1_y_pml )( i, j, k ) * pow(sigma_y_p[j],2) * dA_over_dy / pow(kappa_y_p[j],4) ;
                    // time operation on u3 : Be carefull, u3 has to be considered like an envelop * a carrier wave
                    ( *u3_np1_y_pml )( i, j, k ) = ( ( *u3_np1_y_pml )( i, j, k ) - ( *u3_nm1_y_pml )( i, j, k )*( 1. + 0.5*dt*( i1*k0 + alpha_y_p[j]+sigma_y_p[j]/kappa_y_p[j] ) ) / dt ) * dt / ( 0.5*dt*(i1*k0 + alpha_y_p[j]+sigma_y_p[j]/kappa_y_p[j]) - 1. ) ;
                    // 2. update u2
                    ( *u2_np1_y_pml )( i, j, k ) = (2.*sigma_prime_y_p[j]*kappa_y_p[j]+pow(kappa_y_p[j],2)*alpha_prime_y_p[j]-3.*kappa_prime_y_p[j]*sigma_y_p[j])*dA_over_dy ;
                    ( *u2_np1_y_pml )( i, j, k ) = ( *u2_np1_y_pml )( i, j, k ) + sigma_y_p[j]*kappa_y_p[j]*d2A_over_dy2 ;
                    ( *u2_np1_y_pml )( i, j, k ) = ( *u2_np1_y_pml )( i, j, k ) * sigma_y_p[j] ;
                    ( *u2_np1_y_pml )( i, j, k ) = ( *u2_np1_y_pml )( i, j, k ) - pow(kappa_y_p[j],3)*0.5*( ( *u3_np1_y_pml )( i, j, k ) + ( *u3_nm1_y_pml )( i, j, k ) ) ;
                    ( *u2_np1_y_pml )( i, j, k ) = ( *u2_np1_y_pml )( i, j, k ) / pow(kappa_y_p[j],4) ;
                    // time operation on u2 : Be carefull, u2 has to be considered like an envelop * a carrier wave
                    ( *u2_np1_y_pml )( i, j, k ) = ( ( *u2_np1_y_pml )( i, j, k ) - ( *u2_nm1_y_pml )( i, j, k )*( 1. + 0.5*dt*( i1*k0 + alpha_y_p[j]+sigma_y_p[j]/kappa_y_p[j] ) ) / dt ) * dt / ( 0.5*dt*(i1*k0 + alpha_y_p[j]+sigma_y_p[j]/kappa_y_p[j]) - 1. ) ;
                    // 3. update u1
                    ( *u1_np1_y_pml )( i, j, k ) = ( sigma_prime_y_p[j]*kappa_y_p[j] - 3*kappa_prime_y_p[j]*sigma_y_p[j] ) * dA_over_dy ;
                    ( *u1_np1_y_pml )( i, j, k ) = ( *u1_np1_y_pml )( i, j, k ) + 2.*sigma_y_p[j]*kappa_y_p[j]*d2A_over_dy2 ;
                    ( *u1_np1_y_pml )( i, j, k ) = ( *u1_np1_y_pml )( i, j, k ) - pow(kappa_y_p[j],3)*0.5*( ( *u2_np1_y_pml )( i, j, k ) + ( *u2_nm1_y_pml )( i, j, k ) ) ;
                    ( *u1_np1_y_pml )( i, j, k ) = ( *u1_np1_y_pml )( i, j, k ) / pow(kappa_y_p[j],4) ;
                    // time operation on u1 : Be carefull, u1 has to be considered like an envelop * a carrier wave
                    ( *u1_np1_y_pml )( i, j, k ) = ( ( *u1_np1_y_pml )( i, j, k ) - ( *u1_nm1_y_pml )( i, j, k )*( 1. + 0.5*dt*( i1*k0 + alpha_y_p[j]+sigma_y_p[j]/kappa_y_p[j] ) ) / dt ) *dt / ( 0.5*dt*(i1*k0 + alpha_y_p[j]+sigma_y_p[j]/kappa_y_p[j]) - 1. ) ;
                    // ----
                    // Envelop udpate with correction/source terms
                    // ----
                    // 4.a update A : Correction/source terms
                    source_term_x = ( kappa_x_p[i] - pow(kappa_x_p[i],3) )*d2A_over_dx2 ;
                    source_term_x = source_term_x - kappa_prime_x_p[i]*dA_over_dx ;
                    source_term_x = source_term_x - pow(kappa_x_p[i],3)*0.5*( ( *u1_np1_x_pml )( i, j, k ) + ( *u1_nm1_x_pml )( i, j, k ) ) ;
                    source_term_x = dt*dt*source_term_x / pow(kappa_x_p[i],3) ;
                    // ----
                    source_term_y = ( kappa_y_p[j] - pow(kappa_y_p[j],3) )*d2A_over_dy2 ;
                    source_term_y = source_term_y - kappa_prime_y_p[j]*dA_over_dy ;
                    source_term_y = source_term_y - pow(kappa_y_p[j],3)*0.5*( ( *u1_np1_y_pml )( i, j, k ) + ( *u1_nm1_y_pml )( i, j, k ) ) ;
                    source_term_y = dt*dt*source_term_y / pow(kappa_y_p[j],3) ;
                    // ----
                    ( *A_np1_pml )( i, j, k ) = 1.*source_term_x + 1.*source_term_y ;
                    // ( *A_np1_pml )( i, j, k ) = 0;
                    // 4.b standard envelope FDTD
                    ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + dt*dt*d2A_over_dz2 ;
                    ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + dt*dt*d2A_over_dy2 ;
                    ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + dt*dt*d2A_over_dx2 ;
                    ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + dt*dt*k0*k0*( *A_n_pml )( i, j, k ) ;
                    ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) - (1.+i1*k0*dt) * ( *A_nm1_pml )( i, j, k ) ;
                    ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + 2.*( *A_n_pml )( i, j, k ) ;
                    ( *A_np1_pml )( i, j, k ) = ( ( 1.+i1*k0*dt) / (1.+k0*k0*dt*dt) )*( *A_np1_pml )( i, j, k );
                } // end y loop
            } // end x loop
        }

        for( unsigned int i=0 ; i<nx_p ; i++ ) { // x loop
            for( unsigned int j=0 ; j < ny_p ; j++ ) { // y loop
                for( unsigned int k=0 ; k < nz_p ; k++ ) { // z loop
                    // X-PML-ADE
                    ( *u3_nm1_x_pml )( i, j, k )        = 1.*( *u3_np1_x_pml )( i, j, k );
                    ( *u2_nm1_x_pml )( i, j, k )        = 1.*( *u2_np1_x_pml )( i, j, k );
                    ( *u1_nm1_x_pml )( i, j, k )        = 1.*( *u1_np1_x_pml )( i, j, k );
                    // Y-PML-ADE
                    ( *u3_nm1_y_pml )( i, j, k )        = 1.*( *u3_np1_y_pml )( i, j, k );
                    ( *u2_nm1_y_pml )( i, j, k )        = 1.*( *u2_np1_y_pml )( i, j, k );
                    ( *u1_nm1_y_pml )( i, j, k )        = 1.*( *u1_np1_y_pml )( i, j, k );
                    // A-field
                    ( *A_nm1_pml )( i, j, k )       = 1.*( *A_n_pml )( i, j, k );
                    ( *A_n_pml )( i, j, k )         = 1.*( *A_np1_pml )( i, j, k );
                } // end z loop
            } // end y loop
        } // end x loop
    }
    else if (iDim == 2) {
        // A (p,p,p) Remind that in PML, there no current
        for( unsigned int i=1 ; i<nx_p-1; i++ ) { // x loop
            for( unsigned int j=1 ; j < ny_p-1 ; j++ ) { // y loop
                for( unsigned int k=solvermin ; k < solvermax ; k++ ) { // z loop
                    // ====
                    // STD Solver for propagation in vacuum
                    // ====
                    // ( *A_np1_pml )( i, j, k ) = - (1+i1*k0*dt) * ( *A_nm1_pml )( i, j, k ) ;
                    // ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + 2. * ( *A_n_pml )( i, j, k ) ;
                    // ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + 2.*i1*k0*dt*dt*( ( *A_n_pml )( i+1, j, k )-( *A_n_pml )( i-1, j, k ) )/(2.*dx) ;
                    // ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + dt*dt*( ( *A_n_pml )( i-1, j, k )-2.*( *A_n_pml )( i, j, k )+( *A_n_pml )( i+1, j, k ) )/(dx*dx) ;
                    // ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + dt*dt*( ( *A_n_pml )( i, j-1, k )-2.*( *A_n_pml )( i, j, k )+( *A_n_pml )( i, j+1, k ) )/(dy*dy) ;
                    // ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + dt*dt*( ( *A_n_pml )( i, j, k-1 )-2.*( *A_n_pml )( i, j, k )+( *A_n_pml )( i, j, k+1 ) )/(dz*dz) ;
                    // ( *A_np1_pml )( i, j, k ) = ( ( 1+i1*k0*dt) / (1+k0*k0*dt*dt) )*( *A_np1_pml )( i, j, k );
                    // ====
                    // ADE-CFS-PML
                    // ====
                    // 1. update u3
                    // 2. update u2
                    // 3. update u1
                    // 4. use u1 like a source terme in the standard enveloppe scheme
                    // ====
                    // Envelope hypothesis with the carrier wave which propagate along x :
                    // ----
                    // dA/dx = dA/dx + ik0 A
                    std::complex<double> dA_over_dx_fdtd = ( ( *A_n_pml )( i+1, j, k )-( *A_n_pml )( i-1, j, k ) )/(2.*dx) ;
                    std::complex<double> dA_over_dx = dA_over_dx_fdtd
                                                      + i1*k0*( *A_n_pml )( i, j, k ) ;
                    // d2A/dx^2 = d2A/dx^2 + 2ik0 dA/dx - k0^2 A
                    std::complex<double> d2A_over_dx2_fdtd = ( ( *A_n_pml )( i-1, j, k )-2.*( *A_n_pml )( i, j, k )+( *A_n_pml )( i+1, j, k ) )/(dx*dx) ;
                    std::complex<double> d2A_over_dx2 = d2A_over_dx2_fdtd
                                                        + 2.*i1*k0*( ( *A_n_pml )( i+1, j, k )-( *A_n_pml )( i-1, j, k ) )/(2.*dx)
                                                        - k0*k0*( *A_n_pml )( i, j, k ) ;
                    // dA/dy = dA/dy
                    std::complex<double> dA_over_dy = ( ( *A_n_pml )( i, j+1, k )-( *A_n_pml )( i, j-1, k ) )/(2.*dy) ;
                    // d2A/dy^2 = d2A/dy^2
                    std::complex<double> d2A_over_dy2 = ( ( *A_n_pml )( i, j-1, k )-2.*( *A_n_pml )( i, j, k )+( *A_n_pml )( i, j+1, k ) )/(dy*dy) ;
                    // dA/dz = dA/dz
                    std::complex<double> dA_over_dz = ( ( *A_n_pml )( i, j, k+1 )-( *A_n_pml )( i, j, k-1 ) )/(2.*dz) ;
                    // d2A/dz^2 = d2A/dz^2
                    std::complex<double> d2A_over_dz2 = ( ( *A_n_pml )( i, j, k-1 )-2.*( *A_n_pml )( i, j, k )+( *A_n_pml )( i, j, k+1 ) )/(dz*dz) ;
                    // ====
                    // ADE update
                    // ----
                    // La modification des opérateurs différentiels est independant les uns des autres
                    // Il ammène a résoudres plusieurs ADE et simplement d'additionner le resultat dans les regions ou les coefficients ne sont pas nul
                    // Ainsi il faut resoudre u3_x -> u2_x -> u1_x -> source_terme_x
                    // Puis                   u3_y -> u2_y -> u1_y -> source_terme_y
                    // Enfin                  Equations diff avec les termes sources !
                    // ----
                    // 1. update u3
                    ( *u3_np1_x_pml )( i, j, k ) = -kappa_prime_x_p[i]*sigma_x_p[i] ;
                    ( *u3_np1_x_pml )( i, j, k ) = ( *u3_np1_x_pml )( i, j, k ) + sigma_prime_x_p[i]*kappa_x_p[i] ;
                    ( *u3_np1_x_pml )( i, j, k ) = ( *u3_np1_x_pml )( i, j, k ) + alpha_prime_x_p[i]*pow(kappa_x_p[i],2) ;
                    ( *u3_np1_x_pml )( i, j, k ) = ( *u3_np1_x_pml )( i, j, k ) * pow(sigma_x_p[i],2) * dA_over_dx / pow(kappa_x_p[i],4) ;
                    // time operation on u3 : Be carefull, u3 has to be considered like an envelop * a carrier wave
                    ( *u3_np1_x_pml )( i, j, k ) = ( ( *u3_np1_x_pml )( i, j, k ) - ( *u3_nm1_x_pml )( i, j, k )*( 1. + 0.5*dt*( i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) * dt / ( 0.5*dt*(i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    // 2. update u2
                    ( *u2_np1_x_pml )( i, j, k ) = (2.*sigma_prime_x_p[i]*kappa_x_p[i]+pow(kappa_x_p[i],2)*alpha_prime_x_p[i]-3.*kappa_prime_x_p[i]*sigma_x_p[i])*dA_over_dx ;
                    ( *u2_np1_x_pml )( i, j, k ) = ( *u2_np1_x_pml )( i, j, k ) + sigma_x_p[i]*kappa_x_p[i]*d2A_over_dx2 ;
                    ( *u2_np1_x_pml )( i, j, k ) = ( *u2_np1_x_pml )( i, j, k ) * sigma_x_p[i] ;
                    ( *u2_np1_x_pml )( i, j, k ) = ( *u2_np1_x_pml )( i, j, k ) - pow(kappa_x_p[i],3)*0.5*( ( *u3_np1_x_pml )( i, j, k ) + ( *u3_nm1_x_pml )( i, j, k ) ) ;
                    ( *u2_np1_x_pml )( i, j, k ) = ( *u2_np1_x_pml )( i, j, k ) / pow(kappa_x_p[i],4) ;
                    // time operation on u2 : Be carefull, u2 has to be considered like an envelop * a carrier wave
                    ( *u2_np1_x_pml )( i, j, k ) = ( ( *u2_np1_x_pml )( i, j, k ) - ( *u2_nm1_x_pml )( i, j, k )*( 1. + 0.5*dt*( i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) * dt / ( 0.5*dt*(i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    // 3. update u1
                    ( *u1_np1_x_pml )( i, j, k ) = ( sigma_prime_x_p[i]*kappa_x_p[i] - 3*kappa_prime_x_p[i]*sigma_x_p[i] ) * dA_over_dx ;
                    ( *u1_np1_x_pml )( i, j, k ) = ( *u1_np1_x_pml )( i, j, k ) + 2.*sigma_x_p[i]*kappa_x_p[i]*d2A_over_dx2 ;
                    ( *u1_np1_x_pml )( i, j, k ) = ( *u1_np1_x_pml )( i, j, k ) - pow(kappa_x_p[i],3)*0.5*( ( *u2_np1_x_pml )( i, j, k ) + ( *u2_nm1_x_pml )( i, j, k ) ) ;
                    ( *u1_np1_x_pml )( i, j, k ) = ( *u1_np1_x_pml )( i, j, k ) / pow(kappa_x_p[i],4) ;
                    // time operation on u1 : Be carefull, u1 has to be considered like an envelop * a carrier wave
                    ( *u1_np1_x_pml )( i, j, k ) = ( ( *u1_np1_x_pml )( i, j, k ) - ( *u1_nm1_x_pml )( i, j, k )*( 1. + 0.5*dt*( i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) *dt / ( 0.5*dt*(i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    // 1. update u3
                    ( *u3_np1_y_pml )( i, j, k ) = -kappa_prime_y_p[j]*sigma_y_p[j] ;
                    ( *u3_np1_y_pml )( i, j, k ) = ( *u3_np1_y_pml )( i, j, k ) + sigma_prime_y_p[j]*kappa_y_p[j] ;
                    ( *u3_np1_y_pml )( i, j, k ) = ( *u3_np1_y_pml )( i, j, k ) + alpha_prime_y_p[j]*pow(kappa_y_p[j],2) ;
                    ( *u3_np1_y_pml )( i, j, k ) = ( *u3_np1_y_pml )( i, j, k ) * pow(sigma_y_p[j],2) * dA_over_dy / pow(kappa_y_p[j],4) ;
                    // time operation on u3 : Be carefull, u3 has to be considered like an envelop * a carrier wave
                    ( *u3_np1_y_pml )( i, j, k ) = ( ( *u3_np1_y_pml )( i, j, k ) - ( *u3_nm1_y_pml )( i, j, k )*( 1. + 0.5*dt*( i1*k0 + alpha_y_p[j]+sigma_y_p[j]/kappa_y_p[j] ) ) / dt ) * dt / ( 0.5*dt*(i1*k0 + alpha_y_p[j]+sigma_y_p[j]/kappa_y_p[j]) - 1. ) ;
                    // 2. update u2
                    ( *u2_np1_y_pml )( i, j, k ) = (2.*sigma_prime_y_p[j]*kappa_y_p[j]+pow(kappa_y_p[j],2)*alpha_prime_y_p[j]-3.*kappa_prime_y_p[j]*sigma_y_p[j])*dA_over_dy ;
                    ( *u2_np1_y_pml )( i, j, k ) = ( *u2_np1_y_pml )( i, j, k ) + sigma_y_p[j]*kappa_y_p[j]*d2A_over_dy2 ;
                    ( *u2_np1_y_pml )( i, j, k ) = ( *u2_np1_y_pml )( i, j, k ) * sigma_y_p[j] ;
                    ( *u2_np1_y_pml )( i, j, k ) = ( *u2_np1_y_pml )( i, j, k ) - pow(kappa_y_p[j],3)*0.5*( ( *u3_np1_y_pml )( i, j, k ) + ( *u3_nm1_y_pml )( i, j, k ) ) ;
                    ( *u2_np1_y_pml )( i, j, k ) = ( *u2_np1_y_pml )( i, j, k ) / pow(kappa_y_p[j],4) ;
                    // time operation on u2 : Be carefull, u2 has to be considered like an envelop * a carrier wave
                    ( *u2_np1_y_pml )( i, j, k ) = ( ( *u2_np1_y_pml )( i, j, k ) - ( *u2_nm1_y_pml )( i, j, k )*( 1. + 0.5*dt*( i1*k0 + alpha_y_p[j]+sigma_y_p[j]/kappa_y_p[j] ) ) / dt ) * dt / ( 0.5*dt*(i1*k0 + alpha_y_p[j]+sigma_y_p[j]/kappa_y_p[j]) - 1. ) ;
                    // 3. update u1
                    ( *u1_np1_y_pml )( i, j, k ) = ( sigma_prime_y_p[j]*kappa_y_p[j] - 3*kappa_prime_y_p[j]*sigma_y_p[j] ) * dA_over_dy ;
                    ( *u1_np1_y_pml )( i, j, k ) = ( *u1_np1_y_pml )( i, j, k ) + 2.*sigma_y_p[j]*kappa_y_p[j]*d2A_over_dy2 ;
                    ( *u1_np1_y_pml )( i, j, k ) = ( *u1_np1_y_pml )( i, j, k ) - pow(kappa_y_p[j],3)*0.5*( ( *u2_np1_y_pml )( i, j, k ) + ( *u2_nm1_y_pml )( i, j, k ) ) ;
                    ( *u1_np1_y_pml )( i, j, k ) = ( *u1_np1_y_pml )( i, j, k ) / pow(kappa_y_p[j],4) ;
                    // time operation on u1 : Be carefull, u1 has to be considered like an envelop * a carrier wave
                    ( *u1_np1_y_pml )( i, j, k ) = ( ( *u1_np1_y_pml )( i, j, k ) - ( *u1_nm1_y_pml )( i, j, k )*( 1. + 0.5*dt*( i1*k0 + alpha_y_p[j]+sigma_y_p[j]/kappa_y_p[j] ) ) / dt ) *dt / ( 0.5*dt*(i1*k0 + alpha_y_p[j]+sigma_y_p[j]/kappa_y_p[j]) - 1. ) ;
                    // 1. update u3
                    ( *u3_np1_z_pml )( i, j, k ) = -kappa_prime_z_p[k]*sigma_z_p[k] ;
                    ( *u3_np1_z_pml )( i, j, k ) = ( *u3_np1_z_pml )( i, j, k ) + sigma_prime_z_p[k]*kappa_z_p[k] ;
                    ( *u3_np1_z_pml )( i, j, k ) = ( *u3_np1_z_pml )( i, j, k ) + alpha_prime_z_p[k]*pow(kappa_z_p[k],2) ;
                    ( *u3_np1_z_pml )( i, j, k ) = ( *u3_np1_z_pml )( i, j, k ) * pow(sigma_z_p[k],2) * dA_over_dz / pow(kappa_z_p[k],4) ;
                    // time operation on u3 : Be carefull, u3 has to be considered like an envelop * a carrier wave
                    ( *u3_np1_z_pml )( i, j, k ) = ( ( *u3_np1_z_pml )( i, j, k ) - ( *u3_nm1_z_pml )( i, j, k )*( 1. + 0.5*dt*( i1*k0 + alpha_z_p[k]+sigma_z_p[k]/kappa_z_p[k] ) ) / dt ) * dt / ( 0.5*dt*(i1*k0 + alpha_z_p[k]+sigma_z_p[k]/kappa_z_p[k]) - 1. ) ;
                    // 2. update u2
                    ( *u2_np1_z_pml )( i, j, k ) = (2.*sigma_prime_z_p[k]*kappa_z_p[k]+pow(kappa_z_p[k],2)*alpha_prime_z_p[k]-3.*kappa_prime_z_p[k]*sigma_z_p[k])*dA_over_dz ;
                    ( *u2_np1_z_pml )( i, j, k ) = ( *u2_np1_z_pml )( i, j, k ) + sigma_z_p[k]*kappa_z_p[k]*d2A_over_dz2 ;
                    ( *u2_np1_z_pml )( i, j, k ) = ( *u2_np1_z_pml )( i, j, k ) * sigma_z_p[k] ;
                    ( *u2_np1_z_pml )( i, j, k ) = ( *u2_np1_z_pml )( i, j, k ) - pow(kappa_z_p[k],3)*0.5*( ( *u3_np1_z_pml )( i, j, k ) + ( *u3_nm1_z_pml )( i, j, k ) ) ;
                    ( *u2_np1_z_pml )( i, j, k ) = ( *u2_np1_z_pml )( i, j, k ) / pow(kappa_z_p[k],4) ;
                    // time operation on u2 : Be carefull, u2 has to be considered like an envelop * a carrier wave
                    ( *u2_np1_z_pml )( i, j, k ) = ( ( *u2_np1_z_pml )( i, j, k ) - ( *u2_nm1_z_pml )( i, j, k )*( 1. + 0.5*dt*( i1*k0 + alpha_z_p[k]+sigma_z_p[k]/kappa_z_p[k] ) ) / dt ) * dt / ( 0.5*dt*(i1*k0 + alpha_z_p[k]+sigma_z_p[k]/kappa_z_p[k]) - 1. ) ;
                    // 3. update u1
                    ( *u1_np1_z_pml )( i, j, k ) = ( sigma_prime_z_p[k]*kappa_z_p[k] - 3*kappa_prime_z_p[k]*sigma_z_p[k] ) * dA_over_dz ;
                    ( *u1_np1_z_pml )( i, j, k ) = ( *u1_np1_z_pml )( i, j, k ) + 2.*sigma_z_p[k]*kappa_z_p[k]*d2A_over_dz2 ;
                    ( *u1_np1_z_pml )( i, j, k ) = ( *u1_np1_z_pml )( i, j, k ) - pow(kappa_z_p[k],3)*0.5*( ( *u2_np1_z_pml )( i, j, k ) + ( *u2_nm1_z_pml )( i, j, k ) ) ;
                    ( *u1_np1_z_pml )( i, j, k ) = ( *u1_np1_z_pml )( i, j, k ) / pow(kappa_z_p[k],4) ;
                    // time operation on u1 : Be carefull, u1 has to be considered like an envelop * a carrier wave
                    ( *u1_np1_z_pml )( i, j, k ) = ( ( *u1_np1_z_pml )( i, j, k ) - ( *u1_nm1_z_pml )( i, j, k )*( 1. + 0.5*dt*( i1*k0 + alpha_z_p[k]+sigma_z_p[k]/kappa_z_p[k] ) ) / dt ) *dt / ( 0.5*dt*(i1*k0 + alpha_z_p[k]+sigma_z_p[k]/kappa_z_p[k]) - 1. ) ;
                    // ----
                    // Envelop udpate with correction/source terms
                    // ----
                    // 4.a update A : Correction/source terms
                    source_term_x = ( kappa_x_p[i] - pow(kappa_x_p[i],3) )*d2A_over_dx2 ;
                    source_term_x = source_term_x - kappa_prime_x_p[i]*dA_over_dx ;
                    source_term_x = source_term_x - pow(kappa_x_p[i],3)*0.5*( ( *u1_np1_x_pml )( i, j, k ) + ( *u1_nm1_x_pml )( i, j, k ) ) ;
                    source_term_x = dt*dt*source_term_x / pow(kappa_x_p[i],3) ;
                    // ----
                    source_term_y = ( kappa_y_p[j] - pow(kappa_y_p[j],3) )*d2A_over_dy2 ;
                    source_term_y = source_term_y - kappa_prime_y_p[j]*dA_over_dy ;
                    source_term_y = source_term_y - pow(kappa_y_p[j],3)*0.5*( ( *u1_np1_y_pml )( i, j, k ) + ( *u1_nm1_y_pml )( i, j, k ) ) ;
                    source_term_y = dt*dt*source_term_y / pow(kappa_y_p[j],3) ;
                    // ----
                    source_term_z = ( kappa_z_p[k] - pow(kappa_z_p[k],3) )*d2A_over_dz2 ;
                    source_term_z = source_term_z - kappa_prime_z_p[k]*dA_over_dz ;
                    source_term_z = source_term_z - pow(kappa_z_p[k],3)*0.5*( ( *u1_np1_z_pml )( i, j, k ) + ( *u1_nm1_z_pml )( i, j, k ) ) ;
                    source_term_z = dt*dt*source_term_z / pow(kappa_z_p[k],3) ;
                    // ----
                    ( *A_np1_pml )( i, j, k ) = 1.*source_term_x + 1.*source_term_y + 1.*source_term_z ;
                    // ( *A_np1_pml )( i, j, k ) = 0;
                    // 4.b standard envelope FDTD
                    ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + dt*dt*d2A_over_dz2 ;
                    ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + dt*dt*d2A_over_dy2 ;
                    ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + dt*dt*d2A_over_dx2 ;
                    ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + dt*dt*k0*k0*( *A_n_pml )( i, j, k ) ;
                    ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) - (1.+i1*k0*dt) * ( *A_nm1_pml )( i, j, k ) ;
                    ( *A_np1_pml )( i, j, k ) = ( *A_np1_pml )( i, j, k ) + 2.*( *A_n_pml )( i, j, k ) ;
                    ( *A_np1_pml )( i, j, k ) = ( ( 1.+i1*k0*dt) / (1.+k0*k0*dt*dt) )*( *A_np1_pml )( i, j, k );
                } // end z loop
            } // end y loop
        } // end x loop

        for( unsigned int i=0 ; i<nx_p ; i++ ) { // x loop
            for( unsigned int j=0 ; j < ny_p ; j++ ) { // y loop
                for( unsigned int k=0 ; k < nz_p ; k++ ) { // z loop
                    // X-PML-ADE
                    ( *u3_nm1_x_pml )( i, j, k )        = 1.*( *u3_np1_x_pml )( i, j, k );
                    ( *u2_nm1_x_pml )( i, j, k )        = 1.*( *u2_np1_x_pml )( i, j, k );
                    ( *u1_nm1_x_pml )( i, j, k )        = 1.*( *u1_np1_x_pml )( i, j, k );
                    // Y-PML-ADE
                    ( *u3_nm1_y_pml )( i, j, k )        = 1.*( *u3_np1_y_pml )( i, j, k );
                    ( *u2_nm1_y_pml )( i, j, k )        = 1.*( *u2_np1_y_pml )( i, j, k );
                    ( *u1_nm1_y_pml )( i, j, k )        = 1.*( *u1_np1_y_pml )( i, j, k );
                    // z-PML-ADE
                    ( *u3_nm1_z_pml )( i, j, k )        = 1.*( *u3_np1_z_pml )( i, j, k );
                    ( *u2_nm1_z_pml )( i, j, k )        = 1.*( *u2_np1_z_pml )( i, j, k );
                    ( *u1_nm1_z_pml )( i, j, k )        = 1.*( *u1_np1_z_pml )( i, j, k );
                    // A-field
                    ( *A_nm1_pml )( i, j, k )       = 1.*( *A_n_pml )( i, j, k );
                    ( *A_n_pml )( i, j, k )         = 1.*( *A_np1_pml )( i, j, k );
                } // end z loop
            } // end y loop
        } // end x loop
    }
}
