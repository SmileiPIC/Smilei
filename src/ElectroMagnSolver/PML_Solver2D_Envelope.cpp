#include "PML_Solver2D_Envelope.h"
#include "ElectroMagn.h"
#include "LaserEnvelope.h"
#include "EnvelopeBC2D_PML.h"
#include "Field2D.h"
#include "cField2D.h"
#include "Patch.h"
#include <complex>
#include <algorithm>

PML_Solver2D_Envelope::PML_Solver2D_Envelope( Params &params )
    : Solver2D( params )
{
}

PML_Solver2D_Envelope::~PML_Solver2D_Envelope()
{
}

void PML_Solver2D_Envelope::operator()( ElectroMagn *fields )
{
    ERROR( "This is not a solver for the main domain" );

    //cField2D *A2D_n        = static_cast<cField2D *>( envelope->A_ );   // the envelope at timestep n
    //cField2D *A2D_nm1      = static_cast<cField2D *>( envelope->A0_ );  // the envelope at timestep n-1
}

void PML_Solver2D_Envelope::setDomainSizeAndCoefficients( int iDim, int min_or_max, int ncells_pml_domain, int startpml, int* ncells_pml_min, int* ncells_pml_max, Patch* patch )
{
    if ( iDim == 0 ) {
        nx_p = ncells_pml_domain;
        nx_d = ncells_pml_domain+1;
    }
    else if ( iDim == 1 ) {
        ny_p = ncells_pml_domain;
        ny_d = ncells_pml_domain+1;
        nx_p += ncells_pml_min[0]-1*(patch->isXmin()) + ncells_pml_max[0]-1*(patch->isXmax());
        nx_d += ncells_pml_min[0] + ncells_pml_max[0];
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

    // While min and max of each dim have got the same ncells_pml, coefficient are the same
    // for min and max PML

    double alpha_cx = 1.01 ; // Try to use a more practical timestep !
    double alpha_cy = 0.10 ; // 0.25 ; // 0.8 for transverse is ok

    // Best parameter
    // for kappa_max = 1.2 ; sigma_max = 2 (in reality 1.59) in order to have dt = 0.8 dt_cfl
    // for kappa_max = 1.0 ; sigma_max = 2 (in reality 1.59) dt = 0.7 dt_cfl in order to not diverge (magical timestep approx)
    // It's difficult for now to stabilize PML more with kappa and it depend on resolution
    // For now, without changing kappa, we fix sigma_max = 1 and p = 2 and 10 cells pml for dx=2 and 20 cells pml for dx=1
    // alpha is constant and > sigma_max
    // for dx = 2 dt <= 0.8
    // for dx = 1 dt <= 0.9
    // For Yee : dt <= dx / sqrt(2) / sqrt(1 + sigma**2 * dx**2/8)
    // for dx = 2, result in the box seems to be not converge (impulse is enlarging a lot due to numerical dispersion)
    // for dx = 1, impulse seems to be stable
    // The size of one PML layer for the test case is 10/128 = 8% of the transverse size
    // For a production case at dx=1 PML layer is 20 cells and box size is 400-600 cells 5%-3% of the transverse simulation size and around 1% for the longitudinal size
    // For the moment we have an absorption around -60db -50db so a 0.001 factor which is good for PML
    // For LWFA simulation dy is larger than dx. With the same number of PML cells, there will be a greater absorption for transverse field

    // Have to do :
    // [] PML parameters have to be modify in the namelist
    // [] Force alpha to be always greater than sigma max with respect of the discretization

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
            alpha_x_p[i] = 0. ;
            kappa_prime_x_p[i] = 0. ;
            sigma_prime_x_p[i] = 0. ;
            alpha_prime_x_p[i] = 0. ;
        }
        // Params for other cells (PML Media) when i>=3
        kappa_x_max = 1.0 ;
        sigma_x_max = 1.36 ; // 1.16 for 20 cells ; // 1.36 for 10 cells ;
        alpha_x_max = 0.0 ;
        power_pml_kappa_x = 3.;
        power_pml_sigma_x = 2.;
        power_pml_alpha_x = 1.;
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
        /*
        Here is to check all the PML Parameters !
        Have to be comment or suppress at the end
        */ 
        // std::cout << "iDim/MinMax : " << iDim <<"/"<< min_or_max << std::endl ;
        // std::cout << "kappa : " << std::endl << "[";
        // for(double x: kappa_x_p) {
        //     std::cout <<  x << ", " ;
        // }
        // std::cout << "]" << std::endl;
        // std::cout << "sigma : " << std::endl << "[";
        // for(double x: sigma_x_p) {
        //     std::cout << x << ", ";
        // }
        // std::cout << "]" << std::endl;
        // std::cout << "alpha : " << std::endl << "[";
        // for(double x: alpha_x_p) {
        //     std::cout << x << ", ";
        // }
        // std::cout << "]" << std::endl;
        // std::cout << "kappa_prime : " << std::endl << "[";
        // for(double x: kappa_prime_x_p) {
        //     std::cout << x << ", ";
        // }
        // std::cout << "]" << std::endl;
        // std::cout << "sigma_prime : " << std::endl << "[";
        // for(double x: sigma_prime_x_p) {
        //     std::cout << x << ", ";
        // }
        // std::cout << "]" << std::endl;
        // std::cout << "alpha_prime : " << std::endl << "[";
        // for(double x: alpha_prime_x_p) {
        //     std::cout << x << ", ";
        // }
        // std::cout << "]" << std::endl;
        // Y-direction
        for ( int j=0 ; j<ny_p ; j++ ) {
            kappa_y_p[j] = 1. ;
            sigma_y_p[j] = 0. ;
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
            kappa_x_max = 1. ;
            sigma_x_max = 1.36 ; // 1.16 for 20 cells ; // 1.36 for 10 cells ;
            alpha_x_max = 0.0 ;
            power_pml_kappa_x = 3.;
            power_pml_sigma_x = 2.;
            power_pml_alpha_x = 1.;
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
            kappa_x_max = 1. ;
            sigma_x_max = 1.36 ; // 1.16 for 20 cells ; // 1.36 for 10 cells ;
            alpha_x_max = 0.0 ;
            power_pml_kappa_x = 3.;
            power_pml_sigma_x = 2. ;
            power_pml_alpha_x = 1.;
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
        /*
        Here is to check all the PML Parameters !
        Have to be comment or suppress at the end
        */ 
        // std::cout << "iDim/MinMax (x along y): " << iDim <<"/"<< min_or_max << std::endl ;
        // std::cout << "kappa : " << std::endl << "[";
        // for(double x: kappa_x_p) {
        //     std::cout <<  x << ", " ;
        // }
        // std::cout << "]" << std::endl;
        // std::cout << "sigma : " << std::endl << "[";
        // for(double x: sigma_x_p) {
        //     std::cout << x << ", ";
        // }
        // std::cout << "]" << std::endl;
        // std::cout << "alpha : " << std::endl << "[";
        // for(double x: alpha_x_p) {
        //     std::cout << x << ", ";
        // }
        // std::cout << "]" << std::endl;
        // std::cout << "kappa_prime : " << std::endl << "[";
        // for(double x: kappa_prime_x_p) {
        //     std::cout << x << ", ";
        // }
        // std::cout << "]" << std::endl;
        // std::cout << "sigma_prime : " << std::endl << "[";
        // for(double x: sigma_prime_x_p) {
        //     std::cout << x << ", ";
        // }
        // std::cout << "]" << std::endl;
        // std::cout << "alpha_prime : " << std::endl << "[";
        // for(double x: alpha_prime_x_p) {
        //     std::cout << x << ", ";
        // }
        // std::cout << "]" << std::endl;
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
        kappa_y_max = 1. ;
        sigma_y_max = 4. ; // 2.32 ; // 2.32 for 20 cells ; // 2.72 for 10 cells ;
        alpha_y_max = 0.0 ;
        power_pml_kappa_y = 3.;
        power_pml_sigma_y = 2.;
        power_pml_alpha_y = 1.;
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
        /*
        Here is to check all the PML Parameters !
        Have to be comment or suppress at the end
        */ 
        // std::cout << "iDim/MinMax : " << iDim <<"/"<< min_or_max << std::endl ;
        // std::cout << "kappa : " << std::endl << "[";
        // for(double x: kappa_y_p) {
        //     std::cout <<  x << ", " ;
        // }
        // std::cout << "]" << std::endl;
        // std::cout << "sigma : " << std::endl << "[";
        // for(double x: sigma_y_p) {
        //     std::cout << x << ", ";
        // }
        // std::cout << "]" << std::endl;
        // std::cout << "alpha : " << std::endl << "[";
        // for(double x: alpha_y_p) {
        //     std::cout << x << ", ";
        // }
        // std::cout << "]" << std::endl;
        // std::cout << "kappa_prime : " << std::endl << "[";
        // for(double x: kappa_prime_y_p) {
        //     std::cout << x << ", ";
        // }
        // std::cout << "]" << std::endl;
        // std::cout << "sigma_prime : " << std::endl << "[";
        // for(double x: sigma_prime_y_p) {
        //     std::cout << x << ", ";
        // }
        // std::cout << "]" << std::endl;
        // std::cout << "alpha_prime : " << std::endl << "[";
        // for(double x: alpha_prime_y_p) {
        //     std::cout << x << ", ";
        // }
        // std::cout << "]" << std::endl;
    }
}

void PML_Solver2D_Envelope::compute_A_from_G( LaserEnvelope *envelope, int iDim, int min_or_max, int solvermin, int solvermax )
{
    EnvelopeBC2D_PML* pml_fields = static_cast<EnvelopeBC2D_PML*>( envelope->EnvBoundCond[iDim*2+min_or_max] );

    cField2D* A2D_nm1_pml = NULL;
    cField2D* u1_nm1_x_pml = NULL;
    cField2D* u2_nm1_x_pml = NULL;
    cField2D* u3_nm1_x_pml = NULL;
    cField2D* u1_nm1_y_pml = NULL;
    cField2D* u2_nm1_y_pml = NULL;
    cField2D* u3_nm1_y_pml = NULL;

    cField2D* A2D_n_pml = NULL;

    cField2D* A2D_np1_pml = NULL;
    cField2D* u1_np1_x_pml = NULL;
    cField2D* u2_np1_x_pml = NULL;
    cField2D* u3_np1_x_pml = NULL;
    cField2D* u1_np1_y_pml = NULL;
    cField2D* u2_np1_y_pml = NULL;
    cField2D* u3_np1_y_pml = NULL;

    A2D_nm1_pml = pml_fields->A2D_nm1_;
    u1_nm1_x_pml = pml_fields->u1_nm1_x_;
    u2_nm1_x_pml = pml_fields->u2_nm1_x_;
    u3_nm1_x_pml = pml_fields->u3_nm1_x_;
    u1_nm1_y_pml = pml_fields->u1_nm1_y_;
    u2_nm1_y_pml = pml_fields->u2_nm1_y_;
    u3_nm1_y_pml = pml_fields->u3_nm1_y_;

    A2D_n_pml = pml_fields->A2D_n_;

    A2D_np1_pml = pml_fields->A2D_np1_;
    u1_np1_x_pml = pml_fields->u1_np1_x_;
    u2_np1_x_pml = pml_fields->u2_np1_x_;
    u3_np1_x_pml = pml_fields->u3_np1_x_;
    u1_np1_y_pml = pml_fields->u1_np1_y_;
    u2_np1_y_pml = pml_fields->u2_np1_y_;
    u3_np1_y_pml = pml_fields->u3_np1_y_;

    // Auxiliary Quantities
    std::complex<double> i1 = std::complex<double>( 0., 1. ); // imaginary unit
    double k0 = 1.; // laser wavenumber
    std::complex<double> source_term_x ;
    std::complex<double> source_term_y ;

    if (iDim == 0) {
        double c_yx_kappa = 1.00 ;
        double c_yx_sigma = 0.00 ;
        double c_yx_alpha = 0.00 ;
        // A (p,p,p) Remind that in PML, there no current
        for( unsigned int k=0 ; k<1 ; k++ ) {
            //// explicit solver
            for( unsigned int i=solvermin ; i<solvermax; i++ ) { // x loop
                for( unsigned int j=1 ; j < ny_p-1 ; j++ ) { // y loop
                    // ====
                    // STD Solver for propagation in vacuum
                    // ====
                    // ( *A2D_np1_pml )( i, j ) = - (1+i1*k0*dt) * ( *A2D_nm1_pml )( i, j ) ;
                    // ( *A2D_np1_pml )( i, j ) = ( *A2D_np1_pml )( i, j ) + 2. * ( *A2D_n_pml )( i, j ) ;
                    // ( *A2D_np1_pml )( i, j ) = ( *A2D_np1_pml )( i, j ) + 2.*i1*k0*dt*dt*( ( *A2D_n_pml )( i+1, j )-( *A2D_n_pml )( i-1, j ) )/(2.*dx) ;
                    // ( *A2D_np1_pml )( i, j ) = ( *A2D_np1_pml )( i, j ) + dt*dt*( ( *A2D_n_pml )( i-1, j )-2.*( *A2D_n_pml )( i, j )+( *A2D_n_pml )( i+1, j ) )/(dx*dx) ;
                    // ( *A2D_np1_pml )( i, j ) = ( *A2D_np1_pml )( i, j ) + dt*dt*( ( *A2D_n_pml )( i, j-1 )-2.*( *A2D_n_pml )( i, j )+( *A2D_n_pml )( i, j+1 ) )/(dy*dy) ;
                    // ( *A2D_np1_pml )( i, j ) = ( ( 1+i1*k0*dt) / (1+k0*k0*dt*dt) )*( *A2D_np1_pml )( i, j );
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
                    std::complex<double> dA_over_dx_fdtd = ( ( *A2D_n_pml )( i+1, j )-( *A2D_n_pml )( i-1, j ) )/(2.*dx) ;
                    std::complex<double> dA_over_dx = dA_over_dx_fdtd
                                                      + i1*k0*( *A2D_n_pml )( i, j ) ;
                    // d2A/dx^2 = d2A/dx^2 + 2ik0 dA/dx - k0^2 A
                    std::complex<double> d2A_over_dx2_fdtd = ( ( *A2D_n_pml )( i-1, j )-2.*( *A2D_n_pml )( i, j )+( *A2D_n_pml )( i+1, j ) )/(dx*dx) ;
                    std::complex<double> d2A_over_dx2 = d2A_over_dx2_fdtd
                                                        + 2.*i1*k0*( ( *A2D_n_pml )( i+1, j )-( *A2D_n_pml )( i-1, j ) )/(2.*dx)
                                                        - k0*k0*( *A2D_n_pml )( i, j ) ;
                    // d2A/dy^2 = d2A/dy^2
                    std::complex<double> d2A_over_dy2 = ( ( *A2D_n_pml )( i, j-1 )-2.*( *A2D_n_pml )( i, j )+( *A2D_n_pml )( i, j+1 ) )/(dy*dy) ;
                    // ====
                    // ADE update
                    // ----
                    // M-PML Block
                    // 1. update u2
                    // ( *u2_np1_y_pml )( i, j ) = pow(c_yx_sigma*sigma_x_p[i],2)/pow(c_yx_kappa*kappa_x_p[i],2)*d2A_over_dy2 ;
                    // // time operation on u2 : Be carefull, u2 has to be considered like an envelop * a carrier wave
                    // ( *u2_np1_y_pml )( i, j ) = ( ( *u2_np1_y_pml )( i, j ) - ( *u2_nm1_y_pml )( i, j )*( 1. + 0.5*dt*( i1*k0 + c_yx_alpha*alpha_x_p[i]+c_yx_sigma*sigma_x_p[i]/(c_yx_kappa*kappa_x_p[i]) ) ) / dt ) * dt / ( 0.5*dt*( i1*k0 + c_yx_alpha*alpha_x_p[i]+c_yx_sigma*sigma_x_p[i]/(c_yx_kappa*kappa_x_p[i])) - 1. ) ;
                    // // 2. update u1
                    // ( *u1_np1_y_pml )( i, j ) = 2.*c_yx_sigma*sigma_x_p[i]*d2A_over_dy2 ;
                    // ( *u1_np1_y_pml )( i, j ) = ( *u1_np1_y_pml )( i, j ) - pow(c_yx_kappa*kappa_x_p[i],2)*0.5*( ( *u2_np1_y_pml )( i, j ) + ( *u2_nm1_y_pml )( i, j ) ) ;
                    // ( *u1_np1_y_pml )( i, j ) = ( *u1_np1_y_pml )( i, j ) / pow(c_yx_kappa*kappa_x_p[i],2) ;
                    // // time operation on u1 : Be carefull, u1 has to be considered like an envelop * a carrier wave
                    // ( *u1_np1_y_pml )( i, j ) = ( ( *u1_np1_y_pml )( i, j ) - ( *u1_nm1_y_pml )( i, j )*( 1. + 0.5*dt*( i1*k0 + c_yx_alpha*alpha_x_p[i]+c_yx_sigma*sigma_x_p[i]/(c_yx_kappa*kappa_x_p[i]) ) ) / dt ) * dt / ( 0.5*dt*(i1*k0 + c_yx_alpha*alpha_x_p[i]+c_yx_sigma*sigma_x_p[i]/(c_yx_kappa*kappa_x_p[i])) - 1. ) ;
                    // CFS-PML Block
                    // 1. update u3
                    ( *u3_np1_x_pml )( i, j ) = -kappa_prime_x_p[i]*sigma_x_p[i] ;
                    ( *u3_np1_x_pml )( i, j ) = ( *u3_np1_x_pml )( i, j ) + sigma_prime_x_p[i]*kappa_x_p[i] ;
                    ( *u3_np1_x_pml )( i, j ) = ( *u3_np1_x_pml )( i, j ) + alpha_prime_x_p[i]*pow(kappa_x_p[i],2) ;
                    ( *u3_np1_x_pml )( i, j ) = ( *u3_np1_x_pml )( i, j ) * pow(sigma_x_p[i],2) * dA_over_dx / pow(kappa_x_p[i],4) ;
                    // time operation on u3 : Be carefull, u3 has to be considered like an envelop * a carrier wave
                    ( *u3_np1_x_pml )( i, j ) = ( ( *u3_np1_x_pml )( i, j ) - ( *u3_nm1_x_pml )( i, j )*( 1. + 0.5*dt*( i1*k0 + alpha_x_p[i] + sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) * dt / ( 0.5*dt*(i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    //( *u3_np1_x_pml )( i, j ) = ( ( *u3_np1_x_pml )( i, j ) - ( *u3_nm1_x_pml )( i, j )*( 1. + 1.0*dt*( i1*k0 + alpha_x_p[i] + sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) * dt / ( 1.0*dt*(i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    //( *u3_np1_x_pml )( i, j ) = ( ( *u3_np1_x_pml )( i, j ) - ( *u3_nm1_x_pml )( i, j )*( 1. + 0.5*dt*( alpha_x_p[i] + sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) * dt / ( 0.5*dt*( alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    // 2. update u2
                    ( *u2_np1_x_pml )( i, j ) = (2.*sigma_prime_x_p[i]*kappa_x_p[i]+pow(kappa_x_p[i],2)*alpha_prime_x_p[i]-3.*kappa_prime_x_p[i]*sigma_x_p[i])*dA_over_dx ;
                    ( *u2_np1_x_pml )( i, j ) = ( *u2_np1_x_pml )( i, j ) + sigma_x_p[i]*kappa_x_p[i]*d2A_over_dx2 ;
                    ( *u2_np1_x_pml )( i, j ) = ( *u2_np1_x_pml )( i, j ) * sigma_x_p[i] ;
                    ( *u2_np1_x_pml )( i, j ) = ( *u2_np1_x_pml )( i, j ) - pow(kappa_x_p[i],3)*0.5*( ( *u3_np1_x_pml )( i, j ) + ( *u3_nm1_x_pml )( i, j ) ) ;
                    //( *u2_np1_x_pml )( i, j ) = ( *u2_np1_x_pml )( i, j ) - pow(kappa_x_p[i],3)*( *u3_np1_x_pml )( i, j ) ;
                    ( *u2_np1_x_pml )( i, j ) = ( *u2_np1_x_pml )( i, j ) / pow(kappa_x_p[i],4) ;
                    // time operation on u2 : Be carefull, u2 has to be considered like an envelop * a carrier wave
                    ( *u2_np1_x_pml )( i, j ) = ( ( *u2_np1_x_pml )( i, j ) - ( *u2_nm1_x_pml )( i, j )*( 1. + 0.5*dt*( i1*k0 + alpha_x_p[i] + sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) * dt / ( 0.5*dt*(i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    //( *u2_np1_x_pml )( i, j ) = ( ( *u2_np1_x_pml )( i, j ) - ( *u2_nm1_x_pml )( i, j )*( 1. + 1.0*dt*( i1*k0 + alpha_x_p[i] + sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) * dt / ( 1.0*dt*(i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    //( *u2_np1_x_pml )( i, j ) = ( ( *u2_np1_x_pml )( i, j ) - ( *u2_nm1_x_pml )( i, j )*( 1. + 0.5*dt*( alpha_x_p[i] + sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) * dt / ( 0.5*dt*( alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    // 3. update u1
                    ( *u1_np1_x_pml )( i, j ) = ( sigma_prime_x_p[i]*kappa_x_p[i] - 3*kappa_prime_x_p[i]*sigma_x_p[i] ) * dA_over_dx ;
                    ( *u1_np1_x_pml )( i, j ) = ( *u1_np1_x_pml )( i, j ) + 2.*sigma_x_p[i]*kappa_x_p[i]*d2A_over_dx2 ;
                    ( *u1_np1_x_pml )( i, j ) = ( *u1_np1_x_pml )( i, j ) - pow(kappa_x_p[i],3)*0.5*( ( *u2_np1_x_pml )( i, j ) + ( *u2_nm1_x_pml )( i, j ) ) ;
                    //( *u1_np1_x_pml )( i, j ) = ( *u1_np1_x_pml )( i, j ) - pow(kappa_x_p[i],3)*( *u2_np1_x_pml )( i, j ) ;
                    ( *u1_np1_x_pml )( i, j ) = ( *u1_np1_x_pml )( i, j ) / pow(kappa_x_p[i],4) ;
                    // time operation on u1 : Be carefull, u1 has to be considered like an envelop * a carrier wave
                    ( *u1_np1_x_pml )( i, j ) = ( ( *u1_np1_x_pml )( i, j ) - ( *u1_nm1_x_pml )( i, j )*( 1. + 0.5*dt*( i1*k0 + alpha_x_p[i] + sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) * dt / ( 0.5*dt*(i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    //( *u1_np1_x_pml )( i, j ) = ( ( *u1_np1_x_pml )( i, j ) - ( *u1_nm1_x_pml )( i, j )*( 1. + 1.0*dt*( i1*k0 + alpha_x_p[i] + sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) * dt / ( 1.0*dt*(i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    //( *u1_np1_x_pml )( i, j ) = ( ( *u1_np1_x_pml )( i, j ) - ( *u1_nm1_x_pml )( i, j )*( 1. + 0.5*dt*( alpha_x_p[i] + sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) * dt / ( 0.5*dt*( alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    // ----
                    // Envelop udpate with correction/source terms
                    // ----
                    // 4.a update A : Correction/source terms
                    source_term_x = ( kappa_x_p[i] - pow(kappa_x_p[i],3) )*d2A_over_dx2 ;
                    source_term_x = source_term_x - kappa_prime_x_p[i]*dA_over_dx ;
                    source_term_x = source_term_x - pow(kappa_x_p[i],3)*0.5*( ( *u1_np1_x_pml )( i, j ) + ( *u1_nm1_x_pml )( i, j ) ) ;
                    //source_term_x = source_term_x - pow(kappa_x_p[i],3)*( *u1_np1_x_pml )( i, j ) ;
                    source_term_x = dt*dt*source_term_x / pow(kappa_x_p[i],3) ;
                    // ----
                    // source_term_y = ( 1. - pow(c_yx_kappa*kappa_x_p[i],2) )*d2A_over_dy2 ;
                    // source_term_y = source_term_y - pow(c_yx_kappa*kappa_x_p[i],2)*0.5*( ( *u1_np1_y_pml )( i, j ) + ( *u1_nm1_y_pml )( i, j ) ) ;
                    // source_term_y = dt*dt*source_term_y / pow(c_yx_kappa*kappa_x_p[i],2) ;
                    // ----
                    ( *A2D_np1_pml )( i, j ) = 1.*source_term_x ; // + 1.*source_term_y ;
                    //( *A2D_np1_pml )( i, j ) = 0;
                    // 4.b standard envelope FDTD 
                    ( *A2D_np1_pml )( i, j ) = ( *A2D_np1_pml )( i, j ) + dt*dt*d2A_over_dy2 ;
                    ( *A2D_np1_pml )( i, j ) = ( *A2D_np1_pml )( i, j ) + dt*dt*d2A_over_dx2 ;
                    ( *A2D_np1_pml )( i, j ) = ( *A2D_np1_pml )( i, j ) + dt*dt*k0*k0*( *A2D_n_pml )( i, j ) ;
                    ( *A2D_np1_pml )( i, j ) = ( *A2D_np1_pml )( i, j ) - (1.+i1*k0*dt) * ( *A2D_nm1_pml )( i, j ) ;
                    ( *A2D_np1_pml )( i, j ) = ( *A2D_np1_pml )( i, j ) + 2.*( *A2D_n_pml )( i, j ) ;
                    ( *A2D_np1_pml )( i, j ) = ( ( 1.+i1*k0*dt) / (1.+k0*k0*dt*dt) )*( *A2D_np1_pml )( i, j );
                } // end y loop
            } // end x loop

            for( unsigned int i=0 ; i<nx_p ; i++ ) { // x loop
                for( unsigned int j=0 ; j < ny_p ; j++ ) { // y loop
                    // final back-substitution
                    // Auxillary Variable
                    // X-PML-ADE
                    ( *u3_nm1_x_pml )( i, j )        = 1.*( *u3_np1_x_pml )( i, j );
                    ( *u2_nm1_x_pml )( i, j )        = 1.*( *u2_np1_x_pml )( i, j );
                    ( *u1_nm1_x_pml )( i, j )        = 1.*( *u1_np1_x_pml )( i, j );
                    // Y-PML-ADE
                    ( *u2_nm1_y_pml )( i, j )        = 1.*( *u2_np1_y_pml )( i, j );
                    ( *u1_nm1_y_pml )( i, j )        = 1.*( *u1_np1_y_pml )( i, j );
                    // A-field
                    ( *A2D_nm1_pml )( i, j )       = 1.*( *A2D_n_pml )( i, j );
                    ( *A2D_n_pml )( i, j )         = 1.*( *A2D_np1_pml )( i, j );
                } // end y loop
            } // end x loop
        }
    }

    else if (iDim == 1) {
        // A (p,p,p) Remind that in PML, there no current
        for( unsigned int k=0 ; k<1 ; k++ ) {
            //// explicit solver
            for( unsigned int i=1 ; i<nx_p-1; i++ ) { // x loop
                for( unsigned int j=solvermin ; j < solvermax ; j++ ) { // y loop
                    // ====
                    // STD Solver for propagation in vacuum
                    // ====
                    // ( *A2D_np1_pml )( i, j ) = - (1+i1*k0*dt) * ( *A2D_nm1_pml )( i, j ) ;
                    // ( *A2D_np1_pml )( i, j ) = ( *A2D_np1_pml )( i, j ) + 2. * ( *A2D_n_pml )( i, j ) ;
                    // ( *A2D_np1_pml )( i, j ) = ( *A2D_np1_pml )( i, j ) + 2.*i1*k0*dt*dt*( ( *A2D_n_pml )( i+1, j )-( *A2D_n_pml )( i-1, j ) )/(2.*dx) ;
                    // ( *A2D_np1_pml )( i, j ) = ( *A2D_np1_pml )( i, j ) + dt*dt*( ( *A2D_n_pml )( i-1, j )-2.*( *A2D_n_pml )( i, j )+( *A2D_n_pml )( i+1, j ) )/(dx*dx) ;
                    // ( *A2D_np1_pml )( i, j ) = ( *A2D_np1_pml )( i, j ) + dt*dt*( ( *A2D_n_pml )( i, j-1 )-2.*( *A2D_n_pml )( i, j )+( *A2D_n_pml )( i, j+1 ) )/(dy*dy) ;
                    // ( *A2D_np1_pml )( i, j ) = ( ( 1+i1*k0*dt) / (1+k0*k0*dt*dt) )*( *A2D_np1_pml )( i, j );
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
                    std::complex<double> dA_over_dx_fdtd = ( ( *A2D_n_pml )( i+1, j )-( *A2D_n_pml )( i-1, j ) )/(2.*dx) ;
                    std::complex<double> dA_over_dx = dA_over_dx_fdtd
                                                      + i1*k0*( *A2D_n_pml )( i, j ) ;
                    // d2A/dx^2 = d2A/dx^2 + 2ik0 dA/dx - k0^2 A
                    std::complex<double> d2A_over_dx2_fdtd = ( ( *A2D_n_pml )( i-1, j )-2.*( *A2D_n_pml )( i, j )+( *A2D_n_pml )( i+1, j ) )/(dx*dx) ;
                    std::complex<double> d2A_over_dx2 = d2A_over_dx2_fdtd
                                                        + 2.*i1*k0*( ( *A2D_n_pml )( i+1, j )-( *A2D_n_pml )( i-1, j ) )/(2.*dx)
                                                        - k0*k0*( *A2D_n_pml )( i, j ) ;
                    // dA/dy = dA/dy
                    std::complex<double> dA_over_dy = ( ( *A2D_n_pml )( i, j+1 )-( *A2D_n_pml )( i, j-1 ) )/(2.*dy) ;
                    // d2A/dy^2 = d2A/dy^2
                    std::complex<double> d2A_over_dy2 = ( ( *A2D_n_pml )( i, j-1 )-2.*( *A2D_n_pml )( i, j )+( *A2D_n_pml )( i, j+1 ) )/(dy*dy) ;
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
                    ( *u3_np1_x_pml )( i, j ) = -kappa_prime_x_p[i]*sigma_x_p[i] ;
                    ( *u3_np1_x_pml )( i, j ) = ( *u3_np1_x_pml )( i, j ) + sigma_prime_x_p[i]*kappa_x_p[i] ;
                    ( *u3_np1_x_pml )( i, j ) = ( *u3_np1_x_pml )( i, j ) + alpha_prime_x_p[i]*pow(kappa_x_p[i],2) ;
                    ( *u3_np1_x_pml )( i, j ) = ( *u3_np1_x_pml )( i, j ) * pow(sigma_x_p[i],2) * dA_over_dx / pow(kappa_x_p[i],4) ;
                    // time operation on u3 : Be carefull, u3 has to be considered like an envelop * a carrier wave
                    ( *u3_np1_x_pml )( i, j ) = ( ( *u3_np1_x_pml )( i, j ) - ( *u3_nm1_x_pml )( i, j )*( 1. + 0.5*dt*( i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) * dt / ( 0.5*dt*(i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    // 2. update u2
                    ( *u2_np1_x_pml )( i, j ) = (2.*sigma_prime_x_p[i]*kappa_x_p[i]+pow(kappa_x_p[i],2)*alpha_prime_x_p[i]-3.*kappa_prime_x_p[i]*sigma_x_p[i])*dA_over_dx ;
                    ( *u2_np1_x_pml )( i, j ) = ( *u2_np1_x_pml )( i, j ) + sigma_x_p[i]*kappa_x_p[i]*d2A_over_dx2 ;
                    ( *u2_np1_x_pml )( i, j ) = ( *u2_np1_x_pml )( i, j ) * sigma_x_p[i] ;
                    ( *u2_np1_x_pml )( i, j ) = ( *u2_np1_x_pml )( i, j ) - pow(kappa_x_p[i],3)*0.5*( ( *u3_np1_x_pml )( i, j ) + ( *u3_nm1_x_pml )( i, j ) ) ;
                    ( *u2_np1_x_pml )( i, j ) = ( *u2_np1_x_pml )( i, j ) / pow(kappa_x_p[i],4) ;
                    // time operation on u2 : Be carefull, u2 has to be considered like an envelop * a carrier wave
                    ( *u2_np1_x_pml )( i, j ) = ( ( *u2_np1_x_pml )( i, j ) - ( *u2_nm1_x_pml )( i, j )*( 1. + 0.5*dt*( i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) * dt / ( 0.5*dt*(i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    // 3. update u1
                    ( *u1_np1_x_pml )( i, j ) = ( sigma_prime_x_p[i]*kappa_x_p[i] - 3*kappa_prime_x_p[i]*sigma_x_p[i] ) * dA_over_dx ;
                    ( *u1_np1_x_pml )( i, j ) = ( *u1_np1_x_pml )( i, j ) + 2.*sigma_x_p[i]*kappa_x_p[i]*d2A_over_dx2 ;
                    ( *u1_np1_x_pml )( i, j ) = ( *u1_np1_x_pml )( i, j ) - pow(kappa_x_p[i],3)*0.5*( ( *u2_np1_x_pml )( i, j ) + ( *u2_nm1_x_pml )( i, j ) ) ;
                    ( *u1_np1_x_pml )( i, j ) = ( *u1_np1_x_pml )( i, j ) / pow(kappa_x_p[i],4) ;
                    // time operation on u1 : Be carefull, u1 has to be considered like an envelop * a carrier wave
                    ( *u1_np1_x_pml )( i, j ) = ( ( *u1_np1_x_pml )( i, j ) - ( *u1_nm1_x_pml )( i, j )*( 1. + 0.5*dt*( i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i] ) ) / dt ) *dt / ( 0.5*dt*(i1*k0 + alpha_x_p[i]+sigma_x_p[i]/kappa_x_p[i]) - 1. ) ;
                    // 1. update u3
                    ( *u3_np1_y_pml )( i, j ) = -kappa_prime_y_p[j]*sigma_y_p[j] ;
                    ( *u3_np1_y_pml )( i, j ) = ( *u3_np1_y_pml )( i, j ) + sigma_prime_y_p[j]*kappa_y_p[j] ;
                    ( *u3_np1_y_pml )( i, j ) = ( *u3_np1_y_pml )( i, j ) + alpha_prime_y_p[j]*pow(kappa_y_p[j],2) ;
                    ( *u3_np1_y_pml )( i, j ) = ( *u3_np1_y_pml )( i, j ) * pow(sigma_y_p[j],2) * dA_over_dy / pow(kappa_y_p[j],4) ;
                    // time operation on u3 : Be carefull, u3 has to be considered like an envelop * a carrier wave
                    ( *u3_np1_y_pml )( i, j ) = ( ( *u3_np1_y_pml )( i, j ) - ( *u3_nm1_y_pml )( i, j )*( 1. + 0.5*dt*( i1*k0 + alpha_y_p[j]+sigma_y_p[j]/kappa_y_p[j] ) ) / dt ) * dt / ( 0.5*dt*(i1*k0 + alpha_y_p[j]+sigma_y_p[j]/kappa_y_p[j]) - 1. ) ;
                    // 2. update u2
                    ( *u2_np1_y_pml )( i, j ) = (2.*sigma_prime_y_p[j]*kappa_y_p[j]+pow(kappa_y_p[j],2)*alpha_prime_y_p[j]-3.*kappa_prime_y_p[j]*sigma_y_p[j])*dA_over_dy ;
                    ( *u2_np1_y_pml )( i, j ) = ( *u2_np1_y_pml )( i, j ) + sigma_y_p[j]*kappa_y_p[j]*d2A_over_dy2 ;
                    ( *u2_np1_y_pml )( i, j ) = ( *u2_np1_y_pml )( i, j ) * sigma_y_p[j] ;
                    ( *u2_np1_y_pml )( i, j ) = ( *u2_np1_y_pml )( i, j ) - pow(kappa_y_p[j],3)*0.5*( ( *u3_np1_y_pml )( i, j ) + ( *u3_nm1_y_pml )( i, j ) ) ;
                    ( *u2_np1_y_pml )( i, j ) = ( *u2_np1_y_pml )( i, j ) / pow(kappa_y_p[j],4) ;
                    // time operation on u2 : Be carefull, u2 has to be considered like an envelop * a carrier wave
                    ( *u2_np1_y_pml )( i, j ) = ( ( *u2_np1_y_pml )( i, j ) - ( *u2_nm1_y_pml )( i, j )*( 1. + 0.5*dt*( i1*k0 + alpha_y_p[j]+sigma_y_p[j]/kappa_y_p[j] ) ) / dt ) * dt / ( 0.5*dt*(i1*k0 + alpha_y_p[j]+sigma_y_p[j]/kappa_y_p[j]) - 1. ) ;
                    // 3. update u1
                    ( *u1_np1_y_pml )( i, j ) = ( sigma_prime_y_p[j]*kappa_y_p[j] - 3*kappa_prime_y_p[j]*sigma_y_p[j] ) * dA_over_dy ;
                    ( *u1_np1_y_pml )( i, j ) = ( *u1_np1_y_pml )( i, j ) + 2.*sigma_y_p[j]*kappa_y_p[j]*d2A_over_dy2 ;
                    ( *u1_np1_y_pml )( i, j ) = ( *u1_np1_y_pml )( i, j ) - pow(kappa_y_p[j],3)*0.5*( ( *u2_np1_y_pml )( i, j ) + ( *u2_nm1_y_pml )( i, j ) ) ;
                    ( *u1_np1_y_pml )( i, j ) = ( *u1_np1_y_pml )( i, j ) / pow(kappa_y_p[j],4) ;
                    // time operation on u1 : Be carefull, u1 has to be considered like an envelop * a carrier wave
                    ( *u1_np1_y_pml )( i, j ) = ( ( *u1_np1_y_pml )( i, j ) - ( *u1_nm1_y_pml )( i, j )*( 1. + 0.5*dt*( i1*k0 + alpha_y_p[j]+sigma_y_p[j]/kappa_y_p[j] ) ) / dt ) *dt / ( 0.5*dt*(i1*k0 + alpha_y_p[j]+sigma_y_p[j]/kappa_y_p[j]) - 1. ) ;
                    // ----
                    // Envelop udpate with correction/source terms
                    // ----
                    // 4.a update A : Correction/source terms
                    source_term_x = ( kappa_x_p[i] - pow(kappa_x_p[i],3) )*d2A_over_dx2 ;
                    source_term_x = source_term_x - kappa_prime_x_p[i]*dA_over_dx ;
                    source_term_x = source_term_x - pow(kappa_x_p[i],3)*0.5*( ( *u1_np1_x_pml )( i, j ) + ( *u1_nm1_x_pml )( i, j ) ) ;
                    source_term_x = dt*dt*source_term_x / pow(kappa_x_p[i],3) ;
                    // ----
                    source_term_y = ( kappa_y_p[j] - pow(kappa_y_p[j],3) )*d2A_over_dy2 ;
                    source_term_y = source_term_y - kappa_prime_y_p[j]*dA_over_dy ;
                    source_term_y = source_term_y - pow(kappa_y_p[j],3)*0.5*( ( *u1_np1_y_pml )( i, j ) + ( *u1_nm1_y_pml )( i, j ) ) ;
                    source_term_y = dt*dt*source_term_y / pow(kappa_y_p[j],3) ;
                    // ----
                    ( *A2D_np1_pml )( i, j ) = 1*source_term_x + 1.*source_term_y ;
                    // ( *A2D_np1_pml )( i, j ) = 0;
                    // 4.b standard envelope FDTD 
                    ( *A2D_np1_pml )( i, j ) = ( *A2D_np1_pml )( i, j ) + dt*dt*d2A_over_dy2 ;
                    ( *A2D_np1_pml )( i, j ) = ( *A2D_np1_pml )( i, j ) + dt*dt*d2A_over_dx2 ;
                    ( *A2D_np1_pml )( i, j ) = ( *A2D_np1_pml )( i, j ) + dt*dt*k0*k0*( *A2D_n_pml )( i, j ) ;
                    ( *A2D_np1_pml )( i, j ) = ( *A2D_np1_pml )( i, j ) - (1.+i1*k0*dt) * ( *A2D_nm1_pml )( i, j ) ;
                    ( *A2D_np1_pml )( i, j ) = ( *A2D_np1_pml )( i, j ) + 2.*( *A2D_n_pml )( i, j ) ;
                    ( *A2D_np1_pml )( i, j ) = ( ( 1.+i1*k0*dt) / (1.+k0*k0*dt*dt) )*( *A2D_np1_pml )( i, j );
                } // end y loop
            } // end x loop

            for( unsigned int i=0 ; i<nx_p ; i++ ) { // x loop
                for( unsigned int j=0 ; j < ny_p ; j++ ) { // y loop
                    // X-PML-ADE
                    ( *u3_nm1_x_pml )( i, j )        = 1.*( *u3_np1_x_pml )( i, j );
                    ( *u2_nm1_x_pml )( i, j )        = 1.*( *u2_np1_x_pml )( i, j );
                    ( *u1_nm1_x_pml )( i, j )        = 1.*( *u1_np1_x_pml )( i, j );
                    // Y-PML-ADE
                    ( *u3_nm1_y_pml )( i, j )        = 1.*( *u3_np1_y_pml )( i, j );
                    ( *u2_nm1_y_pml )( i, j )        = 1.*( *u2_np1_y_pml )( i, j );
                    ( *u1_nm1_y_pml )( i, j )        = 1.*( *u1_np1_y_pml )( i, j );
                    // A-field
                    ( *A2D_nm1_pml )( i, j )       = 1.*( *A2D_n_pml )( i, j );
                    ( *A2D_n_pml )( i, j )         = 1.*( *A2D_np1_pml )( i, j );
                } // end y loop
            } // end x loop
        }
    }
}
