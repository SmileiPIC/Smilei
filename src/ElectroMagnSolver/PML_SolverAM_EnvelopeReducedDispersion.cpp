#include "PML_SolverAM_EnvelopeReducedDispersion.h"
#include "ElectroMagnAM.h"
#include "EnvelopeBCAM_PML.h"
#include "LaserEnvelope.h"
#include "cField2D.h"
#include <complex>
#include "dcomplex.h"
#include <algorithm>
#include "Patch.h"

PML_SolverAM_EnvelopeReducedDispersion::PML_SolverAM_EnvelopeReducedDispersion( Params &params )
    : SolverAM( params )
{
    // delta = 0. ;
    // For now, there is stability issue with Reduced-Dispersion Stencil. So delta = 1.
    // All other parameter are the same th standard stencil
    // delta = (1.-(dt/dl)*(dt/dl))/3. ;
    delta = (1.-(dt/dl)*(dt/dl))/3. ; delta *= 1.0 ;

    // Parameters fixed for :
    // dx = 1. and dy=4*dx
    // dt = 0.8*dx
    // 20 cells PMLs

    /*
     * If you decrease dx (/2), then you have to multiply sigma_max by the same factor (*2) in order to keep the same efficiency
     * If you use step smaller than dx=1, then you can use greater sigma without losing stability
     * If you use step smaller than dx=1, then you can use dt -> dx
     */

    // Absorption efficiency in x for 1 pass in PML A -> A/100

    //Define here the value of coefficient kappa_l_max, power_kappa_l, sigma_l_max, power_sigma_l
    // Vaccum
    // kappa_l_max = 1.0 ;
    // sigma_l_max = 0.0 ;
    // alpha_l_max = 0.0 ;
    // alpha_cl    = 0.0 ;
    // power_pml_kappa_l = 1.;
    // power_pml_sigma_l = 1.;
    // power_pml_alpha_l = 1.;
    // Abs
    kappa_cl          = params.envelope_pml_kappa_parameters[0][0] ;
    kappa_l_max       = params.envelope_pml_kappa_parameters[0][1] ;
    power_pml_kappa_l = params.envelope_pml_kappa_parameters[0][2] ;
    sigma_l_max       = params.envelope_pml_sigma_parameters[0][0] ;
    power_pml_sigma_l = params.envelope_pml_sigma_parameters[0][1] ;
    alpha_cl          = params.envelope_pml_alpha_parameters[0][0] ;
    alpha_l_max       = params.envelope_pml_alpha_parameters[0][1] ;
    power_pml_alpha_l = params.envelope_pml_alpha_parameters[0][2] ;
    
    /*
    For the moment x-PML are force to be vacuum, with Dirichlet condition at x-max and x-min
    Due to stability issues. There is no PML for x for now.
    */

    //Define here the value of coefficient kappa_y_max, power_kappa_y, sigma_y_max, power_sigma_y
    // Vaccum
    // kappa_r_max = 1. ;
    // sigma_r_max = 0. ;
    // alpha_r_max = 0.0 ;
    // alpha_cr    = 0.0 ;
    // power_pml_kappa_r = 1.;
    // power_pml_sigma_r = 1.;
    // power_pml_alpha_r = 1.;
    // Abs
    kappa_cr          = params.envelope_pml_kappa_parameters[1][0] ; // Default 1 int order to have box/PML interface
    kappa_r_max       = params.envelope_pml_kappa_parameters[1][1] ; // Increase kappa can increase stability and absorbtion
    power_pml_kappa_r = params.envelope_pml_kappa_parameters[1][2] ; // 2,3,4... 3 by default
    sigma_r_max       = params.envelope_pml_sigma_parameters[1][0] ; // 10 for 20 cells dx = 1 ; 20 for 20 cells dx = 0.5
    power_pml_sigma_r = params.envelope_pml_sigma_parameters[1][1] ; // 2 is a good choice
    alpha_cr          = params.envelope_pml_alpha_parameters[1][0] ; // 0.05 seems sufficient to ensure stability for "long-time"
    alpha_r_max       = params.envelope_pml_alpha_parameters[1][1] ; // Same as alpha_cr to ensure alpha is constant across PML 
    power_pml_alpha_r = params.envelope_pml_alpha_parameters[1][2] ; // No importance
}

PML_SolverAM_EnvelopeReducedDispersion::~PML_SolverAM_EnvelopeReducedDispersion(){
}

void PML_SolverAM_EnvelopeReducedDispersion::operator()( ElectroMagn * )
{
    ERROR( "This is not a solver for the main domain" );

    //cField2D *A_n        = static_cast<cField2D *>( envelope->A_ );   // the envelope at timestep n
    //cField2D *A_nm1      = static_cast<cField2D *>( envelope->A0_ );  // the envelope at timestep n-1
}

void PML_SolverAM_EnvelopeReducedDispersion::setDomainSizeAndCoefficients( int iDim, int min_or_max, std::vector<unsigned int> dimPrim, int ncells_pml_domain, int startpml, int* ncells_pml_min, int* ncells_pml_max, Patch* patch )
{
    const unsigned int nl_p = dimPrim[0];
    const unsigned int nr_p = dimPrim[1];
    
    if( iDim == 0 ) {
        j_glob_pml = patch->getCellStartingGlobalIndex( 1 );
    } else if( iDim == 1 ) {
        j_glob_pml = patch->getCellStartingGlobalIndex( 1 )+patch->size_[1] + patch->oversize[1] - 1; // For norder=4
    }

    isYmin = (patch->isYmin());

    //PML Coeffs Kappa,Sigma ...
    //Primal
    kappa_l_p.resize( nl_p );
    sigma_l_p.resize( nl_p );
    alpha_l_p.resize( nl_p );
    kappa_prime_l_p.resize( nl_p );
    sigma_prime_l_p.resize( nl_p );
    alpha_prime_l_p.resize( nl_p );

    kappa_r_p.resize( nr_p );
    sigma_r_p.resize( nr_p );
    alpha_r_p.resize( nr_p );
    kappa_prime_r_p.resize( nr_p );
    sigma_prime_r_p.resize( nr_p );
    alpha_prime_r_p.resize( nr_p );
    integrate_kappa_r_p.resize( nr_p ) ;
    integrate_sigma_r_p.resize( nr_p ) ;
    integrate_alpha_r_p.resize( nr_p ) ;

    // While min and max of each dim have got the same ncells_pml, coefficient are the same
    // for min and max PML

    rmax = patch->getDomainLocalMax( 1 ) ;
    // std::cout << rmax << std::endl;
    //std::cout << '('<< iDim << ',' <<j_glob_pml*dr << ')' << std::endl;
    r0 = rmax + (patch->oversize[1] + 1. )*dr ;

    if ( iDim == 0 ) {
        // 3 cells (oversize) are vaccum so the PML media begin at r0 which is :
        // Eventually the size of PML media is :
        length_r_pml = 0 ;
        length_l_pml = (ncells_pml_domain-startpml+0.5)*dl ;
        // Primal grid
        // Longitudinal
        // Params for first cell of PML-patch (vacuum) i = 0,1,2
        for ( int i=0 ; i<startpml ; i++ ) {
            // Coeffs for the first cell
            kappa_l_p[i] = 1. ;
            sigma_l_p[i] = 0. ;
            alpha_l_p[i] = 0. ;
            kappa_prime_l_p[i] = 0. ;
            sigma_prime_l_p[i] = 0. ;
            alpha_prime_l_p[i] = 0. ;
        }
        // Params for other cells (PML Media) when i>=3
        for ( int i=startpml; i< (int) nl_p ; i++ ) {
            // Parameters
            kappa_l_p[i] = kappa_cl + (kappa_l_max - kappa_cl) * pow( (i-startpml)*dl , power_pml_kappa_l ) / pow( length_l_pml , power_pml_kappa_l ) ;
            sigma_l_p[i] = sigma_l_max * pow( (i-startpml)*dl , power_pml_sigma_l ) / pow( length_l_pml , power_pml_sigma_l ) ;
            alpha_l_p[i] = alpha_cl + (alpha_l_max-alpha_cl) * pow(1. - (i-startpml)*dl/length_l_pml  , power_pml_alpha_l );
            // Derivatives
            kappa_prime_l_p[i] = (kappa_l_max - kappa_cl) * power_pml_kappa_l * pow( (i-startpml)*dl , power_pml_kappa_l-1 ) / pow( length_l_pml , power_pml_kappa_l ) ;
            sigma_prime_l_p[i] = sigma_l_max * power_pml_sigma_l * pow( (i-startpml)*dl , power_pml_sigma_l-1 ) / pow( length_l_pml , power_pml_sigma_l ) ;
            alpha_prime_l_p[i] = - power_pml_alpha_l*(alpha_l_max-alpha_cl)/length_l_pml* pow(1. - (i-startpml)*dl/length_l_pml , power_pml_alpha_l-1 );
        }
        if (min_or_max==0) {
            std::reverse(kappa_l_p.begin(), kappa_l_p.end());
            std::reverse(sigma_l_p.begin(), sigma_l_p.end());
            std::reverse(alpha_l_p.begin(), alpha_l_p.end());
            std::reverse(kappa_prime_l_p.begin(), kappa_prime_l_p.end());
            std::reverse(sigma_prime_l_p.begin(), sigma_prime_l_p.end());
            std::reverse(alpha_prime_l_p.begin(), alpha_prime_l_p.end());
            for (unsigned int i=0 ; i<nl_p ; i++){
                // Due to SMILEI convention for propagating wave
                // kappa_l_p[i] *= +1;
                // sigma_l_p[i] *= -1;
                // alpha_l_p[i] *= -1;
                kappa_l_p[i] = +1;
                sigma_l_p[i] = 0;
                alpha_l_p[i] = 0;
                // Due to SMILEI convention for propagating wave
                // For the envelope it's not a good solution for the min value !
                kappa_prime_l_p[i] = 0;
                sigma_prime_l_p[i] = 0;
                alpha_prime_l_p[i] = 0;
                // kappa_prime_l_p[i] *= -1.;
                // sigma_prime_l_p[i] *= +1.;
                // alpha_prime_l_p[i] *= +1.;
            }
        }
        if (min_or_max==1) {
            for (unsigned int i=0 ; i<nl_p ; i++){
                // Due to SMILEI convention for propagating wave
                // kappa_l_p[i] *= +1;
                // sigma_l_p[i] *= -1;
                // alpha_l_p[i] *= -1;
                kappa_l_p[i] = 1;
                sigma_l_p[i] = 0;
                alpha_l_p[i] = 0;
                // Due to SMILEI convention for propagating wave
                kappa_prime_l_p[i] = 0;
                sigma_prime_l_p[i] = 0;
                alpha_prime_l_p[i] = 0;
                // kappa_prime_l_p[i] *= +1;
                // sigma_prime_l_p[i] *= -1;
                // alpha_prime_l_p[i] *= -1;
            }
        }
        // Radial direction
        for ( unsigned int j=0 ; j<nr_p ; j++ ) {
            kappa_r_p[j] = 1. ;
            sigma_r_p[j] = 0. ;
            alpha_r_p[j] = 0. ;
            kappa_prime_r_p[j] = 0. ;
            sigma_prime_r_p[j] = 0. ;
            alpha_prime_r_p[j] = 0. ;
            integrate_kappa_r_p[j] = 1*( rmax + j*dr - r0  - 1.*dr) ;
            integrate_sigma_r_p[j] = 0. ;
            integrate_alpha_r_p[j] = 0. ;
        }
    }
    if ( iDim == 1 ) {
        // 3 cells are vaccum so the PML media begin at r0 which is :
        // Eventually the size of PML media is :
        length_r_pml = (ncells_pml_domain-startpml+0.5)*dr ;
        length_l_pml_lmax = (ncells_pml_max[0]+0.5)*dl ;
        length_l_pml_lmin = (ncells_pml_min[0]+0.5)*dl ;
        for ( unsigned int i=0 ; i<nl_p ; i++ ) {
            kappa_l_p[i] = 1. ;
            sigma_l_p[i] = 0. ;
            alpha_l_p[i] = 0. ;
            kappa_prime_l_p[i] = 0. ;
            sigma_prime_l_p[i] = 0. ;
            alpha_prime_l_p[i] = 0. ;
        }
        if (ncells_pml_min[0] != 0 ){
            for ( int i=0 ; i<ncells_pml_min[0] ; i++ ) {
                // Parameters
                kappa_l_p[i] = kappa_cl + (kappa_l_max - kappa_cl) * pow( ( ncells_pml_min[0] - 1 - i )*dl , power_pml_kappa_l ) / pow( length_l_pml_lmin , power_pml_kappa_l ) ;
                sigma_l_p[i] = sigma_l_max * pow( ( ncells_pml_min[0] - 1 - i )*dl , power_pml_sigma_l ) / pow( length_l_pml_lmin , power_pml_sigma_l ) ;
                alpha_l_p[i] = alpha_cl + (alpha_l_max-alpha_cl) * pow (1. - ( ncells_pml_min[0] - 1 - i )*dl/length_l_pml_lmin , power_pml_alpha_l );
                // Derivatives
                kappa_prime_l_p[i] = (kappa_l_max - kappa_cl) * power_pml_kappa_l * pow( ( ncells_pml_min[0] - 1 - i )*dl , power_pml_kappa_l-1 ) / pow( length_l_pml_lmin , power_pml_kappa_l ) ;
                sigma_prime_l_p[i] = sigma_l_max * power_pml_sigma_l * pow( ( ncells_pml_min[0] - 1 - i )*dl , power_pml_sigma_l-1 ) / pow( length_l_pml_lmin , power_pml_sigma_l ) ;
                alpha_prime_l_p[i] = -(alpha_l_max-alpha_cl) * power_pml_alpha_l / length_l_pml_lmin * pow( 1. - ( ncells_pml_min[0] - 1 - i )*dl/length_l_pml_lmin , power_pml_alpha_l-1 ) ;
                // Convention Envelop Smilei
                // kappa_l_p[i] *= +1 ;
                // sigma_l_p[i] *= -1 ;
                // alpha_l_p[i] *= -1 ;
                kappa_l_p[i] = +1 ;
                sigma_l_p[i] = 0 ;
                alpha_l_p[i] = 0 ;
                // kappa_prime_l_p[i] *= -1 ;
                // sigma_prime_l_p[i] *= +1 ;
                // alpha_prime_l_p[i] *= +1 ;
                kappa_prime_l_p[i] *= 0 ;
                sigma_prime_l_p[i] *= 0 ;
                alpha_prime_l_p[i] *= 0 ;
            }
        }
        if (ncells_pml_max[0] != 0 ){
            for ( int i=(nl_p-1)-(ncells_pml_max[0]-1) ; i< (int) nl_p ; i++ ) {
                // Parameters
                kappa_l_p[i] = kappa_cl + (kappa_l_max - kappa_cl) * pow( ( i - ( (nl_p-1)-(ncells_pml_max[0]-1) ) )*dl , power_pml_kappa_l ) / pow( length_l_pml_lmax , power_pml_kappa_l ) ;
                sigma_l_p[i] = sigma_l_max * pow( (i - ( (nl_p-1)-(ncells_pml_max[0]-1) ) )*dl , power_pml_sigma_l ) / pow( length_l_pml_lmax, power_pml_sigma_l ) ;
                alpha_l_p[i] = alpha_cl + (alpha_l_max-alpha_cl) * pow(1. - ( i - ( (nl_p-1)-(ncells_pml_max[0]-1) ) )*dl/length_l_pml_lmax , power_pml_alpha_l );
                // Derivatives
                kappa_prime_l_p[i] = (kappa_l_max - kappa_cl) * power_pml_kappa_l * pow( ( i - ( (nl_p-1)-(ncells_pml_max[0]-1) ) )*dl , power_pml_kappa_l-1 ) / pow( length_l_pml_lmax , power_pml_kappa_l ) ;
                sigma_prime_l_p[i] = sigma_l_max * power_pml_sigma_l * pow( ( i - ( (nl_p-1)-(ncells_pml_max[0]-1) ) )*dl , power_pml_sigma_l-1 ) / pow( length_l_pml_lmax , power_pml_sigma_l ) ;
                alpha_prime_l_p[i] = -(alpha_l_max-alpha_cl) * power_pml_alpha_l / length_l_pml_lmax * pow( 1. - ( i - ( (nl_p-1)-(ncells_pml_max[0]-1) ) )*dl / length_l_pml_lmax , power_pml_alpha_l-1 ) ;
                // Convention Envelop Smilei
                // kappa_l_p[i] *= +1 ;
                // sigma_l_p[i] *= -1 ;
                // alpha_l_p[i] *= -1 ;
                kappa_l_p[i] = +1 ;
                sigma_l_p[i] = 0 ;
                alpha_l_p[i] = 0 ;
                // kappa_prime_l_p[i] *= +1 ;
                // sigma_prime_l_p[i] *= -1 ;
                // alpha_prime_l_p[i] *= -1 ;
                kappa_prime_l_p[i] = 0 ;
                sigma_prime_l_p[i] = 0 ;
                alpha_prime_l_p[i] = 0 ;
            }
        }
        // R-direction
        for ( int j=0 ; j<startpml ; j++ ) {
            // Coeffs for the first cell
            kappa_r_p[j] = kappa_cr ;
            sigma_r_p[j] = 0. ;
            alpha_r_p[j] = 0. ;
            kappa_prime_r_p[j] = 0. ;
            sigma_prime_r_p[j] = 0. ;
            alpha_prime_r_p[j] = 0. ;
            integrate_kappa_r_p[j] = kappa_cr*( rmax + j*dr - r0 - 1.*dr ) ;
            integrate_sigma_r_p[j] = 0. ;
            integrate_alpha_r_p[j] = 0. ;
        }
        // Params for other cells (PML Media) when i>=3
        for ( int j=startpml; j< (int) nr_p ; j++ ) {
            // Parameters
            kappa_r_p[j] = kappa_cr + (kappa_r_max - kappa_cr) * pow( (j-startpml)*dr , power_pml_kappa_r ) / pow( length_r_pml , power_pml_kappa_r ) ;
            sigma_r_p[j] = sigma_r_max * pow( (j-startpml)*dr , power_pml_sigma_r ) / pow( length_r_pml , power_pml_sigma_r ) ;
            alpha_r_p[j] = alpha_cr + (alpha_r_max-alpha_cr) * pow(1. - (j-startpml)*dr / length_r_pml , power_pml_alpha_r );
            // Derivatives
            kappa_prime_r_p[j] = (kappa_r_max-kappa_cr) * power_pml_kappa_r * pow( (j-startpml)*dr , power_pml_kappa_r-1 ) / pow( length_r_pml , power_pml_kappa_r ) ;
            sigma_prime_r_p[j] = sigma_r_max * power_pml_sigma_r * pow( (j-startpml)*dr , power_pml_sigma_r-1 ) / pow( length_r_pml , power_pml_sigma_r ) ;
            alpha_prime_r_p[j] = -(alpha_r_max-alpha_cr) * power_pml_alpha_r / length_r_pml * pow( 1. - (j-startpml)*dr / length_r_pml , power_pml_alpha_r-1 ) ;
            // Integrates
            integrate_kappa_r_p[j] = kappa_cr *( rmax + j*dr - r0 - 1.*dr ) + (kappa_r_max - kappa_cr ) / pow( length_r_pml , power_pml_kappa_r ) * pow( (j-startpml)*dr , power_pml_kappa_r+1 ) / (power_pml_kappa_r+1) ;
            integrate_sigma_r_p[j] = sigma_r_max / pow( length_r_pml , power_pml_sigma_r ) * pow( (j-startpml)*dr , power_pml_sigma_r+1 ) / ( power_pml_sigma_r+1 ) ;
            integrate_alpha_r_p[j] = 1*alpha_r_p[j] ;
        }
        if (min_or_max==0) {
            std::reverse(kappa_r_p.begin(), kappa_r_p.end());
            std::reverse(sigma_r_p.begin(), sigma_r_p.end());
            std::reverse(alpha_r_p.begin(), alpha_r_p.end());
            std::reverse(kappa_prime_r_p.begin(), kappa_prime_r_p.end());
            std::reverse(sigma_prime_r_p.begin(), sigma_prime_r_p.end());
            std::reverse(alpha_prime_r_p.begin(), alpha_prime_r_p.end());
            for (unsigned int j=0 ; j<nr_p ; j++){
                // Due to SMILEI convention for propagating wave
                kappa_r_p[j] *= +1;
                sigma_r_p[j] *= -1;
                alpha_r_p[j] *= -1;
                // Due to SMILEI convention for propagating wave
                kappa_prime_r_p[j] *= -1.;
                sigma_prime_r_p[j] *= +1.;
                alpha_prime_r_p[j] *= +1.;
                // kappa_prime_r_p[j] *= +1.;
                // sigma_prime_r_p[j] *= -1.;
                // alpha_prime_r_p[j] *= -1.;
                // Due to SMILEI convention for propagating wave
                integrate_kappa_r_p[j] *= -1;
                integrate_sigma_r_p[j] *= +1;
                integrate_alpha_r_p[j] *= +1;
            }
        }
        if (min_or_max==1) {
            for (unsigned int j=0 ; j<nr_p ; j++){
                // Due to SMILEI convention for propagating wave
                kappa_r_p[j] *= +1;
                sigma_r_p[j] *= -1;
                alpha_r_p[j] *= -1;
                // Due to SMILEI convention for propagating wave
                kappa_prime_r_p[j] *= +1;
                sigma_prime_r_p[j] *= -1;
                alpha_prime_r_p[j] *= -1;
                // Due to SMILEI convention for propagating wave
                integrate_kappa_r_p[j] *= +1;
                integrate_sigma_r_p[j] *= -1;
                integrate_alpha_r_p[j] *= -1;
            }
        }
    }
}

void PML_SolverAM_EnvelopeReducedDispersion::compute_A_from_G( LaserEnvelope *envelope, int iDim, int min_or_max, std::vector<unsigned int> dimPrim, unsigned int solvermin, unsigned int solvermax )
{
    const unsigned int nl_p = dimPrim[0];
    const unsigned int nr_p = dimPrim[1];
    EnvelopeBCAM_PML* pml_fields = static_cast<EnvelopeBCAM_PML*>( envelope->EnvBoundCond[iDim*2+min_or_max] );

    cField2D* A_nm1_pml = NULL;
    cField2D* G_nm1_pml = NULL;
    cField2D* u1_nm1_l_pml = NULL;
    cField2D* u2_nm1_l_pml = NULL;
    cField2D* u3_nm1_l_pml = NULL;
    cField2D* u1_nm1_r_pml = NULL;
    cField2D* u2_nm1_r_pml = NULL;
    cField2D* u3_nm1_r_pml = NULL;

    cField2D* A_n_pml = NULL;
    cField2D* G_n_pml = NULL;
    Field2D* Chi_n_pml = NULL;

    cField2D* A_np1_pml = NULL;
    cField2D* G_np1_pml = NULL;
    cField2D* u1_np1_l_pml = NULL;
    cField2D* u2_np1_l_pml = NULL;
    cField2D* u3_np1_l_pml = NULL;
    cField2D* u1_np1_r_pml = NULL;
    cField2D* u2_np1_r_pml = NULL;
    cField2D* u3_np1_r_pml = NULL;

    A_nm1_pml = pml_fields->A_nm1_;
    G_nm1_pml = pml_fields->G_nm1_;
    u1_nm1_l_pml = pml_fields->u1_nm1_l_;
    u2_nm1_l_pml = pml_fields->u2_nm1_l_;
    u3_nm1_l_pml = pml_fields->u3_nm1_l_;
    u1_nm1_r_pml = pml_fields->u1_nm1_r_;
    u2_nm1_r_pml = pml_fields->u2_nm1_r_;
    u3_nm1_r_pml = pml_fields->u3_nm1_r_;

    A_n_pml = pml_fields->A_n_;
    G_n_pml = pml_fields->G_n_;
    Chi_n_pml = pml_fields->Chi_;

    A_np1_pml = pml_fields->A_np1_;
    G_np1_pml = pml_fields->G_np1_;
    u1_np1_l_pml = pml_fields->u1_np1_l_;
    u2_np1_l_pml = pml_fields->u2_np1_l_;
    u3_np1_l_pml = pml_fields->u3_np1_l_;
    u1_np1_r_pml = pml_fields->u1_np1_r_;
    u2_np1_r_pml = pml_fields->u2_np1_r_;
    u3_np1_r_pml = pml_fields->u3_np1_r_;

    // Auxiliary Quantities
    std::complex<double> i1 = std::complex<double>( 0., 1. ); // imaginary unit
    double k0 = 1.; // laser wavenumber
    std::complex<double> source_term_x ;
    std::complex<double> source_term_y ;

    if (iDim == 0) {
        for( unsigned int k=0 ; k<1 ; k++ ) {
            // explicit solver
            for( unsigned int i=solvermin ; i<solvermax; i++ ) { // x loop
                for( unsigned int j=std::max(3*isYmin,1) ; j < nr_p-1 ; j++ ) { // y loop
                    // dA/dx = dA/dx + ik0 A
                    // r dA/dx = r dA/dx + ik0 rA <=> dG/dx = dG/dx + ik0 G
                    std::complex<double> dG_over_dx_fdtd = (1.+delta)*( ( *G_n_pml )( i+1, j )-( *G_n_pml )( i-1, j ) )/(2.*dl) - delta*( ( *G_n_pml )( i+2, j )-( *G_n_pml )( i-2, j ) )/(4.*dl) ;
                    // std::complex<double> dG_over_dx = dG_over_dx_fdtd + i1*k0*( *G_n_pml )( i, j ) ;
                    // d2A/dx^2 = d2A/dx^2 + 2ik0 dA/dx - k0^2 A
                    // r d2A/dx^2 = r d2A/dx^2 + r 2ik0 dA/dx - r k0^2 A <=> d2G/dx^2 = d2G/dx^2 + 2ik0 dG/dx - k0^2 G
                    std::complex<double> d2G_over_dx2_fdtd = (1.+delta)*( ( *G_n_pml )( i-1, j )-2.*( *G_n_pml )( i, j )+( *G_n_pml )( i+1, j ) )/(dl*dl)-delta*( ( *G_n_pml )( i-2, j )-2.*( *G_n_pml )( i, j )+( *G_n_pml )( i+2, j ) )/(4.*dl*dl) ;
                    std::complex<double> d2G_over_dx2 = d2G_over_dx2_fdtd
                                                        + 2.*i1*k0*dG_over_dx_fdtd
                                                        - k0*k0*( *G_n_pml )( i, j ) ;
                    // dA/dy = dA/dy
                    std::complex<double> dA_over_dy = ( ( *A_n_pml )( i, j+1 )-( *A_n_pml )( i, j-1 ) )/(2.*dr) ;
                    // d2G/dy^2 = d2G/dy^2
                    std::complex<double> d2G_over_dy2 = ( ( *G_n_pml )( i, j-1 )-2.*( *G_n_pml )( i, j )+( *G_n_pml )( i, j+1 ) )/(dr*dr) ;
                    // ====
                    // STD Solver for propagation in vacuum
                    // ====
                    // ( *A_np1_pml )( i, j ) = 0. ;
                    // ( *A_np1_pml )( i, j ) += d2A_over_dy2;
                    // ( *A_np1_pml )( i, j ) += dA_over_dy/( (double) ( j_glob_pml+j )*dr );
                    // ( *A_np1_pml )( i, j ) += d2A_over_dx2_fdtd;
                    // ( *A_np1_pml )( i, j ) += 2.*i1*k0*dA_over_dx_fdtd;
                    // ( *A_np1_pml )( i, j ) = ( *A_np1_pml )( i, j )*dt*dt;
                    // ( *A_np1_pml )( i, j ) += 2.*( *A_n_pml )( i, j )-(1.+i1*k0*dt)*( *A_nm1_pml )( i, j );
                    // ( *A_np1_pml )( i, j ) = ( *A_np1_pml )( i, j )*(1.+i1*k0*dt)/(1.+k0*k0*dt*dt);
                    // ====
                    // 1. update u3
                    ( *u3_np1_l_pml )( i, j ) = +kappa_prime_l_p[i]*sigma_l_p[i] ;
                    ( *u3_np1_l_pml )( i, j ) = ( *u3_np1_l_pml )( i, j ) - sigma_prime_l_p[i]*kappa_l_p[i] ;
                    ( *u3_np1_l_pml )( i, j ) = ( *u3_np1_l_pml )( i, j ) - alpha_prime_l_p[i]*kappa_l_p[i]*kappa_l_p[i] ;
                    ( *u3_np1_l_pml )( i, j ) = ( *u3_np1_l_pml )( i, j ) * -1. *  sigma_l_p[i]*sigma_l_p[i] * dG_over_dx_fdtd / (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]) ;
                    // time operation on u3 : Be carefull, u3 has to be considered like an envelop * a carrier wave
                    ( *u3_np1_l_pml )( i, j ) = (2.*dt)/(2.-dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )*( *u3_np1_l_pml )( i, j ) ;
                    ( *u3_np1_l_pml )( i, j ) = ( *u3_np1_l_pml )( i, j ) + ( 2.+dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )/( 2.-dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )*( *u3_nm1_l_pml )( i, j ) ;
                    // 2. update u2
                    ( *u2_np1_l_pml )( i, j ) = -1*(2.*sigma_prime_l_p[i]*kappa_l_p[i]+kappa_l_p[i]*kappa_l_p[i]*alpha_prime_l_p[i]-3.*kappa_prime_l_p[i]*sigma_l_p[i])*dG_over_dx_fdtd ;
                    ( *u2_np1_l_pml )( i, j ) = ( *u2_np1_l_pml )( i, j ) - sigma_l_p[i]*kappa_l_p[i]*d2G_over_dx2_fdtd ;
                    ( *u2_np1_l_pml )( i, j ) = ( *u2_np1_l_pml )( i, j ) * sigma_l_p[i] ;
                    ( *u2_np1_l_pml )( i, j ) = ( *u2_np1_l_pml )( i, j ) - (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i])*0.5*( ( *u3_np1_l_pml )( i, j ) + ( *u3_nm1_l_pml )( i, j ) ) ;
                    ( *u2_np1_l_pml )( i, j ) = ( *u2_np1_l_pml )( i, j ) / (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]) ;
                    // time operation on u2 : Be carefull, u2 has to be considered like an envelop * a carrier wave
                    ( *u2_np1_l_pml )( i, j ) = (2.*dt)/(2.-dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )*( *u2_np1_l_pml )( i, j ) ;
                    ( *u2_np1_l_pml )( i, j ) = ( *u2_np1_l_pml )( i, j ) + ( 2.+dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )/( 2.-dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )*( *u2_nm1_l_pml )( i, j ) ;
                    // 3. update u1
                    ( *u1_np1_l_pml )( i, j ) = -1.*( 3*kappa_prime_l_p[i]*sigma_l_p[i]  - sigma_prime_l_p[i]*kappa_l_p[i] ) * dG_over_dx_fdtd ;
                    ( *u1_np1_l_pml )( i, j ) = ( *u1_np1_l_pml )( i, j ) + 2.*sigma_l_p[i]*kappa_l_p[i]*d2G_over_dx2 ;
                    ( *u1_np1_l_pml )( i, j ) = ( *u1_np1_l_pml )( i, j ) + 2.*i1*k0*sigma_l_p[i]*kappa_l_p[i]*kappa_l_p[i] * dG_over_dx_fdtd ;
                    ( *u1_np1_l_pml )( i, j ) = ( *u1_np1_l_pml )( i, j ) - (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i])*0.5*( ( *u2_np1_l_pml )( i, j ) + ( *u2_nm1_l_pml )( i, j ) ) ;
                    ( *u1_np1_l_pml )( i, j ) = ( *u1_np1_l_pml )( i, j ) / (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]) ;
                    // time operation on u1 : Be carefull, u1 has to be considered like an envelop * a carrier wave
                    ( *u1_np1_l_pml )( i, j ) = (2.*dt)/(2.-dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )*( *u1_np1_l_pml )( i, j ) ;
                    ( *u1_np1_l_pml )( i, j ) = ( *u1_np1_l_pml )( i, j ) + ( 2.+dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )/( 2.-dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )*( *u1_nm1_l_pml )( i, j ) ;
                    // ----
                    // Envelop udpate with correction/source terms
                    // ----
                    // 4.a update A : Correction/source terms
                    source_term_x = ( kappa_l_p[i] - (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]) )*d2G_over_dx2_fdtd ;
                    source_term_x = source_term_x - kappa_prime_l_p[i]*dG_over_dx_fdtd ;
                    source_term_x = source_term_x + ( 2.*i1*k0*kappa_l_p[i]*kappa_l_p[i] - 2.*i1*k0*(kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]) ) * dG_over_dx_fdtd;
                    source_term_x = source_term_x - (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i])*0.5*( ( *u1_np1_l_pml )( i, j ) + ( *u1_nm1_l_pml )( i, j ) ) ;
                    source_term_x = dt*dt*source_term_x / (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]) ;
                    // // Test ADE Scheme
                    // ( *G_np1_pml )( i, j ) = 0. ;
                    ( *G_np1_pml )( i, j ) = 1.*source_term_x - dt*dt*( *G_n_pml )( i, j )*( *Chi_n_pml )(i, j) ;
                    // 4.b Envelope FDTD with intermediate variable
                    ( *G_np1_pml )( i, j ) = ( *G_np1_pml )( i, j ) + dt*dt*d2G_over_dy2 ;
                    ( *G_np1_pml )( i, j ) = ( *G_np1_pml )( i, j ) - dt*dt*dA_over_dy ;
                    ( *G_np1_pml )( i, j ) = ( *G_np1_pml )( i, j ) + dt*dt*d2G_over_dx2 ;
                    ( *G_np1_pml )( i, j ) = ( *G_np1_pml )( i, j ) + dt*dt*k0*k0*( *G_n_pml )( i, j ) ;
                    ( *G_np1_pml )( i, j ) = ( *G_np1_pml )( i, j ) - (1.+i1*k0*dt) * ( *G_nm1_pml )( i, j ) ;
                    ( *G_np1_pml )( i, j ) = ( *G_np1_pml )( i, j ) + 2.*( *G_n_pml )( i, j ) ;
                    ( *G_np1_pml )( i, j ) = ( ( 1.+i1*k0*dt) / (1.+k0*k0*dt*dt) )*( *G_np1_pml )( i, j );
                    // ----
                    ( *A_np1_pml )( i, j ) = ( *G_np1_pml )( i, j ) / ( (double) ( j_glob_pml+j )*dr ) ;
                } // end y loop
            } // end x loop

            if (isYmin){
                // For the moment, STD AM FDTD is use on the axis without G variable (No PML in this region for sure)
                // Maybe it could be interesting to find boundary axis condition for G propagation equation
                for( unsigned int i=solvermin ; i<solvermax; i++ ) {
                    unsigned int j = 2; // j_p = 2 corresponds to r=0
                    std::complex<double> dA_over_dx_fdtd = (1.+delta)*( ( *A_n_pml )( i+1, j )-( *A_n_pml )( i-1, j ) )/(2.*dl) - delta*( ( *A_n_pml )( i+2, j )-( *A_n_pml )( i-2, j ) )/(4.*dl) ;
                    // std::complex<double> dA_over_dx = dA_over_dx_fdtd + i1*k0*( *A_n_pml )( i, j ) ;
                    // d2A/dx^2 = d2A/dx^2 + 2ik0 dA/dx - k0^2 A
                    // r d2A/dx^2 = r d2A/dx^2 + r 2ik0 dA/dx - r k0^2 A <=> d2G/dx^2 = d2G/dx^2 + 2ik0 dG/dx - k0^2 G
                    std::complex<double> d2A_over_dx2_fdtd = (1.+delta)*( ( *A_n_pml )( i-1, j )-2.*( *A_n_pml )( i, j )+( *A_n_pml )( i+1, j ) )/(dl*dl)-delta*( ( *A_n_pml )( i-2, j )-2.*( *A_n_pml )( i, j )+( *A_n_pml )( i+2, j ) )/(4.*dl*dl) ;
                    std::complex<double> d2A_over_dx2 = d2A_over_dx2_fdtd
                                                        + 2.*i1*k0*dA_over_dx_fdtd
                                                        - k0*k0*( *A_n_pml )( i, j ) ;
                    // dA/dy = dA/dy
                    // d2G/dy^2 = d2G/dy^2
                    // ====
                    // STD Solver for propagation in vacuum
                    // ====
                    // ( *A_np1_pml )( i, j ) = 0. ;
                    // ( *A_np1_pml )( i, j ) += 4.*( ( *A_n_pml )( i, j+1 )-( *A_n_pml )( i, j ) )/(dr*dr) ;
                    // ( *A_np1_pml )( i, j ) += d2A_over_dx2_fdtd;
                    // ( *A_np1_pml )( i, j ) += 2.*i1*k0*dA_over_dx_fdtd;
                    // ( *A_np1_pml )( i, j ) = ( *A_np1_pml )( i, j )*dt*dt;
                    // ( *A_np1_pml )( i, j ) += 2.*( *A_n_pml )( i, j )-(1.+i1*k0*dt)*( *A_nm1_pml )( i, j );
                    // ( *A_np1_pml )( i, j ) = ( *A_np1_pml )( i, j )*(1.+i1*k0*dt)/(1.+k0*k0*dt*dt);
                    // // Update G on the axis
                    // ( *G_np1_pml )( i, j ) = 0. ;
                    // ====
                    // STD Solver for propagation in vacuum with or without PML source
                    // ====
                    // 1. update u3
                    ( *u3_np1_l_pml )( i, j ) = +kappa_prime_l_p[i]*sigma_l_p[i] ;
                    ( *u3_np1_l_pml )( i, j ) = ( *u3_np1_l_pml )( i, j ) - sigma_prime_l_p[i]*kappa_l_p[i] ;
                    ( *u3_np1_l_pml )( i, j ) = ( *u3_np1_l_pml )( i, j ) - alpha_prime_l_p[i]*kappa_l_p[i]*kappa_l_p[i] ;
                    ( *u3_np1_l_pml )( i, j ) = ( *u3_np1_l_pml )( i, j ) * -1. *  sigma_l_p[i]*sigma_l_p[i] * dA_over_dx_fdtd / (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]) ;
                    // time operation on u3 : Be carefull, u3 has to be considered like an envelop * a carrier wave
                    ( *u3_np1_l_pml )( i, j ) = (2.*dt)/(2.-dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )*( *u3_np1_l_pml )( i, j ) ;
                    ( *u3_np1_l_pml )( i, j ) = ( *u3_np1_l_pml )( i, j ) + ( 2.+dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )/( 2.-dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )*( *u3_nm1_l_pml )( i, j ) ;
                    // 2. update u2
                    ( *u2_np1_l_pml )( i, j ) = -1*(2.*sigma_prime_l_p[i]*kappa_l_p[i]+kappa_l_p[i]*kappa_l_p[i]*alpha_prime_l_p[i]-3.*kappa_prime_l_p[i]*sigma_l_p[i])*dA_over_dx_fdtd ;
                    ( *u2_np1_l_pml )( i, j ) = ( *u2_np1_l_pml )( i, j ) - sigma_l_p[i]*kappa_l_p[i]*d2A_over_dx2_fdtd ;
                    ( *u2_np1_l_pml )( i, j ) = ( *u2_np1_l_pml )( i, j ) * sigma_l_p[i] ;
                    ( *u2_np1_l_pml )( i, j ) = ( *u2_np1_l_pml )( i, j ) - (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i])*0.5*( ( *u3_np1_l_pml )( i, j ) + ( *u3_nm1_l_pml )( i, j ) ) ;
                    ( *u2_np1_l_pml )( i, j ) = ( *u2_np1_l_pml )( i, j ) / (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]) ;
                    // time operation on u2 : Be carefull, u2 has to be considered like an envelop * a carrier wave
                    ( *u2_np1_l_pml )( i, j ) = (2.*dt)/(2.-dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )*( *u2_np1_l_pml )( i, j ) ;
                    ( *u2_np1_l_pml )( i, j ) = ( *u2_np1_l_pml )( i, j ) + ( 2.+dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )/( 2.-dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )*( *u2_nm1_l_pml )( i, j ) ;
                    // 3. update u1
                    ( *u1_np1_l_pml )( i, j ) = -1.*( 3*kappa_prime_l_p[i]*sigma_l_p[i]  - sigma_prime_l_p[i]*kappa_l_p[i] ) * dA_over_dx_fdtd ;
                    ( *u1_np1_l_pml )( i, j ) = ( *u1_np1_l_pml )( i, j ) + 2.*sigma_l_p[i]*kappa_l_p[i]*d2A_over_dx2 ;
                    ( *u1_np1_l_pml )( i, j ) = ( *u1_np1_l_pml )( i, j ) + 2.*i1*k0*sigma_l_p[i]*kappa_l_p[i]*kappa_l_p[i] * dA_over_dx_fdtd ;
                    ( *u1_np1_l_pml )( i, j ) = ( *u1_np1_l_pml )( i, j ) - (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i])*0.5*( ( *u2_np1_l_pml )( i, j ) + ( *u2_nm1_l_pml )( i, j ) ) ;
                    ( *u1_np1_l_pml )( i, j ) = ( *u1_np1_l_pml )( i, j ) / (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]) ;
                    // time operation on u1 : Be carefull, u1 has to be considered like an envelop * a carrier wave
                    ( *u1_np1_l_pml )( i, j ) = (2.*dt)/(2.-dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )*( *u1_np1_l_pml )( i, j ) ;
                    ( *u1_np1_l_pml )( i, j ) = ( *u1_np1_l_pml )( i, j ) + ( 2.+dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )/( 2.-dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )*( *u1_nm1_l_pml )( i, j ) ;
                    // ----
                    // Envelop udpate with correction/source terms
                    // ----
                    // 4.a update A : Correction/source terms
                    source_term_x = ( kappa_l_p[i] - (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]) )*d2A_over_dx2_fdtd ;
                    source_term_x = source_term_x - kappa_prime_l_p[i]*dA_over_dx_fdtd ;
                    source_term_x = source_term_x + ( 2.*i1*k0*kappa_l_p[i]*kappa_l_p[i] - 2.*i1*k0*(kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]) ) * dA_over_dx_fdtd;
                    source_term_x = source_term_x - (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i])*0.5*( ( *u1_np1_l_pml )( i, j ) + ( *u1_nm1_l_pml )( i, j ) ) ;
                    source_term_x = dt*dt*source_term_x / (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]) ;
                    // ( *A_np1_pml )( i, j ) = 0. ;
                    // 4.b Envelope FDTD with intermediate variable
                    ( *A_np1_pml )( i, j ) = 1.*source_term_x - dt*dt*( *A_n_pml )( i, j )*( *Chi_n_pml )(i, j) ;
                    ( *A_np1_pml )( i, j ) += 4.*( ( *A_n_pml )( i, j+1 )-( *A_n_pml )( i, j ) )/(dr*dr) ;
                    ( *A_np1_pml )( i, j ) += d2A_over_dx2_fdtd;
                    ( *A_np1_pml )( i, j ) += 2.*i1*k0*dA_over_dx_fdtd;
                    ( *A_np1_pml )( i, j ) = ( *A_np1_pml )( i, j )*dt*dt;
                    ( *A_np1_pml )( i, j ) += 2.*( *A_n_pml )( i, j )-(1.+i1*k0*dt)*( *A_nm1_pml )( i, j );
                    ( *A_np1_pml )( i, j ) = ( *A_np1_pml )( i, j )*(1.+i1*k0*dt)/(1.+k0*k0*dt*dt);
                    // Update G on the axis
                    ( *G_np1_pml )( i, j ) = 0. ;
                }
            }

            for( unsigned int i=0 ; i<nl_p ; i++ ) { // x loop
                for( unsigned int j=0 ; j < nr_p ; j++ ) { // y loop
                    // X-PML-ADE
                    ( *u3_nm1_l_pml )( i, j )        = 1.*( *u3_np1_l_pml )( i, j );
                    ( *u2_nm1_l_pml )( i, j )        = 1.*( *u2_np1_l_pml )( i, j );
                    ( *u1_nm1_l_pml )( i, j )        = 1.*( *u1_np1_l_pml )( i, j );
                    // // Y-PML-ADE
                    // ( *u3_nm1_y_pml )( i, j )        = 1.*( *u3_np1_y_pml )( i, j );
                    // ( *u2_nm1_y_pml )( i, j )        = 1.*( *u2_np1_y_pml )( i, j );
                    // ( *u1_nm1_y_pml )( i, j )        = 1.*( *u1_np1_y_pml )( i, j );
                    // G-field
                    ( *G_nm1_pml )( i, j )       = 1.*( *G_n_pml )( i, j );
                    ( *G_n_pml )( i, j )         = 1.*( *G_np1_pml )( i, j );
                    // A-field
                    ( *A_nm1_pml )( i, j )       = 1.*( *A_n_pml )( i, j );
                    ( *A_n_pml )( i, j )         = 1.*( *A_np1_pml )( i, j );
                } // end y loop
            } // end x loop
        }
    }

    else if (iDim == 1) {
        // A (p,p,p) Remind that in PML, there no current
        for( unsigned int k=0 ; k<1 ; k++ ) {
            // explicit solver
            for( unsigned int i=2 ; i<nl_p-2; i++ ) { // x loop
                for( unsigned int j=solvermin ; j < solvermax ; j++ ) { // y loop
                    std::complex<double> dG_over_dx_fdtd = (1.+delta)*( ( *G_n_pml )( i+1, j )-( *G_n_pml )( i-1, j ) )/(2.*dl) - delta*( ( *G_n_pml )( i+2, j )-( *G_n_pml )( i-2, j ) )/(4.*dl) ;
                    // std::complex<double> dG_over_dx = dG_over_dx_fdtd + i1*k0*( *G_n_pml )( i, j ) ;
                    // d2A/dx^2 = d2A/dx^2 + 2ik0 dA/dx - k0^2 A
                    // r d2A/dx^2 = r d2A/dx^2 + r 2ik0 dA/dx - r k0^2 A <=> d2G/dx^2 = d2G/dx^2 + 2ik0 dG/dx - k0^2 G
                    std::complex<double> d2G_over_dx2_fdtd = (1.+delta)*( ( *G_n_pml )( i-1, j )-2.*( *G_n_pml )( i, j )+( *G_n_pml )( i+1, j ) )/(dl*dl)-delta*( ( *G_n_pml )( i-2, j )-2.*( *G_n_pml )( i, j )+( *G_n_pml )( i+2, j ) )/(4.*dl*dl) ;
                    std::complex<double> d2G_over_dx2 = d2G_over_dx2_fdtd
                                                        + 2.*i1*k0*dG_over_dx_fdtd
                                                        - k0*k0*( *G_n_pml )( i, j ) ;
                    // dA/dy = dA/dy
                    std::complex<double> dA_over_dy = ( ( *A_n_pml )( i, j+1 )-( *A_n_pml )( i, j-1 ) )/(2.*dr) ;
                    // dG/dy = dG/dy
                    std::complex<double> dG_over_dy = ( ( *G_n_pml )( i, j+1 )-( *G_n_pml )( i, j-1 ) )/(2.*dr) ;
                    // d2G/dy^2 = d2G/dy^2
                    std::complex<double> d2G_over_dy2 = ( ( *G_n_pml )( i, j-1 )-2.*( *G_n_pml )( i, j )+( *G_n_pml )( i, j+1 ) )/(dr*dr) ;
                    // ====
                    // STD Solver for propagation in vacuum
                    // ====
                    // ( *A_np1_pml )( i, j ) = 0. ;
                    // ( *A_np1_pml )( i, j ) += d2A_over_dy2;
                    // ( *A_np1_pml )( i, j ) += dA_over_dy/( (double) ( j_glob_pml+j )*dr );
                    // ( *A_np1_pml )( i, j ) += d2A_over_dx2_fdtd;
                    // ( *A_np1_pml )( i, j ) += 2.*i1*k0*dA_over_dx_fdtd;
                    // ( *A_np1_pml )( i, j ) = ( *A_np1_pml )( i, j )*dt*dt;
                    // ( *A_np1_pml )( i, j ) += 2.*( *A_n_pml )( i, j )-(1.+i1*k0*dt)*( *A_nm1_pml )( i, j );
                    // ( *A_np1_pml )( i, j ) = ( *A_np1_pml )( i, j )*(1.+i1*k0*dt)/(1.+k0*k0*dt*dt);
                    // ====
                    // ADE CFS M-PML
                    // ====
                    // ADE Scheme for X-PML (depend only on x (i index))
                    // Present on iDim = 1 in order to be able to treat corner
                    // In order to take account the M-PML method sigma(x) -> sigma (x,y) <=> sigma_l(x) -> sigma_l(x)+sigma_l(y) = sigma_l(x)+0.1*sigma_r(y)
                    // Derivative are on x direction, so sigma_l_prime is unchanged
                    // ====
                    // ADE Scheme for X-PML (depend only on x (i index))
                    // ====
                    // 1. update u3
                    ( *u3_np1_l_pml )( i, j ) = +kappa_prime_l_p[i]*sigma_l_p[i] ;
                    ( *u3_np1_l_pml )( i, j ) = ( *u3_np1_l_pml )( i, j ) - sigma_prime_l_p[i]*kappa_l_p[i] ;
                    ( *u3_np1_l_pml )( i, j ) = ( *u3_np1_l_pml )( i, j ) - alpha_prime_l_p[i]*kappa_l_p[i]*kappa_l_p[i] ;
                    ( *u3_np1_l_pml )( i, j ) = - ( *u3_np1_l_pml )( i, j ) *  sigma_l_p[i]*sigma_l_p[i] * dG_over_dx_fdtd / (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]) ;
                    // time operation on u3 : Be carefull, u3 has to be considered like an envelop * a carrier wave
                    ( *u3_np1_l_pml )( i, j ) = (2.*dt)/(2.-dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )*( *u3_np1_l_pml )( i, j ) ;
                    ( *u3_np1_l_pml )( i, j ) = ( *u3_np1_l_pml )( i, j ) + ( 2.+dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )/( 2.-dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )*( *u3_nm1_l_pml )( i, j ) ;
                    // 2. update u2
                    ( *u2_np1_l_pml )( i, j ) = -1*(2.*sigma_prime_l_p[i]*kappa_l_p[i]+kappa_l_p[i]*kappa_l_p[i]*alpha_prime_l_p[i]-3.*kappa_prime_l_p[i]*sigma_l_p[i])*dG_over_dx_fdtd ;
                    ( *u2_np1_l_pml )( i, j ) = ( *u2_np1_l_pml )( i, j ) - sigma_l_p[i]*kappa_l_p[i]*d2G_over_dx2_fdtd ;
                    ( *u2_np1_l_pml )( i, j ) = ( *u2_np1_l_pml )( i, j ) * sigma_l_p[i] ;
                    ( *u2_np1_l_pml )( i, j ) = ( *u2_np1_l_pml )( i, j ) - (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i])*0.5*( ( *u3_np1_l_pml )( i, j ) + ( *u3_nm1_l_pml )( i, j ) ) ;
                    ( *u2_np1_l_pml )( i, j ) = ( *u2_np1_l_pml )( i, j ) / (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]) ;
                    // time operation on u2 : Be carefull, u2 has to be considered like an envelop * a carrier wave
                    ( *u2_np1_l_pml )( i, j ) = (2.*dt)/(2.-dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )*( *u2_np1_l_pml )( i, j ) ;
                    ( *u2_np1_l_pml )( i, j ) = ( *u2_np1_l_pml )( i, j ) + ( 2.+dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )/( 2.-dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )*( *u2_nm1_l_pml )( i, j ) ;
                    // 3. update u1
                    ( *u1_np1_l_pml )( i, j ) = -1.*( 3*kappa_prime_l_p[i]*sigma_l_p[i]  - sigma_prime_l_p[i]*kappa_l_p[i] ) * dG_over_dx_fdtd ;
                    ( *u1_np1_l_pml )( i, j ) = ( *u1_np1_l_pml )( i, j ) + 2.*sigma_l_p[i]*kappa_l_p[i]*d2G_over_dx2_fdtd ;
                    ( *u1_np1_l_pml )( i, j ) = ( *u1_np1_l_pml )( i, j ) + 2.*i1*k0*sigma_l_p[i]*kappa_l_p[i]*kappa_l_p[i] * dG_over_dx_fdtd ;
                    ( *u1_np1_l_pml )( i, j ) = ( *u1_np1_l_pml )( i, j ) - (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i])*0.5*( ( *u2_np1_l_pml )( i, j ) + ( *u2_nm1_l_pml )( i, j ) ) ;
                    ( *u1_np1_l_pml )( i, j ) = ( *u1_np1_l_pml )( i, j ) / (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]) ;
                    // time operation on u1 : Be carefull, u1 has to be considered like an envelop * a carrier wave
                    ( *u1_np1_l_pml )( i, j ) = (2.*dt)/(2.-dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )*( *u1_np1_l_pml )( i, j ) ;
                    ( *u1_np1_l_pml )( i, j ) = ( *u1_np1_l_pml )( i, j ) + ( 2.+dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )/( 2.-dt*(i1*k0 + alpha_l_p[i]+sigma_l_p[i]/kappa_l_p[i] ) )*( *u1_nm1_l_pml )( i, j ) ;
                    // ====
                    // ADE Scheme for Y-PML (depend only on y (j index))
                    // ====
                    // 1. update u3
                    ( *u3_np1_r_pml )( i, j ) = +kappa_prime_r_p[j]*sigma_r_p[j] ;
                    ( *u3_np1_r_pml )( i, j ) = ( *u3_np1_r_pml )( i, j ) - sigma_prime_r_p[j]*kappa_r_p[j] ;
                    ( *u3_np1_r_pml )( i, j ) = ( *u3_np1_r_pml )( i, j ) - alpha_prime_r_p[j]*kappa_r_p[j]*kappa_r_p[j] ;
                    ( *u3_np1_r_pml )( i, j ) = - ( *u3_np1_r_pml )( i, j ) * pow(sigma_r_p[j],2) * dG_over_dy / (kappa_r_p[j]*kappa_r_p[j]*kappa_r_p[j]*kappa_r_p[j]) ;
                    // time operation on u3 : Be carefull, u3 has to be considered like an envelop * a carrier wave
                    ( *u3_np1_r_pml )( i, j ) = (2.*dt)/(2.-dt*(i1*k0 + alpha_r_p[j]+sigma_r_p[j]/kappa_r_p[j] ) )*( *u3_np1_r_pml )( i, j ) ;
                    ( *u3_np1_r_pml )( i, j ) = ( *u3_np1_r_pml )( i, j ) + ( 2.+dt*(i1*k0 + alpha_r_p[j]+sigma_r_p[j]/kappa_r_p[j] ) )/( 2.-dt*(i1*k0 + alpha_r_p[j]+sigma_r_p[j]/kappa_r_p[j] ) )*( *u3_nm1_r_pml )( i, j ) ;
                    // 2. update u2
                    ( *u2_np1_r_pml )( i, j ) = -1.*(2.*sigma_prime_r_p[j]*kappa_r_p[j]+kappa_r_p[j]*kappa_r_p[j]*alpha_prime_r_p[j]-3.*kappa_prime_r_p[j]*sigma_r_p[j])*dG_over_dy ;
                    ( *u2_np1_r_pml )( i, j ) = ( *u2_np1_r_pml )( i, j ) - sigma_r_p[j]*kappa_r_p[j]*d2G_over_dy2 ;
                    ( *u2_np1_r_pml )( i, j ) = ( *u2_np1_r_pml )( i, j ) * sigma_r_p[j] ;
                    ( *u2_np1_r_pml )( i, j ) = ( *u2_np1_r_pml )( i, j ) - (kappa_r_p[j]*kappa_r_p[j]*kappa_r_p[j])*0.5*( ( *u3_np1_r_pml )( i, j ) + ( *u3_nm1_r_pml )( i, j ) ) ;
                    ( *u2_np1_r_pml )( i, j ) = ( *u2_np1_r_pml )( i, j ) / (kappa_r_p[j]*kappa_r_p[j]*kappa_r_p[j]*kappa_r_p[j]) ;
                    // time operation on u2 : Be carefull, u2 has to be considered like an envelop * a carrier wave
                    ( *u2_np1_r_pml )( i, j ) = (2.*dt)/(2.-dt*(i1*k0 + alpha_r_p[j]+sigma_r_p[j]/kappa_r_p[j] ) )*( *u2_np1_r_pml )( i, j ) ;
                    ( *u2_np1_r_pml )( i, j ) = ( *u2_np1_r_pml )( i, j ) + ( 2.+dt*(i1*k0 + alpha_r_p[j]+sigma_r_p[j]/kappa_r_p[j] ) )/( 2.-dt*(i1*k0 + alpha_r_p[j]+sigma_r_p[j]/kappa_r_p[j] ) )*( *u2_nm1_r_pml )( i, j ) ;
                    // 3. update u1
                    ( *u1_np1_r_pml )( i, j ) = -1.*( 3*kappa_prime_r_p[j]*sigma_r_p[j] - sigma_prime_r_p[j]*kappa_r_p[j] ) * dG_over_dy ;
                    ( *u1_np1_r_pml )( i, j ) = ( *u1_np1_r_pml )( i, j ) + 2.*sigma_r_p[j]*kappa_r_p[j]*d2G_over_dy2 ;
                    ( *u1_np1_r_pml )( i, j ) = ( *u1_np1_r_pml )( i, j ) - sigma_r_p[j]*kappa_r_p[j]*kappa_r_p[j] * dA_over_dy ;
                    ( *u1_np1_r_pml )( i, j ) = ( *u1_np1_r_pml )( i, j ) - (kappa_r_p[j]*kappa_r_p[j]*kappa_r_p[j])*0.5*( ( *u2_np1_r_pml )( i, j ) + ( *u2_nm1_r_pml )( i, j ) ) ;
                    ( *u1_np1_r_pml )( i, j ) = ( *u1_np1_r_pml )( i, j ) / (kappa_r_p[j]*kappa_r_p[j]*kappa_r_p[j]*kappa_r_p[j]) ;
                    // time operation on u1 : Be carefull, u1 has to be considered like an envelop * a carrier wave
                    ( *u1_np1_r_pml )( i, j ) = (2.*dt)/(2.-dt*(i1*k0 + alpha_r_p[j]+sigma_r_p[j]/kappa_r_p[j] ) )*( *u1_np1_r_pml )( i, j ) ;
                    ( *u1_np1_r_pml )( i, j ) = ( *u1_np1_r_pml )( i, j ) + ( 2.+dt*(i1*k0 + alpha_r_p[j]+sigma_r_p[j]/kappa_r_p[j] ) )/( 2.-dt*(i1*k0 + alpha_r_p[j]+sigma_r_p[j]/kappa_r_p[j] ) )*( *u1_nm1_r_pml )( i, j ) ;
                    // ----
                    // Envelop udpate with correction/source terms
                    // ----
                    // 4.a update A : Correction/source terms
                    source_term_x = ( kappa_l_p[i] - (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]) )*d2G_over_dx2_fdtd ;
                    source_term_x = source_term_x - kappa_prime_l_p[i]*dG_over_dx_fdtd ;
                    source_term_x = source_term_x + ( 2.*i1*k0*kappa_l_p[i]*kappa_l_p[i] - 2.*i1*k0*(kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]) ) * dG_over_dx_fdtd;
                    source_term_x = source_term_x + (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i])*0.5*( ( *u1_np1_l_pml )( i, j ) + ( *u1_nm1_l_pml )( i, j ) ) ;
                    source_term_x = dt*dt*source_term_x / (kappa_l_p[i]*kappa_l_p[i]*kappa_l_p[i]) ;

                    source_term_y = ( kappa_r_p[j] - (kappa_r_p[j]*kappa_r_p[j]*kappa_r_p[j]) )*d2G_over_dy2 ;
                    source_term_y = source_term_y - kappa_prime_r_p[j]*dG_over_dy ;
                    source_term_y = source_term_y + ( (kappa_r_p[j]*kappa_r_p[j]*kappa_r_p[j]) - kappa_r_p[j]*kappa_r_p[j] )*dA_over_dy ;
                    source_term_y = source_term_y + (kappa_r_p[j]*kappa_r_p[j]*kappa_r_p[j])*0.5*( ( *u1_np1_r_pml )( i, j ) + ( *u1_nm1_r_pml )( i, j ) ) ;
                    source_term_y = dt*dt*source_term_y / (kappa_r_p[j]*kappa_r_p[j]*kappa_r_p[j]) ;
                    // ----
                    // Test ADE Scheme
                    // ( *G_np1_pml )( i, j ) = 0 ; // No decay
                    // ( *G_np1_pml )( i, j ) = 1.*source_term_y ; // Only y decay
                    ( *G_np1_pml )( i, j ) = 1.*source_term_x + 1.*source_term_y - dt*dt*( *G_n_pml )( i, j )*( *Chi_n_pml )(i, j) ;
                    // 4.b Envelope FDTD with intermediate variable
                    ( *G_np1_pml )( i, j ) = ( *G_np1_pml )( i, j ) + dt*dt*d2G_over_dy2 ;
                    ( *G_np1_pml )( i, j ) = ( *G_np1_pml )( i, j ) - dt*dt*dA_over_dy ;
                    ( *G_np1_pml )( i, j ) = ( *G_np1_pml )( i, j ) + dt*dt*d2G_over_dx2 ;
                    ( *G_np1_pml )( i, j ) = ( *G_np1_pml )( i, j ) + dt*dt*k0*k0*( *G_n_pml )( i, j ) ;
                    ( *G_np1_pml )( i, j ) = ( *G_np1_pml )( i, j ) - (1.+i1*k0*dt) * ( *G_nm1_pml )( i, j ) ;
                    ( *G_np1_pml )( i, j ) = ( *G_np1_pml )( i, j ) + 2.*( *G_n_pml )( i, j ) ;
                    ( *G_np1_pml )( i, j ) = ( ( 1.+i1*k0*dt) / (1.+k0*k0*dt*dt) )*( *G_np1_pml )( i, j );
                    // ----
                    // ( *A_np1_pml )( i, j ) = ( *G_np1_pml )( i, j ) / ( (double) ( j_glob_pml+j )*dr ) ;
                    // ( *A_np1_pml )( i, j ) = (
                    //         ( ( (double) ( j_glob_pml+j )*dr ) + dt*( i1*k0*( (double) ( j_glob_pml+j )*dr) ) )*( *A_nm1_pml )( i, j )
                    //         + ( *G_np1_pml )( i, j ) * ( 1. - dt*( i1*k0 ) )
                    //         - ( *G_nm1_pml )( i, j ) * ( 1. + dt*( i1*k0 ) )
                    //     ) / ( ( (double) ( j_glob_pml+j )*dr ) - dt*( i1*k0*( (double) ( j_glob_pml+j )*dr ) ) ) ;
                    // std::cout << ( (double) ( j_glob_pml+j )*dr ) - (r0 + integrate_kappa_r_p[j]) << std::endl;
                    // ( *A_np1_pml )( i, j ) = (
                    //         ( (r0 + integrate_kappa_r_p[j]) + dt*( i1*k0*(r0 + integrate_kappa_r_p[j]) ) )*( *A_nm1_pml )( i, j )
                    //         + ( *G_np1_pml )( i, j ) * ( 1. - dt*( i1*k0 ) )
                    //         - ( *G_nm1_pml )( i, j ) * ( 1. + dt*( i1*k0 ) )
                    //     ) / ( (r0 + integrate_kappa_r_p[j]) - dt*( i1*k0*(r0 + integrate_kappa_r_p[j]) ) ) ;
                    // std::cout << ( (double) ( j_glob_pml+j )*dr ) - (r0 + integrate_kappa_r_p[j]) << std::endl;
                    ( *A_np1_pml )( i, j ) = (
                            ( (r0 + integrate_kappa_r_p[j]) + dt*( alpha_r_p[j]*(r0 + integrate_kappa_r_p[j]) + integrate_sigma_r_p[j] + i1*k0*(r0 + integrate_kappa_r_p[j]) ) )*( *A_nm1_pml )( i, j )
                            + ( *G_np1_pml )( i, j ) * ( 1. - dt*( alpha_r_p[j] + i1*k0 ) )
                            - ( *G_nm1_pml )( i, j ) * ( 1. + dt*( alpha_r_p[j] + i1*k0 ) )
                        ) / ( (r0 + integrate_kappa_r_p[j]) - dt*( alpha_r_p[j]*(r0 + integrate_kappa_r_p[j]) + integrate_sigma_r_p[j] + i1*k0*(r0 + integrate_kappa_r_p[j]) ) ) ;
                } // end y loop
            } // end x loop

            // for( unsigned int i=2 ; i<nl_p-2; i++ ) { // x loop
            //     for( unsigned int j=solvermin ; j < solvermax ; j++ ) { // y loop
            //         // An_f = An + nu/2.*(1.*Anm1-2.*An+1.*Anp1)
            //         // With advanced scheme is instable
            //         ( *G_n_pml )( i, j ) += (0.02/2.)*(( *G_nm1_pml )( i, j ) - 2.*( *G_n_pml )( i, j ) + ( *G_np1_pml )( i, j ) );
            //         ( *A_n_pml )( i, j ) += (0.02/2.)*(( *A_nm1_pml )( i, j ) - 2.*( *A_n_pml )( i, j ) + ( *A_np1_pml )( i, j ) );
            //      }
            // }

            for( unsigned int i=0 ; i<nl_p ; i++ ) { // x loop
                for( unsigned int j=0 ; j < nr_p ; j++ ) { // y loop
                    // X-PML-ADE
                    ( *u3_nm1_l_pml )( i, j )        = 1.*( *u3_np1_l_pml )( i, j );
                    ( *u2_nm1_l_pml )( i, j )        = 1.*( *u2_np1_l_pml )( i, j );
                    ( *u1_nm1_l_pml )( i, j )        = 1.*( *u1_np1_l_pml )( i, j );
                    // Y-PML-ADE
                    ( *u3_nm1_r_pml )( i, j )        = 1.*( *u3_np1_r_pml )( i, j );
                    ( *u2_nm1_r_pml )( i, j )        = 1.*( *u2_np1_r_pml )( i, j );
                    ( *u1_nm1_r_pml )( i, j )        = 1.*( *u1_np1_r_pml )( i, j );
                    // G-field
                    ( *G_nm1_pml )( i, j )       = 1.*( *G_n_pml )( i, j );
                    ( *G_n_pml )( i, j )         = 1.*( *G_np1_pml )( i, j );
                    // A-field
                    ( *A_nm1_pml )( i, j )       = 1.*( *A_n_pml )( i, j );
                    ( *A_n_pml )( i, j )         = 1.*( *A_np1_pml )( i, j );
                } // end y loop
            } // end x loop
        }
    }
}
