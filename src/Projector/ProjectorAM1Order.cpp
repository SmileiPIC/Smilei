#include "ProjectorAM1Order.h"

#include <cmath>
#include <iostream>
#include <complex>
#include "dcomplex.h"
#include "ElectroMagnAM.h"
#include "cField2D.h"
#include "Particles.h"
#include "Tools.h"
#include "Patch.h"
#include "PatchAM.h"

using namespace std;

//ProjectorAM1Order is written for exclusive use of spectral solvers with non staggered grids, primal size grids with half cell length shift along r.

// ---------------------------------------------------------------------------------------------------------------------
// Constructor for ProjectorAM1Order
// ---------------------------------------------------------------------------------------------------------------------
ProjectorAM1Order::ProjectorAM1Order( Params &params, Patch *patch ) : ProjectorAM( params, patch )
{
    dt = params.timestep;
    dr = params.cell_length[1];
    dl_inv_   = 1.0/params.cell_length[0];
    dl_ov_dt  = params.cell_length[0] / params.timestep;
    dr_ov_dt  = params.cell_length[1] / params.timestep;
    dr_inv_   = 1.0 / dr;
    one_ov_dt  = 1.0 / params.timestep;
    Nmode=params.nmodes;
    i_domain_begin = patch->getCellStartingGlobalIndex( 0 );
    j_domain_begin = patch->getCellStartingGlobalIndex( 1 );

    oversizeR = params.oversize[1]; 
    nprimr = params.n_space[1] + 2*params.oversize[1] + 1;
    npriml = params.n_space[0] + 2*params.oversize[0] + 1;

    invR = &((static_cast<PatchAM *>( patch )->invR)[0]);
    invRd = &((static_cast<PatchAM *>( patch )->invRd)[0]);

    if (!params.is_pxr)
        staggered = true;
    else
        staggered = false;
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for ProjectorAM1Order
// ---------------------------------------------------------------------------------------------------------------------
ProjectorAM1Order::~ProjectorAM1Order()
{
}

// ---------------------------------------------------------------------------------------------------------------------
//! Non charge conserving projector for diags, frozen species or spectral solver - mode >= 0
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM1Order::basicForComplex( complex<double> *rhoj, Particles &particles, unsigned int ipart, unsigned int type, int imode )
{
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    
    int iloc, nr( nprimr );
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) );
    
    if( type > 0 ) { //if current density
        charge_weight *= 1./sqrt( 1.0 + particles.momentum( 0, ipart )*particles.momentum( 0, ipart )
                                  + particles.momentum( 1, ipart )*particles.momentum( 1, ipart )
                                  + particles.momentum( 2, ipart )*particles.momentum( 2, ipart ) );
        if( type == 1 ) { //if Jl
            charge_weight *= particles.momentum( 0, ipart );
        } else if( type == 2 ) { //if Jr
            charge_weight *= ( particles.momentum( 1, ipart )*particles.position( 1, ipart ) + particles.momentum( 2, ipart )*particles.position( 2, ipart ) ) / r ;
        } else { //if Jt
            charge_weight *= ( -particles.momentum( 1, ipart )*particles.position( 2, ipart ) + particles.momentum( 2, ipart )*particles.position( 1, ipart ) ) / r ;
        }
    }
    
    complex<double> e_theta = ( particles.position( 1, ipart ) + Icpx*particles.position( 2, ipart ) )/r;
    complex<double> C_m = 1.;
    if( imode > 0 ) {
        C_m = 2.;
    }
    for( unsigned int i=0; i<( unsigned int )imode; i++ ) {
        C_m *= e_theta;
    }
    
    double xpn, rpn;
    double delta;
    double Sl1[2], Sr1[2];
    
    
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = (particles.position( 0, ipart ) ) * dl_inv_ ;
    int ip = floor( xpn );
    delta  = xpn - ( double )ip;
    Sl1[0] = 1. - delta ;
    Sl1[1] = delta;

    rpn = r * dr_inv_ - 0.5 ; //-0.5 because of cells being shifted by dr/2
    int jp = floor( rpn );
    delta  = rpn - ( double )jp;
    Sr1[0] = 1. - delta;
    Sr1[1] = delta;
   
    if (rpn < 0.){ // If particle is between 0 and dr/2.
        jp = 0;
        Sr1[1] = 0.; 
        if(imode%2 == 0){ 
            Sr1[0] = 1.;
        } else {
            Sr1[0] = 2.*delta - 1.; 
        }
    }
 
    ip -= i_domain_begin ;
    jp -= j_domain_begin ;
    // j_domain_begin should always be zero in spectral since no paralellization along r.
    
    for( unsigned int i=0 ; i<2 ; i++ ) {
        iloc = ( i+ip )*nr+jp;
        for( unsigned int j=0 ; j<2 ; j++ ) {
            rhoj [iloc+j] += C_m*charge_weight* Sl1[i]*Sr1[j] * invR[j+jp];
        }
    }//i
} // END Project for diags local current densities

// ---------------------------------------------------------------------------------------------------------------------
//! Project local currents for all modes, not charge conserving
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM1Order::currents( ElectroMagnAM *emAM, Particles &particles, unsigned int ipart, double invgf, bool diag_flag, int ispec)
{

    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    int iloc, jloc, linindex;
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    
    // variable declaration
    double xpn, ypn;
    double delta;

    double  Sl1[2], Sr1[2],Sl1d[2], Sr1d[2];
    complex<double> e_theta, C_m = 1.; 
    complex<double> *Jl, *Jr, *Jt, *rho;
    
    double rp = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) );
    double crt_p = charge_weight * ( particles.momentum( 2, ipart )*particles.position( 1, ipart )-particles.momentum( 1, ipart )*particles.position( 2, ipart ) )/( rp )*invgf;
    double crl_p = charge_weight * ( particles.momentum( 0, ipart )) *invgf;
    double crr_p = charge_weight * ( particles.momentum( 1, ipart )*particles.position( 1, ipart ) + particles.momentum( 2, ipart )*particles.position( 2, ipart ))/rp*invgf;
    e_theta = ( particles.position( 1, ipart ) + Icpx*particles.position( 2, ipart ) )/rp;
    
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position( 0, ipart ) * dl_inv_ ;
    int ip = floor( xpn );
    delta  = xpn - ( double )ip;
    Sl1[0] = 1. - delta;
    Sl1[1] = delta;
    
    ypn = rp *dr_inv_ -0.5 ; //-0.5 because grid is shifted by dr/2
    int jp = floor( ypn );
    delta  = ypn - ( double )jp;
    Sr1[0] = 1. - delta;
    Sr1[1] = delta;

    if (ypn < 0.){ // If particle is between 0 and dr/2.
        jp = 0;
        Sr1[1] = 0.; // Only account for deposition above axis. Symetry is handled in interpolation.
    }


    ip  -= i_domain_begin ;
    jp  -= j_domain_begin ;

    double *invR_local  = &(invR[jp]);

    for( unsigned int imode=0; imode<( unsigned int )Nmode; imode++ ) {
        if( imode == 1 ) {
            C_m = 2.;
        }
        if( imode > 0 ) {
            C_m *= e_theta;
        }

        if (ypn < 0.){ // If particle is between 0 and dr/2.
            if(imode%2 == 0){ 
                Sr1[0] = 1.;
            } else {
                Sr1[0] = 2.*delta - 1.; 
            }
        }
        
        if (!diag_flag){
            Jl =  &( *emAM->Jl_[imode] )( 0 );
            Jr =  &( *emAM->Jr_[imode] )( 0 );
            Jt =  &( *emAM->Jt_[imode] )( 0 );
            rho = &( *emAM->rho_AM_[imode] )( 0 ) ; // In spectral, always project density
        } else {
            unsigned int n_species = emAM->Jl_.size() / Nmode;
            unsigned int ifield = imode*n_species+ispec;
            Jl  = emAM->Jl_s    [ifield] ? &( * ( emAM->Jl_s    [ifield] ) )( 0 ) : &( *emAM->Jl_    [imode] )( 0 ) ;
            Jr  = emAM->Jr_s    [ifield] ? &( * ( emAM->Jr_s    [ifield] ) )( 0 ) : &( *emAM->Jr_    [imode] )( 0 ) ;
            Jt  = emAM->Jt_s    [ifield] ? &( * ( emAM->Jt_s    [ifield] ) )( 0 ) : &( *emAM->Jt_    [imode] )( 0 ) ;
            rho = emAM->rho_AM_s[ifield] ? &( * ( emAM->rho_AM_s[ifield] ) )( 0 ) : &( *emAM->rho_AM_[imode] )( 0 ) ;

        }

        // Rho
        for( unsigned int i=0 ; i<2 ; i++ ) {
            iloc = ( i+ip )*nprimr+jp;
            for( unsigned int j=0 ; j<2 ; j++ ) {
                linindex = iloc+j;
                rho [linindex] += C_m*charge_weight* Sl1[i]*Sr1[j]*invR_local[j];
            }
        }//i
        
        // Jl^(d,p)
        for( unsigned int i=0 ; i<2 ; i++ ) {
            iloc = ( i+ip )*nprimr + jp;
            for( unsigned int j=0 ; j<2 ; j++ ) {
                linindex = iloc+j;
                Jl [linindex] += C_m * crl_p* Sl1[i]*Sr1[j]*invR_local[j] ;
            }
        }//i
        // Jr^(p,d)
        for( unsigned int i=0 ; i<2 ; i++ ) {
            iloc = ( i+ip )* nprimr + jp;
            for( unsigned int j=0 ; j<2 ; j++ ) {
                linindex = iloc+j;
                Jr [linindex] += C_m * crr_p* Sl1[i]*Sr1[j]*invR_local[j] ;
            }
        }//i
        // Jt^(p,p)
        for( unsigned int i=0 ; i<2 ; i++ ) {
            iloc = ( i+ip )*nprimr + jp;
            for( unsigned int j=0 ; j<2 ; j++ ) {
                linindex = iloc+j;
                Jt [linindex] += C_m * crt_p* Sl1[i]*Sr1[j]*invR_local[j] ;
            }
        }
    }// end loop on modes
    
} // END Project local current densities (Jl, Jr, Jt)

//------------------------------------//
//Wrapper for projection
void ProjectorAM1Order::currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int icell, int ipart_ref )
{
        
    std::vector<double> *invgf = &( smpi->dynamics_invgf[ithread] );
    
    ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( EMfields );

    for( int ipart=istart ; ipart<iend; ipart++ ) {
        currents( emAM, particles,  ipart, ( *invgf )[ipart], diag_flag, ispec);
    }
}

// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities : ionization NOT DONE YET
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM1Order::ionizationCurrents( Field *Jl, Field *Jr, Field *Jt, Particles &particles, int ipart, LocalFields Jion )
{
}
