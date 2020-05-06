#include "Projector1D4Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field1D.h"
#include "Particles.h"
#include "Tools.h"
#include "Patch.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Projector1D4Order
// ---------------------------------------------------------------------------------------------------------------------
Projector1D4Order::Projector1D4Order( Params &params, Patch *patch )
    : Projector1D( params, patch )
{
    dx_inv_  = 1.0/params.cell_length[0];
    dx_ov_dt = params.cell_length[0] / params.timestep;
    
    //double defined for use in coefficients
    
    index_domain_begin = patch->getCellStartingGlobalIndex( 0 );
    
    DEBUG( "cell_length "<< params.cell_length[0] );
    
}

Projector1D4Order::~Projector1D4Order()
{
}

// ---------------------------------------------------------------------------------------------------------------------
//! Project current densities : main projector
// ---------------------------------------------------------------------------------------------------------------------
void Projector1D4Order::currents( double *Jx, double *Jy, double *Jz, Particles &particles, unsigned int ipart, double invgf, int *iold, double *delta )
{
    // Declare local variables
    int ipo, ip;
    int ip_m_ipo;
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    double xjn, xj_m_xipo, xj_m_xipo2, xj_m_xipo3, xj_m_xipo4, xj_m_xip, xj_m_xip2, xj_m_xip3, xj_m_xip4;
    double crx_p = charge_weight*dx_ov_dt;                // current density for particle moving in the x-direction
    double cry_p = charge_weight*particles.momentum( 1, ipart )*invgf;  // current density in the y-direction of the macroparticle
    double crz_p = charge_weight*particles.momentum( 2, ipart )*invgf;  // current density allow the y-direction of the macroparticle
    double S0[7], S1[7], Wl[7], Wt[7], Jx_p[7];            // arrays used for the Esirkepov projection method
    // Initialize variables
    for( unsigned int i=0; i<7; i++ ) {
        S0[i]=0.;
        S1[i]=0.;
        Wl[i]=0.;
        Wt[i]=0.;
        Jx_p[i]=0.;
    }//i
    
    
    // Locate particle old position on the primal grid
    xj_m_xipo  = *delta;                   // normalized distance to the nearest grid point
    xj_m_xipo2 = xj_m_xipo  * xj_m_xipo;                 // square of the normalized distance to the nearest grid point
    xj_m_xipo3 = xj_m_xipo2 * xj_m_xipo;              // cube of the normalized distance to the nearest grid point
    xj_m_xipo4 = xj_m_xipo3 * xj_m_xipo;              // 4th power of the normalized distance to the nearest grid point
    
    // Locate particle new position on the primal grid
    xjn       = particles.position( 0, ipart ) * dx_inv_;
    ip        = round( xjn );                         // index of the central node
    xj_m_xip  = xjn - ( double )ip;                   // normalized distance to the nearest grid point
    xj_m_xip2 = xj_m_xip  * xj_m_xip;                    // square of the normalized distance to the nearest grid point
    xj_m_xip3 = xj_m_xip2 * xj_m_xip;                 // cube of the normalized distance to the nearest grid point
    xj_m_xip4 = xj_m_xip3 * xj_m_xip;                 // 4th power of the normalized distance to the nearest grid point
    
    
    // coefficients 4th order interpolation on 5 nodes
    
    S0[1] = dble_1_ov_384   - dble_1_ov_48  * xj_m_xipo  + dble_1_ov_16 * xj_m_xipo2 - dble_1_ov_12 * xj_m_xipo3 + dble_1_ov_24 * xj_m_xipo4;
    S0[2] = dble_19_ov_96   - dble_11_ov_24 * xj_m_xipo  + dble_1_ov_4 * xj_m_xipo2  + dble_1_ov_6  * xj_m_xipo3 - dble_1_ov_6  * xj_m_xipo4;
    S0[3] = dble_115_ov_192 - dble_5_ov_8   * xj_m_xipo2 + dble_1_ov_4 * xj_m_xipo4;
    S0[4] = dble_19_ov_96   + dble_11_ov_24 * xj_m_xipo  + dble_1_ov_4 * xj_m_xipo2  - dble_1_ov_6  * xj_m_xipo3 - dble_1_ov_6  * xj_m_xipo4;
    S0[5] = dble_1_ov_384   + dble_1_ov_48  * xj_m_xipo  + dble_1_ov_16 * xj_m_xipo2 + dble_1_ov_12 * xj_m_xipo3 + dble_1_ov_24 * xj_m_xipo4;
    
    // coefficients 2nd order interpolation on 5 nodes
    ipo        = *iold;                          // index of the central node
    ip_m_ipo = ip-ipo-index_domain_begin;
    
    S1[ip_m_ipo+1] = dble_1_ov_384   - dble_1_ov_48  * xj_m_xip  + dble_1_ov_16 * xj_m_xip2 - dble_1_ov_12 * xj_m_xip3 + dble_1_ov_24 * xj_m_xip4;
    S1[ip_m_ipo+2] = dble_19_ov_96   - dble_11_ov_24 * xj_m_xip  + dble_1_ov_4 * xj_m_xip2  + dble_1_ov_6  * xj_m_xip3 - dble_1_ov_6  * xj_m_xip4;
    S1[ip_m_ipo+3] = dble_115_ov_192 - dble_5_ov_8   * xj_m_xip2 + dble_1_ov_4 * xj_m_xip4;
    S1[ip_m_ipo+4] = dble_19_ov_96   + dble_11_ov_24 * xj_m_xip  + dble_1_ov_4 * xj_m_xip2  - dble_1_ov_6  * xj_m_xip3 - dble_1_ov_6  * xj_m_xip4;
    S1[ip_m_ipo+5] = dble_1_ov_384   + dble_1_ov_48  * xj_m_xip  + dble_1_ov_16 * xj_m_xip2 + dble_1_ov_12 * xj_m_xip3 + dble_1_ov_24 * xj_m_xip4;
    
    // coefficients used in the Esirkepov method
    for( unsigned int i=0; i<7; i++ ) {
        Wl[i] = S0[i] - S1[i];           // for longitudinal current (x)
        Wt[i] = 0.5 * ( S0[i] + S1[i] ); // for transverse currents (y,z)
    }//i
    
    // local current created by the particle
    // calculate using the charge conservation equation
    for( unsigned int i=1; i<7; i++ ) {
        Jx_p[i] = Jx_p[i-1] + crx_p * Wl[i-1];
    }
    
    ipo -= 3 ;
    
    // 4th order projection for the total currents & charge density
    // At the 4th order, oversize = 3.
    for( unsigned int i=0; i<7; i++ ) {
        Jx[i  + ipo]  += Jx_p[i];
        Jy[i  + ipo]  += cry_p * Wt[i];
        Jz[i  + ipo]  += crz_p * Wt[i];
    }//i
    
}


// ---------------------------------------------------------------------------------------------------------------------
//!  Project current densities & charge : diagFields timstep
// ---------------------------------------------------------------------------------------------------------------------
void Projector1D4Order::currentsAndDensity( double *Jx, double *Jy, double *Jz, double *rho, Particles &particles, unsigned int ipart, double invgf, int *iold, double *delta )
{
    // Declare local variables
    int ipo, ip;
    int ip_m_ipo;
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    double xjn, xj_m_xipo, xj_m_xipo2, xj_m_xipo3, xj_m_xipo4, xj_m_xip, xj_m_xip2, xj_m_xip3, xj_m_xip4;
    double crx_p = charge_weight*dx_ov_dt;                // current density for particle moving in the x-direction
    double cry_p = charge_weight*particles.momentum( 1, ipart )*invgf;  // current density in the y-direction of the macroparticle
    double crz_p = charge_weight*particles.momentum( 2, ipart )*invgf;  // current density allow the y-direction of the macroparticle
    double S0[7], S1[7], Wl[7], Wt[7], Jx_p[7];            // arrays used for the Esirkepov projection method
    // Initialize variables
    for( unsigned int i=0; i<7; i++ ) {
        S0[i]=0.;
        S1[i]=0.;
        Wl[i]=0.;
        Wt[i]=0.;
        Jx_p[i]=0.;
    }//i
    
    
    // Locate particle old position on the primal grid
    xj_m_xipo  = *delta;                   // normalized distance to the nearest grid point
    xj_m_xipo2 = xj_m_xipo  * xj_m_xipo;                 // square of the normalized distance to the nearest grid point
    xj_m_xipo3 = xj_m_xipo2 * xj_m_xipo;              // cube of the normalized distance to the nearest grid point
    xj_m_xipo4 = xj_m_xipo3 * xj_m_xipo;              // 4th power of the normalized distance to the nearest grid point
    
    // Locate particle new position on the primal grid
    xjn       = particles.position( 0, ipart ) * dx_inv_;
    ip        = round( xjn );                         // index of the central node
    xj_m_xip  = xjn - ( double )ip;                   // normalized distance to the nearest grid point
    xj_m_xip2 = xj_m_xip  * xj_m_xip;                    // square of the normalized distance to the nearest grid point
    xj_m_xip3 = xj_m_xip2 * xj_m_xip;                 // cube of the normalized distance to the nearest grid point
    xj_m_xip4 = xj_m_xip3 * xj_m_xip;                 // 4th power of the normalized distance to the nearest grid point
    
    
    // coefficients 4th order interpolation on 5 nodes
    
    S0[1] = dble_1_ov_384   - dble_1_ov_48  * xj_m_xipo  + dble_1_ov_16 * xj_m_xipo2 - dble_1_ov_12 * xj_m_xipo3 + dble_1_ov_24 * xj_m_xipo4;
    S0[2] = dble_19_ov_96   - dble_11_ov_24 * xj_m_xipo  + dble_1_ov_4 * xj_m_xipo2  + dble_1_ov_6  * xj_m_xipo3 - dble_1_ov_6  * xj_m_xipo4;
    S0[3] = dble_115_ov_192 - dble_5_ov_8   * xj_m_xipo2 + dble_1_ov_4 * xj_m_xipo4;
    S0[4] = dble_19_ov_96   + dble_11_ov_24 * xj_m_xipo  + dble_1_ov_4 * xj_m_xipo2  - dble_1_ov_6  * xj_m_xipo3 - dble_1_ov_6  * xj_m_xipo4;
    S0[5] = dble_1_ov_384   + dble_1_ov_48  * xj_m_xipo  + dble_1_ov_16 * xj_m_xipo2 + dble_1_ov_12 * xj_m_xipo3 + dble_1_ov_24 * xj_m_xipo4;
    
    // coefficients 2nd order interpolation on 5 nodes
    ipo        = *iold;                          // index of the central node
    ip_m_ipo = ip-ipo-index_domain_begin;
    
    S1[ip_m_ipo+1] = dble_1_ov_384   - dble_1_ov_48  * xj_m_xip  + dble_1_ov_16 * xj_m_xip2 - dble_1_ov_12 * xj_m_xip3 + dble_1_ov_24 * xj_m_xip4;
    S1[ip_m_ipo+2] = dble_19_ov_96   - dble_11_ov_24 * xj_m_xip  + dble_1_ov_4 * xj_m_xip2  + dble_1_ov_6  * xj_m_xip3 - dble_1_ov_6  * xj_m_xip4;
    S1[ip_m_ipo+3] = dble_115_ov_192 - dble_5_ov_8   * xj_m_xip2 + dble_1_ov_4 * xj_m_xip4;
    S1[ip_m_ipo+4] = dble_19_ov_96   + dble_11_ov_24 * xj_m_xip  + dble_1_ov_4 * xj_m_xip2  - dble_1_ov_6  * xj_m_xip3 - dble_1_ov_6  * xj_m_xip4;
    S1[ip_m_ipo+5] = dble_1_ov_384   + dble_1_ov_48  * xj_m_xip  + dble_1_ov_16 * xj_m_xip2 + dble_1_ov_12 * xj_m_xip3 + dble_1_ov_24 * xj_m_xip4;
    
    // coefficients used in the Esirkepov method
    for( unsigned int i=0; i<7; i++ ) {
        Wl[i] = S0[i] - S1[i];           // for longitudinal current (x)
        Wt[i] = 0.5 * ( S0[i] + S1[i] ); // for transverse currents (y,z)
    }//i
    
    // local current created by the particle
    // calculate using the charge conservation equation
    for( unsigned int i=1; i<7; i++ ) {
        Jx_p[i] = Jx_p[i-1] + crx_p * Wl[i-1];
    }
    
    ipo -= 3;
    
    // 4th order projection for the total currents & charge density
    // At the 4th order, oversize = 3.
    for( unsigned int i=0; i<7; i++ ) {
        Jx[i  + ipo ]  += Jx_p[i];
        Jy[i  + ipo ]  += cry_p * Wt[i];
        Jz[i  + ipo ]  += crz_p * Wt[i];
        rho[i  + ipo ] += charge_weight * S1[i];
    }//i
    
}//END Project local current densities (sort)


// ---------------------------------------------------------------------------------------------------------------------
//! Project charge : frozen & diagFields timstep
// ---------------------------------------------------------------------------------------------------------------------
void Projector1D4Order::basic( double *rhoj, Particles &particles, unsigned int ipart, unsigned int type )
{

    //Warning : this function is used for frozen species or initialization only and doesn't use the standard scheme.
    //rho type = 0
    //Jx type = 1
    //Jy type = 2
    //Jz type = 3
    
    
    // Declare local variables
    //int ipo, ip, iloc;
    int ip;
    //int ip_m_ipo;
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    double xjn, xj_m_xip, xj_m_xip2, xj_m_xip3, xj_m_xip4;
    double S1[7];            // arrays used for the Esirkepov projection method
    // Initialize variables
    for( unsigned int i=0; i<7; i++ ) {
        S1[i]=0.;
    }//i
    
    if( type > 0 ) {
        charge_weight *= 1./sqrt( 1.0 + particles.momentum( 0, ipart )*particles.momentum( 0, ipart )
                                  + particles.momentum( 1, ipart )*particles.momentum( 1, ipart )
                                  + particles.momentum( 2, ipart )*particles.momentum( 2, ipart ) );
                                  
        if( type == 1 ) {
            charge_weight *= particles.momentum( 0, ipart );
        } else if( type == 2 ) {
            charge_weight *= particles.momentum( 1, ipart );
        } else {
            charge_weight *= particles.momentum( 2, ipart );
        }
    }
    
    // Locate particle new position on the primal grid
    xjn       = particles.position( 0, ipart ) * dx_inv_;
    ip        = round( xjn + 0.5 * ( type==1 ) );                       // index of the central node
    xj_m_xip  = xjn - ( double )ip;                   // normalized distance to the nearest grid point
    xj_m_xip2 = xj_m_xip  * xj_m_xip;                    // square of the normalized distance to the nearest grid point
    xj_m_xip3 = xj_m_xip2 * xj_m_xip;                 // cube of the normalized distance to the nearest grid point
    xj_m_xip4 = xj_m_xip3 * xj_m_xip;                 // 4th power of the normalized distance to the nearest grid point
    
    // coefficients 2nd order interpolation on 5 nodes
    //ip_m_ipo = ip-ipo;
    
    S1[1] = dble_1_ov_384   - dble_1_ov_48  * xj_m_xip  + dble_1_ov_16 * xj_m_xip2 - dble_1_ov_12 * xj_m_xip3 + dble_1_ov_24 * xj_m_xip4;
    S1[2] = dble_19_ov_96   - dble_11_ov_24 * xj_m_xip  + dble_1_ov_4 * xj_m_xip2  + dble_1_ov_6  * xj_m_xip3 - dble_1_ov_6  * xj_m_xip4;
    S1[3] = dble_115_ov_192 - dble_5_ov_8   * xj_m_xip2 + dble_1_ov_4 * xj_m_xip4;
    S1[4] = dble_19_ov_96   + dble_11_ov_24 * xj_m_xip  + dble_1_ov_4 * xj_m_xip2  - dble_1_ov_6  * xj_m_xip3 - dble_1_ov_6  * xj_m_xip4;
    S1[5] = dble_1_ov_384   + dble_1_ov_48  * xj_m_xip  + dble_1_ov_16 * xj_m_xip2 + dble_1_ov_12 * xj_m_xip3 + dble_1_ov_24 * xj_m_xip4;
    
    ip -= index_domain_begin + 3 ;
    
    // 4th order projection for the charge density
    // At the 4th order, oversize = 3.
    for( unsigned int i=0; i<7; i++ ) {
        //iloc = i  + ipo - 3;
        rhoj[i  + ip ] += charge_weight * S1[i];
    }//i
    
}


// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities (ionize)
// ---------------------------------------------------------------------------------------------------------------------
void Projector1D4Order::ionizationCurrents( Field *Jx, Field *Jy, Field *Jz, Particles &particles, int ipart, LocalFields Jion )
{
    Field1D *Jx1D  = static_cast<Field1D *>( Jx );
    Field1D *Jy1D  = static_cast<Field1D *>( Jy );
    Field1D *Jz1D  = static_cast<Field1D *>( Jz );
    
    
    //Declaration of local variables
    int i, im2, im1, ip1, ip2;
    double xjn, xjmxi, xjmxi2, xjmxi3, xjmxi4;
    double cim2, cim1, ci, cip1, cip2;
    
    // weighted currents
    double weight = inv_cell_volume * particles.weight( ipart );
    double Jx_ion = Jion.x * weight;
    double Jy_ion = Jion.y * weight;
    double Jz_ion = Jion.z * weight;
    
    //Locate particle on the grid
    xjn    = particles.position( 0, ipart ) * dx_inv_; // normalized distance to the first node
    
    
    // Compute Jx_ion on the dual grid
    // -------------------------------
    
    i      = round( xjn+0.5 );             // index of the central node
    xjmxi  = xjn - ( double )i + 0.5;      // normalized distance to the nearest grid point
    xjmxi2 = xjmxi*xjmxi;                  // square of the normalized distance to the nearest grid point
    xjmxi3 = xjmxi2*xjmxi;                 // cube
    xjmxi4 = xjmxi2*xjmxi2;                 // fourth-power
    
    i  -= index_domain_begin;
    im2 = i-2;
    im1 = i-1;
    ip1 = i+1;
    ip2 = i+2;

    cim2 = dble_1_ov_384   - dble_1_ov_48  * xjmxi  + dble_1_ov_16 * xjmxi2 - dble_1_ov_12 * xjmxi3 + dble_1_ov_24 * xjmxi4;
    cim1 = dble_19_ov_96   - dble_11_ov_24 * xjmxi  + dble_1_ov_4 * xjmxi2  + dble_1_ov_6  * xjmxi3 - dble_1_ov_6  * xjmxi4;
    ci   = dble_115_ov_192 - dble_5_ov_8   * xjmxi2 + dble_1_ov_4 * xjmxi4;
    cip1 = dble_19_ov_96   + dble_11_ov_24 * xjmxi  + dble_1_ov_4 * xjmxi2  - dble_1_ov_6  * xjmxi3 - dble_1_ov_6  * xjmxi4;
    cip2 = dble_1_ov_384   + dble_1_ov_48  * xjmxi  + dble_1_ov_16 * xjmxi2 + dble_1_ov_12 * xjmxi3 + dble_1_ov_24 * xjmxi4;
    
    // Jx
    ( *Jx1D )( im2 )  += cim2 * Jx_ion;
    ( *Jx1D )( im1 )  += cim1 * Jx_ion;
    ( *Jx1D )( i )    += ci   * Jx_ion;
    ( *Jx1D )( ip1 )  += cip1 * Jx_ion;
    ( *Jx1D )( ip2 )  += cip2 * Jx_ion;
    
    
    // Compute Jy_ion & Jz_ion on the primal grid
    // ------------------------------------------
    
    i      = round( xjn );                 // index of the central node
    xjmxi  = xjn - ( double )i;            // normalized distance to the nearest grid point
    xjmxi2 = xjmxi*xjmxi;                  // square of the normalized distance to the nearest grid point
    
    i  -= index_domain_begin;
    im2 = i-2;
    im1 = i-1;
    ip1 = i+1;
    ip2 = i+2;
    
    cim2 = dble_1_ov_384   - dble_1_ov_48  * xjmxi  + dble_1_ov_16 * xjmxi2 - dble_1_ov_12 * xjmxi3 + dble_1_ov_24 * xjmxi4;
    cim1 = dble_19_ov_96   - dble_11_ov_24 * xjmxi  + dble_1_ov_4 * xjmxi2  + dble_1_ov_6  * xjmxi3 - dble_1_ov_6  * xjmxi4;
    ci   = dble_115_ov_192 - dble_5_ov_8   * xjmxi2 + dble_1_ov_4 * xjmxi4;
    cip1 = dble_19_ov_96   + dble_11_ov_24 * xjmxi  + dble_1_ov_4 * xjmxi2  - dble_1_ov_6  * xjmxi3 - dble_1_ov_6  * xjmxi4;
    cip2 = dble_1_ov_384   + dble_1_ov_48  * xjmxi  + dble_1_ov_16 * xjmxi2 + dble_1_ov_12 * xjmxi3 + dble_1_ov_24 * xjmxi4;
    
    // Jy
    ( *Jy1D )( im2 )  += cim2 * Jy_ion;
    ( *Jy1D )( im1 )  += cim1 * Jy_ion;
    ( *Jy1D )( i )    += ci   * Jy_ion;
    ( *Jy1D )( ip1 )  += cip1 * Jy_ion;
    ( *Jy1D )( ip2 )  += cip2 * Jy_ion;
    
    // Jz
    ( *Jz1D )( im2 )  += cim2 * Jz_ion;
    ( *Jz1D )( im1 )  += cim1 * Jz_ion;
    ( *Jz1D )( i )    += ci   * Jz_ion;
    ( *Jz1D )( ip1 )  += cip1 * Jz_ion;
    ( *Jz1D )( ip2 )  += cip2 * Jz_ion;
    
} // END Project global current densities (ionize)


void Projector1D4Order::currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int icell, int ipart_ref )
{
    std::vector<int> *iold = &( smpi->dynamics_iold[ithread] );
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );
    std::vector<double> *invgf = &( smpi->dynamics_invgf[ithread] );
    
    Jx_  =  &( *EMfields->Jx_ )( 0 );
    Jy_  =  &( *EMfields->Jy_ )( 0 );
    Jz_  =  &( *EMfields->Jz_ )( 0 );
    rho_ =  &( *EMfields->rho_ )( 0 );
    
    // If no field diagnostics this timestep, then the projection is done directly on the total arrays
    if( !diag_flag ) {
        if( !is_spectral ) {
            for( int ipart=istart ; ipart<iend; ipart++ ) {
                currents( Jx_, Jy_, Jz_, particles,  ipart, ( *invgf )[ipart], &( *iold )[ipart], &( *delta )[ipart] );
            }
        } else {
            for( int ipart=istart ; ipart<iend; ipart++ ) {
                currentsAndDensity( Jx_, Jy_, Jz_, rho_, particles,  ipart, ( *invgf )[ipart], &( *iold )[ipart], &( *delta )[ipart] );
            }
        }
        // Otherwise, the projection may apply to the species-specific arrays
    } else {
        double *b_Jx  = EMfields->Jx_s [ispec] ? &( *EMfields->Jx_s [ispec] )( 0 ) : &( *EMfields->Jx_ )( 0 ) ;
        double *b_Jy  = EMfields->Jy_s [ispec] ? &( *EMfields->Jy_s [ispec] )( 0 ) : &( *EMfields->Jy_ )( 0 ) ;
        double *b_Jz  = EMfields->Jz_s [ispec] ? &( *EMfields->Jz_s [ispec] )( 0 ) : &( *EMfields->Jz_ )( 0 ) ;
        double *b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
        for( int ipart=istart ; ipart<iend; ipart++ ) {
            currentsAndDensity( b_Jx, b_Jy, b_Jz, b_rho, particles,  ipart, ( *invgf )[ipart], &( *iold )[ipart], &( *delta )[ipart] );
        }
    }
    
}

// Projector for susceptibility used as source term in envelope equation
void Projector1D4Order::susceptibility( ElectroMagn *EMfields, Particles &particles, double species_mass, SmileiMPI *smpi, int istart, int iend,  int ithread, int ibin, int ipart_ref )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}
