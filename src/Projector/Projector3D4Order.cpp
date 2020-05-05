#include "Projector3D4Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field3D.h"
#include "Particles.h"
#include "Tools.h"
#include "Patch.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Projector3D4Order
// ---------------------------------------------------------------------------------------------------------------------
Projector3D4Order::Projector3D4Order( Params &params, Patch *patch ) : Projector3D( params, patch )
{
    dx_inv_   = 1.0/params.cell_length[0];
    dx_ov_dt  = params.cell_length[0] / params.timestep;
    dy_inv_   = 1.0/params.cell_length[1];
    dy_ov_dt  = params.cell_length[1] / params.timestep;
    dz_inv_   = 1.0/params.cell_length[2];
    dz_ov_dt  = params.cell_length[2] / params.timestep;
    
    i_domain_begin = patch->getCellStartingGlobalIndex( 0 );
    j_domain_begin = patch->getCellStartingGlobalIndex( 1 );
    k_domain_begin = patch->getCellStartingGlobalIndex( 2 );
    
    nprimz = params.n_space[2] + 2*params.oversize[2] + 1;
    nprimy = params.n_space[1] + 2*params.oversize[1] + 1;
    
    DEBUG( "cell_length "<< params.cell_length[0] );
    
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Projector3D4Order
// ---------------------------------------------------------------------------------------------------------------------
Projector3D4Order::~Projector3D4Order()
{
}


// ---------------------------------------------------------------------------------------------------------------------
//! Project local currents (sort)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D4Order::currents( double *Jx, double *Jy, double *Jz, Particles &particles, unsigned int ipart, double invgf, int *iold, double *deltaold )
{
    int nparts = particles.size();
    
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    double crx_p = charge_weight*dx_ov_dt;
    double cry_p = charge_weight*dy_ov_dt;
    double crz_p = charge_weight*dz_ov_dt;
    
    // variable declaration
    double xpn, ypn, zpn;
    double delta, delta2, delta3, delta4;
    // arrays used for the Esirkepov projection method
    double Sx0[7], Sx1[7], Sy0[7], Sy1[7], Sz0[7], Sz1[7], DSx[7], DSy[7], DSz[7];
    double tmpJx[7][7], tmpJy[7][7], tmpJz[7][7];
    
    for( unsigned int i=0; i<7; i++ ) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
        Sz1[i] = 0.;
    }
    for( unsigned int j=0; j<7; j++ )
        for( unsigned int k=0; k<7; k++ ) {
            tmpJx[j][k] = 0.;
        }
    for( unsigned int i=0; i<7; i++ )
        for( unsigned int k=0; k<7; k++ ) {
            tmpJy[i][k] = 0.;
        }
    for( unsigned int i=0; i<7; i++ )
        for( unsigned int j=0; j<7; j++ ) {
            tmpJz[i][j] = 0.;
        }
        
    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------
    
    // locate the particle on the primal grid at former time-step & calculate coeff. S0
    delta = deltaold[0*nparts];
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sx0[0] = 0.;
    Sx0[1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sx0[2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sx0[3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sx0[4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sx0[5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sx0[6] = 0.;
    
    delta = deltaold[1*nparts];
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sy0[0] = 0.;
    Sy0[1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sy0[2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy0[3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sy0[4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy0[5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sy0[6] = 0.;
    
    delta = deltaold[2*nparts];
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sz0[0] = 0.;
    Sz0[1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sz0[2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sz0[3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sz0[4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sz0[5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sz0[6] = 0.;
    
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position( 0, ipart ) * dx_inv_;
    int ip = round( xpn );
    int ipo = iold[0*nparts];
    int ip_m_ipo = ip-ipo-i_domain_begin;
    delta  = xpn - ( double )ip;
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sx1[ip_m_ipo+1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sx1[ip_m_ipo+2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sx1[ip_m_ipo+3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sx1[ip_m_ipo+4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sx1[ip_m_ipo+5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    
    ypn = particles.position( 1, ipart ) * dy_inv_;
    int jp = round( ypn );
    int jpo = iold[1*nparts];
    int jp_m_jpo = jp-jpo-j_domain_begin;
    delta  = ypn - ( double )jp;
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sy1[jp_m_jpo+1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sy1[jp_m_jpo+2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy1[jp_m_jpo+3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sy1[jp_m_jpo+4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy1[jp_m_jpo+5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    
    zpn = particles.position( 2, ipart ) * dz_inv_;
    int kp = round( zpn );
    int kpo = iold[2*nparts];
    int kp_m_kpo = kp-kpo-k_domain_begin;
    delta  = zpn - ( double )kp;
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sz1[kp_m_kpo+1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sz1[kp_m_kpo+2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sz1[kp_m_kpo+3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sz1[kp_m_kpo+4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sz1[kp_m_kpo+5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    
    // computes Esirkepov coefficients
    for( unsigned int i=0; i < 7; i++ ) {
        DSx[i] = Sx1[i] - Sx0[i];
        DSy[i] = Sy1[i] - Sy0[i];
        DSz[i] = Sz1[i] - Sz0[i];
    }
    
    // ---------------------------
    // Calculate the total current
    // ---------------------------
    
    ipo -= 3;   //This minus 3 come from the order 4 scheme, based on a 7 points stencil from -3 to +3.
    // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
    jpo -= 3;
    kpo -= 3;
    
    int iloc, jloc, kloc, linindex;
    
    // Jx^(d,p,p)
    for( unsigned int i=1 ; i<7 ; i++ ) {
        iloc = i+ipo;
        for( unsigned int j=0 ; j<7 ; j++ ) {
            jloc = j+jpo;
            for( unsigned int k=0 ; k<7 ; k++ ) {
                tmpJx[j][k] -= crx_p * DSx[i-1] * ( Sy0[j]*Sz0[k] + 0.5*DSy[j]*Sz0[k] + 0.5*DSz[k]*Sy0[j] + one_third*DSy[j]*DSz[k] );
                kloc = k+kpo;
                linindex = iloc*nprimz*nprimy+jloc*nprimz+kloc;
                Jx [linindex] += tmpJx[j][k]; // iloc = (i+ipo)*nprimy;
            }
        }
    }//i
    
    // Jy^(p,d,p)
    for( unsigned int i=0 ; i<7 ; i++ ) {
        iloc = i+ipo;
        for( unsigned int j=1 ; j<7 ; j++ ) {
            jloc = j+jpo;
            for( unsigned int k=0 ; k<7 ; k++ ) {
                tmpJy[i][k] -= cry_p * DSy[j-1] * ( Sz0[k]*Sx0[i] + 0.5*DSz[k]*Sx0[i] + 0.5*DSx[i]*Sz0[k] + one_third*DSz[k]*DSx[i] );
                kloc = k+kpo;
                linindex = iloc*nprimz*( nprimy+1 )+jloc*nprimz+kloc;
                Jy [linindex] += tmpJy[i][k]; //
            }
        }
    }//i
    
    // Jz^(p,p,d)
    for( unsigned int i=0 ; i<7 ; i++ ) {
        iloc = i+ipo;
        for( unsigned int j=0 ; j<7 ; j++ ) {
            jloc = j+jpo;
            for( unsigned int k=1 ; k<7 ; k++ ) {
                tmpJz[i][j] -= crz_p * DSz[k-1] * ( Sx0[i]*Sy0[j] + 0.5*DSx[i]*Sy0[j] + 0.5*DSy[j]*Sx0[i] + one_third*DSx[i]*DSy[j] );
                kloc = k+kpo;
                linindex = iloc*( nprimz+1 )*nprimy+jloc*( nprimz+1 )+kloc;
                Jz [linindex] += tmpJz[i][j]; //
            }
        }
    }//i
    
    
} // END Project local current densities (Jx, Jy, Jz, sort)


// ---------------------------------------------------------------------------------------------------------------------
//! Project local current densities (sort)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D4Order::currentsAndDensity( double *Jx, double *Jy, double *Jz, double *rho, Particles &particles, unsigned int ipart, double invgf, int *iold, double *deltaold )
{
    int nparts = particles.size();
    
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    double crx_p = charge_weight*dx_ov_dt;
    double cry_p = charge_weight*dy_ov_dt;
    double crz_p = charge_weight*dz_ov_dt;
    
    // variable declaration
    double xpn, ypn, zpn;
    double delta, delta2, delta3, delta4;
    // arrays used for the Esirkepov projection method
    double Sx0[7], Sx1[7], Sy0[7], Sy1[7], Sz0[7], Sz1[7], DSx[7], DSy[7], DSz[7];
    double tmpJx[7][7], tmpJy[7][7], tmpJz[7][7];
    
    for( unsigned int i=0; i<7; i++ ) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
        Sz1[i] = 0.;
    }
    
    for( unsigned int j=0; j<7; j++ )
        for( unsigned int k=0; k<7; k++ ) {
            tmpJx[j][k] = 0.;
        }
    for( unsigned int i=0; i<7; i++ )
        for( unsigned int k=0; k<7; k++ ) {
            tmpJy[i][k] = 0.;
        }
    for( unsigned int i=0; i<7; i++ )
        for( unsigned int j=0; j<7; j++ ) {
            tmpJz[i][j] = 0.;
        }
    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------
    
    // locate the particle on the primal grid at former time-step & calculate coeff. S0
    delta = deltaold[0*nparts];
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sx0[0] = 0.;
    Sx0[1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sx0[2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sx0[3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sx0[4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sx0[5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sx0[6] = 0.;
    
    delta = deltaold[1*nparts];
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sy0[0] = 0.;
    Sy0[1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sy0[2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy0[3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sy0[4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy0[5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sy0[6] = 0.;
    
    delta = deltaold[2*nparts];
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sz0[0] = 0.;
    Sz0[1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sz0[2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sz0[3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sz0[4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sz0[5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sz0[6] = 0.;
    
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position( 0, ipart ) * dx_inv_;
    int ip = round( xpn );
    int ipo = iold[0*nparts];
    int ip_m_ipo = ip-ipo-i_domain_begin;
    delta  = xpn - ( double )ip;
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sx1[ip_m_ipo+1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sx1[ip_m_ipo+2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sx1[ip_m_ipo+3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sx1[ip_m_ipo+4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sx1[ip_m_ipo+5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    
    ypn = particles.position( 1, ipart ) * dy_inv_;
    int jp = round( ypn );
    int jpo = iold[1*nparts];
    int jp_m_jpo = jp-jpo-j_domain_begin;
    delta  = ypn - ( double )jp;
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sy1[jp_m_jpo+1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sy1[jp_m_jpo+2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy1[jp_m_jpo+3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sy1[jp_m_jpo+4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy1[jp_m_jpo+5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    
    zpn = particles.position( 2, ipart ) * dz_inv_;
    int kp = round( zpn );
    int kpo = iold[2*nparts];
    int kp_m_kpo = kp-kpo-k_domain_begin;
    delta  = zpn - ( double )kp;
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sz1[kp_m_kpo+1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sz1[kp_m_kpo+2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sz1[kp_m_kpo+3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sz1[kp_m_kpo+4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sz1[kp_m_kpo+5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    
    // computes Esirkepov coefficients
    for( unsigned int i=0; i < 7; i++ ) {
        DSx[i] = Sx1[i] - Sx0[i];
        DSy[i] = Sy1[i] - Sy0[i];
        DSz[i] = Sz1[i] - Sz0[i];
    }
    
    // ---------------------------
    // Calculate the total current
    // ---------------------------
    
    ipo -= 3;   //This minus 3 come from the order 3 scheme, based on a 7 points stencil from -3 to +3.
    // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
    jpo -= 3;
    kpo -= 3;
    
    int iloc, jloc, kloc, linindex;
    
    // Jx^(d,p,p)
    for( unsigned int i=1 ; i<7 ; i++ ) {
        iloc = i+ipo;
        for( unsigned int j=0 ; j<7 ; j++ ) {
            jloc = j+jpo;
            for( unsigned int k=0 ; k<7 ; k++ ) {
                tmpJx[j][k] -= crx_p * DSx[i-1] * ( Sy0[j]*Sz0[k] + 0.5*DSy[j]*Sz0[k] + 0.5*DSz[k]*Sy0[j] + one_third*DSy[j]*DSz[k] );
                kloc = k+kpo;
                linindex = iloc*nprimz*nprimy+jloc*nprimz+kloc;
                Jx [linindex] += tmpJx[j][k]; // iloc = (i+ipo)*nprimy;
            }
        }
    }//i
    
    // Jy^(p,d,p)
    for( unsigned int i=0 ; i<7 ; i++ ) {
        iloc = i+ipo;
        for( unsigned int j=1 ; j<7 ; j++ ) {
            jloc = j+jpo;
            for( unsigned int k=0 ; k<7 ; k++ ) {
                tmpJy[i][k] -= cry_p * DSy[j-1] * ( Sz0[k]*Sx0[i] + 0.5*DSz[k]*Sx0[i] + 0.5*DSx[i]*Sz0[k] + one_third*DSz[k]*DSx[i] );
                kloc = k+kpo;
                linindex = iloc*nprimz*( nprimy+1 )+jloc*nprimz+kloc;
                Jy [linindex] += tmpJy[i][k]; //
            }
        }
    }//i
    
    // Jz^(p,p,d)
    for( unsigned int i=0 ; i<7 ; i++ ) {
        iloc = i+ipo;
        for( unsigned int j=0 ; j<7 ; j++ ) {
            jloc = j+jpo;
            for( unsigned int k=1 ; k<7 ; k++ ) {
                tmpJz[i][j] -= crz_p * DSz[k-1] * ( Sx0[i]*Sy0[j] + 0.5*DSx[i]*Sy0[j] + 0.5*DSy[j]*Sx0[i] + one_third*DSx[i]*DSy[j] );
                kloc = k+kpo;
                linindex = iloc*( nprimz+1 )*nprimy+jloc*( nprimz+1 )+kloc;
                Jz [linindex] += tmpJz[i][j]; //
            }
        }
    }//i
    
    // Rho^(p,p,p)
    for( unsigned int i=0 ; i<7 ; i++ ) {
        iloc = i+ipo;
        for( unsigned int j=0 ; j<7 ; j++ ) {
            jloc = j+jpo;
            for( unsigned int k=0 ; k<7 ; k++ ) {
                kloc = k+kpo;
                linindex = iloc*nprimz*nprimy+jloc*nprimz+kloc;
                rho[linindex] += charge_weight * Sx1[i]*Sy1[j]*Sz1[k];
            }
        }
    }//i
    
} // END Project local densities (Jx, Jy, Jz, rho, sort)


// ---------------------------------------------------------------------------------------------------------------------
//! Project local densities only (Frozen species)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D4Order::basic( double *rhoj, Particles &particles, unsigned int ipart, unsigned int type )
{
    //Warning : this function is used for frozen species or initialization only and doesn't use the standard scheme.
    //rho type = 0
    //Jx type = 1
    //Jy type = 2
    //Jz type = 3
    
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    
    int iloc, jloc;
    int ny( nprimy ), nz( nprimz ), nyz;
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    
    if( type > 0 ) {
        charge_weight *= 1./sqrt( 1.0 + particles.momentum( 0, ipart )*particles.momentum( 0, ipart )
                                  + particles.momentum( 1, ipart )*particles.momentum( 1, ipart )
                                  + particles.momentum( 2, ipart )*particles.momentum( 2, ipart ) );
                                  
        if( type == 1 ) {
            charge_weight *= particles.momentum( 0, ipart );
        } else if( type == 2 ) {
            charge_weight *= particles.momentum( 1, ipart );
            ny ++;
        } else {
            charge_weight *= particles.momentum( 2, ipart );
            nz ++;
        }
    }
    nyz = ny*nz;
    
    // variable declaration
    double xpn, ypn, zpn;
    double delta, delta2, delta3, delta4;
    double Sx1[7], Sy1[7], Sz1[7]; // arrays used for the Esirkepov projection method
    
// Initialize all current-related arrays to zero
    for( unsigned int i=0; i<7; i++ ) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
        Sz1[i] = 0.;
    }
    
    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------
    
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position( 0, ipart ) * dx_inv_;
    int ip = round( xpn+ 0.5*( type==1 ) );
    delta  = xpn - ( double )ip;
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sx1[1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sx1[2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sx1[3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sx1[4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sx1[5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    
    ypn = particles.position( 1, ipart ) * dy_inv_;
    int jp = round( ypn+ 0.5*( type==2 ) );
    delta  = ypn - ( double )jp;
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sy1[1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sy1[2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy1[3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sy1[4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy1[5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    
    zpn = particles.position( 2, ipart ) * dz_inv_;
    int kp = round( zpn+ 0.5*( type==3 ) );
    delta  = zpn - ( double )kp;
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sz1[1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sz1[2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sz1[3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sz1[4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sz1[5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    
    // ---------------------------
    // Calculate the total charge
    // ---------------------------
    ip -= i_domain_begin + 3;
    jp -= j_domain_begin + 3;
    kp -= k_domain_begin + 3;
    
    for( unsigned int i=0 ; i<7 ; i++ ) {
        iloc = ( i+ip )*nyz;
        for( unsigned int j=0 ; j<7 ; j++ ) {
            jloc = ( jp+j )*nz;
            for( unsigned int k=0 ; k<7 ; k++ ) {
                rhoj[iloc+jloc+kp+k] += charge_weight * Sx1[i]*Sy1[j]*Sz1[k];
            }
        }
    }//i
    
} // END Project local current densities (Frozen species)

// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities (ionize)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D4Order::ionizationCurrents( Field *Jx, Field *Jy, Field *Jz, Particles &particles, int ipart, LocalFields Jion )
{
    Field3D *Jx3D  = static_cast<Field3D *>( Jx );
    Field3D *Jy3D  = static_cast<Field3D *>( Jy );
    Field3D *Jz3D  = static_cast<Field3D *>( Jz );
    
    
    //Declaration of local variables
    int ip, id, jp, jd, kp, kd;
    double xpn, xpmxip, xpmxip2, xpmxip3, xpmxip4, xpmxid, xpmxid2, xpmxid3, xpmxid4;
    double ypn, ypmyjp, ypmyjp2, ypmyjp3, ypmyjp4, ypmyjd, ypmyjd2, ypmyjd3, ypmyjd4;
    double zpn, zpmzkp, zpmzkp2, zpmzkp3, zpmzkp4, zpmzkd, zpmzkd2, zpmzkd3, zpmzkd4;
    double Sxp[5], Sxd[5], Syp[5], Syd[5], Szp[5], Szd[5];
    
    // weighted currents
    double weight = inv_cell_volume * particles.weight( ipart );
    double Jx_ion = Jion.x * weight;
    double Jy_ion = Jion.y * weight;
    double Jz_ion = Jion.z * weight;
    
    //Locate particle on the grid
    xpn    = particles.position( 0, ipart ) * dx_inv_; // normalized distance to the first node
    ypn    = particles.position( 1, ipart ) * dy_inv_; // normalized distance to the first node
    zpn    = particles.position( 2, ipart ) * dz_inv_; // normalized distance to the first node
    
    // x-primal index
    ip      = round( xpn );                  // x-index of the central node
    xpmxip  = xpn - ( double )ip;            // normalized distance to the nearest grid point
    xpmxip2 = xpmxip*xpmxip;                 // square of the normalized distance to the nearest grid point
    xpmxip3 = xpmxip2*xpmxip;                // cube
    xpmxip4 = xpmxip2*xpmxip2;               // fourth-power
    
    // x-dual index
    id      = round( xpn+0.5 );              // x-index of the central node
    xpmxid  = xpn - ( double )id + 0.5;      // normalized distance to the nearest grid point
    xpmxid2 = xpmxid*xpmxid;                 // square of the normalized distance to the nearest grid point
    xpmxid3 = xpmxid2*xpmxid;                // cube
    xpmxid4 = xpmxid2*xpmxid2;               // fourth-power
    
    // y-primal index
    jp      = round( ypn );                  // y-index of the central node
    ypmyjp  = ypn - ( double )jp;            // normalized distance to the nearest grid point
    ypmyjp2 = ypmyjp*ypmyjp;                 // square of the normalized distance to the nearest grid point
    ypmyjp3 = ypmyjp2*ypmyjp;                // cube
    ypmyjp4 = ypmyjp2*ypmyjp2;               // fourth-power
    
    // y-dual index
    jd      = round( ypn+0.5 );              // y-index of the central node
    ypmyjd  = ypn - ( double )jd + 0.5;      // normalized distance to the nearest grid point
    ypmyjd2 = ypmyjd*ypmyjd;                 // square of the normalized distance to the nearest grid point
    ypmyjd3 = ypmyjd2*ypmyjd;                // cube
    ypmyjd4 = ypmyjd2*ypmyjd2;               // fourth-power
    
    // z-primal index
    kp      = round( zpn );                  // z-index of the central node
    zpmzkp  = zpn - ( double )kp;            // normalized distance to the nearest grid point
    zpmzkp2 = zpmzkp*zpmzkp;                 // square of the normalized distance to the nearest grid point
    zpmzkp3 = zpmzkp2*zpmzkp;                // cube
    zpmzkp4 = zpmzkp2*zpmzkp2;               // fourth-power
    
    // z-dual index
    kd      = round( zpn+0.5 );              // z-index of the central node
    zpmzkd  = zpn - ( double )kd + 0.5;      // normalized distance to the nearest grid point
    zpmzkd2 = zpmzkd*zpmzkd;                 // square of the normalized distance to the nearest grid point
    zpmzkd3 = zpmzkd2*zpmzkd;                // square
    zpmzkd4 = zpmzkd2*zpmzkd2;               // square
    
    Sxp[0] = dble_1_ov_384   - dble_1_ov_48  * xpmxip  + dble_1_ov_16 * xpmxip2 - dble_1_ov_12 * xpmxip3 + dble_1_ov_24 * xpmxip4;
    Sxp[1] = dble_19_ov_96   - dble_11_ov_24 * xpmxip  + dble_1_ov_4  * xpmxip2 + dble_1_ov_6  * xpmxip3 - dble_1_ov_6  * xpmxip4;
    Sxp[2] = dble_115_ov_192 - dble_5_ov_8   * xpmxip2 + dble_1_ov_4  * xpmxip4;
    Sxp[3] = dble_19_ov_96   + dble_11_ov_24 * xpmxip  + dble_1_ov_4  * xpmxip2 - dble_1_ov_6  * xpmxip3 - dble_1_ov_6  * xpmxip4;
    Sxp[4] = dble_1_ov_384   + dble_1_ov_48  * xpmxip  + dble_1_ov_16 * xpmxip2 + dble_1_ov_12 * xpmxip3 + dble_1_ov_24 * xpmxip4;
    
    Sxd[0] = dble_1_ov_384   - dble_1_ov_48  * xpmxid  + dble_1_ov_16 * xpmxid2 - dble_1_ov_12 * xpmxid3 + dble_1_ov_24 * xpmxid4;
    Sxd[1] = dble_19_ov_96   - dble_11_ov_24 * xpmxid  + dble_1_ov_4  * xpmxid2 + dble_1_ov_6  * xpmxid3 - dble_1_ov_6  * xpmxid4;
    Sxd[2] = dble_115_ov_192 - dble_5_ov_8   * xpmxid2 + dble_1_ov_4  * xpmxid4;
    Sxd[3] = dble_19_ov_96   + dble_11_ov_24 * xpmxid  + dble_1_ov_4  * xpmxid2 - dble_1_ov_6  * xpmxid3 - dble_1_ov_6  * xpmxid4;
    Sxd[4] = dble_1_ov_384   + dble_1_ov_48  * xpmxid  + dble_1_ov_16 * xpmxid2 + dble_1_ov_12 * xpmxid3 + dble_1_ov_24 * xpmxid4;
    
    Syp[0] = dble_1_ov_384   - dble_1_ov_48  * ypmyjp  + dble_1_ov_16 * ypmyjp2 - dble_1_ov_12 * ypmyjp3 + dble_1_ov_24 * ypmyjp4;
    Syp[1] = dble_19_ov_96   - dble_11_ov_24 * ypmyjp  + dble_1_ov_4  * ypmyjp2 + dble_1_ov_6  * ypmyjp3 - dble_1_ov_6  * ypmyjp4;
    Syp[2] = dble_115_ov_192 - dble_5_ov_8   * ypmyjp2 + dble_1_ov_4  * ypmyjp4;
    Syp[3] = dble_19_ov_96   + dble_11_ov_24 * ypmyjp  + dble_1_ov_4  * ypmyjp2 - dble_1_ov_6  * ypmyjp3 - dble_1_ov_6  * ypmyjp4;
    Syp[4] = dble_1_ov_384   + dble_1_ov_48  * ypmyjp  + dble_1_ov_16 * ypmyjp2 + dble_1_ov_12 * ypmyjp3 + dble_1_ov_24 * ypmyjp4;
    
    Syd[0] = dble_1_ov_384   - dble_1_ov_48  * ypmyjd  + dble_1_ov_16 * ypmyjd2 - dble_1_ov_12 * ypmyjd3 + dble_1_ov_24 * ypmyjd4;
    Syd[1] = dble_19_ov_96   - dble_11_ov_24 * ypmyjd  + dble_1_ov_4  * ypmyjd2 + dble_1_ov_6  * ypmyjd3 - dble_1_ov_6  * ypmyjd4;
    Syd[2] = dble_115_ov_192 - dble_5_ov_8   * ypmyjd2 + dble_1_ov_4  * ypmyjd4;
    Syd[3] = dble_19_ov_96   + dble_11_ov_24 * ypmyjd  + dble_1_ov_4  * ypmyjd2 - dble_1_ov_6  * ypmyjd3 - dble_1_ov_6  * ypmyjd4;
    Syd[4] = dble_1_ov_384   + dble_1_ov_48  * ypmyjd  + dble_1_ov_16 * ypmyjd2 + dble_1_ov_12 * ypmyjd3 + dble_1_ov_24 * ypmyjd4;
    
    Szp[0] = dble_1_ov_384   - dble_1_ov_48  * zpmzkp  + dble_1_ov_16 * zpmzkp2 - dble_1_ov_12 * zpmzkp3 + dble_1_ov_24 * zpmzkp4;
    Szp[1] = dble_19_ov_96   - dble_11_ov_24 * zpmzkp  + dble_1_ov_4  * zpmzkp2 + dble_1_ov_6  * zpmzkp3 - dble_1_ov_6  * zpmzkp4;
    Szp[2] = dble_115_ov_192 - dble_5_ov_8   * zpmzkp2 + dble_1_ov_4  * zpmzkp4;
    Szp[3] = dble_19_ov_96   + dble_11_ov_24 * zpmzkp  + dble_1_ov_4  * zpmzkp2 - dble_1_ov_6  * zpmzkp3 - dble_1_ov_6  * zpmzkp4;
    Szp[4] = dble_1_ov_384   + dble_1_ov_48  * zpmzkp  + dble_1_ov_16 * zpmzkp2 + dble_1_ov_12 * zpmzkp3 + dble_1_ov_24 * zpmzkp4;
    
    Szd[0] = dble_1_ov_384   - dble_1_ov_48  * zpmzkd  + dble_1_ov_16 * zpmzkd2 - dble_1_ov_12 * zpmzkd3 + dble_1_ov_24 * zpmzkd4;
    Szd[1] = dble_19_ov_96   - dble_11_ov_24 * zpmzkd  + dble_1_ov_4  * zpmzkd2 + dble_1_ov_6  * zpmzkd3 - dble_1_ov_6  * zpmzkd4;
    Szd[2] = dble_115_ov_192 - dble_5_ov_8   * zpmzkd2 + dble_1_ov_4  * zpmzkd4;
    Szd[3] = dble_19_ov_96   + dble_11_ov_24 * zpmzkd  + dble_1_ov_4  * zpmzkd2 - dble_1_ov_6  * zpmzkd3 - dble_1_ov_6  * zpmzkd4;
    Szd[4] = dble_1_ov_384   + dble_1_ov_48  * zpmzkd  + dble_1_ov_16 * zpmzkd2 + dble_1_ov_12 * zpmzkd3 + dble_1_ov_24 * zpmzkd4;
    
    ip  -= i_domain_begin;
    id  -= i_domain_begin;
    jp  -= j_domain_begin;
    jd  -= j_domain_begin;
    kp  -= k_domain_begin;
    kd  -= k_domain_begin;
    
    for( unsigned int i=0 ; i<5 ; i++ ) {
        int iploc=ip+i-2;
        int idloc=id+i-2;
        for( unsigned int j=0 ; j<5 ; j++ ) {
            int jploc=jp+j-2;
            int jdloc=jd+j-2;
            for( unsigned int k=0 ; k<5 ; k++ ) {
                int kploc=kp+k-2;
                int kdloc=kd+k-2;
                // Jx^(d,p,p)
                ( *Jx3D )( idloc, jploc, kploc ) += Jx_ion * Sxd[i]*Syp[j]*Szp[k];
                // Jy^(p,d,p)
                ( *Jy3D )( iploc, jdloc, kploc ) += Jy_ion * Sxp[i]*Syd[j]*Szp[k];
                // Jz^(p,p,d)
                ( *Jz3D )( iploc, jploc, kdloc ) += Jz_ion * Sxp[i]*Syp[j]*Szd[k];
            }//k
        }//j
    }//i
    
    
} // END Project global current densities (ionize)

//Wrapper for projection
void Projector3D4Order::currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int icell, int ipart_ref )
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
void Projector3D4Order::susceptibility( ElectroMagn *EMfields, Particles &particles, double species_mass, SmileiMPI *smpi, int istart, int iend,  int ithread, int icell, int ipart_ref )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}

