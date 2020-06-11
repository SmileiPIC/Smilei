#include "Projector2D4Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field2D.h"
#include "Particles.h"
#include "Tools.h"
#include "Patch.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Projector2D4Order
// ---------------------------------------------------------------------------------------------------------------------
Projector2D4Order::Projector2D4Order( Params &params, Patch *patch ) : Projector2D( params, patch )
{
    dx_inv_   = 1.0/params.cell_length[0];
    dx_ov_dt  = params.cell_length[0] / params.timestep;
    dy_inv_   = 1.0/params.cell_length[1];
    dy_ov_dt  = params.cell_length[1] / params.timestep;
    
    i_domain_begin = patch->getCellStartingGlobalIndex( 0 );
    j_domain_begin = patch->getCellStartingGlobalIndex( 1 );
    
    nprimy = params.n_space[1] + 2*params.oversize[1] + 1;
    
    DEBUG( "cell_length "<< params.cell_length[0] );
    
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Projector2D4Order
// ---------------------------------------------------------------------------------------------------------------------
Projector2D4Order::~Projector2D4Order()
{
}


// ---------------------------------------------------------------------------------------------------------------------
//! Project current densities : main projector
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D4Order::currents( double *Jx, double *Jy, double *Jz, Particles &particles, unsigned int ipart, double invgf, int *iold, double *deltaold )
{
    int nparts = particles.size();
    
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    
    int iloc;
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    double crx_p = charge_weight*dx_ov_dt;
    double cry_p = charge_weight*dy_ov_dt;
    double crz_p = charge_weight*one_third*particles.momentum( 2, ipart )*invgf;
    
    // variable declaration
    double xpn, ypn;
    double delta, delta2, delta3, delta4;
    // arrays used for the Esirkepov projection method
    double  Sx0[7], Sx1[7], Sy0[7], Sy1[7], DSx[7], DSy[7], tmpJx[7];
    
    for( unsigned int i=0; i<7; i++ ) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
        tmpJx[i] = 0.;
    }
    Sx0[0] = 0.;
    Sx0[6] = 0.;
    Sy0[0] = 0.;
    Sy0[6] = 0.;
    
    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------
    
    // locate the particle on the primal grid at former time-step & calculate coeff. S0
    delta = deltaold[0*nparts];
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    
    Sx0[1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sx0[2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sx0[3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sx0[4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sx0[5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    
    delta = deltaold[1*nparts];
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    
    Sy0[1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sy0[2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy0[3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sy0[4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy0[5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    
    
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
    
    for( unsigned int i=0; i < 7; i++ ) {
        DSx[i] = Sx1[i] - Sx0[i];
        DSy[i] = Sy1[i] - Sy0[i];
    }
    
    // calculate Esirkepov coeff. Wx, Wy, Wz when used
    double tmp, tmp2, tmp3, tmpY;
    //Do not compute useless weights.
    // ------------------------------------------------
    // Local current created by the particle
    // calculate using the charge conservation equation
    // ------------------------------------------------
    
    // ---------------------------
    // Calculate the total current
    // ---------------------------
    ipo -= 3; //This minus 3 come from the order 4 scheme, based on a 7 points stencil from -3 to +3.
    jpo -= 3;
    // i =0
    {
        iloc = ipo*nprimy+jpo;
        tmp2 = 0.5*Sx1[0];
        tmp3 =     Sx1[0];
        Jz[iloc]  += crz_p * ( Sy1[0]*tmp3 );
        tmp = 0;
        tmpY = Sx0[0] + 0.5*DSx[0];
        for( unsigned int j=1 ; j<7 ; j++ ) {
            tmp -= cry_p * DSy[j-1] * tmpY;
            Jy[iloc+j+ipo]  += tmp; //Because size of Jy in Y is nprimy+1.
            Jz[iloc+j]  += crz_p * ( Sy0[j]*tmp2 + Sy1[j]*tmp3 );
        }
    }//i
    
    for( unsigned int i=1 ; i<7 ; i++ ) {
        iloc = ( i+ipo )*nprimy+jpo;
        tmpJx[0] -= crx_p *  DSx[i-1] * ( 0.5*DSy[0] );
        Jx[iloc]  += tmpJx[0];
        tmp2 = 0.5*Sx1[i] + Sx0[i];
        tmp3 = 0.5*Sx0[i] + Sx1[i];
        Jz[iloc]  += crz_p * ( Sy1[0]*tmp3 );
        tmp = 0;
        tmpY = Sx0[i] + 0.5*DSx[i];
        for( unsigned int j=1 ; j<7 ; j++ ) {
            tmpJx[j] -= crx_p * DSx[i-1] * ( Sy0[j] + 0.5*DSy[j] );
            Jx[iloc+j]  += tmpJx[j];
            tmp -= cry_p * DSy[j-1] * tmpY;
            Jy[iloc+j+i+ipo]  += tmp; //Because size of Jy in Y is nprimy+1.
            Jz[iloc+j]  += crz_p * ( Sy0[j]*tmp2 + Sy1[j]*tmp3 );
        }
    }//i
    
}


// ---------------------------------------------------------------------------------------------------------------------
//! Project current densities & charge : diagFields timstep
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D4Order::currentsAndDensity( double *Jx, double *Jy, double *Jz, double *rho, Particles &particles, unsigned int ipart, double invgf, int *iold, double *deltaold )
{
    int nparts = particles.size();
    
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    
    int iloc;
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    double crx_p = charge_weight*dx_ov_dt;
    double cry_p = charge_weight*dy_ov_dt;
    double crz_p = charge_weight*one_third*particles.momentum( 2, ipart )*invgf;
    
    // variable declaration
    double xpn, ypn;
    double delta, delta2, delta3, delta4;
    // arrays used for the Esirkepov projection method
    double  Sx0[7], Sx1[7], Sy0[7], Sy1[7], DSx[7], DSy[7], tmpJx[7];
    
    for( unsigned int i=0; i<7; i++ ) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
        tmpJx[i] = 0.;
    }
    Sx0[0] = 0.;
    Sx0[6] = 0.;
    Sy0[0] = 0.;
    Sy0[6] = 0.;
    
    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------
    
    // locate the particle on the primal grid at former time-step & calculate coeff. S0
    delta = deltaold[0*nparts];
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    
    Sx0[1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sx0[2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sx0[3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sx0[4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sx0[5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    
    delta = deltaold[1*nparts];
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    
    Sy0[1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sy0[2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy0[3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sy0[4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy0[5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    
    
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
    
    for( unsigned int i=0; i < 7; i++ ) {
        DSx[i] = Sx1[i] - Sx0[i];
        DSy[i] = Sy1[i] - Sy0[i];
    }
    
    // calculate Esirkepov coeff. Wx, Wy, Wz when used
    double tmp, tmp2, tmp3, tmpY;
    //Do not compute useless weights.
    // ------------------------------------------------
    // Local current created by the particle
    // calculate using the charge conservation equation
    // ------------------------------------------------
    
    // ---------------------------
    // Calculate the total current
    // ---------------------------
    ipo -= 3; //This minus 3 come from the order 4 scheme, based on a 7 points stencil from -3 to +3.
    jpo -= 3;
    // i =0
    {
        iloc = ipo*nprimy+jpo;
        tmp2 = 0.5*Sx1[0];
        tmp3 =     Sx1[0];
        Jz[iloc]  += crz_p * ( Sy1[0]*tmp3 );
        rho[iloc] += charge_weight * Sx1[0]*Sy1[0];
        tmp = 0;
        tmpY = Sx0[0] + 0.5*DSx[0];
        for( unsigned int j=1 ; j<7 ; j++ ) {
            tmp -= cry_p * DSy[j-1] * tmpY;
            Jy[iloc+j+ipo]  += tmp; //Because size of Jy in Y is nprimy+1.
            Jz[iloc+j]  += crz_p * ( Sy0[j]*tmp2 + Sy1[j]*tmp3 );
            rho[iloc+j] += charge_weight * Sx1[0]*Sy1[j];
        }
    }//i
    
    for( unsigned int i=1 ; i<7 ; i++ ) {
        iloc = ( i+ipo )*nprimy+jpo;
        tmpJx[0] -= crx_p *  DSx[i-1] * ( 0.5*DSy[0] );
        Jx[iloc]  += tmpJx[0];
        tmp2 = 0.5*Sx1[i] + Sx0[i];
        tmp3 = 0.5*Sx0[i] + Sx1[i];
        Jz[iloc]  += crz_p * ( Sy1[0]*tmp3 );
        rho[iloc] += charge_weight * Sx1[i]*Sy1[0];
        tmp = 0;
        tmpY = Sx0[i] + 0.5*DSx[i];
        for( unsigned int j=1 ; j<7 ; j++ ) {
            tmpJx[j] -= crx_p * DSx[i-1] * ( Sy0[j] + 0.5*DSy[j] );
            Jx[iloc+j]  += tmpJx[j];
            tmp -= cry_p * DSy[j-1] * tmpY;
            Jy[iloc+j+i+ipo]  += tmp; //Because size of Jy in Y is nprimy+1.
            Jz[iloc+j]  += crz_p * ( Sy0[j]*tmp2 + Sy1[j]*tmp3 );
            rho[iloc+j] += charge_weight * Sx1[i]*Sy1[j];
        }
    }//i
    
    
}


// ---------------------------------------------------------------------------------------------------------------------
//! Project charge : frozen & diagFields timstep
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D4Order::basic( double *rhoj, Particles &particles, unsigned int ipart, unsigned int type )
{
    //Warning : this function is used for frozen species or initialization only and doesn't use the standard scheme.
    //rho type = 0
    //Jx type = 1
    //Jy type = 2
    //Jz type = 3
    
    int iloc;
    int ny( nprimy );
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
        }
    }
    
    // variable declaration
    double xpn, ypn;
    double delta, delta2, delta3, delta4;
    // arrays used for the Esirkepov projection method
    double  Sx1[7], Sy1[7];
    
    for( unsigned int i=0; i<7; i++ ) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
    }
    
    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position( 0, ipart ) * dx_inv_;
    int ip        = round( xpn + 0.5 * ( type==1 ) );                       // index of the central node
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
    int jp = round( ypn + 0.5*( type==2 ) );
    delta  = ypn - ( double )jp;
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    
    Sy1[1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sy1[2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy1[3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sy1[4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy1[5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    
    // ---------------------------
    // Calculate the total current
    // ---------------------------
    ip -= i_domain_begin + 3;
    jp -= j_domain_begin + 3;
    
    for( unsigned int i=0 ; i<7 ; i++ ) {
        iloc = ( i+ip )*ny+jp;
        for( unsigned int j=0 ; j<7 ; j++ ) {
            rhoj[iloc+j] += charge_weight * Sx1[i]*Sy1[j];
        }
    }//i
}


// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities : ionization
// ---------------------------------------------------------------------------------------------------------------------
void  Projector2D4Order::ionizationCurrents( Field *Jx, Field *Jy, Field *Jz, Particles &particles, int ipart, LocalFields Jion )
{
    Field2D *Jx2D  = static_cast<Field2D *>( Jx );
    Field2D *Jy2D  = static_cast<Field2D *>( Jy );
    Field2D *Jz2D  = static_cast<Field2D *>( Jz );
    
    
    //Declaration of local variables
    int ip, id, jp, jd;
    double xpn, xpmxip, xpmxip2, xpmxip3, xpmxip4, xpmxid, xpmxid2, xpmxid3, xpmxid4;
    double ypn, ypmyjp, ypmyjp2, ypmyjp3, ypmyjp4, ypmyjd, ypmyjd2, ypmyjd3, ypmyjd4;
    double Sxp[5], Sxd[5], Syp[5], Syd[5];
    
    // weighted currents
    double weight = inv_cell_volume * particles.weight( ipart );
    double Jx_ion = Jion.x * weight;
    double Jy_ion = Jion.y * weight;
    double Jz_ion = Jion.z * weight;
    
    //Locate particle on the grid
    xpn    = particles.position( 0, ipart ) * dx_inv_; // normalized distance to the first node
    ypn    = particles.position( 1, ipart ) * dy_inv_; // normalized distance to the first node
    
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

    ip  -= i_domain_begin;
    id  -= i_domain_begin;
    jp  -= j_domain_begin;
    jd  -= j_domain_begin;

    for (unsigned int i=0 ; i<5 ; i++) {
        int iploc=ip+i-2;
        int idloc=id+i-2;
        for (unsigned int j=0 ; j<5 ; j++) {
            int jploc=jp+j-2;
            int jdloc=jd+j-2;
            // Jx^(d,p)
            (*Jx2D)(idloc,jploc) += Jx_ion * Sxd[i]*Syp[j];
            // Jy^(p,d)
            (*Jy2D)(iploc,jdloc) += Jy_ion * Sxp[i]*Syd[j];
            // Jz^(p,p)
            (*Jz2D)(iploc,jploc) += Jz_ion * Sxp[i]*Syp[j];
        }
    }
}


// ---------------------------------------------------------------------------------------------------------------------
//! Wrapper for projection
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D4Order::currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int icell, int ipart_ref )
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
void Projector2D4Order::susceptibility( ElectroMagn *EMfields, Particles &particles, double species_mass, SmileiMPI *smpi, int istart, int iend,  int ithread, int icell, int ipart_ref )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}
