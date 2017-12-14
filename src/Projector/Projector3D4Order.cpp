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
Projector3D4Order::Projector3D4Order (Params& params, Patch* patch) : Projector3D(params, patch)
{
    dx_inv_   = 1.0/params.cell_length[0];
    dx_ov_dt  = params.cell_length[0] / params.timestep;
    dy_inv_   = 1.0/params.cell_length[1];
    dy_ov_dt  = params.cell_length[1] / params.timestep;
    dz_inv_   = 1.0/params.cell_length[2];
    dz_ov_dt  = params.cell_length[2] / params.timestep;
    
    one_third = 1.0/3.0;

    //double defined for use in coefficients
    dble_1_ov_384 = 1.0/384.0;
    dble_1_ov_48 = 1.0/48.0;
    dble_1_ov_16 = 1.0/16.0;
    dble_1_ov_12 = 1.0/12.0;
    dble_1_ov_24 = 1.0/24.0;
    dble_19_ov_96 = 19.0/96.0;
    dble_11_ov_24 = 11.0/24.0;
    dble_1_ov_4 = 1.0/4.0;
    dble_1_ov_6 = 1.0/6.0;
    dble_115_ov_192 = 115.0/192.0;
    dble_5_ov_8 = 5.0/8.0;

    i_domain_begin = patch->getCellStartingGlobalIndex(0);
    j_domain_begin = patch->getCellStartingGlobalIndex(1);
    k_domain_begin = patch->getCellStartingGlobalIndex(2);

    DEBUG("cell_length "<< params.cell_length[0]);

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
void Projector3D4Order::operator() (double* Jx, double* Jy, double* Jz, Particles &particles, unsigned int ipart, double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaold)
{
    int nparts = particles.size();

    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
    double crx_p = charge_weight*dx_ov_dt;
    double cry_p = charge_weight*dy_ov_dt;
    double crz_p = charge_weight*dz_ov_dt;
    
    // variable declaration
    double xpn, ypn, zpn;
    double delta, delta2, delta3, delta4;
    // arrays used for the Esirkepov projection method
    double Sx0[7], Sx1[7], Sy0[7], Sy1[7], Sz0[7], Sz1[7], DSx[7], DSy[7], DSz[7];
    double tmpJx[7][7], tmpJy[7][7], tmpJz[7][7];
    
    for (unsigned int i=0; i<7; i++) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
        Sz1[i] = 0.;
    }
    for (unsigned int j=0; j<7; j++)
        for (unsigned int k=0; k<7; k++)
            tmpJx[j][k] = 0.;
    for (unsigned int i=0; i<7; i++)
        for (unsigned int k=0; k<7; k++)
            tmpJy[i][k] = 0.;
    for (unsigned int i=0; i<7; i++)
        for (unsigned int j=0; j<7; j++)
            tmpJz[i][j] = 0.;
    
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
    xpn = particles.position(0, ipart) * dx_inv_;
    int ip = round(xpn);
    int ipo = iold[0*nparts];
    int ip_m_ipo = ip-ipo-i_domain_begin;
    delta  = xpn - (double)ip;
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sx1[ip_m_ipo+1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sx1[ip_m_ipo+2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sx1[ip_m_ipo+3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sx1[ip_m_ipo+4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sx1[ip_m_ipo+5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    
    ypn = particles.position(1, ipart) * dy_inv_;
    int jp = round(ypn);
    int jpo = iold[1*nparts];
    int jp_m_jpo = jp-jpo-j_domain_begin;
    delta  = ypn - (double)jp;
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sy1[jp_m_jpo+1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sy1[jp_m_jpo+2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy1[jp_m_jpo+3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sy1[jp_m_jpo+4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy1[jp_m_jpo+5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    
    zpn = particles.position(2, ipart) * dz_inv_;
    int kp = round(zpn);
    int kpo = iold[2*nparts];
    int kp_m_kpo = kp-kpo-k_domain_begin;
    delta  = zpn - (double)kp;
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sz1[kp_m_kpo+1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sz1[kp_m_kpo+2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sz1[kp_m_kpo+3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sz1[kp_m_kpo+4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sz1[kp_m_kpo+5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    
    // computes Esirkepov coefficients
    for (unsigned int i=0; i < 7; i++) {
        DSx[i] = Sx1[i] - Sx0[i];
        DSy[i] = Sy1[i] - Sy0[i];
        DSz[i] = Sz1[i] - Sz0[i];
    }
    
    // ---------------------------
    // Calculate the total current
    // ---------------------------
    
    ipo -= bin+3;   //This minus 3 come from the order 4 scheme, based on a 7 points stencil from -3 to +3.
    // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
    jpo -= 3;
    kpo -= 3;
    
    int iloc, jloc, kloc, linindex;

    // Jx^(d,p,p)
    for (unsigned int i=1 ; i<7 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<7 ; j++) {
            jloc = j+jpo;
            for (unsigned int k=0 ; k<7 ; k++) {
                tmpJx[j][k] -= crx_p * DSx[i-1] * (Sy0[j]*Sz0[k] + 0.5*DSy[j]*Sz0[k] + 0.5*DSz[k]*Sy0[j] + one_third*DSy[j]*DSz[k]);
                kloc = k+kpo;
                linindex = iloc*b_dim[2]*b_dim[1]+jloc*b_dim[2]+kloc;
                Jx [linindex] += tmpJx[j][k]; // iloc = (i+ipo)*b_dim[1];
            }
        }
    }//i
    
    // Jy^(p,d,p)
    for (unsigned int i=0 ; i<7 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=1 ; j<7 ; j++) {
            jloc = j+jpo;
            for (unsigned int k=0 ; k<7 ; k++) {
                tmpJy[i][k] -= cry_p * DSy[j-1] * (Sz0[k]*Sx0[i] + 0.5*DSz[k]*Sx0[i] + 0.5*DSx[i]*Sz0[k] + one_third*DSz[k]*DSx[i]);
                kloc = k+kpo;
                linindex = iloc*b_dim[2]*(b_dim[1]+1)+jloc*b_dim[2]+kloc;
                Jy [linindex] += tmpJy[i][k]; //
            }
        }
    }//i
    
    // Jz^(p,p,d)
    for (unsigned int i=0 ; i<7 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<7 ; j++) {
            jloc = j+jpo;
            for (unsigned int k=1 ; k<7 ; k++) {
                tmpJz[i][j] -= crz_p * DSz[k-1] * (Sx0[i]*Sy0[j] + 0.5*DSx[i]*Sy0[j] + 0.5*DSy[j]*Sx0[i] + one_third*DSx[i]*DSy[j]);
                kloc = k+kpo;
                linindex = iloc*(b_dim[2]+1)*b_dim[1]+jloc*(b_dim[2]+1)+kloc;
                Jz [linindex] += tmpJz[i][j]; //
            }
        }
    }//i
    
    
} // END Project local current densities (Jx, Jy, Jz, sort)


// ---------------------------------------------------------------------------------------------------------------------
//! Project local current densities (sort)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D4Order::operator() (double* Jx, double* Jy, double* Jz, double* rho, Particles &particles, unsigned int ipart, double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaold)
{
    int nparts = particles.size();

    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
    double crx_p = charge_weight*dx_ov_dt;
    double cry_p = charge_weight*dy_ov_dt;
    double crz_p = charge_weight*dz_ov_dt;
    
    // variable declaration
    double xpn, ypn, zpn;
    double delta, delta2, delta3, delta4;
    // arrays used for the Esirkepov projection method
    double Sx0[7], Sx1[7], Sy0[7], Sy1[7], Sz0[7], Sz1[7], DSx[7], DSy[7], DSz[7];
    double tmpJx[7][7], tmpJy[7][7], tmpJz[7][7];
    
    for (unsigned int i=0; i<7; i++) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
        Sz1[i] = 0.;
    }
    
    for (unsigned int j=0; j<7; j++)
        for (unsigned int k=0; k<7; k++)
            tmpJx[j][k] = 0.;
    for (unsigned int i=0; i<7; i++)
        for (unsigned int k=0; k<7; k++)
            tmpJy[i][k] = 0.;
    for (unsigned int i=0; i<7; i++)
        for (unsigned int j=0; j<7; j++)
            tmpJz[i][j] = 0.;
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
    xpn = particles.position(0, ipart) * dx_inv_;
    int ip = round(xpn);
    int ipo = iold[0*nparts];
    int ip_m_ipo = ip-ipo-i_domain_begin;
    delta  = xpn - (double)ip;
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sx1[ip_m_ipo+1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sx1[ip_m_ipo+2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sx1[ip_m_ipo+3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sx1[ip_m_ipo+4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sx1[ip_m_ipo+5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    
    ypn = particles.position(1, ipart) * dy_inv_;
    int jp = round(ypn);
    int jpo = iold[1*nparts];
    int jp_m_jpo = jp-jpo-j_domain_begin;
    delta  = ypn - (double)jp;
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sy1[jp_m_jpo+1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sy1[jp_m_jpo+2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy1[jp_m_jpo+3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sy1[jp_m_jpo+4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy1[jp_m_jpo+5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    
    zpn = particles.position(2, ipart) * dz_inv_;
    int kp = round(zpn);
    int kpo = iold[2*nparts];
    int kp_m_kpo = kp-kpo-k_domain_begin;
    delta  = zpn - (double)kp;
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sz1[kp_m_kpo+1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sz1[kp_m_kpo+2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sz1[kp_m_kpo+3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sz1[kp_m_kpo+4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sz1[kp_m_kpo+5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    
    // computes Esirkepov coefficients
    for (unsigned int i=0; i < 7; i++) {
        DSx[i] = Sx1[i] - Sx0[i];
        DSy[i] = Sy1[i] - Sy0[i];
        DSz[i] = Sz1[i] - Sz0[i];
    }
        
    // ---------------------------
    // Calculate the total current
    // ---------------------------

    ipo -= bin+3;   //This minus 3 come from the order 3 scheme, based on a 7 points stencil from -3 to +3.
                    // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
    jpo -= 3;
    kpo -= 3;
    
    int iloc, jloc, kloc, linindex;
    
    // Jx^(d,p,p)
    for (unsigned int i=1 ; i<7 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<7 ; j++) {
            jloc = j+jpo;
            for (unsigned int k=0 ; k<7 ; k++) {
                tmpJx[j][k] -= crx_p * DSx[i-1] * (Sy0[j]*Sz0[k] + 0.5*DSy[j]*Sz0[k] + 0.5*DSz[k]*Sy0[j] + one_third*DSy[j]*DSz[k]);
                kloc = k+kpo;
                linindex = iloc*b_dim[2]*b_dim[1]+jloc*b_dim[2]+kloc;
                Jx [linindex] += tmpJx[j][k]; // iloc = (i+ipo)*b_dim[1];
            }
        }
    }//i
    
    // Jy^(p,d,p)
    for (unsigned int i=0 ; i<7 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=1 ; j<7 ; j++) {
            jloc = j+jpo;
            for (unsigned int k=0 ; k<7 ; k++) {
                tmpJy[i][k] -= cry_p * DSy[j-1] * (Sz0[k]*Sx0[i] + 0.5*DSz[k]*Sx0[i] + 0.5*DSx[i]*Sz0[k] + one_third*DSz[k]*DSx[i]);
                kloc = k+kpo;
                linindex = iloc*b_dim[2]*(b_dim[1]+1)+jloc*b_dim[2]+kloc;
                Jy [linindex] += tmpJy[i][k]; //
            }
        }
    }//i
    
    // Jz^(p,p,d)
    for (unsigned int i=0 ; i<7 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<7 ; j++) {
            jloc = j+jpo;
            for (unsigned int k=1 ; k<7 ; k++) {
                tmpJz[i][j] -= crz_p * DSz[k-1] * (Sx0[i]*Sy0[j] + 0.5*DSx[i]*Sy0[j] + 0.5*DSy[j]*Sx0[i] + one_third*DSx[i]*DSy[j]);
                kloc = k+kpo;
                linindex = iloc*(b_dim[2]+1)*b_dim[1]+jloc*(b_dim[2]+1)+kloc;
                Jz [linindex] += tmpJz[i][j]; //
            }
        }
    }//i

    // Rho^(p,p,p)
    for (unsigned int i=0 ; i<7 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<7 ; j++) {
            jloc = j+jpo;
            for (unsigned int k=0 ; k<7 ; k++) {
                kloc = k+kpo;
                linindex = iloc*b_dim[2]*b_dim[1]+jloc*b_dim[2]+kloc;
                rho[linindex] += charge_weight * Sx1[i]*Sy1[j]*Sz1[k];
            }
        }
    }//i

} // END Project local densities (Jx, Jy, Jz, rho, sort)


// ---------------------------------------------------------------------------------------------------------------------
//! Project local densities only (Frozen species)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D4Order::operator() (double* rho, Particles &particles, unsigned int ipart, unsigned int bin, std::vector<unsigned int> &b_dim)
{
    //Warning : this function is used for frozen species only. It is assumed that position = position_old !!!

    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------

    int iloc,jloc;
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);

    // variable declaration
    double xpn, ypn, zpn;
    double delta, delta2, delta3, delta4;
    double Sx1[7], Sy1[7], Sz1[7]; // arrays used for the Esirkepov projection method

// Initialize all current-related arrays to zero
    for (unsigned int i=0; i<7; i++) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
        Sz1[i] = 0.;
    }

    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------

    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position(0, ipart) * dx_inv_;
    int ip = round(xpn);
    delta  = xpn - (double)ip;
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sx1[1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sx1[2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sx1[3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sx1[4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sx1[5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;

    ypn = particles.position(1, ipart) * dy_inv_;
    int jp = round(ypn);
    delta  = ypn - (double)jp;
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;
    Sy1[1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sy1[2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy1[3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sy1[4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy1[5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;

    zpn = particles.position(2, ipart) * dz_inv_;
    int kp = round(zpn);
    delta  = zpn - (double)kp;
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
    ip -= i_domain_begin + bin +3;
    jp -= j_domain_begin + 3;
    kp -= k_domain_begin + 3;

    for (unsigned int i=0 ; i<7 ; i++) {
        iloc = (i+ip)*b_dim[2]*b_dim[1];
        for (unsigned int j=0 ; j<7 ; j++) {
            jloc = (jp+j)*b_dim[2];
            for (unsigned int k=0 ; k<7 ; k++) {
                rho[iloc+jloc+kp+k] += charge_weight * Sx1[i]*Sy1[j]*Sz1[k];
            }
        }
    }//i

} // END Project local current densities (Frozen species)

// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities (ionize)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D4Order::operator() (Field* Jx, Field* Jy, Field* Jz, Particles &particles, int ipart, LocalFields Jion)
{
    ERROR("Projection of ionization current not yet defined for 3D 4nd order");

} // END Project global current densities (ionize)

//Wrapper for projection
void Projector3D4Order::operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread, int ibin, int clrw, bool diag_flag, bool is_spectral, std::vector<unsigned int> &b_dim, int ispec)
{
    std::vector<int> *iold = &(smpi->dynamics_iold[ithread]);
    std::vector<double> *delta = &(smpi->dynamics_deltaold[ithread]);
    std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);
    
    int dim1 = EMfields->dimPrim[1];
    int dim2 = EMfields->dimPrim[2];

    // If no field diagnostics this timestep, then the projection is done directly on the total arrays
    if (!diag_flag){ 
        if (!is_spectral) {
            double* b_Jx =  &(*EMfields->Jx_ )(ibin*clrw* dim1   * dim2   );
            double* b_Jy =  &(*EMfields->Jy_ )(ibin*clrw*(dim1+1)* dim2   );
            double* b_Jz =  &(*EMfields->Jz_ )(ibin*clrw* dim1   *(dim2+1));
            for ( int ipart=istart ; ipart<iend; ipart++ )
                (*this)(b_Jx , b_Jy , b_Jz , particles,  ipart, (*invgf)[ipart], ibin*clrw, b_dim, &(*iold)[ipart], &(*delta)[ipart]);
        }
        else {
            double* b_Jx =  &(*EMfields->Jx_ )(ibin*clrw* dim1   * dim2   );
            double* b_Jy =  &(*EMfields->Jy_ )(ibin*clrw*(dim1+1)* dim2   );
            double* b_Jz =  &(*EMfields->Jz_ )(ibin*clrw* dim1   *(dim2+1));
            double* b_rho=  &(*EMfields->rho_)(ibin*clrw* dim1   * dim2   );
            for ( int ipart=istart ; ipart<iend; ipart++ )
                (*this)(b_Jx , b_Jy , b_Jz , b_rho , particles,  ipart, (*invgf)[ipart], ibin*clrw, b_dim, &(*iold)[ipart], &(*delta)[ipart]);
        }           
        // Otherwise, the projection may apply to the species-specific arrays
    } else {
        double* b_Jx  = EMfields->Jx_s [ispec] ? &(*EMfields->Jx_s [ispec])(ibin*clrw* dim1   *dim2) : &(*EMfields->Jx_ )(ibin*clrw* dim1   *dim2) ;
        double* b_Jy  = EMfields->Jy_s [ispec] ? &(*EMfields->Jy_s [ispec])(ibin*clrw*(dim1+1)*dim2) : &(*EMfields->Jy_ )(ibin*clrw*(dim1+1)*dim2) ;
        double* b_Jz  = EMfields->Jz_s [ispec] ? &(*EMfields->Jz_s [ispec])(ibin*clrw*dim1*(dim2+1)) : &(*EMfields->Jz_ )(ibin*clrw*dim1*(dim2+1)) ;
        double* b_rho = EMfields->rho_s[ispec] ? &(*EMfields->rho_s[ispec])(ibin*clrw* dim1   *dim2) : &(*EMfields->rho_)(ibin*clrw* dim1   *dim2) ;
        for ( int ipart=istart ; ipart<iend; ipart++ )
            (*this)(b_Jx , b_Jy , b_Jz ,b_rho, particles,  ipart, (*invgf)[ipart], ibin*clrw, b_dim, &(*iold)[ipart], &(*delta)[ipart]);
    }

}
