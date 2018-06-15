#include "Projector3D2Order_susceptibility.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field3D.h"
#include "Particles.h"
#include "Tools.h"
#include "Patch.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Projector3D2Order_susceptibility
// ---------------------------------------------------------------------------------------------------------------------
Projector3D2Order_susceptibility::Projector3D2Order_susceptibility (Params& params, Patch* patch) : Projector3D(params, patch)
{
    dx_inv_   = 1.0/params.cell_length[0];
    dy_inv_   = 1.0/params.cell_length[1];
    dz_inv_   = 1.0/params.cell_length[2];
  
    one_third = 1.0/3.0;

    i_domain_begin = patch->getCellStartingGlobalIndex(0);
    j_domain_begin = patch->getCellStartingGlobalIndex(1);
    k_domain_begin = patch->getCellStartingGlobalIndex(2);

    dt             = params.timestep;
    dts2           = params.timestep/2.;
    dts4           = params.timestep/4.;

    DEBUG("cell_length "<< params.cell_length[0]);

}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Projector3D2Order
// ---------------------------------------------------------------------------------------------------------------------
Projector3D2Order_susceptibility::~Projector3D2Order_susceptibility()
{
}


// ---------------------------------------------------------------------------------------------------------------------
//! Project local currents (sort)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D2Order_susceptibility::operator() (double* Jx, double* Jy, double* Jz, Particles &particles, unsigned int ipart, double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaold)
{
    // int nparts = particles.size();
    // 
    // // -------------------------------------
    // // Variable declaration & initialization
    // // -------------------------------------
    // 
    // // (x,y,z) components of the current density for the macro-particle
    // double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
    // double crx_p = charge_weight*dx_ov_dt;
    // double cry_p = charge_weight*dy_ov_dt;
    // double crz_p = charge_weight*dz_ov_dt;
    // 
    // // variable declaration
    // double xpn, ypn, zpn;
    // double delta, delta2;
    // // arrays used for the Esirkepov projection method
    // double Sx0[5], Sx1[5], Sy0[5], Sy1[5], Sz0[5], Sz1[5], DSx[5], DSy[5], DSz[5];
    // double tmpJx[5][5], tmpJy[5][5], tmpJz[5][5];
    // 
    // for (unsigned int i=0; i<5; i++) {
    //     Sx1[i] = 0.;
    //     Sy1[i] = 0.;
    //     Sz1[i] = 0.;
    // }
    // for (unsigned int j=0; j<5; j++)
    //     for (unsigned int k=0; k<5; k++)
    //         tmpJx[j][k] = 0.;
    // for (unsigned int i=0; i<5; i++)
    //     for (unsigned int k=0; k<5; k++)
    //         tmpJy[i][k] = 0.;
    // for (unsigned int i=0; i<5; i++)
    //     for (unsigned int j=0; j<5; j++)
    //         tmpJz[i][j] = 0.;
    // 
    // // --------------------------------------------------------
    // // Locate particles & Calculate Esirkepov coef. S, DS and W
    // // --------------------------------------------------------
    // 
    // // locate the particle on the primal grid at former time-step & calculate coeff. S0
    // delta = deltaold[0*nparts];
    // delta2 = delta*delta;
    // Sx0[0] = 0.;
    // Sx0[1] = 0.5 * (delta2-delta+0.25);
    // Sx0[2] = 0.75-delta2;
    // Sx0[3] = 0.5 * (delta2+delta+0.25);
    // Sx0[4] = 0.;
    // 
    // delta = deltaold[1*nparts];
    // delta2 = delta*delta;
    // Sy0[0] = 0.;
    // Sy0[1] = 0.5 * (delta2-delta+0.25);
    // Sy0[2] = 0.75-delta2;
    // Sy0[3] = 0.5 * (delta2+delta+0.25);
    // Sy0[4] = 0.;
    // 
    // delta = deltaold[2*nparts];
    // delta2 = delta*delta;
    // Sz0[0] = 0.;
    // Sz0[1] = 0.5 * (delta2-delta+0.25);
    // Sz0[2] = 0.75-delta2;
    // Sz0[3] = 0.5 * (delta2+delta+0.25);
    // Sz0[4] = 0.;
    // 
    // // locate the particle on the primal grid at current time-step & calculate coeff. S1
    // xpn = particles.position(0, ipart) * dx_inv_;
    // int ip = round(xpn);
    // int ipo = iold[0*nparts];
    // int ip_m_ipo = ip-ipo-i_domain_begin;
    // delta  = xpn - (double)ip;
    // delta2 = delta*delta;
    // Sx1[ip_m_ipo+1] = 0.5 * (delta2-delta+0.25);
    // Sx1[ip_m_ipo+2] = 0.75-delta2;
    // Sx1[ip_m_ipo+3] = 0.5 * (delta2+delta+0.25);
    // 
    // ypn = particles.position(1, ipart) * dy_inv_;
    // int jp = round(ypn);
    // int jpo = iold[1*nparts];
    // int jp_m_jpo = jp-jpo-j_domain_begin;
    // delta  = ypn - (double)jp;
    // delta2 = delta*delta;
    // Sy1[jp_m_jpo+1] = 0.5 * (delta2-delta+0.25);
    // Sy1[jp_m_jpo+2] = 0.75-delta2;
    // Sy1[jp_m_jpo+3] = 0.5 * (delta2+delta+0.25);
    // 
    // zpn = particles.position(2, ipart) * dz_inv_;
    // int kp = round(zpn);
    // int kpo = iold[2*nparts];
    // int kp_m_kpo = kp-kpo-k_domain_begin;
    // delta  = zpn - (double)kp;
    // delta2 = delta*delta;
    // Sz1[kp_m_kpo+1] = 0.5 * (delta2-delta+0.25);
    // Sz1[kp_m_kpo+2] = 0.75-delta2;
    // Sz1[kp_m_kpo+3] = 0.5 * (delta2+delta+0.25);
    // 
    // // computes Esirkepov coefficients
    // for (unsigned int i=0; i < 5; i++) {
    //     DSx[i] = Sx1[i] - Sx0[i];
    //     DSy[i] = Sy1[i] - Sy0[i];
    //     DSz[i] = Sz1[i] - Sz0[i];
    // }
    // 
    // // ---------------------------
    // // Calculate the total current
    // // ---------------------------
    // 
    // ipo -= (int)bin+2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
    // // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
    // jpo -= 2;
    // kpo -= 2;
    // 
    // int linindex, linindex_x, linindex_y;
    // double tmp, tmp2;
    // double vtmp[5];
    // 
    // // Jx^(d,p,p)
    // int  z_size = b_dim[2];
    // int yz_size = b_dim[2]*b_dim[1];
    // int linindex0 = ipo*yz_size+jpo*z_size+kpo;
    // tmp = 0.;
    // linindex = linindex0;
    // tmp2 = crx_p * (one_third*Sy1[0]*Sz1[0]);
    // for (int i=1 ; i<5 ; i++) {
    //     tmp -= DSx[i-1] * tmp2;
    //     linindex += yz_size;
    //     Jx [linindex] += tmp; // iloc = (i+ipo)*b_dim[1];
    // }//i
    // for ( unsigned int i=0 ; i<5 ; i++) vtmp[i] = 0.;
    // linindex_x = linindex0;
    // for (int k=1 ; k<5 ; k++) {
    //     linindex_x += 1;
    //     linindex    = linindex_x;
    //     tmp = crx_p * (0.5*Sy1[0]*Sz0[k] + one_third*Sy1[0]*DSz[k]);
    //     for (int i=1 ; i<5 ; i++) {
    //         vtmp[k] -= DSx[i-1] * tmp;
    //         linindex += yz_size;
    //         Jx [linindex] += vtmp[k]; // iloc = (i+ipo)*b_dim[1];
    //     }
    // }//i
    // for ( unsigned int i=0 ; i<5 ; i++) vtmp[i] = 0.;
    // linindex_x = linindex0;
    // for (int j=1 ; j<5 ; j++) {
    //     linindex_x += z_size;
    //     linindex    = linindex_x;
    //     tmp = crx_p * (0.5*Sz1[0]*Sy0[j] + one_third*DSy[j]*Sz1[0]);
    //     for (int i=1 ; i<5 ; i++) {
    //         vtmp[j] -= DSx[i-1] * tmp;
    //         linindex += yz_size;
    //         Jx [linindex] += vtmp[j]; // iloc = (i+ipo)*b_dim[1];
    //     }
    // }//i
    // linindex_x = linindex0;
    // for (int j=1 ; j<5 ; j++) {
    //     linindex_x += z_size;
    //     linindex_y  = linindex_x;
    //     for (int k=1 ; k<5 ; k++) {
    //         linindex_y += 1;
    //         linindex    = linindex_y;
    //         tmp = crx_p * (Sy0[j]*Sz0[k] + 0.5*DSy[j]*Sz0[k] + 0.5*DSz[k]*Sy0[j] + one_third*DSy[j]*DSz[k]);
    //         for (int i=1 ; i<5 ; i++) {
    //             tmpJx[j][k] -= DSx[i-1] * tmp;
    //             linindex += yz_size;
    //             Jx [linindex] += tmpJx[j][k]; // iloc = (i+ipo)*b_dim[1];
    //         }
    //     }
    // }//i
    // 
    // 
    // // Jy^(p,d,p)
    // yz_size = b_dim[2]*(b_dim[1]+1);
    // linindex0 = ipo*yz_size+jpo*z_size+kpo;
    // tmp = 0.;
    // linindex = linindex0;
    // tmp2 = cry_p * (one_third*Sz1[0]*Sx1[0]);
    // for (int j=1 ; j<5 ; j++) {
    //     tmp -= DSy[j-1] * tmp2;
    //     linindex += z_size;
    //     Jy [linindex] += tmp; //
    // }//i
    // for ( unsigned int i=0 ; i<5 ; i++) vtmp[i] = 0.;
    // linindex_x = linindex0;
    // for (int k=1 ; k<5 ; k++) {
    //     linindex_x += 1;
    //     linindex    = linindex_x;
    //     tmp  = cry_p * (0.5*Sx1[0]*Sz0[k] + one_third*DSz[k]*Sx1[0]);
    //     for (int j=1 ; j<5 ; j++) {
    //         vtmp[k] -= DSy[j-1] * tmp;
    //         linindex += z_size;
    //         Jy [linindex] += vtmp[k]; //
    //     }
    // }
    // for ( unsigned int i=0 ; i<5 ; i++) vtmp[i] = 0.;
    // linindex_x = linindex0;
    // for (int i=1 ; i<5 ; i++) {
    //     linindex_x += yz_size;
    //     linindex    = linindex_x;
    //     tmp = cry_p * (0.5*Sz1[0]*Sx0[i] + one_third*Sz1[0]*DSx[i]); 
    //     for (int j=1 ; j<5 ; j++) {
    //         vtmp[i] -= DSy[j-1] * tmp;
    //         linindex += z_size;
    //         Jy [linindex] += vtmp[i]; //
    //     }
    // }//i
    // linindex_x = linindex0;
    // for (int i=1 ; i<5 ; i++) {
    //     linindex_x += yz_size;
    //     linindex_y  = linindex_x;
    //     for (int k=1 ; k<5 ; k++) {
    //         linindex_y += 1;
    //         linindex    = linindex_y;
    //         tmp = cry_p * (Sz0[k]*Sx0[i] + 0.5*DSz[k]*Sx0[i] + 0.5*DSx[i]*Sz0[k] + one_third*DSz[k]*DSx[i]);
    //         for (int j=1 ; j<5 ; j++) {
    //             tmpJy[i][k] -= DSy[j-1] * tmp;
    //             linindex +=z_size;
    //             Jy [linindex] += tmpJy[i][k]; //
    //         }
    //     }
    // }//i
    // 
    // // Jz^(p,p,d)
    // z_size =  b_dim[2]+1;
    // yz_size = (b_dim[2]+1)*b_dim[1];
    // linindex0 = ipo*yz_size+jpo*z_size+kpo;
    // tmp = 0.;
    // linindex = linindex0;
    // tmp2 = crz_p * (one_third*Sx1[0]*Sy1[0]);
    // for (int k=1 ; k<5 ; k++) {
    //     tmp -= DSz[k-1] * tmp2;
    //     linindex += 1;
    //     Jz [linindex] += tmp; //
    // }//i
    // for ( unsigned int i=0 ; i<5 ; i++) vtmp[i] = 0.;
    // linindex_x = linindex0;
    // for (int j=1 ; j<5 ; j++) {
    //     linindex_x += z_size;
    //     linindex    = linindex_x; 
    //     tmp = crz_p * (0.5*Sx1[0]*Sy0[j] + one_third*Sx1[0]*DSy[j]);
    //     for (int k=1 ; k<5 ; k++) {
    //         vtmp[j] -= DSz[k-1] * tmp;
    //         linindex += 1;
    //         Jz [linindex] += vtmp[j]; //
    //      }
    // }//i
    // for ( unsigned int i=0 ; i<5 ; i++) vtmp[i] = 0.;
    // linindex_x = linindex0;
    // for (int i=1 ; i<5 ; i++) {
    //     linindex_x += yz_size;
    //     linindex    = linindex_x;
    //     tmp = crz_p * (0.5*Sy1[0]*Sx0[i] + one_third*DSx[i]*Sy1[0]);
    //     for (int k=1 ; k<5 ; k++) {
    //         vtmp[i] -= DSz[k-1] * tmp;
    //         linindex += 1;
    //         Jz [linindex] += vtmp[i]; //
    //     }
    // }//i
    // linindex_x = linindex0;
    // for (int i=1 ; i<5 ; i++) {
    //     linindex_x += yz_size;
    //     linindex_y  = linindex_x;
    //     for (int j=1 ; j<5 ; j++) {
    //         linindex_y += z_size;
    //         linindex    = linindex_y;
    //         tmp = crz_p*(Sx0[i]*Sy0[j] + 0.5*DSx[i]*Sy0[j] + 0.5*DSy[j]*Sx0[i] + one_third*DSx[i]*DSy[j]);
    //         for (int k=1 ; k<5 ; k++) {
    //             tmpJz[i][j] -= DSz[k-1] * tmp;
    //             linindex += 1;
    //             Jz [linindex] += tmpJz[i][j]; //
    //         }
    //     }
    // }//i
    // 
    
} // END Project local current densities (Jx, Jy, Jz, sort)


// ---------------------------------------------------------------------------------------------------------------------
//! Project local current densities (sort)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D2Order_susceptibility::operator() (double* Jx, double* Jy, double* Jz, double* rho, Particles &particles, unsigned int ipart, double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaold)
{
    // int nparts = particles.size();
    // 
    // // -------------------------------------
    // // Variable declaration & initialization
    // // -------------------------------------
    // 
    // // (x,y,z) components of the current density for the macro-particle
    // double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
    // double crx_p = charge_weight*dx_ov_dt;
    // double cry_p = charge_weight*dy_ov_dt;
    // double crz_p = charge_weight*dz_ov_dt;
    // 
    // // variable declaration
    // double xpn, ypn, zpn;
    // double delta, delta2;
    // // arrays used for the Esirkepov projection method
    // double Sx0[5], Sx1[5], Sy0[5], Sy1[5], Sz0[5], Sz1[5], DSx[5], DSy[5], DSz[5];
    // double tmpJx[5][5], tmpJy[5][5], tmpJz[5][5];
    // 
    // for (unsigned int i=0; i<5; i++) {
    //     Sx1[i] = 0.;
    //     Sy1[i] = 0.;
    //     Sz1[i] = 0.;
    // }
    // 
    // for (unsigned int j=0; j<5; j++)
    //     for (unsigned int k=0; k<5; k++)
    //         tmpJx[j][k] = 0.;
    // for (unsigned int i=0; i<5; i++)
    //     for (unsigned int k=0; k<5; k++)
    //         tmpJy[i][k] = 0.;
    // for (unsigned int i=0; i<5; i++)
    //     for (unsigned int j=0; j<5; j++)
    //         tmpJz[i][j] = 0.;
    // // --------------------------------------------------------
    // // Locate particles & Calculate Esirkepov coef. S, DS and W
    // // --------------------------------------------------------
    // 
    // // locate the particle on the primal grid at former time-step & calculate coeff. S0
    // delta = deltaold[0*nparts];
    // delta2 = delta*delta;
    // Sx0[0] = 0.;
    // Sx0[1] = 0.5 * (delta2-delta+0.25);
    // Sx0[2] = 0.75-delta2;
    // Sx0[3] = 0.5 * (delta2+delta+0.25);
    // Sx0[4] = 0.;
    // 
    // delta = deltaold[1*nparts];
    // delta2 = delta*delta;
    // Sy0[0] = 0.;
    // Sy0[1] = 0.5 * (delta2-delta+0.25);
    // Sy0[2] = 0.75-delta2;
    // Sy0[3] = 0.5 * (delta2+delta+0.25);
    // Sy0[4] = 0.;
    // 
    // delta = deltaold[2*nparts];
    // delta2 = delta*delta;
    // Sz0[0] = 0.;
    // Sz0[1] = 0.5 * (delta2-delta+0.25);
    // Sz0[2] = 0.75-delta2;
    // Sz0[3] = 0.5 * (delta2+delta+0.25);
    // Sz0[4] = 0.;
    // 
    // // locate the particle on the primal grid at current time-step & calculate coeff. S1
    // xpn = particles.position(0, ipart) * dx_inv_;
    // int ip = round(xpn);
    // int ipo = iold[0*nparts];
    // int ip_m_ipo = ip-ipo-i_domain_begin;
    // delta  = xpn - (double)ip;
    // delta2 = delta*delta;
    // Sx1[ip_m_ipo+1] = 0.5 * (delta2-delta+0.25);
    // Sx1[ip_m_ipo+2] = 0.75-delta2;
    // Sx1[ip_m_ipo+3] = 0.5 * (delta2+delta+0.25);
    // 
    // ypn = particles.position(1, ipart) * dy_inv_;
    // int jp = round(ypn);
    // int jpo = iold[1*nparts];
    // int jp_m_jpo = jp-jpo-j_domain_begin;
    // delta  = ypn - (double)jp;
    // delta2 = delta*delta;
    // Sy1[jp_m_jpo+1] = 0.5 * (delta2-delta+0.25);
    // Sy1[jp_m_jpo+2] = 0.75-delta2;
    // Sy1[jp_m_jpo+3] = 0.5 * (delta2+delta+0.25);
    // 
    // zpn = particles.position(2, ipart) * dz_inv_;
    // int kp = round(zpn);
    // int kpo = iold[2*nparts];
    // int kp_m_kpo = kp-kpo-k_domain_begin;
    // delta  = zpn - (double)kp;
    // delta2 = delta*delta;
    // Sz1[kp_m_kpo+1] = 0.5 * (delta2-delta+0.25);
    // Sz1[kp_m_kpo+2] = 0.75-delta2;
    // Sz1[kp_m_kpo+3] = 0.5 * (delta2+delta+0.25);
    // 
    // // computes Esirkepov coefficients
    // for (unsigned int i=0; i < 5; i++) {
    //     DSx[i] = Sx1[i] - Sx0[i];
    //     DSy[i] = Sy1[i] - Sy0[i];
    //     DSz[i] = Sz1[i] - Sz0[i];
    // }
    // 
    // // ---------------------------
    // // Calculate the total current
    // // ---------------------------
    // 
    // ipo -= bin+2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
    //                 // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
    // jpo -= 2;
    // kpo -= 2;
    // 
    // int iloc, jloc, kloc, linindex;
    // 
    // // Jx^(d,p,p)
    // for (unsigned int i=1 ; i<5 ; i++) {
    //     iloc = i+ipo;
    //     for (unsigned int j=0 ; j<5 ; j++) {
    //         jloc = j+jpo;
    //         for (unsigned int k=0 ; k<5 ; k++) {
    //             tmpJx[j][k] -= crx_p * DSx[i-1] * (Sy0[j]*Sz0[k] + 0.5*DSy[j]*Sz0[k] + 0.5*DSz[k]*Sy0[j] + one_third*DSy[j]*DSz[k]);
    //             kloc = k+kpo;
    //             linindex = iloc*b_dim[2]*b_dim[1]+jloc*b_dim[2]+kloc;
    //             Jx [linindex] += tmpJx[j][k]; // iloc = (i+ipo)*b_dim[1];
    //         }
    //     }
    // }//i
    // 
    // // Jy^(p,d,p)
    // for (unsigned int i=0 ; i<5 ; i++) {
    //     iloc = i+ipo;
    //     for (unsigned int j=1 ; j<5 ; j++) {
    //         jloc = j+jpo;
    //         for (unsigned int k=0 ; k<5 ; k++) {
    //             tmpJy[i][k] -= cry_p * DSy[j-1] * (Sz0[k]*Sx0[i] + 0.5*DSz[k]*Sx0[i] + 0.5*DSx[i]*Sz0[k] + one_third*DSz[k]*DSx[i]);
    //             kloc = k+kpo;
    //             linindex = iloc*b_dim[2]*(b_dim[1]+1)+jloc*b_dim[2]+kloc;
    //             Jy [linindex] += tmpJy[i][k]; //
    //         }
    //     }
    // }//i
    // 
    // // Jz^(p,p,d)
    // for (unsigned int i=0 ; i<5 ; i++) {
    //     iloc = i+ipo;
    //     for (unsigned int j=0 ; j<5 ; j++) {
    //         jloc = j+jpo;
    //         for (unsigned int k=1 ; k<5 ; k++) {
    //             tmpJz[i][j] -= crz_p * DSz[k-1] * (Sx0[i]*Sy0[j] + 0.5*DSx[i]*Sy0[j] + 0.5*DSy[j]*Sx0[i] + one_third*DSx[i]*DSy[j]);
    //             kloc = k+kpo;
    //             linindex = iloc*(b_dim[2]+1)*b_dim[1]+jloc*(b_dim[2]+1)+kloc;
    //             Jz [linindex] += tmpJz[i][j]; //
    //         }
    //     }
    // }//i
    // 
    // // Rho^(p,p,p)
    // for (unsigned int i=0 ; i<5 ; i++) {
    //     iloc = i+ipo;
    //     for (unsigned int j=0 ; j<5 ; j++) {
    //         jloc = j+jpo;
    //         for (unsigned int k=0 ; k<5 ; k++) {
    //             kloc = k+kpo;
    //             linindex = iloc*b_dim[2]*b_dim[1]+jloc*b_dim[2]+kloc;
    //             rho[linindex] += charge_weight * Sx1[i]*Sy1[j]*Sz1[k];
    //         }
    //     }
    // }//i

} // END Project local densities (Jx, Jy, Jz, rho, sort)


// ---------------------------------------------------------------------------------------------------------------------
//! Project local densities only (Frozen species)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D2Order_susceptibility::operator() (double* rho, Particles &particles, unsigned int ipart, unsigned int bin, std::vector<unsigned int> &b_dim)
{
//     //Warning : this function is used for frozen species only. It is assumed that position = position_old !!!
// 
//     // -------------------------------------
//     // Variable declaration & initialization
//     // -------------------------------------
// 
//     int iloc,jloc;
//     // (x,y,z) components of the current density for the macro-particle
//     double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
// 
//     // variable declaration
//     double xpn, ypn, zpn;
//     double delta, delta2;
//     double Sx1[5], Sy1[5], Sz1[5]; // arrays used for the Esirkepov projection method
// 
// // Initialize all current-related arrays to zero
//     for (unsigned int i=0; i<5; i++) {
//         Sx1[i] = 0.;
//         Sy1[i] = 0.;
//         Sz1[i] = 0.;
//     }
// 
//     // --------------------------------------------------------
//     // Locate particles & Calculate Esirkepov coef. S, DS and W
//     // --------------------------------------------------------
// 
//     // locate the particle on the primal grid at current time-step & calculate coeff. S1
//     xpn = particles.position(0, ipart) * dx_inv_;
//     int ip = round(xpn);
//     delta  = xpn - (double)ip;
//     delta2 = delta*delta;
//     Sx1[1] = 0.5 * (delta2-delta+0.25);
//     Sx1[2] = 0.75-delta2;
//     Sx1[3] = 0.5 * (delta2+delta+0.25);
// 
//     ypn = particles.position(1, ipart) * dy_inv_;
//     int jp = round(ypn);
//     delta  = ypn - (double)jp;
//     delta2 = delta*delta;
//     Sy1[1] = 0.5 * (delta2-delta+0.25);
//     Sy1[2] = 0.75-delta2;
//     Sy1[3] = 0.5 * (delta2+delta+0.25);
// 
//     zpn = particles.position(2, ipart) * dz_inv_;
//     int kp = round(zpn);
//     delta  = zpn - (double)kp;
//     delta2 = delta*delta;
//     Sz1[1] = 0.5 * (delta2-delta+0.25);
//     Sz1[2] = 0.75-delta2;
//     Sz1[3] = 0.5 * (delta2+delta+0.25);
// 
//     // ---------------------------
//     // Calculate the total charge
//     // ---------------------------
//     ip -= i_domain_begin + bin +2;
//     jp -= j_domain_begin + 2;
//     kp -= k_domain_begin + 2;
// 
//     for (unsigned int i=0 ; i<5 ; i++) {
//         iloc = (i+ip)*b_dim[2]*b_dim[1];
//         for (unsigned int j=0 ; j<5 ; j++) {
//             jloc = (jp+j)*b_dim[2];
//             for (unsigned int k=0 ; k<5 ; k++) {
//                 rho[iloc+jloc+kp+k] += charge_weight * Sx1[i]*Sy1[j]*Sz1[k];
//             }
//         }
//     }//i

} // END Project local current densities (Frozen species)

// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities (ionize)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D2Order_susceptibility::operator() (Field* Jx, Field* Jy, Field* Jz, Particles &particles, int ipart, LocalFields Jion)
{
    // Field3D* Jx3D  = static_cast<Field3D*>(Jx);
    // Field3D* Jy3D  = static_cast<Field3D*>(Jy);
    // Field3D* Jz3D  = static_cast<Field3D*>(Jz);
    // 
    // 
    // //Declaration of local variables
    // int ip, id, jp, jd, kp, kd;
    // double xpn, xpmxip, xpmxip2, xpmxid, xpmxid2;
    // double ypn, ypmyjp, ypmyjp2, ypmyjd, ypmyjd2;
    // double zpn, zpmzkp, zpmzkp2, zpmzkd, zpmzkd2;
    // double Sxp[3], Sxd[3], Syp[3], Syd[3], Szp[3], Szd[3];
    // 
    // // weighted currents
    // double Jx_ion = Jion.x * particles.weight(ipart);
    // double Jy_ion = Jion.y * particles.weight(ipart);
    // double Jz_ion = Jion.z * particles.weight(ipart);
    // 
    // //Locate particle on the grid
    // xpn    = particles.position(0, ipart) * dx_inv_;  // normalized distance to the first node
    // ypn    = particles.position(1, ipart) * dy_inv_;  // normalized distance to the first node
    // zpn    = particles.position(1, ipart) * dz_inv_;  // normalized distance to the first node
    // 
    // // x-primal index
    // ip      = round(xpn);                    // x-index of the central node
    // xpmxip  = xpn - (double)ip;              // normalized distance to the nearest grid point
    // xpmxip2 = xpmxip*xpmxip;                 // square of the normalized distance to the nearest grid point
    // 
    // // x-dual index
    // id      = round(xpn+0.5);                // x-index of the central node
    // xpmxid  = xpn - (double)id + 0.5;        // normalized distance to the nearest grid point
    // xpmxid2 = xpmxid*xpmxid;                 // square of the normalized distance to the nearest grid point
    // 
    // // y-primal index
    // jp      = round(ypn);                    // y-index of the central node
    // ypmyjp  = ypn - (double)jp;              // normalized distance to the nearest grid point
    // ypmyjp2 = ypmyjp*ypmyjp;                 // square of the normalized distance to the nearest grid point
    // 
    // // y-dual index
    // jd      = round(ypn+0.5);                // y-index of the central node
    // ypmyjd  = ypn - (double)jd + 0.5;        // normalized distance to the nearest grid point
    // ypmyjd2 = ypmyjd*ypmyjd;                 // square of the normalized distance to the nearest grid point
    // 
    // // z-primal index
    // kp      = round(zpn);                    // z-index of the central node
    // zpmzkp  = zpn - (double)kp;              // normalized distance to the nearest grid point
    // zpmzkp2 = zpmzkp*zpmzkp;                 // square of the normalized distance to the nearest grid point
    // 
    // // z-dual index
    // kd      = round(zpn+0.5);                // z-index of the central node
    // zpmzkd  = zpn - (double)kd + 0.5;        // normalized distance to the nearest grid point
    // zpmzkd2 = zpmzkd*zpmzkd;                 // square of the normalized distance to the nearest grid point
    // 
    // Sxp[0] = 0.5 * (xpmxip2-xpmxip+0.25);
    // Sxp[1] = (0.75-xpmxip2);
    // Sxp[2] = 0.5 * (xpmxip2+xpmxip+0.25);
    // 
    // Sxd[0] = 0.5 * (xpmxid2-xpmxid+0.25);
    // Sxd[1] = (0.75-xpmxid2);
    // Sxd[2] = 0.5 * (xpmxid2+xpmxid+0.25);
    // 
    // Syp[0] = 0.5 * (ypmyjp2-ypmyjp+0.25);
    // Syp[1] = (0.75-ypmyjp2);
    // Syp[2] = 0.5 * (ypmyjp2+ypmyjp+0.25);
    // 
    // Syd[0] = 0.5 * (ypmyjd2-ypmyjd+0.25);
    // Syd[1] = (0.75-ypmyjd2);
    // Syd[2] = 0.5 * (ypmyjd2+ypmyjd+0.25);
    // 
    // Szp[0] = 0.5 * (zpmzkp2-zpmzkp+0.25);
    // Szp[1] = (0.75-zpmzkp2);
    // Szp[2] = 0.5 * (zpmzkp2+zpmzkp+0.25);
    // 
    // Szd[0] = 0.5 * (zpmzkd2-zpmzkd+0.25);
    // Szd[1] = (0.75-zpmzkd2);
    // Szd[2] = 0.5 * (zpmzkd2+zpmzkd+0.25);
    // 
    // ip  -= i_domain_begin;
    // id  -= i_domain_begin;
    // jp  -= j_domain_begin;
    // jd  -= j_domain_begin;
    // kp  -= k_domain_begin;
    // kd  -= k_domain_begin;
    // 
    // for (unsigned int i=0 ; i<3 ; i++) {
    //     int iploc=ip+i-1;
    //     int idloc=id+i-1;
    //     for (unsigned int j=0 ; j<3 ; j++) {
    //         int jploc=jp+j-1;
    //         int jdloc=jd+j-1;
    //         for (unsigned int k=0 ; k<3 ; k++) {
    //             int kploc=kp+k-1;
    //             int kdloc=kd+k-1;
    //             // Jx^(d,p,p)
    //             (*Jx3D)(idloc,jploc,kploc) += Jx_ion * Sxd[i]*Syp[j]*Szp[k];
    //             // Jy^(p,d,p)
    //             (*Jy3D)(iploc,jdloc,kploc) += Jy_ion * Sxp[i]*Syd[j]*Szp[k];
    //             // Jz^(p,p,d)
    //             (*Jz3D)(iploc,jploc,kdloc) += Jz_ion * Sxp[i]*Syp[j]*Szd[k];
    //         }//k
    //     }//j
    // }//i
    // 


} // END Project global current densities (ionize)

//Wrapper for projection
void Projector3D2Order_susceptibility::operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread, int ibin, int clrw, bool diag_flag, bool is_spectral, std::vector<unsigned int> &b_dim, int ispec)
{
    // std::vector<int> *iold = &(smpi->dynamics_iold[ithread]);
    // std::vector<double> *delta = &(smpi->dynamics_deltaold[ithread]);
    // std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);
    // 
    // int dim1 = EMfields->dimPrim[1];
    // int dim2 = EMfields->dimPrim[2];
    // 
    // // If no field diagnostics this timestep, then the projection is done directly on the total arrays
    // if (!diag_flag){ 
    //     if (!is_spectral) {
    //         double* b_Jx =  &(*EMfields->Jx_ )(ibin*clrw* dim1   * dim2   );
    //         double* b_Jy =  &(*EMfields->Jy_ )(ibin*clrw*(dim1+1)* dim2   );
    //         double* b_Jz =  &(*EMfields->Jz_ )(ibin*clrw* dim1   *(dim2+1));
    //         for ( int ipart=istart ; ipart<iend; ipart++ )
    //             (*this)(b_Jx , b_Jy , b_Jz , particles,  ipart, (*invgf)[ipart], ibin*clrw, b_dim, &(*iold)[ipart], &(*delta)[ipart]);
    //     }
    //     else {
    //         double* b_Jx =  &(*EMfields->Jx_ )(ibin*clrw* dim1   * dim2   );
    //         double* b_Jy =  &(*EMfields->Jy_ )(ibin*clrw*(dim1+1)* dim2   );
    //         double* b_Jz =  &(*EMfields->Jz_ )(ibin*clrw* dim1   *(dim2+1));
    //         double* b_rho=  &(*EMfields->rho_)(ibin*clrw* dim1   * dim2   );
    //         for ( int ipart=istart ; ipart<iend; ipart++ )
    //             (*this)(b_Jx , b_Jy , b_Jz , b_rho , particles,  ipart, (*invgf)[ipart], ibin*clrw, b_dim, &(*iold)[ipart], &(*delta)[ipart]);
    //     }
    // // Otherwise, the projection may apply to the species-specific arrays
    // } else {
    //     double* b_Jx  = EMfields->Jx_s [ispec] ? &(*EMfields->Jx_s [ispec])(ibin*clrw* dim1   *dim2) : &(*EMfields->Jx_ )(ibin*clrw* dim1   *dim2) ;
    //     double* b_Jy  = EMfields->Jy_s [ispec] ? &(*EMfields->Jy_s [ispec])(ibin*clrw*(dim1+1)*dim2) : &(*EMfields->Jy_ )(ibin*clrw*(dim1+1)*dim2) ;
    //     double* b_Jz  = EMfields->Jz_s [ispec] ? &(*EMfields->Jz_s [ispec])(ibin*clrw*dim1*(dim2+1)) : &(*EMfields->Jz_ )(ibin*clrw*dim1*(dim2+1)) ;
    //     double* b_rho = EMfields->rho_s[ispec] ? &(*EMfields->rho_s[ispec])(ibin*clrw* dim1   *dim2) : &(*EMfields->rho_)(ibin*clrw* dim1   *dim2) ;
    //     for ( int ipart=istart ; ipart<iend; ipart++ )
    //         (*this)(b_Jx , b_Jy , b_Jz ,b_rho, particles,  ipart, (*invgf)[ipart], ibin*clrw, b_dim, &(*iold)[ipart], &(*delta)[ipart]);
    // }

}


// Wrapper for projector of susceptibility










// Projector for susceptibility used as source term in envelope equation
                                                     
void Projector3D2Order_susceptibility::project_susceptibility(double* Chi_envelope, Particles &particles, unsigned int ipart, unsigned int bin, std::vector<unsigned int> &b_dim, SmileiMPI* smpi, int ithread, double species_mass)
{
        std::vector<double> *Epart       = &(smpi->dynamics_Epart[ithread]);
        std::vector<double> *Phipart     = &(smpi->dynamics_PHIpart[ithread]);
        std::vector<double> *GradPhipart = &(smpi->dynamics_GradPHIpart[ithread]);

        int iloc,jloc;
     
        double momentum[3];
        
        double gamma_ponderomotive,gamma0,gamma0_sq;
        double charge_over_mass_dts2,charge_sq_over_mass_dts4,charge_sq_over_mass_sq;
        double pxsm, pysm, pzsm;
        double one_over_mass=1./species_mass;

        int nparts = particles.size();
        double* Ex       = &( (*Epart)[0*nparts] );
        double* Ey       = &( (*Epart)[1*nparts] );
        double* Ez       = &( (*Epart)[2*nparts] );
        double* Phi      = &( (*Phipart)[0*nparts] );
        double* GradPhix = &( (*GradPhipart)[0*nparts] );
        double* GradPhiy = &( (*GradPhipart)[1*nparts] );
        double* GradPhiz = &( (*GradPhipart)[2*nparts] );

    
    
        charge_over_mass_dts2    = (double)(particles.charge(ipart))*dts2*one_over_mass;
        // ! ponderomotive force is proportional to charge squared and the field is divided by 4 instead of 2
        charge_sq_over_mass_dts4 = (double)(particles.charge(ipart))*(double)(particles.charge(ipart))*dts4*one_over_mass;      
        // (charge over mass)^2
        charge_sq_over_mass_sq   = (double)(particles.charge(ipart))*(double)(particles.charge(ipart))*one_over_mass*one_over_mass;

        for ( int i = 0 ; i<3 ; i++ )
            momentum[i] = particles.momentum(i,ipart);

        // compute initial ponderomotive gamma 
        gamma0_sq = 1. + momentum[0]*momentum[0]+ momentum[1]*momentum[1] + momentum[2]*momentum[2] + *(Phi+ipart)*charge_sq_over_mass_sq ;
        gamma0    = sqrt(gamma0_sq) ;

        // ( electric field + ponderomotive force for ponderomotive gamma advance ) scalar multiplied by momentum
        pxsm = (gamma0 * charge_over_mass_dts2*(*(Ex+ipart)) - charge_sq_over_mass_dts4*(*(GradPhix+ipart)) ) * momentum[0] / gamma0_sq;
        pysm = (gamma0 * charge_over_mass_dts2*(*(Ey+ipart)) - charge_sq_over_mass_dts4*(*(GradPhiy+ipart)) ) * momentum[1] / gamma0_sq;
        pzsm = (gamma0 * charge_over_mass_dts2*(*(Ez+ipart)) - charge_sq_over_mass_dts4*(*(GradPhiz+ipart)) ) * momentum[2] / gamma0_sq;
        
        // update of gamma ponderomotive 
        gamma_ponderomotive = gamma0 + (pxsm+pysm+pzsm)*0.5 ;

        // susceptibility for the macro-particle
        double charge_weight = (double)(particles.charge(ipart))*(double)(particles.charge(ipart))*particles.weight(ipart)*one_over_mass/gamma_ponderomotive; 

        // variable declaration
        double xpn, ypn, zpn;
        double delta, delta2;
        double Sx1[5], Sy1[5], Sz1[5]; // arrays used for the Esirkepov projection method

    // Initialize all current-related arrays to zero
        for (unsigned int i=0; i<5; i++) {
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
        Sx1[1] = 0.5 * (delta2-delta+0.25);
        Sx1[2] = 0.75-delta2;
        Sx1[3] = 0.5 * (delta2+delta+0.25);

        ypn = particles.position(1, ipart) * dy_inv_;
        int jp = round(ypn);
        delta  = ypn - (double)jp;
        delta2 = delta*delta;
        Sy1[1] = 0.5 * (delta2-delta+0.25);
        Sy1[2] = 0.75-delta2;
        Sy1[3] = 0.5 * (delta2+delta+0.25);

        zpn = particles.position(2, ipart) * dz_inv_;
        int kp = round(zpn);
        delta  = zpn - (double)kp;
        delta2 = delta*delta;
        Sz1[1] = 0.5 * (delta2-delta+0.25);
        Sz1[2] = 0.75-delta2;
        Sz1[3] = 0.5 * (delta2+delta+0.25);

        // ---------------------------
        // Calculate the total charge
        // ---------------------------
        ip -= i_domain_begin + bin +2;
        jp -= j_domain_begin + 2;
        kp -= k_domain_begin + 2;

        for (unsigned int i=0 ; i<5 ; i++) { // i loop
            iloc = (i+ip)*b_dim[2]*b_dim[1];
            for (unsigned int j=0 ; j<5 ; j++) { // j loop
                jloc = (jp+j)*b_dim[2];
                for (unsigned int k=0 ; k<5 ; k++) { // k loop
                    Chi_envelope[iloc+jloc+kp+k] += charge_weight * Sx1[i]*Sy1[j]*Sz1[k];
                } // end k loop
            } // end j loop
        } // end i loop



}
