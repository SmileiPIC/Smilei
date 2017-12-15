#include "Projector3D2OrderV.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field3D.h"
#include "Particles.h"
#include "Tools.h"
#include "Patch.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Projector3D2OrderV
// ---------------------------------------------------------------------------------------------------------------------
Projector3D2OrderV::Projector3D2OrderV (Params& params, Patch* patch) : Projector3D(params, patch)
{
    dx_inv_   = 1.0/params.cell_length[0];
    dx_ov_dt  = params.cell_length[0] / params.timestep;
    dy_inv_   = 1.0/params.cell_length[1];
    dy_ov_dt  = params.cell_length[1] / params.timestep;
    dz_inv_   = 1.0/params.cell_length[2];
    dz_ov_dt  = params.cell_length[2] / params.timestep;

    one_third = 1.0/3.0;

    i_domain_begin = patch->getCellStartingGlobalIndex(0);
    j_domain_begin = patch->getCellStartingGlobalIndex(1);
    k_domain_begin = patch->getCellStartingGlobalIndex(2);

    nprimy = params.n_space[1] + 1;
    nprimz = params.n_space[2] + 1;
    oversize[0] = params.oversize[0];
    oversize[1] = params.oversize[1];
    oversize[2] = params.oversize[2];
    dq_inv[0] = dx_inv_;
    dq_inv[1] = dy_inv_;
    dq_inv[2] = dz_inv_;


    DEBUG("cell_length "<< params.cell_length[0]);

}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Projector3D2OrderV
// ---------------------------------------------------------------------------------------------------------------------
Projector3D2OrderV::~Projector3D2OrderV()
{
}

// ---------------------------------------------------------------------------------------------------------------------
//!  Project current densities & charge : diagFields timstep (not vectorized)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D2OrderV::operator() (double* Jx, double* Jy, double* Jz, double* rho, Particles &particles, unsigned int ipart, double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaold)
{
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
    double delta, delta2;
    // arrays used for the Esirkepov projection method
    double Sx0[5], Sx1[5], Sy0[5], Sy1[5], Sz0[5], Sz1[5], DSx[5], DSy[5], DSz[5];
    double tmpJx[5][5], tmpJy[5][5], tmpJz[5][5];
    
    for (unsigned int i=0; i<5; i++) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
        Sz1[i] = 0.;
    }
    
    for (unsigned int j=0; j<5; j++)
        for (unsigned int k=0; k<5; k++)
            tmpJx[j][k] = 0.;
    for (unsigned int i=0; i<5; i++)
        for (unsigned int k=0; k<5; k++)
            tmpJy[i][k] = 0.;
    for (unsigned int i=0; i<5; i++)
        for (unsigned int j=0; j<5; j++)
            tmpJz[i][j] = 0.;

    int nparts = particles.size();
    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------
    
    // locate the particle on the primal grid at former time-step & calculate coeff. S0
    delta = *deltaold;
    delta2 = delta*delta;
    Sx0[0] = 0.;
    Sx0[1] = 0.5 * (delta2-delta+0.25);
    Sx0[2] = 0.75-delta2;
    Sx0[3] = 0.5 * (delta2+delta+0.25);
    Sx0[4] = 0.;
    
    delta = *(deltaold+nparts);
    delta2 = delta*delta;
    Sy0[0] = 0.;
    Sy0[1] = 0.5 * (delta2-delta+0.25);
    Sy0[2] = 0.75-delta2;
    Sy0[3] = 0.5 * (delta2+delta+0.25);
    Sy0[4] = 0.;
    
    delta = *(deltaold+2*nparts);
    delta2 = delta*delta;
    Sz0[0] = 0.;
    Sz0[1] = 0.5 * (delta2-delta+0.25);
    Sz0[2] = 0.75-delta2;
    Sz0[3] = 0.5 * (delta2+delta+0.25);
    Sz0[4] = 0.;
    
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position(0, ipart) * dx_inv_;
    int ip = round(xpn);
    int ipo = iold[0];
    int ip_m_ipo = ip-ipo-i_domain_begin;
    delta  = xpn - (double)ip;
    delta2 = delta*delta;
    Sx1[ip_m_ipo+1] = 0.5 * (delta2-delta+0.25);
    Sx1[ip_m_ipo+2] = 0.75-delta2;
    Sx1[ip_m_ipo+3] = 0.5 * (delta2+delta+0.25);
    
    ypn = particles.position(1, ipart) * dy_inv_;
    int jp = round(ypn);
    int jpo = iold[1];
    int jp_m_jpo = jp-jpo-j_domain_begin;
    delta  = ypn - (double)jp;
    delta2 = delta*delta;
    Sy1[jp_m_jpo+1] = 0.5 * (delta2-delta+0.25);
    Sy1[jp_m_jpo+2] = 0.75-delta2;
    Sy1[jp_m_jpo+3] = 0.5 * (delta2+delta+0.25);
    
    zpn = particles.position(2, ipart) * dz_inv_;
    int kp = round(zpn);
    int kpo = iold[2];
    int kp_m_kpo = kp-kpo-k_domain_begin;
    delta  = zpn - (double)kp;
    delta2 = delta*delta;
    Sz1[kp_m_kpo+1] = 0.5 * (delta2-delta+0.25);
    Sz1[kp_m_kpo+2] = 0.75-delta2;
    Sz1[kp_m_kpo+3] = 0.5 * (delta2+delta+0.25);
    
    // computes Esirkepov coefficients
    for (unsigned int i=0; i < 5; i++) {
        DSx[i] = Sx1[i] - Sx0[i];
        DSy[i] = Sy1[i] - Sy0[i];
        DSz[i] = Sz1[i] - Sz0[i];
    }
        
    // ---------------------------
    // Calculate the total current
    // ---------------------------

    ipo -= bin+2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
                    // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
    jpo -= 2;
    kpo -= 2;
    
    int iloc, jloc, kloc, linindex;
    
    // Jx^(d,p,p)
    for (unsigned int i=1 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            for (unsigned int k=0 ; k<5 ; k++) {
                tmpJx[j][k] -= crx_p * DSx[i-1] * (Sy0[j]*Sz0[k] + 0.5*DSy[j]*Sz0[k] + 0.5*DSz[k]*Sy0[j] + one_third*DSy[j]*DSz[k]);
                kloc = k+kpo;
                linindex = iloc*b_dim[2]*b_dim[1]+jloc*b_dim[2]+kloc;
                Jx [linindex] += tmpJx[j][k]; // iloc = (i+ipo)*b_dim[1];
            }
        }
    }//i
    
    // Jy^(p,d,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=1 ; j<5 ; j++) {
            jloc = j+jpo;
            for (unsigned int k=0 ; k<5 ; k++) {
                tmpJy[i][k] -= cry_p * DSy[j-1] * (Sz0[k]*Sx0[i] + 0.5*DSz[k]*Sx0[i] + 0.5*DSx[i]*Sz0[k] + one_third*DSz[k]*DSx[i]);
                kloc = k+kpo;
                linindex = iloc*b_dim[2]*(b_dim[1]+1)+jloc*b_dim[2]+kloc;
                Jy [linindex] += tmpJy[i][k]; //
            }
        }
    }//i
    
    // Jz^(p,p,d)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            for (unsigned int k=1 ; k<5 ; k++) {
                tmpJz[i][j] -= crz_p * DSz[k-1] * (Sx0[i]*Sy0[j] + 0.5*DSx[i]*Sy0[j] + 0.5*DSy[j]*Sx0[i] + one_third*DSx[i]*DSy[j]);
                kloc = k+kpo;
                linindex = iloc*(b_dim[2]+1)*b_dim[1]+jloc*(b_dim[2]+1)+kloc;
                Jz [linindex] += tmpJz[i][j]; //
            }
        }
    }//i

    // Rho^(p,p,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            for (unsigned int k=0 ; k<5 ; k++) {
                kloc = k+kpo;
                linindex = iloc*b_dim[2]*b_dim[1]+jloc*b_dim[2]+kloc;
                rho[linindex] += charge_weight * Sx1[i]*Sy1[j]*Sz1[k];
            }
        }
    }//i

} // END Project local current densities at dag timestep.


// ---------------------------------------------------------------------------------------------------------------------
//! Project charge : frozen & diagFields timstep (not vectorized)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D2OrderV::operator() (double* rho, Particles &particles, unsigned int ipart, unsigned int bin, std::vector<unsigned int> &b_dim)
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

    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = (i+ip)*b_dim[2]*b_dim[1];
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = (jp+j)*b_dim[2];
            for (unsigned int k=0 ; k<5 ; k++) {
                rho[iloc+jloc+kp+k] += charge_weight * Sx1[i]*Sy1[j]*Sz1[k];
            }
        }
    }//i

} // END Project local current densities frozen.


// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities : ionization
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D2OrderV::operator() (Field* Jx, Field* Jy, Field* Jz, Particles &particles, int ipart, LocalFields Jion)
{
    Field3D* Jx3D  = static_cast<Field3D*>(Jx);
    Field3D* Jy3D  = static_cast<Field3D*>(Jy);
    Field3D* Jz3D  = static_cast<Field3D*>(Jz);
    
    
    //Declaration of local variables
    int ip, id, jp, jd, kp, kd;
    double xpn, xpmxip, xpmxip2, xpmxid, xpmxid2;
    double ypn, ypmyjp, ypmyjp2, ypmyjd, ypmyjd2;
    double zpn, zpmzkp, zpmzkp2, zpmzkd, zpmzkd2;
    double Sxp[3], Sxd[3], Syp[3], Syd[3], Szp[3], Szd[3];
    
    // weighted currents
    double Jx_ion = Jion.x * particles.weight(ipart);
    double Jy_ion = Jion.y * particles.weight(ipart);
    double Jz_ion = Jion.z * particles.weight(ipart);
    
    //Locate particle on the grid
    xpn    = particles.position(0, ipart) * dx_inv_;  // normalized distance to the first node
    ypn    = particles.position(1, ipart) * dy_inv_;  // normalized distance to the first node
    zpn    = particles.position(1, ipart) * dz_inv_;  // normalized distance to the first node
    
    // x-primal index
    ip      = round(xpn);                    // x-index of the central node
    xpmxip  = xpn - (double)ip;              // normalized distance to the nearest grid point
    xpmxip2 = xpmxip*xpmxip;                 // square of the normalized distance to the nearest grid point
    
    // x-dual index
    id      = round(xpn+0.5);                // x-index of the central node
    xpmxid  = xpn - (double)id + 0.5;        // normalized distance to the nearest grid point
    xpmxid2 = xpmxid*xpmxid;                 // square of the normalized distance to the nearest grid point
    
    // y-primal index
    jp      = round(ypn);                    // y-index of the central node
    ypmyjp  = ypn - (double)jp;              // normalized distance to the nearest grid point
    ypmyjp2 = ypmyjp*ypmyjp;                 // square of the normalized distance to the nearest grid point
    
    // y-dual index
    jd      = round(ypn+0.5);                // y-index of the central node
    ypmyjd  = ypn - (double)jd + 0.5;        // normalized distance to the nearest grid point
    ypmyjd2 = ypmyjd*ypmyjd;                 // square of the normalized distance to the nearest grid point
    
    // z-primal index
    kp      = round(zpn);                    // z-index of the central node
    zpmzkp  = zpn - (double)kp;              // normalized distance to the nearest grid point
    zpmzkp2 = zpmzkp*zpmzkp;                 // square of the normalized distance to the nearest grid point
    
    // z-dual index
    kd      = round(zpn+0.5);                // z-index of the central node
    zpmzkd  = zpn - (double)kd + 0.5;        // normalized distance to the nearest grid point
    zpmzkd2 = zpmzkd*zpmzkd;                 // square of the normalized distance to the nearest grid point
    
    Sxp[0] = 0.5 * (xpmxip2-xpmxip+0.25);
    Sxp[1] = (0.75-xpmxip2);
    Sxp[2] = 0.5 * (xpmxip2+xpmxip+0.25);
    
    Sxd[0] = 0.5 * (xpmxid2-xpmxid+0.25);
    Sxd[1] = (0.75-xpmxid2);
    Sxd[2] = 0.5 * (xpmxid2+xpmxid+0.25);
    
    Syp[0] = 0.5 * (ypmyjp2-ypmyjp+0.25);
    Syp[1] = (0.75-ypmyjp2);
    Syp[2] = 0.5 * (ypmyjp2+ypmyjp+0.25);
    
    Syd[0] = 0.5 * (ypmyjd2-ypmyjd+0.25);
    Syd[1] = (0.75-ypmyjd2);
    Syd[2] = 0.5 * (ypmyjd2+ypmyjd+0.25);
    
    Szp[0] = 0.5 * (zpmzkp2-zpmzkp+0.25);
    Szp[1] = (0.75-zpmzkp2);
    Szp[2] = 0.5 * (zpmzkp2+zpmzkp+0.25);
    
    Szd[0] = 0.5 * (zpmzkd2-zpmzkd+0.25);
    Szd[1] = (0.75-zpmzkd2);
    Szd[2] = 0.5 * (zpmzkd2+zpmzkd+0.25);
    
    ip  -= i_domain_begin;
    id  -= i_domain_begin;
    jp  -= j_domain_begin;
    jd  -= j_domain_begin;
    kp  -= k_domain_begin;
    kd  -= k_domain_begin;
    
    for (unsigned int i=0 ; i<3 ; i++) {
        int iploc=ip+i-1;
        int idloc=id+i-1;
        for (unsigned int j=0 ; j<3 ; j++) {
            int jploc=jp+j-1;
            int jdloc=jd+j-1;
            for (unsigned int k=0 ; k<3 ; k++) {
                int kploc=kp+k-1;
                int kdloc=kd+k-1;
                // Jx^(d,p,p)
                (*Jx3D)(idloc,jploc,kploc) += Jx_ion * Sxd[i]*Syp[j]*Szp[k];
                // Jy^(p,d,p)
                (*Jy3D)(iploc,jdloc,kploc) += Jy_ion * Sxp[i]*Syd[j]*Szp[k];
                // Jz^(p,p,d)
                (*Jz3D)(iploc,jploc,kdloc) += Jz_ion * Sxp[i]*Syp[j]*Szd[k];
            }//k
        }//j
    }//i
    
} // END Project global current densities (ionize)


// ---------------------------------------------------------------------------------------------------------------------
//! Project current densities : main projector vectorized
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D2OrderV::operator() (double* Jx, double* Jy, double* Jz, Particles &particles, unsigned int istart, unsigned int iend, std::vector<double> *invgf, std::vector<unsigned int> &b_dim, int* iold, double *deltaold)
{
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------

    int npart_total = particles.size();
    int ipo = iold[0];
    int jpo = iold[1];
    int kpo = iold[2];
    int ipom2 = ipo-2;
    int jpom2 = jpo-2;
    int kpom2 = kpo-2;

    int vecSize = 8;
    int bsize = 5*5*5*vecSize;

    double bJx[bsize] __attribute__((aligned(64)));

    double Sx0_buff_vect[32] __attribute__((aligned(64)));
    double Sy0_buff_vect[32] __attribute__((aligned(64)));
    double Sz0_buff_vect[32] __attribute__((aligned(64)));
    double DSx[40] __attribute__((aligned(64)));
    double DSy[40] __attribute__((aligned(64)));
    double DSz[40] __attribute__((aligned(64)));
    double charge_weight[8] __attribute__((aligned(64)));

    // Closest multiple of 8 higher or equal than npart = iend-istart.
 
    int cell_nparts( (int)iend-(int)istart );

    for (int ivect=0 ; ivect < cell_nparts; ivect += vecSize ){

        int np_computed = min(cell_nparts-ivect,vecSize);

        /*if (np_computed != 8) {
            #pragma omp simd
            for (unsigned int i=0; i<200; i++)
                bJx[i] = 0.;
        }*/

        #pragma omp simd
        for (int ipart=0 ; ipart<np_computed; ipart++ ){

            // locate the particle on the primal grid at former time-step & calculate coeff. S0
            //                            X                                 //
            double delta = deltaold[ivect+ipart+istart];
            double delta2 = delta*delta;

            Sx0_buff_vect[          ipart] = 0.5 * (delta2-delta+0.25);
            Sx0_buff_vect[  vecSize+ipart] = 0.75-delta2;
            Sx0_buff_vect[2*vecSize+ipart] = 0.5 * (delta2+delta+0.25);
            Sx0_buff_vect[3*vecSize+ipart] = 0.;

            //                            Y                                 //
            delta = deltaold[ivect+ipart+istart+npart_total];
            delta2 = delta*delta;

            Sy0_buff_vect[          ipart] = 0.5 * (delta2-delta+0.25);
            Sy0_buff_vect[  vecSize+ipart] = 0.75-delta2;
            Sy0_buff_vect[2*vecSize+ipart] = 0.5 * (delta2+delta+0.25);
            Sy0_buff_vect[3*vecSize+ipart] = 0.;

            //                            Z                                 //
            delta = deltaold[ivect+ipart+istart+2*npart_total];
            delta2 = delta*delta;

            Sz0_buff_vect[          ipart] = 0.5 * (delta2-delta+0.25);
            Sz0_buff_vect[  vecSize+ipart] = 0.75-delta2;
            Sz0_buff_vect[2*vecSize+ipart] = 0.5 * (delta2+delta+0.25);
            Sz0_buff_vect[3*vecSize+ipart] = 0.;


            // locate the particle on the primal grid at current time-step & calculate coeff. S1
            //                            X                                 //
            double pos = particles.position(0, ivect+ipart+istart) * dx_inv_;
            int cell = round(pos);
            int cell_shift = cell-ipo-i_domain_begin;
            delta  = pos - (double)cell;
            delta2 = delta*delta;
            double deltam =  0.5 * (delta2-delta+0.25);
            double deltap =  0.5 * (delta2+delta+0.25);
            delta2 = 0.75 - delta2;
            double m1 = (cell_shift == -1);
            double c0 = (cell_shift ==  0);
            double p1 = (cell_shift ==  1);
            DSx [          ipart] = m1 * deltam                             ;
            DSx [  vecSize+ipart] = c0 * deltam + m1 * delta2               -  Sx0_buff_vect[          ipart];
            DSx [2*vecSize+ipart] = p1 * deltam + c0 * delta2 + m1* deltap  -  Sx0_buff_vect[  vecSize+ipart];
            DSx [3*vecSize+ipart] =               p1 * delta2 + c0* deltap  -  Sx0_buff_vect[2*vecSize+ipart];
            DSx [4*vecSize+ipart] =                             p1* deltap  ;
            //                            Y                                 //
            pos = particles.position(1, ivect+ipart+istart) * dy_inv_;
            cell = round(pos);
            cell_shift = cell-jpo-j_domain_begin;
            delta  = pos - (double)cell;
            delta2 = delta*delta;
            deltam =  0.5 * (delta2-delta+0.25);
            deltap =  0.5 * (delta2+delta+0.25);
	    delta2 = 0.75 - delta2;
            m1 = (cell_shift == -1);
            c0 = (cell_shift ==  0);
            p1 = (cell_shift ==  1);
            DSy [          ipart] = m1 * deltam                            ;
            DSy [  vecSize+ipart] = c0 * deltam + m1 * delta2              -  Sy0_buff_vect[          ipart]                 ;
            DSy [2*vecSize+ipart] = p1 * deltam + c0 * delta2 + m1* deltap -  Sy0_buff_vect[  vecSize+ipart] ;
            DSy [3*vecSize+ipart] =               p1 * delta2 + c0* deltap -  Sy0_buff_vect[2*vecSize+ipart] ;
            DSy [4*vecSize+ipart] =                             p1* deltap  ;
            //                            Z                                 //
            pos = particles.position(2, ivect+ipart+istart) * dz_inv_;
            cell = round(pos);
            cell_shift = cell-kpo-k_domain_begin;
            delta  = pos - (double)cell;
            delta2 = delta*delta;
            deltam =  0.5 * (delta2-delta+0.25);
            deltap =  0.5 * (delta2+delta+0.25);
	    delta2 = 0.75 - delta2;
            m1 = (cell_shift == -1);
            c0 = (cell_shift ==  0);
            p1 = (cell_shift ==  1);
            DSz [          ipart] = m1 * deltam                                                            ;
            DSz [  vecSize+ipart] = c0 * deltam + m1 * delta2              -  Sz0_buff_vect[          ipart]                 ;
            DSz [2*vecSize+ipart] = p1 * deltam + c0 * delta2 + m1* deltap -  Sz0_buff_vect[  vecSize+ipart] ;
            DSz [3*vecSize+ipart] =               p1 * delta2 + c0* deltap -  Sz0_buff_vect[2*vecSize+ipart] ;
            DSz [4*vecSize+ipart] =                             p1* deltap  ;


            charge_weight[ipart] = (double)(particles.charge(ivect+istart+ipart))*particles.weight(ivect+istart+ipart);
        }

            // Jx^(d,p,p)
        #pragma omp simd
        for (int ipart=0 ; ipart<np_computed; ipart++ ){

            //optrpt complains about the following loop but not unrolling it actually seems to give better result.
            double crx_p = charge_weight[ipart]*dx_ov_dt;

            double sum[5];
            sum[0] = 0.;
            for (unsigned int k=1 ; k<5 ; k++) {
                sum[k] = sum[k-1]-DSx[(k-1)*vecSize+ipart];
            }

            double tmp(  crx_p * ( one_third*DSy[ipart]*DSz[ipart] ) );
            for (unsigned int i=1 ; i<5 ; i++) {
                bJx [ ( (i)*25 )*vecSize+ipart] = sum[i]*tmp;
            }

            for (unsigned int k=1 ; k<5 ; k++) {
                double tmp( crx_p * ( 0.5*DSy[ipart]*Sz0_buff_vect[(k-1)*vecSize+ipart] + one_third*DSy[ipart]*DSz[k*vecSize+ipart] ) );
                int index( ( k )*vecSize+ipart );
                for (unsigned int i=1 ; i<5 ; i++) {
                    bJx [ index+25*(i)*vecSize ] = sum[i]*tmp;
               }
                
            }
            for (unsigned int j=1 ; j<5 ; j++) {
                double tmp( crx_p * ( 0.5*DSz[ipart]*Sy0_buff_vect[(j-1)*vecSize+ipart] + one_third*DSy[j*vecSize+ipart]*DSz[ipart] ) );
                int index( ( j*5 )*vecSize+ipart );
                for (unsigned int i=1 ; i<5 ; i++) {
                    bJx [ index+25*(i)*vecSize ] = sum[i]*tmp;
                }
            }//i
            for ( int j=1 ; j<5 ; j++) {
                for ( int k=1 ; k<5 ; k++) {
                    double tmp( crx_p * ( Sy0_buff_vect[(j-1)*vecSize+ipart]*Sz0_buff_vect[(k-1)*vecSize+ipart] 
                                          + 0.5*DSy[j*vecSize+ipart]*Sz0_buff_vect[(k-1)*vecSize+ipart] 
                                          + 0.5*DSz[k*vecSize+ipart]*Sy0_buff_vect[(j-1)*vecSize+ipart] 
                                          + one_third*DSy[j*vecSize+ipart]*DSz[k*vecSize+ipart] ) );
                    int index( ( j*5 + k )*vecSize+ipart );
                    for ( int i=1 ; i<5 ; i++) {
                        bJx [ index+25*(i)*vecSize ] = sum[i]*tmp;
                   }
                }
            }//i


        } // END ipart (compute coeffs)


        int iloc0 = ipom2*b_dim[1]*b_dim[2]+jpom2*b_dim[2]+kpom2;


        int iloc  = iloc0;
        if (np_computed==8) {
            for (unsigned int i=1 ; i<5 ; i++) {
                iloc += b_dim[1]*b_dim[2];
                for (unsigned int j=0 ; j<5 ; j++) {
                    #pragma omp simd
                    for (unsigned int k=0 ; k<5 ; k++) {
                        double tmpJx = 0.;
                        int ilocal = ((i)*25+j*5+k)*vecSize;
                        #pragma unroll(8)
                        for (int ipart=0 ; ipart<8; ipart++ ){
                            tmpJx += bJx [ilocal+ipart];
                        }
                        Jx[iloc+j*b_dim[2]+k]         += tmpJx;
                    }
                }
            }
        }
        else {
            for (unsigned int i=1 ; i<5 ; i++) {
                iloc += b_dim[1]*b_dim[2];
                for (unsigned int j=0 ; j<5 ; j++) {
                    #pragma omp simd
                    for (unsigned int k=0 ; k<5 ; k++) {
                        double tmpJx = 0.;
                        int ilocal = ((i)*25+j*5+k)*vecSize;
                        for (int ipart=0 ; ipart<np_computed; ipart++ ){
                            tmpJx += bJx [ilocal+ipart];
                        }
                        Jx[iloc+j*b_dim[2]+k]         += tmpJx;
                    }
                }
            }
        }


        // Jy^(p,d,p)

        //for (unsigned int j=0; j<1000; j=j+200)
        //    #pragma omp simd
        //    for (unsigned int i=0; i<40; i++)
        //        bJx[j+i] = 0.;

        #pragma omp simd
        for (int ipart=0 ; ipart<np_computed; ipart++ ){
            //optrpt complains about the following loop but not unrolling it actually seems to give better result.
            double cry_p = charge_weight[ipart]*dy_ov_dt;

            double sum[5];
            sum[0] = 0.;
            for (unsigned int k=1 ; k<5 ; k++) {
                sum[k] = sum[k-1]-DSy[(k-1)*vecSize+ipart];
            }

            double tmp( cry_p *one_third*DSz[ipart]*DSx[ipart] );
            for (unsigned int j=1 ; j<5 ; j++) {
                bJx [((j)*5)*vecSize+ipart] = sum[j]*tmp;
            }
            for (unsigned int k=1 ; k<5 ; k++) {
                double tmp( cry_p * (0.5*DSx[0]*Sz0_buff_vect[(k-1)*vecSize+ipart] + one_third*DSz[k*vecSize+ipart]*DSx[ipart]) );
                int index( ( k )*vecSize+ipart );
                for (unsigned int j=1 ; j<5 ; j++) {
                    bJx [ index+5*j*vecSize ] = sum[j]*tmp;
                }
            }
            for (unsigned int i=1 ; i<5 ; i++) {
                double tmp( cry_p * (0.5*DSz[ipart]*Sx0_buff_vect[(i-1)*vecSize+ipart] + one_third*DSz[0]*DSx[i*vecSize+ipart]) );
                int index( ( i*25 )*vecSize+ipart );
                for (unsigned int j=1 ; j<5 ; j++) {
                    bJx [ index+5*j*vecSize ] = sum[j]*tmp;
                }
            }//i
            for (unsigned int i=1 ; i<5 ; i++) {
                for (unsigned int k=1 ; k<5 ; k++) {
                    double tmp( cry_p * (Sz0_buff_vect[(k-1)*vecSize+ipart]*Sx0_buff_vect[(i-1)*vecSize+ipart]
                                         + 0.5*DSz[k*vecSize+ipart]*Sx0_buff_vect[(i-1)*vecSize+ipart]
                                         + 0.5*DSx[i*vecSize+ipart]*Sz0_buff_vect[(k-1)*vecSize+ipart]
                                         + one_third*DSz[k*vecSize+ipart]*DSx[i*vecSize+ipart]) );
                    int index( ( i*25 + k )*vecSize+ipart );
                    for (unsigned int j=1 ; j<5 ; j++) {
                        bJx [ index+5*j*vecSize ] = sum[j]*tmp;
                   }
                }
            }//i
     

        } // END ipart (compute coeffs)


        iloc = iloc0+ipom2*b_dim[2];
        if (np_computed==8) {
            for (unsigned int i=0 ; i<5 ; i++) {
                for (unsigned int j=1 ; j<5 ; j++) {
                    #pragma omp simd
                    for (unsigned int k=0 ; k<5 ; k++) {
                        double tmpJy = 0.;
                        int ilocal = ((i)*25+j*5+k)*vecSize;
                        #pragma unroll(8)
                        for (int ipart=0 ; ipart<8; ipart++ ){
                            tmpJy += bJx [ilocal+ipart];                  
                        }
                        Jy[iloc+j*b_dim[2]+k] += tmpJy;
                    }
                }
                iloc += (b_dim[1]+1)*b_dim[2];
            }
        }
        else {
            for (unsigned int i=0 ; i<5 ; i++) {
                for (unsigned int j=1 ; j<5 ; j++) {
                    #pragma omp simd
                    for (unsigned int k=0 ; k<5 ; k++) {
                        double tmpJy = 0.;
                        int ilocal = ((i)*25+j*5+k)*vecSize;
                        for (int ipart=0 ; ipart<np_computed; ipart++ ){
                            tmpJy += bJx [ilocal+ipart];                  
                        }
                        Jy[iloc+j*b_dim[2]+k] += tmpJy;
                    }
                }
                iloc += (b_dim[1]+1)*b_dim[2];
            }
        }


        // Jz^(p,p,d)
        //for (unsigned int j=0; j<1000; j=j+40)
        //    #pragma omp simd
        //    for (unsigned int i=0; i<8; i++)
        //        bJx[j+i] = 0.;

        #pragma omp simd
        for (int ipart=0 ; ipart<np_computed; ipart++ ){
            //optrpt complains about the following loop but not unrolling it actually seems to give better result.
            double crz_p = charge_weight[ipart]*dz_ov_dt;

            double sum[5];
            sum[0] = 0.;
            for (unsigned int k=1 ; k<5 ; k++) {
                sum[k] = sum[k-1]-DSz[(k-1)*vecSize+ipart];
            }

            double tmp( crz_p *one_third*DSx[ipart]*DSy[ipart] );
            for (unsigned int k=1 ; k<5 ; k++) {
                bJx[( k )*vecSize+ipart] = sum[k]*tmp;
            }
            for (unsigned int j=1 ; j<5 ; j++) {
                double tmp( crz_p * (0.5*DSx[ipart]*Sy0_buff_vect[(j-1)*vecSize+ipart] + one_third*DSx[ipart]*DSy[j*vecSize+ipart]) );
                int index( ( j*5 )*vecSize+ipart );
                for (unsigned int k=1 ; k<5 ; k++) {
                    bJx [ index+k*vecSize ] = sum[k]*tmp;
                }
            }
            for (unsigned int i=1 ; i<5 ; i++) {
                double tmp( crz_p * (0.5*DSy[ipart]*Sx0_buff_vect[(i-1)*vecSize+ipart] + one_third*DSx[i*vecSize+ipart]*DSy[ipart]) );
                int index( ( i*25 )*vecSize+ipart );
                for (unsigned int k=1 ; k<5 ; k++) {
                    bJx [ index+k*vecSize ] = sum[k]*tmp;
                }
            }//i
            for (unsigned int i=1 ; i<5 ; i++) {
                for (unsigned int j=1 ; j<5 ; j++) {
                    double tmp( crz_p * (Sx0_buff_vect[(i-1)*vecSize+ipart]*Sy0_buff_vect[(j-1)*vecSize+ipart]
                                         + 0.5*DSx[i*vecSize+ipart]*Sy0_buff_vect[(j-1)*vecSize+ipart]
                                         + 0.5*DSy[j*vecSize+ipart]*Sx0_buff_vect[(i-1)*vecSize+ipart]
                                         + one_third*DSx[i*vecSize+ipart]*DSy[j*vecSize+ipart]) );
                    int index( ( i*25 + j*5 )*vecSize+ipart );
                    for (unsigned int k=1 ; k<5 ; k++) {
                        bJx [ index+k*vecSize ] = sum[k]*tmp;
                    }
                }
            }//i


        } // END ipart (compute coeffs)


        iloc = iloc0  + jpom2 +ipom2*b_dim[1];
        if (np_computed==8) {
            for (unsigned int i=0 ; i<5 ; i++) {
                for (unsigned int j=0 ; j<5 ; j++) {
                    #pragma omp simd
                    for (unsigned int k=1 ; k<5 ; k++) {
                        double tmpJz = 0.;
                        int ilocal = ((i)*25+j*5+k)*vecSize;
                        #pragma unroll(8)
                        for (int ipart=0 ; ipart<8; ipart++ ){
                            tmpJz +=  bJx[ilocal+ipart];
                        }
                        Jz [iloc + (j)*(b_dim[2]+1) + k] +=  tmpJz;
                    }
                }
                iloc += b_dim[1]*(b_dim[2]+1);
            }
        }
        else {
            for (unsigned int i=0 ; i<5 ; i++) {
                for (unsigned int j=0 ; j<5 ; j++) {
                    #pragma omp simd
                    for (unsigned int k=1 ; k<5 ; k++) {
                        double tmpJz = 0.;
                        int ilocal = ((i)*25+j*5+k)*vecSize;
                        for (int ipart=0 ; ipart<np_computed; ipart++ ){
                            tmpJz +=  bJx[ilocal+ipart];
                        }
                        Jz [iloc + (j)*(b_dim[2]+1) + k] +=  tmpJz;
                    }
                }
                iloc += b_dim[1]*(b_dim[2]+1);
            }
        }//i
        
    } // ivect
    
} // END Project vectorized


// ---------------------------------------------------------------------------------------------------------------------
//! Wrapper for projection
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D2OrderV::operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread, int scell, int clrw, bool diag_flag, std::vector<unsigned int> &b_dim, int ispec)
{
    if ( istart == iend ) return; //Don't treat empty cells.

    //Independent of cell. Should not be here
    //{   
    std::vector<double> *delta = &(smpi->dynamics_deltaold[ithread]);
    std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);
    //}
    int iold[3];


    iold[0] = scell/(nprimy*nprimz)+oversize[0];

    iold[1] = ( (scell%(nprimy*nprimz)) / nprimz )+oversize[1];
    iold[2] = ( (scell%(nprimy*nprimz)) % nprimz )+oversize[2];
   
    
    // If no field diagnostics this timestep, then the projection is done directly on the total arrays
    if (!diag_flag){ 
        double* b_Jx =  &(*EMfields->Jx_ )(0);
        double* b_Jy =  &(*EMfields->Jy_ )(0);
        double* b_Jz =  &(*EMfields->Jz_ )(0);
        (*this)(b_Jx , b_Jy , b_Jz , particles,  istart, iend, invgf, b_dim, iold, &(*delta)[0]);
            
    // Otherwise, the projection may apply to the species-specific arrays
    } else {
        int ibin = 0; // Trick to make it compatible for the moment.
        int dim1 = EMfields->dimPrim[1];
        int dim2 = EMfields->dimPrim[2];
        double* b_Jx  = EMfields->Jx_s [ispec] ? &(*EMfields->Jx_s [ispec])(ibin*clrw* dim1   *dim2) : &(*EMfields->Jx_ )(ibin*clrw* dim1   *dim2) ;
        double* b_Jy  = EMfields->Jy_s [ispec] ? &(*EMfields->Jy_s [ispec])(ibin*clrw*(dim1+1)*dim2) : &(*EMfields->Jy_ )(ibin*clrw*(dim1+1)*dim2) ;
        double* b_Jz  = EMfields->Jz_s [ispec] ? &(*EMfields->Jz_s [ispec])(ibin*clrw*dim1*(dim2+1)) : &(*EMfields->Jz_ )(ibin*clrw*dim1*(dim2+1)) ;
        double* b_rho = EMfields->rho_s[ispec] ? &(*EMfields->rho_s[ispec])(ibin*clrw* dim1   *dim2) : &(*EMfields->rho_)(ibin*clrw* dim1   *dim2) ;
        for (int ipart=istart ; ipart<iend; ipart++ ) {
            //Do not use cells sorting for now : f(ipart) for now, f(istart) laterfor now,
            //(*iold)[ipart       ] = round( particles.position(0, ipart)* dx_inv_ - dt*particles.momentum(0, ipart)*(*invgf)[ipart] * dx_inv_ ) - i_domain_begin ;
            //(*iold)[ipart+nparts] = round( particles.position(1, ipart)* dy_inv_ - dt*particles.momentum(1, ipart)*(*invgf)[ipart] * dy_inv_ ) - j_domain_begin ;
            (*this)(b_Jx , b_Jy , b_Jz ,b_rho, particles,  ipart, (*invgf)[ipart], ibin*clrw, b_dim, iold, &(*delta)[ipart]);
        }
    }
}
