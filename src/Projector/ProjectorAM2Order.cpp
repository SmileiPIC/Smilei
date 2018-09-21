#include "ProjectorAM2Order.h"

#include <cmath>
#include <iostream>
#include <complex>
#include "dcomplex.h"    
#include "ElectroMagnAM.h"
#include "cField2D.h"
#include "Particles.h"
#include "Tools.h"
#include "Patch.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for ProjectorAM2Order
// ---------------------------------------------------------------------------------------------------------------------
ProjectorAM2Order::ProjectorAM2Order (Params& params, Patch* patch) : ProjectorAM(params, patch)
{
    dl_inv_   = 1.0/params.cell_length[0];
    dl_ov_dt  = params.cell_length[0] / params.timestep;
    dr_inv_   = 1.0/params.cell_length[1];
    dt = params.timestep;
    dr_ov_dt  = params.cell_length[1] / params.timestep;
    Nmode=params.nmodes; 
    one_third = 1.0/3.0;
    dr = params.cell_length[1];
    i_domain_begin = patch->getCellStartingGlobalIndex(0);
    j_domain_begin = patch->getCellStartingGlobalIndex(1);
    n_species = patch->vecSpecies.size();

    nprimr = params.n_space[1] + 2*params.oversize[1] + 1;

    DEBUG("cell_length "<< params.cell_length[0]);

}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for ProjectorAM2Order
// ---------------------------------------------------------------------------------------------------------------------
ProjectorAM2Order::~ProjectorAM2Order()
{
}


// ---------------------------------------------------------------------------------------------------------------------
//! Project local currents for mode=0
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM2Order::operator() (complex<double>* Jl, complex<double>* Jr, complex<double>* Jt, Particles &particles, unsigned int ipart, double invgf, int* iold, double* deltaold)
{   int nparts= particles.size();
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------   int iloc,
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
    double crl_p = charge_weight*dl_ov_dt;
    double crr_p = charge_weight*dr_ov_dt;

    // variable declaration
    double xpn, ypn, rp;
    double delta, delta2;
    // arrays used for the Esirkepov projection method
    double  Sx0[5], Sx1[5], Sy0[5], Sy1[5], DSx[5], DSy[5];
    double  Wx[5][5], Wy[5][5], Wz[5][5], Jx_p[5][5], Jy_p[5][5], Jz_p[5][5];
    for (unsigned int i=0; i<5; i++) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
    }
    Sx0[0] = 0.;
    Sx0[4] = 0.;
    Sy0[0] = 0.;
    Sy0[4] = 0.;


    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------
    
    // locate the particle on the primal grid at former time-step & calculate coeff. S0
    delta = deltaold[0*nparts];
    delta2 = delta*delta;
    Sx0[1] = 0.5 * (delta2-delta+0.25);
    Sx0[2] = 0.75-delta2;
    Sx0[3] = 0.5 * (delta2+delta+0.25);
    
    delta = deltaold[1*nparts];
    delta2 = delta*delta;
    Sy0[1] = 0.5 * (delta2-delta+0.25);
    Sy0[2] = 0.75-delta2;
    Sy0[3] = 0.5 * (delta2+delta+0.25);
    
    
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position(0, ipart) * dl_inv_;
    int ip = round(xpn);
    int ipo = iold[0*nparts];
    int ip_m_ipo = ip-ipo-i_domain_begin;
    delta  = xpn - (double)ip;
    delta2 = delta*delta;
    Sx1[ip_m_ipo+1] = 0.5 * (delta2-delta+0.25);
    Sx1[ip_m_ipo+2] = 0.75-delta2;
    Sx1[ip_m_ipo+3] = 0.5 * (delta2+delta+0.25);
    rp = sqrt (particles.position(1, ipart)*particles.position(1, ipart)+particles.position(2, ipart)*particles.position(2, ipart));
    ypn =  rp * dr_inv_ ;
    double crt_p= charge_weight*(particles.momentum(2,ipart)*particles.position(1,ipart)-particles.momentum(1,ipart)*particles.position(2,ipart))/(rp)*invgf;
    int jp = round(ypn);
    int jpo = iold[1*nparts];
    int jp_m_jpo = jp-jpo-j_domain_begin;
    delta  = ypn - (double)jp;
    delta2 = delta*delta;
    Sy1[jp_m_jpo+1] = 0.5 * (delta2-delta+0.25);
    Sy1[jp_m_jpo+2] = 0.75-delta2;
    Sy1[jp_m_jpo+3] = 0.5 * (delta2+delta+0.25);
    
    for (unsigned int i=0; i < 5; i++) {
        DSx[i] = Sx1[i] - Sx0[i];
        DSy[i] = Sy1[i] - Sy0[i];
    }
    
    // ------------------------------------------------
    // Local current created by the particle
    // calculate using the charge conservation equation
    // ------------------------------------------------
    

    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
                Wx[i][j] = DSx[i] * (Sy0[j] + 0.5*DSy[j]);
                Wy[i][j] = DSy[j] * (Sx0[i] + 0.5*DSx[i]);
		Wz[i][j] = Sx0[i]*Sy0[j] + 0.5*DSx[i]*Sy0[j]+0.5*Sx0[i]*DSy[j]+one_third*DSx[i]*DSy[j];
        }
    }
    
    // ------------------------------------------------
    // Local current created by the particle
    // calculate using the charge conservation equation
    // ------------------------------------------------
    for (unsigned int j=0 ; j<5 ; j++) Jx_p[0][j]= 0.;
    for (unsigned int i=0 ; i<5 ; i++) Jy_p[i][0]= 0.;

    for (unsigned int i=1 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
            Jx_p[i][j]= Jx_p[i-1][j] - crl_p * Wx[i-1][j];
        }
    }
    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=1 ; j<5 ; j++) {
            Jy_p[i][j] = Jy_p[i][j-1] - crr_p * Wy[i][j-1];
        }
    }
    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
            Jz_p[i][j] =   crt_p  * Wz[i][j];
        }
    }
    // ---------------------------
    // Calculate the total current
    // ---------------------------
    
    ipo -= 2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
    // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
    jpo -= 2;
    
    int iloc, jloc, linindex;
        // Jl^(d,p)
        for (unsigned int i=0 ; i<5 ; i++) {
            iloc = i+ipo;
            for (unsigned int j=0 ; j<5 ; j++) {
                jloc = j+jpo;
                linindex = iloc*nprimr+jloc;
                if (jloc+ j_domain_begin == 0){
                    Jl [linindex] += Jx_p[i][j]*6./dr; 
                }
                else {
                    Jl [linindex] += Jx_p[i][j] /abs((jloc+ j_domain_begin)*dr); 
                }
             }
         }//i
    
         // Jr^(p,d)
        for (unsigned int i=0 ; i<5 ; i++) {
            iloc = i+ipo;
            for (unsigned int j=0 ; j<5 ; j++) {
                jloc = j+jpo;
                linindex = iloc*(nprimr+1)+jloc;
                Jr [linindex] += Jy_p[i][j] /abs((jloc+ j_domain_begin-0.5)*dr); 
             }
        }//i

    
        // Jt^(p,p)
        for (unsigned int i=0 ; i<5 ; i++) {
            iloc = i+ipo;
            for (unsigned int j=0 ; j<5 ; j++) {
                jloc = j+jpo;
                linindex = iloc*nprimr+jloc;
                if (jloc+ j_domain_begin != 0){ //Jt_mode_0 = 0 on axis
                    Jt [linindex] += Jz_p[i][j] /abs((jloc+ j_domain_begin)*dr);
                }
            }
        }//i




} // END Project local current densities (Jl, Jr, Jt, sort)

// ---------------------------------------------------------------------------------------------------------------------
//! Project local currents for m>0
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM2Order::operator() (complex<double>* Jl, complex<double>* Jr, complex<double>* Jt, Particles &particles, unsigned int ipart,double invgf, int* iold, double* deltaold, complex<double>* exp_m_theta_old, int imode)
{   
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------   int iloc,
    int nparts= particles.size();
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
    double crl_p = charge_weight*dl_ov_dt;
    double crr_p = charge_weight*dr_ov_dt;

    // variable declaration
    double xpn, ypn;
    double delta, delta2;
    // arrays used for the Esirkepov projection method
    double  Sx0[5], Sx1[5], Sy0[5], Sy1[5], DSx[5], DSy[5];
    complex<double>  Wx[5][5], Wy[5][5], Wz[5][5], Jx_p[5][5], Jy_p[5][5], Jz_p[5][5];
    complex<double> e_delta, e_delta_m1, e_delta_inv, e_theta,e_theta_old, e_bar, e_bar_m1, C_m;
 
     for (unsigned int i=0; i<5; i++) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
    }
    Sx0[0] = 0.;
    Sx0[4] = 0.;
    Sy0[0] = 0.;
    Sy0[4] = 0.;
    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------
    
    // locate the particle on the primal grid at former time-step & calculate coeff. S0
    delta = deltaold[0*nparts];
    delta2 = delta*delta;
    Sx0[1] = 0.5 * (delta2-delta+0.25);
    Sx0[2] = 0.75-delta2;
    Sx0[3] = 0.5 * (delta2+delta+0.25);
    
    delta = deltaold[1*nparts];
    delta2 = delta*delta;
    Sy0[1] = 0.5 * (delta2-delta+0.25);
    Sy0[2] = 0.75-delta2;
    Sy0[3] = 0.5 * (delta2+delta+0.25);
    //calculate exponential coefficients

    double yp = particles.position(1,ipart);
    double zp = particles.position(2,ipart);
    double rp = sqrt (particles.position(1, ipart)*particles.position(1, ipart)+particles.position(2, ipart)*particles.position(2, ipart));
    e_theta = (yp-Icpx*zp)/rp;
    e_theta_old =exp_m_theta_old[0];
    e_delta = 1.;
    e_bar = 1.;
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position(0, ipart) * dl_inv_;
    int ip = round(xpn);
    int ipo = iold[0*nparts];
    int ip_m_ipo = ip-ipo-i_domain_begin;
    delta  = xpn - (double)ip;
    delta2 = delta*delta;
    Sx1[ip_m_ipo+1] = 0.5 * (delta2-delta+0.25);
    Sx1[ip_m_ipo+2] = 0.75-delta2;
    Sx1[ip_m_ipo+3] = 0.5 * (delta2+delta+0.25);
    
    ypn = rp *dr_inv_ ;
    int jp = round(ypn);
    int jpo = iold[1*nparts];
    int jp_m_jpo = jp-jpo-j_domain_begin;
    delta  = ypn - (double)jp;
    delta2 = delta*delta;
    Sy1[jp_m_jpo+1] = 0.5 * (delta2-delta+0.25);
    Sy1[jp_m_jpo+2] = 0.75-delta2;
    Sy1[jp_m_jpo+3] = 0.5 * (delta2+delta+0.25);

    e_delta_m1 = sqrt(e_theta/e_theta_old);
    e_bar_m1 = sqrt(e_theta*e_theta_old);   
    if (std::real(e_theta)+ std::real(e_theta_old) < 0.){
        if (std::imag(e_theta)*std::imag(e_theta_old) > 0.){
            e_bar_m1 *= -1.;
        } else {
            e_delta_m1 *= -1.;
        }
    }

    for (unsigned int i=0; i<imode; i++){
        e_delta *= e_delta_m1;
        e_bar *= e_bar_m1;   
    }

     e_delta_inv =1./e_delta;
    //defining crt_p 
     complex<double> crt_p = - charge_weight*Icpx/(e_bar*dt*(double)imode)*rp;
    for (unsigned int i=0; i < 5; i++) {
        DSx[i] = Sx1[i] - Sx0[i];
        DSy[i] = Sy1[i] - Sy0[i];
    }
    
    // ------------------------------------------------
    // Local current created by the particle
    // calculate using the charge conservation equation
    // ------------------------------------------------
    

    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
            Wx[i][j] = DSx[i] * (Sy0[j] + 0.5*DSy[j]);
            Wy[i][j] = DSy[j] * (Sx0[i] + 0.5*DSx[i]);
            Wz[i][j] = Sy1[j]*Sx1[i]*(e_delta_inv-1.)-Sy0[j]*Sx0[i]*(e_delta-1.);
            
        }
    }
    
    // ------------------------------------------------
    // Local current created by the particle
    // calculate using the charge conservation equation
    // ------------------------------------------------
    for (unsigned int j=0 ; j<5 ; j++) Jx_p[0][j]= 0.;
    for (unsigned int i=0 ; i<5 ; i++) Jy_p[i][0]= 0.;

    for (unsigned int i=1 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
                Jx_p[i][j]= Jx_p[i-1][j] - crl_p * Wx[i-1][j];
            }
        }
    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=1 ; j<5 ; j++) {
                Jy_p[i][j] = Jy_p[i][j-1] - crr_p  * Wy[i][j-1];
            }
        }
    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
                Jz_p[i][j] = crt_p  * Wz[i][j];
        }
    }

    // ---------------------------
    // Calculate the total current
    // ---------------------------
    
    ipo -= 2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
    // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
    jpo -= 2;
    
    int iloc, jloc, linindex;

    C_m = 1.;
    for (unsigned int i=0; i<imode; i++){
    C_m *= e_theta;
    }
    C_m = 1./C_m; 
    // Jl^(d,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*nprimr+jloc;
            if (jloc+ j_domain_begin != 0){ // Jl_mode_1 = 0 on axis
                Jl [linindex] += C_m * Jx_p[i][j] /abs((jloc+ j_domain_begin)*dr); // iloc = (i+ipo)*nprimr;
            }
        }
    }//i
    
    // Jt^(p,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*nprimr+jloc;
            if (jloc+ j_domain_begin == 0){
                Jt [linindex] += Jz_p[i][j]*6./dr;
            }else{
                Jt [linindex] += Jz_p[i][j] /abs((jloc+ j_domain_begin)*dr);
            }
        }
     }
    // Jr^(p,d)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*(nprimr+1)+jloc;
            Jr [linindex] += C_m * Jy_p[i][j] /abs((jloc+ j_domain_begin-0.5)*dr);
        }
    }//i

} // END Project local current densities (Jl, Jr, Jt, sort)



// ---------------------------------------------------------------------------------------------------------------------
//! Project local currents with diag for mode=0 
// ---------------------------------------------------------------------------------------------------------------------

void ProjectorAM2Order::operator() (complex<double>* Jl, complex<double>* Jr, complex<double>* Jt, complex<double>* rho, Particles &particles, unsigned int ipart, double invgf, int* iold, double* deltaold)
{   // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    int nparts= particles.size();
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
    double crl_p = charge_weight*dl_ov_dt;
    double crr_p = charge_weight*dr_ov_dt;
    // variable declaration
    double xpn, ypn, rp;
    double delta, delta2;
    // arrays used for the Esirkepov projection method
    double  Sx0[5], Sx1[5], Sy0[5], Sy1[5], DSx[5], DSy[5];
    double  Wx[5][5], Wy[5][5], Wz[5][5], Jx_p[5][5], Jy_p[5][5], Jz_p[5][5];
    for (unsigned int i=0; i<5; i++) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
    }
    Sx0[0] = 0.;
    Sx0[4] = 0.;
    Sy0[0] = 0.;
    Sy0[4] = 0.;
    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------
    
    // locate the particle on the primal grid at former time-step & calculate coeff. S0
    delta = deltaold[0*nparts];
    delta2 = delta*delta;
    Sx0[1] = 0.5 * (delta2-delta+0.25);
    Sx0[2] = 0.75-delta2;
    Sx0[3] = 0.5 * (delta2+delta+0.25);
    
    delta = deltaold[1*nparts];
    delta2 = delta*delta;
    Sy0[1] = 0.5 * (delta2-delta+0.25);
    Sy0[2] = 0.75-delta2;
    Sy0[3] = 0.5 * (delta2+delta+0.25);
    
    
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position(0, ipart) * dl_inv_;
    int ip = round(xpn);
    int ipo = iold[0*nparts];
    int ip_m_ipo = ip-ipo-i_domain_begin;
    delta  = xpn - (double)ip;
    delta2 = delta*delta;
    Sx1[ip_m_ipo+1] = 0.5 * (delta2-delta+0.25);
    Sx1[ip_m_ipo+2] = 0.75-delta2;
    Sx1[ip_m_ipo+3] = 0.5 * (delta2+delta+0.25);
    rp = sqrt (particles.position(1, ipart)*particles.position(1, ipart)+particles.position(2, ipart)*particles.position(2, ipart));
    ypn = rp * dr_inv_ ;
    double crt_p= charge_weight*(particles.momentum(2,ipart)*particles.position(1,ipart)-particles.momentum(1,ipart)*particles.position(2,ipart))/(rp)*invgf;
    
    int jp = round(ypn);
    int jpo = iold[1*nparts];
    int jp_m_jpo = jp-jpo-j_domain_begin;
    delta  = ypn - (double)jp;
    delta2 = delta*delta;
    Sy1[jp_m_jpo+1] = 0.5 * (delta2-delta+0.25);
    Sy1[jp_m_jpo+2] = 0.75-delta2;
    Sy1[jp_m_jpo+3] = 0.5 * (delta2+delta+0.25);
    
    for (unsigned int i=0; i < 5; i++) {
        DSx[i] = Sx1[i] - Sx0[i];
        DSy[i] = Sy1[i] - Sy0[i];
    }
    
    // ------------------------------------------------
    // Local current created by the particle
    // calculate using the charge conservation equation
    // ------------------------------------------------
    

    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
                Wx[i][j] = DSx[i] * (Sy0[j] + 0.5*DSy[j]);
                Wy[i][j] = DSy[j] * (Sx0[i] + 0.5*DSx[i]);
		Wz[i][j] = Sx0[i]*Sy0[j] + 0.5*DSx[i]*Sy0[j]+0.5*Sx0[i]*DSy[j]+one_third*DSx[i]*DSy[j];
            }
        }
    
    // ------------------------------------------------
    // Local current created by the particle
    // calculate using the charge conservation equation
    // ------------------------------------------------
    for (unsigned int j=0 ; j<5 ; j++) Jx_p[0][j]= 0.;
    for (unsigned int i=0 ; i<5 ; i++) Jy_p[i][0]= 0.;
        

    for (unsigned int i=1 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
                Jx_p[i][j]= Jx_p[i-1][j] - crl_p * Wx[i-1][j];
            }
        }
    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=1 ; j<5 ; j++) {
                Jy_p[i][j] = Jy_p[i][j-1] - crr_p * Wy[i][j-1];
            }
        }
    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
                Jz_p[i][j] =   crt_p  * Wz[i][j];
            }
        }

    // ---------------------------
    // Calculate the total current
    // ---------------------------
    ipo -= 2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
    // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
    jpo -= 2;
    
    int iloc, jloc, linindex;
    
    // Jl^(d,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*nprimr+jloc;
            if (jloc+ j_domain_begin == 0){
                Jl [linindex] += Jx_p[i][j]*6./dr;
            }
            else {
                Jl [linindex] += Jx_p[i][j] /abs((jloc+ j_domain_begin)*dr); // iloc = (i+ipo)*nprimr;
            }
        }
    }//i
    
    // Jr^(p,d)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*(nprimr+1)+jloc;
            Jr [linindex] += Jy_p[i][j] /abs((jloc+ j_domain_begin-0.5)*dr);
         }
    }//i
    
    // Jt^(p,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*nprimr+jloc;
            if (jloc+ j_domain_begin != 0){ // Jt_mode_0 on axis = 0
                Jt [linindex] += Jz_p[i][j] /abs((jloc+ j_domain_begin)*dr);
            }
        }
    }//i

    // Rho^(p,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;    
            linindex = iloc*nprimr+jloc;
             if (jloc+ j_domain_begin == 0){
                rho [linindex] += charge_weight * Sx1[i]*Sy1[j] * 6./dr;
            }
            else {
                rho [linindex] += charge_weight * Sx1[i]*Sy1[j] /abs((jloc+ j_domain_begin)*dr);
            }
        }
    }//i
} // END Project local densities (Jl, Jr, Jt, rho, sort)

// ---------------------------------------------------------------------------------------------------------------------
//! Project local currents with diag for m>0
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM2Order::operator() (complex<double>* Jl, complex<double>* Jr, complex<double>* Jt, complex<double>* rho, Particles &particles, unsigned int ipart, double invgf, int* iold, double* deltaold,complex<double>* exp_m_theta_old,  int imode)
{   
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------   
    int nparts= particles.size();
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
    double crl_p = charge_weight*dl_ov_dt;
    double crr_p = charge_weight*dr_ov_dt;

    // variable declaration
    double xpn, ypn;
    double delta, delta2;
    // arrays used for the Esirkepov projection method
    double  Sx0[5], Sx1[5], Sy0[5], Sy1[5], DSx[5], DSy[5];
    complex<double>  Wx[5][5], Wy[5][5], Wz[5][5], Jx_p[5][5], Jy_p[5][5], Jz_p[5][5];
    complex<double> e_delta,e_delta_m1, e_delta_inv, e_theta,e_theta_old,e_bar,e_bar_m1, C_m;
 
     for (unsigned int i=0; i<5; i++) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
    }
    Sx0[0] = 0.;
    Sx0[4] = 0.;
    Sy0[0] = 0.;
    Sy0[4] = 0.;
    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------
    
    // locate the particle on the primal grid at former time-step & calculate coeff. S0
    delta = deltaold[0*nparts];
    delta2 = delta*delta;
    Sx0[1] = 0.5 * (delta2-delta+0.25);
    Sx0[2] = 0.75-delta2;
    Sx0[3] = 0.5 * (delta2+delta+0.25);
    
    delta = deltaold[1*nparts];
    delta2 = delta*delta;
    Sy0[1] = 0.5 * (delta2-delta+0.25);
    Sy0[2] = 0.75-delta2;
    Sy0[3] = 0.5 * (delta2+delta+0.25);
    //calculate exponential coefficients

    double yp = particles.position(1,ipart);
    double zp = particles.position(2,ipart);
    double rp = sqrt(yp*yp+zp*zp);
    e_theta = (yp-Icpx*zp)/rp;
    e_theta_old = exp_m_theta_old[0];
    e_delta = 1.;
    e_bar =  1.;
     
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position(0, ipart) * dl_inv_;
    int ip = round(xpn);
    int ipo = iold[0*nparts];
    int ip_m_ipo = ip-ipo-i_domain_begin;
    delta  = xpn - (double)ip;
    delta2 = delta*delta;
    Sx1[ip_m_ipo+1] = 0.5 * (delta2-delta+0.25);
    Sx1[ip_m_ipo+2] = 0.75-delta2;
    Sx1[ip_m_ipo+3] = 0.5 * (delta2+delta+0.25);
    
    ypn = rp *dr_inv_ ;
    int jp = round(ypn);
    int jpo = iold[1*nparts];
    int jp_m_jpo = jp-jpo-j_domain_begin;
    delta  = ypn - (double)jp;
    delta2 = delta*delta;
    Sy1[jp_m_jpo+1] = 0.5 * (delta2-delta+0.25);
    Sy1[jp_m_jpo+2] = 0.75-delta2;
    Sy1[jp_m_jpo+3] = 0.5 * (delta2+delta+0.25);
    
    e_delta_m1 = sqrt(e_theta/e_theta_old);
    e_bar_m1 = sqrt(e_theta*e_theta_old);   
    if (std::real(e_theta)+ std::real(e_theta_old) < 0.){
        if (std::imag(e_theta)*std::imag(e_theta_old) > 0.){
            e_bar_m1 *= -1.;
        } else {
            e_delta_m1 *= -1.;
        }
    }

    for (unsigned int i=0; i<imode; i++){
        e_delta *= e_delta_m1;
        e_bar *= e_bar_m1;   
    }
    e_delta_inv =1./e_delta;
    //defining crt_p 
    complex<double> crt_p = -charge_weight*Icpx/(e_bar*dt*(double)imode)*rp;
    for (unsigned int i=0; i < 5; i++) {
        DSx[i] = Sx1[i] - Sx0[i];
        DSy[i] = Sy1[i] - Sy0[i];
    }
    
    // ------------------------------------------------
    // Local current created by the particle
    // calculate using the charge conservation equation
    // ------------------------------------------------
    

    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
                Wx[i][j] = DSx[i] * (Sy0[j] + 0.5*DSy[j]);
                Wy[i][j] = DSy[j] * (Sx0[i] + 0.5*DSx[i]);
		Wz[i][j] = Sy1[j]*Sx1[i]*(e_delta_inv-1.)-Sy0[j]*Sx0[i]*(e_delta-1.);
        }
    }
    
    // ------------------------------------------------
    // Local current created by the particle
    // calculate using the charge conservation equation
    // ------------------------------------------------
    for (unsigned int j=0 ; j<5 ; j++) Jx_p[0][j]= 0.;
    for (unsigned int i=0 ; i<5 ; i++) Jy_p[i][0]= 0.;

    for (unsigned int i=1 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
                Jx_p[i][j]= Jx_p[i-1][j] - crl_p * Wx[i-1][j];
            }
        }
    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=1 ; j<5 ; j++) {
                Jy_p[i][j] = Jy_p[i][j-1] - crr_p * Wy[i][j-1];
            }
        }
    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
                Jz_p[i][j] = crt_p * Wz[i][j];
        }
    }

    // ---------------------------
    // Calculate the total current
    // ---------------------------
    
    ipo -= 2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
    // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
    jpo -= 2;
    
    int iloc, jloc, linindex;

    C_m = 1.;
    for (unsigned int i=0; i<imode; i++){
    C_m *= e_theta;
    }
    C_m= 1./C_m; 
    // Jl^(d,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*nprimr+jloc;
            if (jloc+ j_domain_begin != 0){ //Jl_mode_1 is 0 on axis
                Jl [linindex] += C_m * Jx_p[i][j] /abs((jloc+ j_domain_begin)*dr);
            }
        }
    }//i
    

    
    // Jt^(p,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*nprimr+jloc;
            if (jloc+ j_domain_begin == 0){
                Jt [linindex] += Jz_p[i][j]*6./dr ; 
            }else {
                Jt [linindex] += Jz_p[i][j] / abs((jloc+ j_domain_begin)*dr)  ; 
            }
        }
    }//i

    // Jr^(p,d)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*(nprimr+1)+jloc;
            Jr [linindex] += C_m * Jy_p[i][j] /abs((jloc+ j_domain_begin-0.5)*dr); //
        }
    }//i

     // Rho^(p,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*nprimr+jloc;
            if (jloc+ j_domain_begin != 0){
                rho [linindex] += C_m*charge_weight* Sx1[i]*Sy1[j]/abs((jloc+ j_domain_begin)*dr); // iloc = (i+ipo)*nprimr;
            }
            else {
                rho [linindex] +=  C_m*charge_weight* Sx1[i]*Sy1[j]*6./dr; // iloc = (i+ipo)*nprimr;
            }
        }
    }//i

 
    
} // END Project local current densities (Jl, Jr, Jt, sort)

void ProjectorAM2Order::operator() (double* rhoj, Particles &particles, unsigned int ipart, unsigned int type, std::vector<unsigned int> &b_dim)
{
// Useless function
}



// ---------------------------------------------------------------------------------------------------------------------
//! Project for diags and frozen species - mode >= 0 
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM2Order::operator() (complex<double>* rhoj, Particles &particles, unsigned int ipart, unsigned int type, std::vector<unsigned int> &b_dim, int imode)
{
    //Warning : this function is not charge conserving.

    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------

    int iloc, nr(nprimr);
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
    double r = sqrt (particles.position(1, ipart)*particles.position(1, ipart)+particles.position(2, ipart)*particles.position(2, ipart));

    if (type > 0) { //if current density
        charge_weight *= 1./sqrt(1.0 + particles.momentum(0,ipart)*particles.momentum(0,ipart)
                                     + particles.momentum(1,ipart)*particles.momentum(1,ipart)
                                     + particles.momentum(2,ipart)*particles.momentum(2,ipart));
        if (type == 1){ //if Jl
            charge_weight *= particles.momentum(0,ipart);
        }
        else if (type == 2){ //if Jr
            charge_weight *= (particles.momentum(1,ipart)*particles.position(1,ipart) + particles.momentum(2,ipart)*particles.position(2,ipart))/r ;
            nr++;
        }
        else { //if Jt
            charge_weight *= (-particles.momentum(1,ipart)*particles.position(2,ipart) + particles.momentum(2,ipart)*particles.position(1,ipart))/r ;
        }
    }

    complex<double> e_theta = ( particles.position(1,ipart) + Icpx*particles.position(2,ipart))/r;
    complex<double> C_m = 1.;
    for (unsigned int i=0; i<imode; i++)
        C_m *= e_theta;

    double xpn, ypn;
    double delta, delta2;
    double Sx1[5], Sy1[5]; 

// Initialize all current-related arrays to zero
    for (unsigned int i=0; i<5; i++) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
    }

    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------

    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position(0, ipart) * dl_inv_;
    int ip = round(xpn + 0.5 * (type==1));
    delta  = xpn - (double)ip;
    delta2 = delta*delta;
    Sx1[1] = 0.5 * (delta2-delta+0.25);
    Sx1[2] = 0.75-delta2;
    Sx1[3] = 0.5 * (delta2+delta+0.25);
    ypn = r * dr_inv_ ;
    int jp = round(ypn + 0.5*(type==2));
    delta  = ypn - (double)jp;
    delta2 = delta*delta;
    Sy1[1] = 0.5 * (delta2-delta+0.25);
    Sy1[2] = 0.75-delta2;
    Sy1[3] = 0.5 * (delta2+delta+0.25);

    // ---------------------------
    // Calculate the total charge
    // ---------------------------
    ip -= i_domain_begin + 2;
    jp -= j_domain_begin + 2;

    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = (i+ip)*nr+jp;
        for (unsigned int j=0 ; j<5 ; j++) {
            if ((type != 2) && (j+jp+j_domain_begin == 0)){
                rhoj [iloc+j] += C_m*charge_weight*6.* Sx1[i]*Sy1[j] /dr; 
            }
            else {
                rhoj [iloc+j] += C_m*charge_weight* Sx1[i]*Sy1[j]/abs((j+jp+ j_domain_begin-0.5*(type==2))*dr); 
                }
            }
    }//i
} // END Project for diags local current densities

// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities : ionization NOT DONE YET
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM2Order::operator() (Field* Jl, Field* Jr, Field* Jt, Particles &particles, int ipart, LocalFields Jion)
{  
    cField2D* Jl3D  = static_cast<cField2D*>(Jl);
    cField2D* Jr3D  = static_cast<cField2D*>(Jr);
    cField2D* Jt3D  = static_cast<cField2D*>(Jt);
    
    
    //Declaration of local variables
    int ip, id, jp, jd;
    double xpn, xpmxip, xpmxip2, xpmxid, xpmxid2;
    double ypn, ypmyjp, ypmyjp2, ypmyjd, ypmyjd2;
    double Sxp[3], Sxd[3], Syp[3], Syd[3];
    
    // weighted currents
    double Jx_ion = Jion.x * particles.weight(ipart);
    double Jy_ion = Jion.y * particles.weight(ipart);
    double Jz_ion = Jion.z * particles.weight(ipart);
    
    //Locate particle on the grid
    xpn    = particles.position(0, ipart) * dl_inv_;  // normalized distance to the first node 
    ypn = sqrt (particles.position(1, ipart)*particles.position(1, ipart)+particles.position(2, ipart)*particles.position(2, ipart))*dr_inv_ ;
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
    
    ip  -= i_domain_begin;
    id  -= i_domain_begin;
    jp  -= j_domain_begin;
    jd  -= j_domain_begin;
    
    for (unsigned int i=0 ; i<3 ; i++) {
        int iploc=ip+i-1;
        int idloc=id+i-1;
        for (unsigned int j=0 ; j<3 ; j++) {
            int jploc=jp+j-1;
            int jdloc=jd+j-1;
            if (jploc+ j_domain_begin ==0){
            // Jx^(d,p)
            (*Jl3D)(idloc,jploc) += Jx_ion*8. /dr * Sxd[i]*Syp[j];
            (*Jr3D)(idloc,jploc) += Jy_ion*8. /dr * Sxp[i]*Syd[j];
            (*Jt3D)(idloc,jploc) += Jz_ion*8. /dr * Sxp[i]*Syp[j];  //A corriger dualite et repliement
            } else {
            (*Jl3D)(idloc,jploc) += Jx_ion /((jploc+ j_domain_begin)*dr) * Sxd[i]*Syp[j];
            (*Jr3D)(idloc,jploc) += Jy_ion /((jploc+ j_domain_begin)*dr) * Sxp[i]*Syd[j];
            (*Jt3D)(idloc,jploc) += Jz_ion /((jploc+ j_domain_begin)*dr) * Sxp[i]*Syp[j];
            }

        }
    }//i


} // END Project global current densities (ionize)

//------------------------------------//
//Wrapper for projection
void ProjectorAM2Order::operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread, int ibin, int clrw, bool diag_flag, bool is_spectral, std::vector<unsigned int> &b_dim, int ispec, int ipart_ref)
{  //std::cout<<"projecting"<<std::endl;
   if (is_spectral)
        ERROR("Not implemented");

    std::vector<int> *iold = &(smpi->dynamics_iold[ithread]);
    std::vector<double> *delta = &(smpi->dynamics_deltaold[ithread]);
    std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);   
    std::vector<std::complex<double>> *exp_m_theta_old = &(smpi->dynamics_thetaold[ithread]);

    ElectroMagnAM* emAM = static_cast<ElectroMagnAM*>( EMfields );

    // If no field diagnostics this timestep, then the projection is done directly on the total arrays
    if (!diag_flag){ 

        for ( unsigned int imode = 0; imode<Nmode;imode++){

        complex< double>* b_Jl =  &(*emAM->Jl_[imode] )(0);
        complex<double>* b_Jr =  &(*emAM->Jr_[imode] )(0);
        complex<double>* b_Jt =  &(*emAM->Jt_[imode] )(0);

        if (imode==0){
            for ( int ipart=istart ; ipart<iend; ipart++ )
                (*this)(b_Jl , b_Jr , b_Jt , particles,  ipart, (*invgf)[ipart], &(*iold)[ipart], &(*delta)[ipart]);
        }
        else{	
            for ( int ipart=istart ; ipart<iend; ipart++ )
                (*this)(b_Jl , b_Jr , b_Jt , particles,  ipart,(*invgf)[ipart], &(*iold)[ipart], &(*delta)[ipart],&(*exp_m_theta_old)[ipart], imode);
            } 
        }       // Otherwise, the projection may apply to the species-specific arrays
    } 
    else {
         //Loop on modes 
        for ( unsigned int imode = 0; imode<Nmode;imode++){

            // Fix for n_species which is not know in constructors : now the projector is inside each species
            n_species = emAM->Jl_.size() / Nmode;

	    int ifield = imode*n_species+ispec;
                complex<double>* b_Jl  = emAM->Jl_s [ifield] ? &(* static_cast<cField2D*>(emAM->Jl_s [ifield]) )(0) : &(*emAM->Jl_[imode] )(0) ;
                complex<double>* b_Jr  = emAM->Jr_s [ifield] ? &(* static_cast<cField2D*>(emAM->Jr_s [ifield]) )(0) : &(*emAM->Jr_[imode] )(0) ;
                complex<double>* b_Jt  = emAM->Jt_s [ifield] ? &(* static_cast<cField2D*>(emAM->Jt_s [ifield]) )(0) : &(*emAM->Jt_[imode] )(0) ;
                complex<double>* b_rho = emAM->rho_AM_s[ifield] ? &(* static_cast<cField2D*>(emAM->rho_AM_s[ifield]) )(0) : &(*emAM->rho_AM_[imode])(0) ;
            if (imode==0){
                for ( int ipart=istart ; ipart<iend; ipart++ )
                   (*this)(b_Jl , b_Jr , b_Jt ,b_rho, particles,  ipart, (*invgf)[ipart], &(*iold)[ipart], &(*delta)[ipart]);
             }
             else{    
                for ( int ipart=istart ; ipart<iend; ipart++ )
                    (*this)(b_Jl , b_Jr , b_Jt ,b_rho, particles,  ipart, (*invgf)[ipart], &(*iold)[ipart], &(*delta)[ipart], &(*exp_m_theta_old)[ipart], imode);
                 }

       }
   }
}
