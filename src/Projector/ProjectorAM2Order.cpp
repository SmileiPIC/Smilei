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
    dt = params.timestep;
    dr = params.cell_length[1];
    dl_inv_   = 1.0/params.cell_length[0];
    dl_ov_dt  = params.cell_length[0] / params.timestep;
    dr_inv_   = 1.0 / dr;
    one_ov_dt  = 1.0 / params.timestep;
    Nmode=params.nmodes; 
    one_third = 1.0/3.0;
    i_domain_begin = patch->getCellStartingGlobalIndex(0);
    j_domain_begin = patch->getCellStartingGlobalIndex(1);
    n_species = patch->vecSpecies.size();

    nprimr = params.n_space[1] + 2*params.oversize[1] + 1;

    rprim.resize(nprimr);
    invV.resize(nprimr);
    invVd.resize(nprimr+1);

    for (int j = 0; j< nprimr; j++){
        rprim[j] = abs((j_domain_begin+j)*dr);
        if (j_domain_begin+j == 0){
            //invV[j] = 6./dr; // Correction de Verboncoeur ordre 1.
            invV[j] = 8./dr;   // No correction.
        } else {
            invV[j] = 1./rprim[j];
        }
    }
    for (int j = 0; j< nprimr+1; j++)
        invVd[j] = 1./abs(j_domain_begin+j-0.5);
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
void ProjectorAM2Order::currents_mode0(complex<double>* Jl, complex<double>* Jr, complex<double>* Jt, Particles &particles, unsigned int ipart, double invgf, int* iold, double* deltaold)
{   int nparts= particles.size();
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------   int iloc,
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
    double crl_p = charge_weight*dl_ov_dt;
    double crr_p = charge_weight*one_ov_dt;

    // variable declaration
    double xpn, ypn, rp;
    double delta, delta2;
    // arrays used for the Esirkepov projection method
    double  Sl0[5], Sl1[5], Sr0[5], Sr1[5], DSl[5], DSr[5];
    double  Wl[5][5], Wr[5][5], Wt[5][5], Jl_p[5][5], Jr_p[5][5], Jt_p[5][5];
    for (unsigned int i=0; i<5; i++) {
        Sl1[i] = 0.;
        Sr1[i] = 0.;
    }
    Sl0[0] = 0.;
    Sl0[4] = 0.;
    Sr0[0] = 0.;
    Sr0[4] = 0.;


    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------
    
    // locate the particle on the primal grid at former time-step & calculate coeff. S0
    delta = deltaold[0*nparts];
    delta2 = delta*delta;
    Sl0[1] = 0.5 * (delta2-delta+0.25);
    Sl0[2] = 0.75-delta2;
    Sl0[3] = 0.5 * (delta2+delta+0.25);
    
    delta = deltaold[1*nparts];
    delta2 = delta*delta;
    Sr0[1] = 0.5 * (delta2-delta+0.25);
    Sr0[2] = 0.75-delta2;
    Sr0[3] = 0.5 * (delta2+delta+0.25);
    
    
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position(0, ipart) * dl_inv_;
    int ip = round(xpn);
    int ipo = iold[0*nparts];
    int ip_m_ipo = ip-ipo-i_domain_begin;
    delta  = xpn - (double)ip;
    delta2 = delta*delta;
    Sl1[ip_m_ipo+1] = 0.5 * (delta2-delta+0.25);
    Sl1[ip_m_ipo+2] = 0.75-delta2;
    Sl1[ip_m_ipo+3] = 0.5 * (delta2+delta+0.25);
    rp = sqrt (particles.position(1, ipart)*particles.position(1, ipart)+particles.position(2, ipart)*particles.position(2, ipart));
    ypn =  rp * dr_inv_ ;
    double crt_p= charge_weight*(particles.momentum(2,ipart)*particles.position(1,ipart)-particles.momentum(1,ipart)*particles.position(2,ipart))/(rp)*invgf;
    int jp = round(ypn);
    int jpo = iold[1*nparts];
    int jp_m_jpo = jp-jpo-j_domain_begin;
    delta  = ypn - (double)jp;
    delta2 = delta*delta;
    Sr1[jp_m_jpo+1] = 0.5 * (delta2-delta+0.25);
    Sr1[jp_m_jpo+2] = 0.75-delta2;
    Sr1[jp_m_jpo+3] = 0.5 * (delta2+delta+0.25);
    
    for (unsigned int i=0; i < 5; i++) {
        DSl[i] = Sl1[i] - Sl0[i];
        DSr[i] = Sr1[i] - Sr0[i];
    }
    
    // ------------------------------------------------
    // Local current created by the particle
    // calculate using the charge conservation equation
    // ------------------------------------------------
    

    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
                Wl[i][j] = DSl[i] * (Sr0[j] + 0.5*DSr[j]);
                Wr[i][j] = DSr[j] * (Sl0[i] + 0.5*DSl[i]);
		//Wt[i][j] = Sl0[i]*Sr0[j] + 0.5*DSl[i]*Sr0[j]+0.5*Sl0[i]*DSr[j]+one_third*DSl[i]*DSr[j];
		Wt[i][j] = 0.5 * (Sl0[i]*Sr0[j] + Sl1[i]*Sr1[j]) ;
        }
    }
    
    ipo -= 2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
    // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
    jpo -= 2;
    int iloc, jloc, linindex;

    // ------------------------------------------------
    // Local current created by the particle
    // calculate using the charge conservation equation
    // ------------------------------------------------
    for (unsigned int j=0 ; j<5 ; j++) Jl_p[0][j]= 0.;
    for (unsigned int i=0 ; i<5 ; i++) Jr_p[i][0]= 0.;

    for (unsigned int i=1 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
            Jl_p[i][j]= Jl_p[i-1][j] - crl_p * Wl[i-1][j];
        }
    }
    //for (unsigned int j=1 ; j<5 ; j++) {
    //    jloc = j+jpo;
    //    double Vd = abs(jloc + j_domain_begin - 1.5) ;
    //    for (unsigned int i=0 ; i<5 ; i++) {
    //        Jr_p[i][j] = (Jr_p[i][j-1] * Vd - crr_p * Wr[i][j-1]) * invVd[jloc] ;
    //    }
    //}

    for (int j=3 ; j>=0 ; j--) {
        jloc = j+jpo+1;
        double Vd = abs(jloc + j_domain_begin + 0.5) ;
        for (unsigned int i=0 ; i<5 ; i++) {
            Jr_p[i][j] = (Jr_p[i][j+1] * Vd + crr_p * Wr[i][j+1]) * invVd[jloc];
        }
    }

    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
            Jt_p[i][j] =   crt_p  * Wt[i][j];
        }
    }
    // ---------------------------
    // Calculate the total current
    // ---------------------------
    
    
        // Jl^(d,p)
        for (unsigned int i=1 ; i<5 ; i++) {
            iloc = i+ipo;
            for (unsigned int j=0 ; j<5 ; j++) {
                jloc = j+jpo;
                linindex = iloc*nprimr+jloc;
                Jl [linindex] += Jl_p[i][j]*invV[jloc]; 
             }
         }//i
    
         // Jr^(p,d)
        for (unsigned int i=0 ; i<5 ; i++) {
            iloc = i+ipo;
            for (unsigned int j=0 ; j<4 ; j++) {
                jloc = j+jpo+1;
                linindex = iloc*(nprimr+1)+jloc;
                Jr [linindex] += Jr_p[i][j]; 
             }
        }//i
    
        // Jt^(p,p)
        for (unsigned int i=0 ; i<5 ; i++) {
            iloc = i+ipo;
            for (unsigned int j=0 ; j<5 ; j++) {
                jloc = j+jpo;
                linindex = iloc*nprimr+jloc;
                Jt [linindex] += Jt_p[i][j] *invV[jloc];
            }
        }//i

} // END Project local current densities (Jl, Jr, Jt, sort)

// ---------------------------------------------------------------------------------------------------------------------
//! Project local currents for m>0
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM2Order::currents(complex<double>* Jl, complex<double>* Jr, complex<double>* Jt, Particles &particles, unsigned int ipart,double invgf, int* iold, double* deltaold, complex<double>* exp_m_theta_old, int imode)
{   
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------   int iloc,
    int nparts= particles.size();
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
    double crl_p = charge_weight*dl_ov_dt;
    double crr_p = charge_weight*one_ov_dt;

    // variable declaration
    double xpn, ypn;
    double delta, delta2;
    // arrays used for the Esirkepov projection method
    double  Sl0[5], Sl1[5], Sr0[5], Sr1[5], DSl[5], DSr[5];
    complex<double>  Wl[5][5], Wr[5][5], Wt[5][5], Jl_p[5][5], Jr_p[5][5], Jt_p[5][5];
    complex<double> e_delta, e_delta_m1, e_delta_inv, e_theta,e_theta_old, e_bar, e_bar_m1, C_m, C_m_old;
 
     for (unsigned int i=0; i<5; i++) {
        Sl1[i] = 0.;
        Sr1[i] = 0.;
    }
    Sl0[0] = 0.;
    Sl0[4] = 0.;
    Sr0[0] = 0.;
    Sr0[4] = 0.;
    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------
    
    // locate the particle on the primal grid at former time-step & calculate coeff. S0
    delta = deltaold[0*nparts];
    delta2 = delta*delta;
    Sl0[1] = 0.5 * (delta2-delta+0.25);
    Sl0[2] = 0.75-delta2;
    Sl0[3] = 0.5 * (delta2+delta+0.25);
    
    delta = deltaold[1*nparts];
    delta2 = delta*delta;
    Sr0[1] = 0.5 * (delta2-delta+0.25);
    Sr0[2] = 0.75-delta2;
    Sr0[3] = 0.5 * (delta2+delta+0.25);
    //calculate exponential coefficients

    double yp = particles.position(1,ipart);
    double zp = particles.position(2,ipart);
    double rp = sqrt (particles.position(1, ipart)*particles.position(1, ipart)+particles.position(2, ipart)*particles.position(2, ipart));
    e_theta = (yp + Icpx*zp)/rp;     //exp(i theta)
    e_theta_old =exp_m_theta_old[0]; //exp(-i theta_old)
    double theta = atan2(zp,yp);
    double theta_old =atan2(-std::imag(e_theta_old), std::real(e_theta_old));
    e_delta = 1.;
    e_bar = 1.;
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position(0, ipart) * dl_inv_;
    int ip = round(xpn);
    int ipo = iold[0*nparts];
    int ip_m_ipo = ip-ipo-i_domain_begin;
    delta  = xpn - (double)ip;
    delta2 = delta*delta;
    Sl1[ip_m_ipo+1] = 0.5 * (delta2-delta+0.25);
    Sl1[ip_m_ipo+2] = 0.75-delta2;
    Sl1[ip_m_ipo+3] = 0.5 * (delta2+delta+0.25);
    
    ypn = rp *dr_inv_ ;
    int jp = round(ypn);
    int jpo = iold[1*nparts];
    int jp_m_jpo = jp-jpo-j_domain_begin;
    delta  = ypn - (double)jp;
    delta2 = delta*delta;
    Sr1[jp_m_jpo+1] = 0.5 * (delta2-delta+0.25);
    Sr1[jp_m_jpo+2] = 0.75-delta2;
    Sr1[jp_m_jpo+3] = 0.5 * (delta2+delta+0.25);

    //e_delta_m1 = sqrt(e_theta/e_theta_old);
    //e_bar_m1 = sqrt(e_theta*e_theta_old);   
    //if (std::real(e_theta)+ std::real(e_theta_old) < 0.){
    //    if (std::imag(e_theta)*std::imag(e_theta_old) > 0.){
    //        e_bar_m1 *= -1.;
    //    } else {
    //        e_delta_m1 *= -1.;
    //    }
    //}

    double dtheta = std::remainder(theta-theta_old, 2*M_PI)/2.; // Otherwise dtheta is overestimated when going from -pi to +pi
    double theta_bar = theta_old+dtheta;
    e_delta_m1 = std::polar(1.0,dtheta);
    e_bar_m1 = std::polar(1.0,theta_bar);

    for (unsigned int i=0; i<(unsigned int)imode; i++){
        e_delta *= e_delta_m1;
        e_bar *= e_bar_m1;   
    }

     e_delta_inv =1./e_delta;
    //defining crt_p 
    //complex<double> crt_p = - charge_weight*Icpx/(e_bar*dt*(double)imode)*2.;
    complex<double> crt_p = charge_weight*Icpx*e_bar / (dt*(double)imode)*2.;
    for (unsigned int i=0; i < 5; i++) {
        DSl[i] = Sl1[i] - Sl0[i];
        DSr[i] = Sr1[i] - Sr0[i];
    }
    
    // ------------------------------------------------
    // Local current created by the particle
    // calculate using the charge conservation equation
    // ------------------------------------------------
    

    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
            Wl[i][j] = DSl[i] * (Sr0[j] + 0.5*DSr[j]);
            Wr[i][j] = DSr[j] * (Sl0[i] + 0.5*DSl[i]);
            Wt[i][j] = Sr1[j]*Sl1[i]*(e_delta_inv-1.)-Sr0[j]*Sl0[i]*(e_delta-1.);
            
        }
    }
    
    // ------------------------------------------------
    // Local current created by the particle
    // calculate using the charge conservation equation
    // ------------------------------------------------
    for (unsigned int j=0 ; j<5 ; j++) Jl_p[0][j]= 0.;
    for (unsigned int i=0 ; i<5 ; i++) Jr_p[i][0]= 0.;

    ipo -= 2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
    // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
    jpo -= 2;
    
    int iloc, jloc, linindex;

    for (unsigned int i=1 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
                Jl_p[i][j]= Jl_p[i-1][j] - crl_p * Wl[i-1][j];
            }
        }
    //for (unsigned int j=1 ; j<5 ; j++) {
    //    jloc = j+jpo;
    //    double Vd = abs(jloc + j_domain_begin - 1.5) ;
    //    for (unsigned int i=0 ; i<5 ; i++) {
    //        Jr_p[i][j] = (Jr_p[i][j-1] * Vd - crr_p * Wr[i][j-1]) * invVd[jloc] ;
    //    }
    //}

    for (int j=3 ; j>=0 ; j--) {
        jloc = j+jpo+1;
        double Vd = abs(jloc + j_domain_begin + 0.5) ;
        for (unsigned int i=0 ; i<5 ; i++) {
            Jr_p[i][j] = (Jr_p[i][j+1] * Vd + crr_p * Wr[i][j+1]) * invVd[jloc];
        }
    }

    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
                Jt_p[i][j] = crt_p  * Wt[i][j];
        }
    }

    // ---------------------------
    // Calculate the total current
    // ---------------------------
    

    C_m = 1.;
    C_m_old = 1.;
    for (unsigned int i=0; i<(unsigned int)imode; i++){
    C_m *= e_theta;
    C_m_old *= e_theta_old;
    }
    C_m = 2. * (C_m + 1./C_m_old)/2. ; //multiply modes > 0 by 2 AND exp(i theta_medium) = ( exp(i theta) + exp(i theta_old) ) /2.
 
    // Jl^(d,p)
    for (unsigned int i=1 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*nprimr+jloc;
            Jl [linindex] += C_m * Jl_p[i][j] * invV[jloc];
        }
    }//i
    
    // Jt^(p,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*nprimr+jloc;
            Jt [linindex] += Jt_p[i][j] * invV[jloc] * rprim[jloc] ;
        }
     }
    // Jr^(p,d)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<4 ; j++) {
            jloc = j+jpo+1;
            linindex = iloc*(nprimr+1)+jloc;
            Jr [linindex] += C_m * Jr_p[i][j] ;
        }
    }//i

} // END Project local current densities (Jl, Jr, Jt, sort)



// ---------------------------------------------------------------------------------------------------------------------
//! Project local currents with diag for mode=0 
// ---------------------------------------------------------------------------------------------------------------------

void ProjectorAM2Order::currentsAndDensity_mode0(complex<double>* Jl, complex<double>* Jr, complex<double>* Jt, complex<double>* rho, Particles &particles, unsigned int ipart, double invgf, int* iold, double* deltaold)
{   // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    int nparts= particles.size();
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
    double crl_p = charge_weight*dl_ov_dt;
    double crr_p = charge_weight*one_ov_dt;
    // variable declaration
    double xpn, ypn, rp;
    double delta, delta2;
    // arrays used for the Esirkepov projection method
    double  Sl0[5], Sl1[5], Sr0[5], Sr1[5], DSl[5], DSr[5];
    double  Wl[5][5], Wr[5][5], Wt[5][5], Jl_p[5][5], Jr_p[5][5], Jt_p[5][5];
    for (unsigned int i=0; i<5; i++) {
        Sl1[i] = 0.;
        Sr1[i] = 0.;
    }
    Sl0[0] = 0.;
    Sl0[4] = 0.;
    Sr0[0] = 0.;
    Sr0[4] = 0.;
    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------
    
    // locate the particle on the primal grid at former time-step & calculate coeff. S0
    delta = deltaold[0*nparts];
    delta2 = delta*delta;
    Sl0[1] = 0.5 * (delta2-delta+0.25);
    Sl0[2] = 0.75-delta2;
    Sl0[3] = 0.5 * (delta2+delta+0.25);
    
    delta = deltaold[1*nparts];
    delta2 = delta*delta;
    Sr0[1] = 0.5 * (delta2-delta+0.25);
    Sr0[2] = 0.75-delta2;
    Sr0[3] = 0.5 * (delta2+delta+0.25);
    
    
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position(0, ipart) * dl_inv_;
    int ip = round(xpn);
    int ipo = iold[0*nparts];
    int ip_m_ipo = ip-ipo-i_domain_begin;
    delta  = xpn - (double)ip;
    delta2 = delta*delta;
    Sl1[ip_m_ipo+1] = 0.5 * (delta2-delta+0.25);
    Sl1[ip_m_ipo+2] = 0.75-delta2;
    Sl1[ip_m_ipo+3] = 0.5 * (delta2+delta+0.25);
    rp = sqrt (particles.position(1, ipart)*particles.position(1, ipart)+particles.position(2, ipart)*particles.position(2, ipart));
    ypn = rp * dr_inv_ ;
    double crt_p= charge_weight*(particles.momentum(2,ipart)*particles.position(1,ipart)-particles.momentum(1,ipart)*particles.position(2,ipart))/(rp)*invgf;
    
    int jp = round(ypn);
    int jpo = iold[1*nparts];
    int jp_m_jpo = jp-jpo-j_domain_begin;
    delta  = ypn - (double)jp;
    delta2 = delta*delta;
    Sr1[jp_m_jpo+1] = 0.5 * (delta2-delta+0.25);
    Sr1[jp_m_jpo+2] = 0.75-delta2;
    Sr1[jp_m_jpo+3] = 0.5 * (delta2+delta+0.25);
    
    for (unsigned int i=0; i < 5; i++) {
        DSl[i] = Sl1[i] - Sl0[i];
        DSr[i] = Sr1[i] - Sr0[i];
    }
    
    // ------------------------------------------------
    // Local current created by the particle
    // calculate using the charge conservation equation
    // ------------------------------------------------
    

    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
                Wl[i][j] = DSl[i] * (Sr0[j] + 0.5*DSr[j]);
                Wr[i][j] = DSr[j] * (Sl0[i] + 0.5*DSl[i]);
		//Wt[i][j] = Sl0[i]*Sr0[j] + 0.5*DSl[i]*Sr0[j]+0.5*Sl0[i]*DSr[j]+one_third*DSl[i]*DSr[j];
		Wt[i][j] = 0.5 * (Sl0[i]*Sr0[j] + Sl1[i]*Sr1[j]) ;
            }
        }
    
    // ------------------------------------------------
    // Local current created by the particle
    // calculate using the charge conservation equation
    // ------------------------------------------------
    for (unsigned int j=0 ; j<5 ; j++) Jl_p[0][j]= 0.;
    for (unsigned int i=0 ; i<5 ; i++) Jr_p[i][0]= 0.;

    ipo -= 2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
    // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
    jpo -= 2;
    
    int iloc, jloc, linindex;
        

    for (unsigned int i=1 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
                Jl_p[i][j]= Jl_p[i-1][j] - crl_p * Wl[i-1][j];
            }
        }
    //for (unsigned int j=1 ; j<5 ; j++) {
    //    jloc = j+jpo;
    //    double Vd = abs(jloc + j_domain_begin - 1.5) ;
    //    for (unsigned int i=0 ; i<5 ; i++) {
    //        Jr_p[i][j] = (Jr_p[i][j-1] * Vd - crr_p * Wr[i][j-1]) * invVd[jloc] ;
    //    }
    //}

    for (int j=3 ; j>=0 ; j--) {
        jloc = j+jpo+1;
        double Vd = abs(jloc + j_domain_begin + 0.5) ;
        for (unsigned int i=0 ; i<5 ; i++) {
            Jr_p[i][j] = (Jr_p[i][j+1] * Vd + crr_p * Wr[i][j+1]) * invVd[jloc];
        }
    }

    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
                Jt_p[i][j] =   crt_p  * Wt[i][j];
            }
        }

    // ---------------------------
    // Calculate the total current
    // ---------------------------
    
    // Jl^(d,p)
    for (unsigned int i=1 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*nprimr+jloc;
            Jl [linindex] += Jl_p[i][j] * invV[jloc];
        }
    }//i
    
    // Jr^(p,d)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<4 ; j++) {
            jloc = j+jpo+1;
            linindex = iloc*(nprimr+1)+jloc;
            Jr [linindex] += Jr_p[i][j] ;
         }
    }//i
    
    // Jt^(p,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*nprimr+jloc;
            Jt [linindex] += Jt_p[i][j] * invV[jloc];
        }
    }//i

    // Rho^(p,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;    
            linindex = iloc*nprimr+jloc;
            rho [linindex] += charge_weight * Sl1[i]*Sr1[j] * invV[jloc];
        }
    }//i
} // END Project local densities (Jl, Jr, Jt, rho, sort)

// ---------------------------------------------------------------------------------------------------------------------
//! Project local currents with diag for m>0
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM2Order::currentsAndDensity(complex<double>* Jl, complex<double>* Jr, complex<double>* Jt, complex<double>* rho, Particles &particles, unsigned int ipart, double invgf, int* iold, double* deltaold,complex<double>* exp_m_theta_old,  int imode)
{   
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------   
    int nparts= particles.size();
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
    double crl_p = charge_weight*dl_ov_dt;
    double crr_p = charge_weight*one_ov_dt;

    // variable declaration
    double xpn, ypn;
    double delta, delta2;
    // arrays used for the Esirkepov projection method
    double  Sl0[5], Sl1[5], Sr0[5], Sr1[5], DSl[5], DSr[5];
    complex<double>  Wl[5][5], Wr[5][5], Wt[5][5], Jl_p[5][5], Jr_p[5][5], Jt_p[5][5];
    complex<double> e_delta,e_delta_m1, e_delta_inv, e_theta,e_theta_old,e_bar,e_bar_m1, C_m, C_m_old;
 
     for (unsigned int i=0; i<5; i++) {
        Sl1[i] = 0.;
        Sr1[i] = 0.;
    }
    Sl0[0] = 0.;
    Sl0[4] = 0.;
    Sr0[0] = 0.;
    Sr0[4] = 0.;
    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------
    
    // locate the particle on the primal grid at former time-step & calculate coeff. S0
    delta = deltaold[0*nparts];
    delta2 = delta*delta;
    Sl0[1] = 0.5 * (delta2-delta+0.25);
    Sl0[2] = 0.75-delta2;
    Sl0[3] = 0.5 * (delta2+delta+0.25);
    
    delta = deltaold[1*nparts];
    delta2 = delta*delta;
    Sr0[1] = 0.5 * (delta2-delta+0.25);
    Sr0[2] = 0.75-delta2;
    Sr0[3] = 0.5 * (delta2+delta+0.25);
    //calculate exponential coefficients

    double yp = particles.position(1,ipart);
    double zp = particles.position(2,ipart);
    double rp = sqrt(yp*yp+zp*zp);
    e_theta = (yp + Icpx*zp)/rp;      //exp(i theta)
    e_theta_old = exp_m_theta_old[0]; //exp(-i theta_old)
    e_delta = 1.;
    e_bar =  1.;
    double theta = atan2(zp,yp);
    double theta_old =atan2(-std::imag(e_theta_old), std::real(e_theta_old));
     
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position(0, ipart) * dl_inv_;
    int ip = round(xpn);
    int ipo = iold[0*nparts];
    int ip_m_ipo = ip-ipo-i_domain_begin;
    delta  = xpn - (double)ip;
    delta2 = delta*delta;
    Sl1[ip_m_ipo+1] = 0.5 * (delta2-delta+0.25);
    Sl1[ip_m_ipo+2] = 0.75-delta2;
    Sl1[ip_m_ipo+3] = 0.5 * (delta2+delta+0.25);
    
    ypn = rp *dr_inv_ ;
    int jp = round(ypn);
    int jpo = iold[1*nparts];
    int jp_m_jpo = jp-jpo-j_domain_begin;
    delta  = ypn - (double)jp;
    delta2 = delta*delta;
    Sr1[jp_m_jpo+1] = 0.5 * (delta2-delta+0.25);
    Sr1[jp_m_jpo+2] = 0.75-delta2;
    Sr1[jp_m_jpo+3] = 0.5 * (delta2+delta+0.25);

    double dtheta = std::remainder(theta-theta_old, 2*M_PI)/2.; // Otherwise dtheta is overestimated when going from -pi to +pi
    double theta_bar = theta_old+dtheta ;
    e_delta_m1 = std::polar(1.0,dtheta);
    e_bar_m1 = std::polar(1.0,theta_bar);

    for (unsigned int i=0; i<(unsigned int)imode; i++){
        e_delta *= e_delta_m1;
        e_bar *= e_bar_m1;   
    }
    e_delta_inv =1./e_delta;
    //defining crt_p 
    //complex<double> crt_p = -charge_weight*Icpx/(e_bar*dt*(double)imode)*2.;
    complex<double> crt_p = charge_weight*Icpx*e_bar / (dt*(double)imode)*2.;
    for (unsigned int i=0; i < 5; i++) {
        DSl[i] = Sl1[i] - Sl0[i];
        DSr[i] = Sr1[i] - Sr0[i];
    }
    
    // ------------------------------------------------
    // Local current created by the particle
    // calculate using the charge conservation equation
    // ------------------------------------------------
    

    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
                Wl[i][j] = DSl[i] * (Sr0[j] + 0.5*DSr[j]);
                Wr[i][j] = DSr[j] * (Sl0[i] + 0.5*DSl[i]);
		Wt[i][j] = Sr1[j]*Sl1[i]*(e_delta_inv-1.)-Sr0[j]*Sl0[i]*(e_delta-1.);
        }
    }
    
    // ------------------------------------------------
    // Local current created by the particle
    // calculate using the charge conservation equation
    // ------------------------------------------------
    for (unsigned int j=0 ; j<5 ; j++) Jl_p[0][j]= 0.;
    for (unsigned int i=0 ; i<5 ; i++) Jr_p[i][0]= 0.;

    ipo -= 2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
    // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
    jpo -= 2;
    
    int iloc, jloc, linindex;

    for (unsigned int i=1 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
                Jl_p[i][j]= Jl_p[i-1][j] - crl_p * Wl[i-1][j];
        }
    }
    //for (unsigned int j=1 ; j<5 ; j++) {
    //    jloc = j+jpo;
    //    double Vd = abs(jloc + j_domain_begin - 1.5) ;
    //    for (unsigned int i=0 ; i<5 ; i++) {
    //        Jr_p[i][j] = (Jr_p[i][j-1] * Vd - crr_p * Wr[i][j-1]) * invVd[jloc] ;
    //    }
    //}

    for (int j=3 ; j>=0 ; j--) {
        jloc = j+jpo+1;
        double Vd = abs(jloc + j_domain_begin + 0.5) ;
        for (unsigned int i=0 ; i<5 ; i++) {
            Jr_p[i][j] = (Jr_p[i][j+1] * Vd + crr_p * Wr[i][j+1]) * invVd[jloc];
        }
    }

    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
                Jt_p[i][j] = crt_p * Wt[i][j];
        }
    }

    // ---------------------------
    // Calculate the total current
    // ---------------------------
    

    C_m = 1.;
    C_m_old = 1.;
    for (unsigned int i=0; i<(unsigned int)imode; i++){
    C_m *= e_theta;
    C_m_old *= e_theta_old;
    }
    C_m = 2. * (C_m + 1./C_m_old)/2. ; //multiply modes > 0 by 2 AND exp(i theta_medium) = ( exp(i theta) + exp(i theta_old) ) /2.
    // Jl^(d,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*nprimr+jloc;
            Jl [linindex] += C_m * Jl_p[i][j] * invV[jloc] ;
        }
    }//i
    
    // Jt^(p,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*nprimr+jloc;
            Jt [linindex] += Jt_p[i][j] * invV[jloc] * rprim[jloc] ; 
        }
    }//i

    // Jr^(p,d)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<4 ; j++) {
            jloc = j+jpo+1;
            linindex = iloc*(nprimr+1)+jloc;
            Jr [linindex] += C_m * Jr_p[i][j]; 
        }
    }//i

     // Rho^(p,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*nprimr+jloc;
            rho [linindex] += C_m*charge_weight* Sl1[i]*Sr1[j] * invV[jloc]; 
        }
    }//i

 
    
} // END Project local current densities (Jl, Jr, Jt, sort)

void ProjectorAM2Order::densityFrozen(double* rhoj, Particles &particles, unsigned int ipart, unsigned int type, std::vector<unsigned int> &b_dim)
{
// Useless function
}



// ---------------------------------------------------------------------------------------------------------------------
//! Project for diags and frozen species - mode >= 0 
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM2Order::densityFrozenComplex(complex<double>* rhoj, Particles &particles, unsigned int ipart, unsigned int type, std::vector<unsigned int> &b_dim, int imode)
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
            charge_weight *= (particles.momentum(1,ipart)*particles.position(1,ipart) + particles.momentum(2,ipart)*particles.position(2,ipart)) * dr_inv_ / r ;
            nr++;
        }
        else { //if Jt
            charge_weight *= (-particles.momentum(1,ipart)*particles.position(2,ipart) + particles.momentum(2,ipart)*particles.position(1,ipart))/r ;
        }
    }

    complex<double> e_theta = ( particles.position(1,ipart) + Icpx*particles.position(2,ipart))/r;
    complex<double> C_m = 1.;
    if (imode > 0) C_m = 2.;
    for (unsigned int i=0; i<(unsigned int)imode; i++)
        C_m *= e_theta;

    double xpn, ypn;
    double delta, delta2;
    double Sl1[5], Sr1[5]; 

// Initialize all current-related arrays to zero
    for (unsigned int i=0; i<5; i++) {
        Sl1[i] = 0.;
        Sr1[i] = 0.;
    }

    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------

    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position(0, ipart) * dl_inv_;
    int ip = round(xpn + 0.5 * (type==1));
    delta  = xpn - (double)ip;
    delta2 = delta*delta;
    Sl1[1] = 0.5 * (delta2-delta+0.25);
    Sl1[2] = 0.75-delta2;
    Sl1[3] = 0.5 * (delta2+delta+0.25);
    ypn = r * dr_inv_ ;
    int jp = round(ypn + 0.5*(type==2));
    delta  = ypn - (double)jp;
    delta2 = delta*delta;
    Sr1[1] = 0.5 * (delta2-delta+0.25);
    Sr1[2] = 0.75-delta2;
    Sr1[3] = 0.5 * (delta2+delta+0.25);

    // ---------------------------
    // Calculate the total charge
    // ---------------------------
    ip -= i_domain_begin + 2;
    jp -= j_domain_begin + 2;

    if (type != 2){
        for (unsigned int i=0 ; i<5 ; i++) {
            iloc = (i+ip)*nr+jp;
            for (unsigned int j=0 ; j<5 ; j++)
                rhoj [iloc+j] += C_m*charge_weight* Sl1[i]*Sr1[j] * invV[j+jp]; 
        }//i
    } else {
        for (unsigned int i=0 ; i<5 ; i++) {
            iloc = (i+ip)*nr+jp;
            for (unsigned int j=0 ; j<5 ; j++)
                rhoj [iloc+j] += C_m*charge_weight* Sl1[i]*Sr1[j] * invVd[j+jp]; 
        }//i
    }
} // END Project for diags local current densities

// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities : ionization NOT DONE YET
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM2Order::ionizationCurrents(Field* Jl, Field* Jr, Field* Jt, Particles &particles, int ipart, LocalFields Jion)
{  
    cField2D* JlAM  = static_cast<cField2D*>(Jl);
    cField2D* JrAM  = static_cast<cField2D*>(Jr);
    cField2D* JtAM  = static_cast<cField2D*>(Jt);
    
    
    //Declaration of local variables
    int ip, id, jp, jd;
    double xpn, xpmxip, xpmxip2, xpmxid, xpmxid2;
    double ypn, ypmyjp, ypmyjp2, ypmyjd, ypmyjd2;
    double Slp[3], Sld[3], Srp[3], Srd[3];
    
    // weighted currents
    double Jl_ion = Jion.x * particles.weight(ipart);
    double Jr_ion = Jion.y * particles.weight(ipart);
    double Jt_ion = Jion.z * particles.weight(ipart);
    
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
    
    Slp[0] = 0.5 * (xpmxip2-xpmxip+0.25);
    Slp[1] = (0.75-xpmxip2);
    Slp[2] = 0.5 * (xpmxip2+xpmxip+0.25);
    
    Sld[0] = 0.5 * (xpmxid2-xpmxid+0.25);
    Sld[1] = (0.75-xpmxid2);
    Sld[2] = 0.5 * (xpmxid2+xpmxid+0.25);
    
    Srp[0] = 0.5 * (ypmyjp2-ypmyjp+0.25);
    Srp[1] = (0.75-ypmyjp2);
    Srp[2] = 0.5 * (ypmyjp2+ypmyjp+0.25);
    
    Srd[0] = 0.5 * (ypmyjd2-ypmyjd+0.25);
    Srd[1] = (0.75-ypmyjd2);
    Srd[2] = 0.5 * (ypmyjd2+ypmyjd+0.25);
    
    ip  -= i_domain_begin;
    id  -= i_domain_begin;
    jp  -= j_domain_begin;
    jd  -= j_domain_begin;
    
    for (unsigned int i=0 ; i<3 ; i++) {
        //int iploc=ip+i-1;
        int idloc=id+i-1;
        for (unsigned int j=0 ; j<3 ; j++) {
            int jploc=jp+j-1;
            //int jdloc=jd+j-1;
            if (jploc+ j_domain_begin ==0){
            // Jl^(d,p)
            (*JlAM)(idloc,jploc) += Jl_ion*8. /dr * Sld[i]*Srp[j];
            (*JrAM)(idloc,jploc) += Jr_ion*8. /dr * Slp[i]*Srd[j];
            (*JtAM)(idloc,jploc) += Jt_ion*8. /dr * Slp[i]*Srp[j];  //A corriger dualite et repliement
            } else {
            (*JlAM)(idloc,jploc) += Jl_ion /((jploc+ j_domain_begin)*dr) * Sld[i]*Srp[j];
            (*JrAM)(idloc,jploc) += Jr_ion /((jploc+ j_domain_begin)*dr) * Slp[i]*Srd[j];
            (*JtAM)(idloc,jploc) += Jt_ion /((jploc+ j_domain_begin)*dr) * Slp[i]*Srp[j];
            }

        }
    }//i


} // END Project global current densities (ionize)

//------------------------------------//
//Wrapper for projection
void ProjectorAM2Order::currentsAndDensityWrapper(ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread, int ibin, int clrw, bool diag_flag, bool is_spectral, std::vector<unsigned int> &b_dim, int ispec, int ipart_ref)
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

        complex<double>* b_Jl =  &(*emAM->Jl_[imode] )(0);
        complex<double>* b_Jr =  &(*emAM->Jr_[imode] )(0);
        complex<double>* b_Jt =  &(*emAM->Jt_[imode] )(0);

        if (imode==0){
            for ( int ipart=istart ; ipart<iend; ipart++ ){
                currents_mode0(b_Jl , b_Jr , b_Jt , particles,  ipart, (*invgf)[ipart], &(*iold)[ipart], &(*delta)[ipart]);
            }
        }
        else{	
            for ( int ipart=istart ; ipart<iend; ipart++ )
                currents(b_Jl , b_Jr , b_Jt , particles,  ipart,(*invgf)[ipart], &(*iold)[ipart], &(*delta)[ipart],&(*exp_m_theta_old)[ipart], imode);
            } 
        }       // Otherwise, the projection may apply to the species-specific arrays
    } 
    else {
         //Loop on modes 
        for ( unsigned int imode = 0; imode<Nmode;imode++){

            // Fix for n_species which is not know in constructors : now the projector is inside each species
            n_species = emAM->Jl_.size() / Nmode;

	    int ifield = imode*n_species+ispec;
                complex<double>* b_Jl  = emAM->Jl_s    [ifield] ? &(* (emAM->Jl_s    [ifield]) )(0) : &(*emAM->Jl_    [imode] )(0) ;
                complex<double>* b_Jr  = emAM->Jr_s    [ifield] ? &(* (emAM->Jr_s    [ifield]) )(0) : &(*emAM->Jr_    [imode] )(0) ;
                complex<double>* b_Jt  = emAM->Jt_s    [ifield] ? &(* (emAM->Jt_s    [ifield]) )(0) : &(*emAM->Jt_    [imode] )(0) ;
                complex<double>* b_rho = emAM->rho_AM_s[ifield] ? &(* (emAM->rho_AM_s[ifield]) )(0) : &(*emAM->rho_AM_[imode] )(0) ;
            if (imode==0){
                for ( int ipart=istart ; ipart<iend; ipart++ )
                   currentsAndDensity_mode0(b_Jl , b_Jr , b_Jt ,b_rho, particles,  ipart, (*invgf)[ipart], &(*iold)[ipart], &(*delta)[ipart]);
             }
             else{    
                for ( int ipart=istart ; ipart<iend; ipart++ )
                    currentsAndDensity(b_Jl , b_Jr , b_Jt ,b_rho, particles,  ipart, (*invgf)[ipart], &(*iold)[ipart], &(*delta)[ipart], &(*exp_m_theta_old)[ipart], imode);
                 }

       }
   }
}
