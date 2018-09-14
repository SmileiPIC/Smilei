#include "ProjectorRZ2Order.h"

#include <cmath>
#include <iostream>
#include <complex>
#include "dcomplex.h"    
#include "ElectroMagn3DRZ.h"
#include "cField2D.h"
#include "Particles.h"
#include "Tools.h"
#include "Patch.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for ProjectorRZ2Order
// ---------------------------------------------------------------------------------------------------------------------
ProjectorRZ2Order::ProjectorRZ2Order (Params& params, Patch* patch) : ProjectorRZ(params, patch)
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
    //n_r_max = params.n_space_global[1];
    DEBUG("cell_length "<< params.cell_length[0]);

}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for ProjectorRZ2Order
// ---------------------------------------------------------------------------------------------------------------------
ProjectorRZ2Order::~ProjectorRZ2Order()
{
}


// ---------------------------------------------------------------------------------------------------------------------
//! Project local currents for mode=0
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorRZ2Order::operator() (complex<double>* Jl, complex<double>* Jr, complex<double>* Jt, Particles &particles, unsigned int ipart, double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaold)
{   int nparts= particles.size();
    //std::cout<<"particle momentum in x in proj"<<particles.momentum(0,ipart)<<std::endl;
    //std::cout<<"particle momentum in y in proj"<<particles.momentum(1,ipart)<<std::endl;
    //std::cout<<"particle momentum in z in proj"<<particles.momentum(2,ipart)<<std::endl; 
    //std::cout<<"particle position in x in proj"<<particles.position(0,ipart)<<std::endl;
    //std::cout<<"particle position in r in proj"<<particles.position(1,ipart)<<std::endl;
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
    
    ipo -= bin+2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
    // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
    jpo -= 2;
    
    int iloc, jloc, linindex;
        // Jl^(d,p)
        for (unsigned int i=0 ; i<5 ; i++) {
            iloc = i+ipo;
            for (unsigned int j=0 ; j<5 ; j++) {
                jloc = j+jpo;
                linindex = iloc*b_dim[1]+jloc;
                if (jloc+ j_domain_begin == 0){
                    Jl [linindex] += Jx_p[i][j]*6./dr; // iloc = (i+ipo)*b_dim[1];
                    //Jl [linindex] =0.;
                }
                else {
                    Jl [linindex] += Jx_p[i][j] /abs((jloc+ j_domain_begin)*dr); // iloc = (i+ipo)*b_dim[1];
                    //Jl [linindex] =0;
                    }
             }
         }//i
    
         // Jr^(p,d)
        for (unsigned int i=0 ; i<5 ; i++) {
            iloc = i+ipo;
            for (unsigned int j=0 ; j<5 ; j++) {
                jloc = j+jpo;
                linindex = iloc*(b_dim[1]+1)+jloc;
                Jr [linindex] += Jy_p[i][j] /abs((jloc+ j_domain_begin-0.5)*dr); //
             }
        }//i

    
        // Jt^(p,p)
        for (unsigned int i=0 ; i<5 ; i++) {
            iloc = i+ipo;
            for (unsigned int j=0 ; j<5 ; j++) {
                jloc = j+jpo;
                linindex = iloc*b_dim[1]+jloc;
                if (jloc+ j_domain_begin == 0){
                    Jt [linindex] = 0.; // iloc = (i+ipo)*b_dim[1];
                }
                else {
                    Jt [linindex] += Jz_p[i][j] /abs((jloc+ j_domain_begin)*dr); // iloc = (i+ipo)*b_dim[1];
                    }
            }
        }//i




} // END Project local current densities (Jl, Jr, Jt, sort)

// ---------------------------------------------------------------------------------------------------------------------
//! Project local currents for m>0
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorRZ2Order::operator() (complex<double>* Jl, complex<double>* Jr, complex<double>* Jt, Particles &particles, unsigned int ipart,double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaold, complex<double>* exp_m_theta_old, int imode)
{   //std::cout<<"particle momentum in x in proj"<<particles.momentum(0,ipart)<<std::endl;
    //std::cout<<"particle momentum in y in proj"<<particles.momentum(1,ipart)<<std::endl;
    //std::cout<<"particle momentum in z in proj"<<particles.momentum(2,ipart)<<std::endl; 
    //if (particles.position(1,ipart)<0.){
    //    std::cout<<"particle position in r in proj"<<particles.position(1,ipart)<<std::endl;
    //    }

    //std::cout<<"particle position in x in proj"<<particles.position(0,ipart)<<std::endl;
    //std::cout<<"particle position in r in proj"<<particles.position(1,ipart)<<std::endl;
    int nparts= particles.size();
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------   int iloc,
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
    double crl_p = charge_weight*dl_ov_dt;
    double crr_p = charge_weight*dr_ov_dt;
   // double crt_p= charge_weight*Icpx*;

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
    //cout << std::setprecision(9) << "y= " << yp << endl;
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
        //cout << std::setprecision(9) <<  " e_theta = " << e_theta << " e_theta_old = " << e_theta_old << " e_delta = " << e_delta << " e_bar= " << e_bar << endl;
    }

     e_delta_inv =1./e_delta;
    //defining crt_p 
     complex<double> crt_p = - charge_weight*Icpx/(e_bar*dt*imode)*rp;
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
                // ?? not sure about this ?
                Jz_p[i][j] = crt_p  * Wz[i][j];
                //cout << std::setprecision(9) << " crt_p = " << crt_p << " Wz = " << Wz[i][j] << " Jz = " << Jz_p[i][j] << endl; 
            }
        }

    // ---------------------------
    // Calculate the total current
    // ---------------------------
    
    ipo -= bin+2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
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
            linindex = iloc*b_dim[1]+jloc;
            if (jloc+ j_domain_begin == 0){
                Jl [linindex] +=  C_m * Jx_p[i][j]*6./dr; // iloc = (i+ipo)*b_dim[1];
                }
            else {
                Jl [linindex] += C_m * Jx_p[i][j] /abs((jloc+ j_domain_begin)*dr); // iloc = (i+ipo)*b_dim[1];
                }
           
        }
    }//i
    
    // Jt^(p,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*b_dim[1]+jloc;
            if (jloc+ j_domain_begin == 0){
                Jt [linindex] += Jz_p[i][j]*6./dr;
            }else{
                Jt [linindex] += Jz_p[i][j] /abs((jloc+ j_domain_begin)*dr);
            }
            //if (jloc+j_domain_begin==0){
            //    Jt[linindex] = - 1./3.* (4.* Icpx * Jr[linindex+1] + Jt[linindex+1]);
            //}
        }
     }
    // Jr^(p,d)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*(b_dim[1]+1)+jloc;
            Jr [linindex] += C_m * Jy_p[i][j] /abs((jloc+ j_domain_begin-0.5)*dr);
            //if (jloc+j_domain_begin ==0){
            //    Jr[linindex]+= 2.*Icpx* Jt[linindex] - Jr[linindex+1];
            //}
        }
    }//i

} // END Project local current densities (Jl, Jr, Jt, sort)



// ---------------------------------------------------------------------------------------------------------------------
//! Project local currents with diag for mode=0 
// ---------------------------------------------------------------------------------------------------------------------

void ProjectorRZ2Order::operator() (complex<double>* Jl, complex<double>* Jr, complex<double>* Jt, complex<double>* rho, Particles &particles, unsigned int ipart, double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaold)
{   //std::cout<<"particle momentum in x in proj"<<particles.momentum(0,ipart)<<std::endl;
    //std::cout<<"particle momentum in y in proj"<<particles.momentum(1,ipart)<<std::endl;
    //std::cout<<"particle momentum in z in proj"<<particles.momentum(2,ipart)<<std::endl; 
    //std::cout<<"particle position in x in proj"<<particles.position(0,ipart)<<std::endl;
    //std::cout<<"particle position in r in proj"<<particles.position(1,ipart)<<std::endl;
    int nparts= particles.size();
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
    ipo -= bin+2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
    // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
    jpo -= 2;
    
    int iloc, jloc, linindex;
    
    // Jl^(d,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*b_dim[1]+jloc;
            if (jloc+ j_domain_begin == 0){
                Jl [linindex] += Jx_p[i][j]*6./dr;
                //Jl [linindex] =0.;
            }
            else {
                //Jl [linindex]=0.;
                Jl [linindex] += Jx_p[i][j] /abs((jloc+ j_domain_begin)*dr); // iloc = (i+ipo)*b_dim[1];
                }
            }
    }//i
    
    // Jr^(p,d)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*(b_dim[1]+1)+jloc;
            Jr [linindex] += Jy_p[i][j] /abs((jloc+ j_domain_begin-0.5)*dr); //
         }
    }//i
    
    // Jt^(p,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*b_dim[1]+jloc;
            if (jloc+ j_domain_begin == 0){
                Jt [linindex] += Jz_p[i][j]*6./dr; // iloc = (i+ipo)*b_dim[1];
                //Jt [linindex] =0.;
                }
            else {
                Jt [linindex] += Jz_p[i][j] /abs((jloc+ j_domain_begin)*dr); // iloc = (i+ipo)*b_dim[1];
                }
            }
    }//i

    // Rho^(p,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;    
            linindex = iloc*b_dim[1]+jloc;
         if (jloc+ j_domain_begin == 0){
                rho [linindex] += charge_weight * Sx1[i]*Sy1[j] * 6./dr;
                }
            else {
                rho [linindex] += charge_weight * Sx1[i]*Sy1[j] /abs((jloc+ j_domain_begin)*dr);
                //rho [linindex+1] += rho [linindex-1];
                //rho [linindex+2] += rho [linindex-2];
                 }
        }
    }//i
} // END Project local densities (Jl, Jr, Jt, rho, sort)

// ---------------------------------------------------------------------------------------------------------------------
//! Project local currents with diag for m>0
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorRZ2Order::operator() (complex<double>* Jl, complex<double>* Jr, complex<double>* Jt, complex<double>* rho, Particles &particles, unsigned int ipart, double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaold,complex<double>* exp_m_theta_old,  int imode)
{   //std::cout<<"particle momentum in x in proj"<<particles.momentum(0,ipart)<<std::endl;
    //std::cout<<"particle momentum in y in proj"<<particles.momentum(1,ipart)<<std::endl;
    //std::cout<<"particle momentum in z in proj"<<particles.momentum(2,ipart)<<std::endl;
    //std::cout<<"particle position in x in proj"<<particles.position(0,ipart)<<std::endl;
    //std::cout<<"particle position in r in proj"<<particles.position(1,ipart)<<std::endl;
    int nparts= particles.size();
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------   int iloc,
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
    double crl_p = charge_weight*dl_ov_dt;
    double crr_p = charge_weight*dr_ov_dt;
   // double crt_p= charge_weight*Icpx*;

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
    complex<double> crt_p = -charge_weight*Icpx/(e_bar*dt*imode)*rp;
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
               // ?? same not sure
                Jz_p[i][j] = crt_p * Wz[i][j];
            }
        }

    // ---------------------------
    // Calculate the total current
    // ---------------------------
    
    ipo -= bin+2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
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
            linindex = iloc*b_dim[1]+jloc;
            if (jloc+ j_domain_begin == 0){
                Jl [linindex] +=  C_m * Jx_p[i][j]*6. /dr ; // iloc = (i+ipo)*b_dim[1];
                }
            else {
                Jl [linindex] += C_m * Jx_p[i][j] /abs((jloc+ j_domain_begin)*dr); // iloc = (i+ipo)*b_dim[1];
                }
            }
    }//i
    

    
    // Jt^(p,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*b_dim[1]+jloc;
            if (jloc+ j_domain_begin == 0){
                Jt [linindex] += Jz_p[i][j]*6./dr ; 
            }else {
                Jt [linindex] += Jz_p[i][j] / abs((jloc+ j_domain_begin)*dr)  ; 
            }
            //if (jloc+j_domain_begin==0){
            //    Jt[linindex] += - 1./3.* (4.* Icpx * Jr[linindex+1] + Jt[linindex+1]);
            //}
       // if (jloc+j_domain_begin == 0){
       //     Jt [linindex+1] += Jt [linindex-1];
       //     Jt [linindex+2] += Jt [linindex-2];
       // }
        //if (jloc+j_domain_begin == 0){
        //    Jt [linindex+1] += Jt [linindex-1];
        //    Jt [linindex+2] += Jt [linindex-2];
        //}
            }
    }//i

    // Jr^(p,d)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*(b_dim[1]+1)+jloc;
            Jr [linindex] += C_m * Jy_p[i][j] /abs((jloc+ j_domain_begin-0.5)*dr); //
            //if (jloc+j_domain_begin ==0){
            //    Jr[linindex]+= 2.*Icpx* Jt[linindex] - Jr[linindex+1];
            //}
        //if (jloc+j_domain_begin == 0){
        //    Jr [linindex+1] += Jr [linindex];
        //    Jr [linindex+2] += Jr [linindex-1];
        //    Jr [linindex+3] += Jr [linindex-2];
        //}
            }
    }//i

     // Rho^(p,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*b_dim[1]+jloc;
            if (jloc+ j_domain_begin != 0){
                rho [linindex] += C_m*charge_weight* Sx1[i]*Sy1[j]/abs((jloc+ j_domain_begin)*dr); // iloc = (i+ipo)*b_dim[1];
                }
            else {
                rho [linindex] +=  C_m*charge_weight* Sx1[i]*Sy1[j]*6./dr; // iloc = (i+ipo)*b_dim[1];
                //rho [linindex+1] += rho [linindex-1];
                //rho [linindex+2] += rho [linindex-2];
                }
        }
    }//i

 
    
} // END Project local current densities (Jl, Jr, Jt, sort)


// ---------------------------------------------------------------------------------------------------------------------
//! Project local densities only (Frozen species) NOT DONE YET
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorRZ2Order::operator() (double* rho, Particles &particles, unsigned int ipart, unsigned int bin, std::vector<unsigned int> &b_dim)
{ 
    //Warning : this function is used for frozen species only. It is assumed that position = position_old !!!

    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------

    int iloc;
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);

    // variable declaration
    double xpn, ypn;
    double delta, delta2;
    double Sx1[5], Sy1[5]; // arrays used for the Esirkepov projection method

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
    int ip = round(xpn);
    delta  = xpn - (double)ip;
    delta2 = delta*delta;
    Sx1[1] = 0.5 * (delta2-delta+0.25);
    Sx1[2] = 0.75-delta2;
    Sx1[3] = 0.5 * (delta2+delta+0.25);
    ypn = sqrt (particles.position(1, ipart)*particles.position(1, ipart)+particles.position(2, ipart)*particles.position(2, ipart))*dr_inv_ ;
    int jp = round(ypn);
    delta  = ypn - (double)jp;
    delta2 = delta*delta;
    Sy1[1] = 0.5 * (delta2-delta+0.25);
    Sy1[2] = 0.75-delta2;
    Sy1[3] = 0.5 * (delta2+delta+0.25);

    // ---------------------------
    // Calculate the total charge
    // ---------------------------
    ip -= i_domain_begin + bin +2;
    jp -= j_domain_begin + 2;

    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = (i+ip)*b_dim[1]+jp;
        for (unsigned int j=0 ; j<5 ; j++) {
            if (j+jp+2* j_domain_begin == 0){
                rho [iloc+j] += charge_weight*6.* Sx1[i]*Sy1[j] /dr; // iloc = (i+ipo)*b_dim[1];
                //rho [iloc+j+1] += rho [iloc+j-1];
                //rho [iloc+j+2] += rho [iloc+j-2];
         }
            else {
                rho [iloc+j] += charge_weight* Sx1[i]*Sy1[j]/abs((j+jp+2* j_domain_begin)*dr); // iloc = (i+ipo)*b_dim[1];
                }
            }
    }//i
} // END Project local current densities (Frozen species)

// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities : ionization NOT DONE YET
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorRZ2Order::operator() (Field* Jl, Field* Jr, Field* Jt, Particles &particles, int ipart, LocalFields Jion)
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
void ProjectorRZ2Order::operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread, int ibin, int clrw, bool diag_flag, bool is_spectral, std::vector<unsigned int> &b_dim, int ispec, int ipart_ref)
{  //std::cout<<"projecting"<<std::endl;
   if (is_spectral)
        ERROR("Not implemented");

    std::vector<int> *iold = &(smpi->dynamics_iold[ithread]);
    std::vector<double> *delta = &(smpi->dynamics_deltaold[ithread]);
    std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);   
    std::vector<std::complex<double>> *exp_m_theta_old = &(smpi->dynamics_thetaold[ithread]);
    int dim1 = EMfields->dimPrim[1];
    //int dim2 = EMfields->dimPrim[2];

    ElectroMagn3DRZ* emRZ = static_cast<ElectroMagn3DRZ*>( EMfields );

    // If no field diagnostics this timestep, then the projection is done directly on the total arrays
    if (!diag_flag){ 

        
        // Loop on modes ?
        for ( unsigned int imode = 0; imode<Nmode;imode++){

            if (imode==0){
                complex< double>* b_Jl =  &(*emRZ->Jl_[imode] )(ibin*clrw* dim1 );
                complex<double>* b_Jr =  &(*emRZ->Jr_[imode] )(ibin*clrw*(dim1+1) );
                complex<double>* b_Jt =  &(*emRZ->Jt_[imode] )(ibin*clrw* dim1 );
                for ( int ipart=istart ; ipart<iend; ipart++ )
                    (*this)(b_Jl , b_Jr , b_Jt , particles,  ipart, (*invgf)[ipart], ibin*clrw, b_dim, &(*iold)[ipart], &(*delta)[ipart]);
	    }
            else{	
                complex<double>* b_Jl =  &(*emRZ->Jl_[imode] )(ibin*clrw* dim1 );
                complex<double>* b_Jr =  &(*emRZ->Jr_[imode] )(ibin*clrw*(dim1+1) );
                complex<double>* b_Jt =  &(*emRZ->Jt_[imode] )(ibin*clrw* dim1 );
                for ( int ipart=istart ; ipart<iend; ipart++ )
                    (*this)(b_Jl , b_Jr , b_Jt , particles,  ipart,(*invgf)[ipart], ibin*clrw, b_dim, &(*iold)[ipart], &(*delta)[ipart],&(*exp_m_theta_old)[ipart], imode);
           } 
         }       // Otherwise, the projection may apply to the species-specific arrays
     } 
     else {
         //Loop on modes 
        for ( unsigned int imode = 0; imode<Nmode;imode++){
	    int ifield = imode*n_species+ispec;
            if (imode==0){
                complex<double>* b_Jl  = emRZ->Jl_s [ifield] ? &(* static_cast<cField2D*>(emRZ->Jl_s [ifield]) )(ibin*clrw* dim1) : &(*emRZ->Jl_[imode] )(ibin*clrw* dim1 ) ;
                complex<double>* b_Jr  = emRZ->Jr_s [ifield] ? &(* static_cast<cField2D*>(emRZ->Jr_s [ifield]) )(ibin*clrw*(dim1+1)) : &(*emRZ->Jr_[imode] )(ibin*clrw*(dim1+1)) ;
                complex<double>* b_Jt  = emRZ->Jt_s [ifield] ? &(* static_cast<cField2D*>(emRZ->Jt_s [ifield]) )(ibin*clrw*dim1) : &(*emRZ->Jt_[imode] )(ibin*clrw*dim1) ;
                complex<double>* b_rho = emRZ->rho_RZ_s[ifield] ? &(* static_cast<cField2D*>(emRZ->rho_RZ_s[ifield]) )(ibin*clrw* dim1) : &(*emRZ->rho_RZ_[imode])(ibin*clrw* dim1 ) ;
                for ( int ipart=istart ; ipart<iend; ipart++ )
                   (*this)(b_Jl , b_Jr , b_Jt ,b_rho, particles,  ipart, (*invgf)[ipart], ibin*clrw, b_dim, &(*iold)[ipart], &(*delta)[ipart]);
}
             else{    
                complex<double>* b_Jl  = emRZ->Jl_s [ifield] ? &(* static_cast<cField2D*>(emRZ->Jl_s [ifield]) )(ibin*clrw* dim1) : &(*emRZ->Jl_[imode] )(ibin*clrw* dim1) ;
                complex<double>* b_Jr  = emRZ->Jr_s [ifield] ? &(* static_cast<cField2D*>(emRZ->Jr_s [ifield]) )(ibin*clrw*(dim1+1)) : &(*emRZ->Jr_[imode] )(ibin*clrw*(dim1+1)) ;
                complex<double>* b_Jt  = emRZ->Jt_s [ifield] ? &(* static_cast<cField2D*>(emRZ->Jt_s [ifield]) )(ibin*clrw*dim1) : &(*emRZ->Jt_[imode])(ibin*clrw*dim1) ;
                complex<double>* b_rho = emRZ->rho_RZ_s[ifield] ? &(* static_cast<cField2D*>(emRZ->rho_RZ_s[ifield]) )(ibin*clrw* dim1 ) : &(*emRZ->rho_RZ_[imode])(ibin*clrw* dim1 ) ;
                for ( int ipart=istart ; ipart<iend; ipart++ )
                    (*this)(b_Jl , b_Jr , b_Jt ,b_rho, particles,  ipart, (*invgf)[ipart], ibin*clrw, b_dim, &(*iold)[ipart], &(*delta)[ipart], &(*exp_m_theta_old)[ipart], imode);
                 }

       }
   }
}
