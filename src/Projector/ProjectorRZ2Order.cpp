#include "ProjectorRZ2Order.h"

#include <cmath>
#include <iostream>

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
    dr_ov_dt  = params.cell_length[1] / params.timestep;
    nmodes=params.nmodes; 
    one_third = 1.0/3.0;

    i_domain_begin = patch->getCellStartingGlobalIndex(0);
    j_domain_begin = patch->getCellStartingGlobalIndex(1);

    DEBUG("cell_length "<< params.cell_length[0]);

}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for ProjectorRZ2Order
// ---------------------------------------------------------------------------------------------------------------------
ProjectorRZ2Order::~ProjectorRZ2Order()
{
}


// ---------------------------------------------------------------------------------------------------------------------
//! Project local currents (sort)
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorRZ2Order::operator() (complex<double>* Jl, complex<double>* Jr, complex<double>* Jt, Particles &particles, unsigned int ipart, double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaold)
{
    int nparts= particles.size();
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------   int iloc,
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
    double crl_p = charge_weight*dl_ov_dt;
    double crr_p = charge_weight*dr_ov_dt;
    double crt_p= charge_weight*particles.momentum(2,ipart)*invgf;

    // variable declaration
    double xpn, ypn;
    double delta, delta2;
    // arrays used for the Esirkepov projection method
    double  Sx0[5], Sx1[5], Sy0[5], Sy1[5], DSx[5], DSy[5], tmpJx[5];
    
/*  for (unsigned int i=0; i<5; i++) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
        // local array to accumulate Jx
        // Jx_p[i][j] = Jx_p[i-1][j] - crx_p * Wx[i-1][j];
        tmpJl[i] = 0.;
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
    
    ypn = sqrt (particles.position(1, ipart)*particles.position(1, ipart)+particles.position(2, ipart)*particles.position(2, ipart))*dr_inv ;
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
    
    // calculate Esirkepov coeff. Wx, Wy, Wz when used
    double tmp, tmp2, tmp3, tmpY;
    //Do not compute useless weights.
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
                Jz_p[i][j][k] = crt_p * Wz[i][j];
            }
        }

    // ---------------------------
    // Calculate the total current
    // ---------------------------
    
    ipo -= bin+2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
    // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
    jpo -= 2;
    
    int iloc, jloc, linindex;
    
    // Jx^(d,p,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*b_dim[1]+jloc;
            Jl [linindex] += Jx_p[i][j]; // iloc = (i+ipo)*b_dim[1];
            }
    }//i
    
    // Jy^(p,d,p)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*(b_dim[1]+1)+jloc*b_dim[2];
            Jr [linindex] += Jy_p[i][j]; //
            }
    }//i
    
    // Jz^(p,p,d)
    for (unsigned int i=0 ; i<5 ; i++) {
        iloc = i+ipo;
        for (unsigned int j=0 ; j<5 ; j++) {
            jloc = j+jpo;
            linindex = iloc*(b_dim[2]+1)*b_dim[1]+jloc*(b_dim[2]+1);
            Jt [linindex] += Jz_p[i][j]; //
            }
    }//i
    
   */ 

    
} // END Project local current densities (Jl, Jr, Jt, sort)


// ---------------------------------------------------------------------------------------------------------------------
//! Project local current densities (sort)
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorRZ2Order::operator() (complex<double>* Jl, complex<double>* Jr, complex<double>* Jt, complex<double>* rho, Particles &particles, unsigned int ipart, double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaold)
{

    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    
    // (x,y,z) components of the current density for the macro-particle
    
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
    double crl_p = charge_weight*dl_ov_dt;
    double crr_p = charge_weight*dr_ov_dt;

   

} // END Project local densities (Jl, Jr, Jt, rho, sort)


// ---------------------------------------------------------------------------------------------------------------------
//! Project local densities only (Frozen species)
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorRZ2Order::operator() (double* rho, Particles &particles, unsigned int ipart, unsigned int bin, std::vector<unsigned int> &b_dim)
{
    //Warning : this function is used for frozen species only. It is assumed that position = position_old !!!

    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------

    
    int iloc,jloc;
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);

   

} // END Project local current densities (Frozen species)

// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities (ionize)
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorRZ2Order::operator() (Field* Jl, Field* Jr, Field* Jt, Particles &particles, int ipart, LocalFields Jion)
{
    cField2D* Jl3D  = static_cast<cField2D*>(Jl);
    cField2D* Jr3D  = static_cast<cField2D*>(Jr);
    cField2D* Jt3D  = static_cast<cField2D*>(Jt);
    
    // weighted currents
    double Jx_ion = Jion.x * particles.weight(ipart);
    double Jy_ion = Jion.y * particles.weight(ipart);
    double Jz_ion = Jion.z * particles.weight(ipart);
    


} // END Project global current densities (ionize)

//Wrapper for projection
void ProjectorRZ2Order::operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread, int ibin, int clrw, bool diag_flag, bool is_spectral, std::vector<unsigned int> &b_dim, int ispec, int ipart_ref)
{
    if (is_spectral)
        ERROR("Not implemented");

    std::vector<int> *iold = &(smpi->dynamics_iold[ithread]);
    std::vector<double> *delta = &(smpi->dynamics_deltaold[ithread]);
    std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);
    
    int dim1 = EMfields->dimPrim[1];
    int dim2 = EMfields->dimPrim[2];

    ElectroMagn3DRZ* emRZ = static_cast<ElectroMagn3DRZ*>( EMfields );

    // If no field diagnostics this timestep, then the projection is done directly on the total arrays
    if (!diag_flag){ 

        
        // Loop on modes ?
        int imode = 0;

        complex<double>* b_Jl =  &(*emRZ->Jl_[imode] )(ibin*clrw* dim1   * dim2   );
        complex<double>* b_Jr =  &(*emRZ->Jr_[imode] )(ibin*clrw*(dim1+1)* dim2   );
        complex<double>* b_Jt =  &(*emRZ->Jt_[imode] )(ibin*clrw* dim1   *(dim2+1));
        for ( int ipart=istart ; ipart<iend; ipart++ )
            (*this)(b_Jl , b_Jr , b_Jt , particles,  ipart, (*invgf)[ipart], ibin*clrw, b_dim, &(*iold)[3*ipart], &(*delta)[3*ipart]);
            
        // Otherwise, the projection may apply to the species-specific arrays
    } else {
        // Loop on modes ?
        int imode = 0;

        complex<double>* b_Jl  = emRZ->Jl_s [ispec] ? &(* static_cast<cField2D*>(emRZ->Jl_s [ispec]) )(ibin*clrw* dim1   *dim2) : &(*emRZ->Jl_[imode] )(ibin*clrw* dim1   *dim2) ;
        complex<double>* b_Jr  = emRZ->Jr_s [ispec] ? &(* static_cast<cField2D*>(emRZ->Jr_s [ispec]) )(ibin*clrw*(dim1+1)*dim2) : &(*emRZ->Jr_[imode] )(ibin*clrw*(dim1+1)*dim2) ;
        complex<double>* b_Jt  = emRZ->Jt_s [ispec] ? &(* static_cast<cField2D*>(emRZ->Jt_s [ispec]) )(ibin*clrw*dim1*(dim2+1)) : &(*emRZ->Jt_[imode] )(ibin*clrw*dim1*(dim2+1)) ;
        complex<double>* b_rho = emRZ->rho_s[ispec] ? &(* static_cast<cField2D*>(emRZ->rho_s[ispec]) )(ibin*clrw* dim1   *dim2) : &(*emRZ->rho_RZ_[imode])(ibin*clrw* dim1   *dim2) ;
        for ( int ipart=istart ; ipart<iend; ipart++ )
            (*this)(b_Jl , b_Jr , b_Jt ,b_rho, particles,  ipart, (*invgf)[ipart], ibin*clrw, b_dim, &(*iold)[3*ipart], &(*delta)[3*ipart]);
    }

}
