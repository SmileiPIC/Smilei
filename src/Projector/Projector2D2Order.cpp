#include "Projector2D2Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field2D.h"
#include "Particles.h"
#include "Tools.h"
#include "SmileiMPI_Cart2D.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Projector2D2Order
// ---------------------------------------------------------------------------------------------------------------------
Projector2D2Order::Projector2D2Order (PicParams& params, SmileiMPI* smpi) : Projector2D(params, smpi)
{
    SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);

    dx_inv_   = 1.0/params.cell_length[0];
    dx_ov_dt  = params.cell_length[0] / params.timestep;
    dy_inv_   = 1.0/params.cell_length[1];
    dy_ov_dt  = params.cell_length[1] / params.timestep;

    one_third = 1.0/3.0;

    i_domain_begin = smpi2D->getCellStartingGlobalIndex(0);
    j_domain_begin = smpi2D->getCellStartingGlobalIndex(1);

    DEBUG("cell_length "<< params.cell_length[0]);

}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Projector2D2Order
// ---------------------------------------------------------------------------------------------------------------------
Projector2D2Order::~Projector2D2Order()
{
}


//! Below, in this order :
//!   Project global current densities (EMfields->Jx_/Jy_/Jz_), not used
//!   Projection by species
//!   Project global current charge
//!   Project local current densities (sort)
//!   Project global current densities (ionize)




// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities (EMfields->Jx_/Jy_/Jz_), not used
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D2Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, double gf)
{

    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------

    // static cast of the currents
    Field2D* Jx2D  = static_cast<Field2D*>(EMfields->Jx_);
    Field2D* Jy2D  = static_cast<Field2D*>(EMfields->Jy_);
    Field2D* Jz2D  = static_cast<Field2D*>(EMfields->Jz_);


    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
    double crx_p = charge_weight*dx_ov_dt;
    double cry_p = charge_weight*dy_ov_dt;
    double crz_p = charge_weight*particles.momentum(2, ipart)/gf;


    // variable declaration
    double xpn, ypn;
    double delta, delta2;
    double Sx0[5], Sx1[5], Sy0[5], Sy1[5], DSx[5], DSy[5]; // arrays used for the Esirkepov projection method
    double Wx[5][5], Wy[5][5], Wz[5][5];                   // idem
    double Jx_p[5][5], Jy_p[5][5];                         // idem


    // Initialize all current-related arrays to zero
    for (unsigned int i=0; i<5; i++) {
        Sx0[i] = 0.;
        Sx1[i] = 0.;
        Sy0[i] = 0.;
        Sy1[i] = 0.;
        DSx[i] = 0.;
        DSy[i] = 0.;
    }
    for (unsigned int i=0; i<5; i++) {
        for (unsigned int j=0; j<5; j++) {
            Wx[i][j]   = 0.;
            Wy[i][j]   = 0.;
            Wz[i][j]   = 0.;
            Jx_p[i][j] = 0.;
            Jy_p[i][j] = 0.;
        }
    }//i


    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------

    // locate the particle on the primal grid at former time-step & calculate coeff. S0
    xpn = particles.position_old(0, ipart) * dx_inv_;
    int ipo = round(xpn);
    delta  = xpn - (double)ipo;
    delta2 = delta*delta;
    Sx0[1] = 0.5 * (delta2-delta+0.25);
    Sx0[2] = 0.75-delta2;
    Sx0[3] = 0.5 * (delta2+delta+0.25);

    ypn = particles.position_old(1, ipart) * dy_inv_;
    int jpo = round(ypn);
    delta  = ypn - (double)jpo;
    delta2 = delta*delta;
    Sy0[1] = 0.5 * (delta2-delta+0.25);
    Sy0[2] = 0.75-delta2;
    Sy0[3] = 0.5 * (delta2+delta+0.25);


    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position(0, ipart) * dx_inv_;
    int ip = round(xpn);
    int ip_m_ipo = ip-ipo;
    delta  = xpn - (double)ip;
    delta2 = delta*delta;
    Sx1[ip_m_ipo+1] = 0.5 * (delta2-delta+0.25);
    Sx1[ip_m_ipo+2] = 0.75-delta2;
    Sx1[ip_m_ipo+3] = 0.5 * (delta2+delta+0.25);

    ypn = particles.position(1, ipart) * dy_inv_;
    int jp = round(ypn);
    int jp_m_jpo = jp-jpo;
    delta  = ypn - (double)jp;
    delta2 = delta*delta;
    Sy1[jp_m_jpo+1] = 0.5 * (delta2-delta+0.25);
    Sy1[jp_m_jpo+2] = 0.75-delta2;
    Sy1[jp_m_jpo+3] = 0.5 * (delta2+delta+0.25);


    // calculate Esirkepov coeff. Wx, Wy, Wz
    for (unsigned int i=0; i < 5; i++) {
        DSx[i] = Sx1[i] - Sx0[i];
        DSy[i] = Sy1[i] - Sy0[i];
    }

    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
            Wx[i][j] = DSx[i] * (Sy0[j] + 0.5*DSy[j]);
            Wy[i][j] = DSy[j] * (Sx0[i] + 0.5*DSx[i]);
            Wz[i][j] = one_third * ( Sx1[i] * (0.5*Sy0[j]+Sy1[j]) + Sx0[i] * (Sy0[j]+0.5*Sy1[j]) );
        }
    }


    // ------------------------------------------------
    // Local current created by the particle
    // calculate using the charge conservation equation
    // ------------------------------------------------
    for (unsigned int i=1 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
            Jx_p[i][j] = Jx_p[i-1][j] - crx_p * Wx[i-1][j];
        }
    }
    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=1 ; j<5 ; j++) {
            Jy_p[i][j] = Jy_p[i][j-1] - cry_p * Wy[i][j-1];
        }
    }


    // ---------------------------
    // Calculate the total current
    // ---------------------------
    ipo -= i_domain_begin;
    ip  -= i_domain_begin;
    jpo -= j_domain_begin;
    jp  -= j_domain_begin;

    for (unsigned int i=0 ; i<5 ; i++) {
        unsigned int iloc = i+ipo-2;

        for (unsigned int j=0 ; j<5 ; j++) {
            unsigned int jloc = j+jpo-2;
            (*Jx2D)(iloc,jloc) += Jx_p[i][j];
            (*Jy2D)(iloc,jloc) += Jy_p[i][j];
            (*Jz2D)(iloc,jloc) += crz_p * Wz[i][j];
        }

    }//i


} // END Project global current densities, not used


// ---------------------------------------------------------------------------------------------------------------------
//!   Projection by species
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D2Order::operator() (Field* Jx, Field* Jy, Field* Jz, Field* rho, Particles &particles, int ipart, double gf)
{

    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------

    // static cast of the currents
    Field2D* Jx2D  = static_cast<Field2D*>(Jx);
    Field2D* Jy2D  = static_cast<Field2D*>(Jy);
    Field2D* Jz2D  = static_cast<Field2D*>(Jz);
    Field2D* rho2D = static_cast<Field2D*>(rho);


    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
    double crx_p = charge_weight*dx_ov_dt;
    double cry_p = charge_weight*dy_ov_dt;
    double crz_p = charge_weight*particles.momentum(2, ipart)/gf;


    // variable declaration
    double xpn, ypn;
    double delta, delta2;
    double Sx0[5], Sx1[5], Sy0[5], Sy1[5], DSx[5], DSy[5]; // arrays used for the Esirkepov projection method
    double Wx[5][5], Wy[5][5], Wz[5][5];                   // idem
    double Jx_p[5][5], Jy_p[5][5];                         // idem


    // Initialize all current-related arrays to zero
    for (unsigned int i=0; i<5; i++) {
        Sx0[i] = 0.;
        Sx1[i] = 0.;
        Sy0[i] = 0.;
        Sy1[i] = 0.;
        DSx[i] = 0.;
        DSy[i] = 0.;
    }
    for (unsigned int i=0; i<5; i++) {
        for (unsigned int j=0; j<5; j++) {
            Wx[i][j]   = 0.;
            Wy[i][j]   = 0.;
            Wz[i][j]   = 0.;
            Jx_p[i][j] = 0.;
            Jy_p[i][j] = 0.;
        }
    }//i


    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------

    // locate the particle on the primal grid at former time-step & calculate coeff. S0
    xpn = particles.position_old(0, ipart) * dx_inv_;
    int ipo = round(xpn);
    delta  = xpn - (double)ipo;
    delta2 = delta*delta;
    Sx0[1] = 0.5 * (delta2-delta+0.25);
    Sx0[2] = 0.75-delta2;
    Sx0[3] = 0.5 * (delta2+delta+0.25);

    ypn = particles.position_old(1, ipart) * dy_inv_;
    int jpo = round(ypn);
    delta  = ypn - (double)jpo;
    delta2 = delta*delta;
    Sy0[1] = 0.5 * (delta2-delta+0.25);
    Sy0[2] = 0.75-delta2;
    Sy0[3] = 0.5 * (delta2+delta+0.25);


    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position(0, ipart) * dx_inv_;
    int ip = round(xpn);
    int ip_m_ipo = ip-ipo;
    delta  = xpn - (double)ip;
    delta2 = delta*delta;
    Sx1[ip_m_ipo+1] = 0.5 * (delta2-delta+0.25);
    Sx1[ip_m_ipo+2] = 0.75-delta2;
    Sx1[ip_m_ipo+3] = 0.5 * (delta2+delta+0.25);

    ypn = particles.position(1, ipart) * dy_inv_;
    int jp = round(ypn);
    int jp_m_jpo = jp-jpo;
    delta  = ypn - (double)jp;
    delta2 = delta*delta;
    Sy1[jp_m_jpo+1] = 0.5 * (delta2-delta+0.25);
    Sy1[jp_m_jpo+2] = 0.75-delta2;
    Sy1[jp_m_jpo+3] = 0.5 * (delta2+delta+0.25);


    // calculate Esirkepov coeff. Wx, Wy, Wz
    for (unsigned int i=0; i < 5; i++) {
        DSx[i] = Sx1[i] - Sx0[i];
        DSy[i] = Sy1[i] - Sy0[i];
    }

    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
            Wx[i][j] = DSx[i] * (Sy0[j] + 0.5*DSy[j]);
            Wy[i][j] = DSy[j] * (Sx0[i] + 0.5*DSx[i]);
            Wz[i][j] = one_third * ( Sx1[i] * (0.5*Sy0[j]+Sy1[j]) + Sx0[i] * (Sy0[j]+0.5*Sy1[j]) );
        }
    }


    // ------------------------------------------------
    // Local current created by the particle
    // calculate using the charge conservation equation
    // ------------------------------------------------
    for (unsigned int i=1 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
            Jx_p[i][j] = Jx_p[i-1][j] - crx_p * Wx[i-1][j];
        }
    }
    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=1 ; j<5 ; j++) {
            Jy_p[i][j] = Jy_p[i][j-1] - cry_p * Wy[i][j-1];
        }
    }


    // ---------------------------
    // Calculate the total current
    // ---------------------------
    ipo -= i_domain_begin;
    ip  -= i_domain_begin;
    jpo -= j_domain_begin;
    jp  -= j_domain_begin;

    for (unsigned int i=0 ; i<5 ; i++) {
        unsigned int iloc = i+ipo-2;

        for (unsigned int j=0 ; j<5 ; j++) {
            unsigned int jloc = j+jpo-2;
            (*Jx2D)(iloc,jloc)  += Jx_p[i][j];
            (*Jy2D)(iloc,jloc)  += Jy_p[i][j];
            (*Jz2D)(iloc,jloc)  += crz_p * Wz[i][j];
            (*rho2D)(iloc,jloc) += charge_weight * Sx1[i]*Sy1[j];
        }

    }//i

}//END Projection by species


// ---------------------------------------------------------------------------------------------------------------------
//! Project global current charge
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D2Order::operator() (Field* rho, Particles &particles, int ipart)
{

    //Static cast of the total charge density
    Field2D* rho2D  = static_cast<Field2D*>(rho);

    //Declaration of local variables
    double delta, delta2;
    double rho_p = particles.charge(ipart)*particles.weight(ipart);   // charge density of the macro-particle
    double Sx[3], Sy[3];             // projection coefficient arrays

    //Locate particle on the primal grid & calculate the projection coefficients
    double       xpn = particles.position(0, ipart) * dx_inv_;  // normalized distance to the first node
    int ic  = round(xpn);                   // index of the central node
    delta  = xpn - (double)ic;                       // normalized distance to the nearest grid point
    delta2 = delta*delta;                            // square of the normalized distance to the nearest grid point
    Sx[0]  = 0.5 * (delta2-delta+0.25);
    Sx[1]  = 0.75-delta2;
    Sx[2]  = 0.5 * (delta2+delta+0.25);

    double       ypn = particles.position(1, ipart) * dy_inv_;  // normalized distance to the first node
    int jc   = round(ypn);                  // index of the central node
    delta  = ypn - (double)jc;                       // normalized distance to the nearest grid point
    delta2 = delta*delta;                            // square of the normalized distance to the nearest grid point
    Sy[0]  = 0.5 * (delta2-delta+0.25);
    Sy[1]  = 0.75-delta2;
    Sy[2]  = 0.5 * (delta2+delta+0.25);

    //cout << "Pos = " << particles.position(0, ipart) << " - i global = " << i << " - i local = " << i-index_domain_begin <<endl;

    int i = ic-i_domain_begin-1; // index of first point for projection in x
    int j = jc-j_domain_begin-1; // index of first point for projection in y

    // 2nd order projection for the total charge density
    for (unsigned int iloc=0 ; iloc<3 ; iloc++) {
        for (unsigned int jloc=0 ; jloc<3 ; jloc++) {
            (*rho2D)(i+iloc,j+jloc) += Sx[iloc]*Sy[jloc]*rho_p;
        }
    }

} // END Project global current charge


// ---------------------------------------------------------------------------------------------------------------------
//! Project local current densities (sort)
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D2Order::operator() (double* Jx, double* Jy, double* Jz, double* rho, Particles &particles, int ipart, double gf, unsigned int bin, unsigned int b_dim1)
{

    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------

    int iloc, jloc;
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
    double crx_p = charge_weight*dx_ov_dt;
    double cry_p = charge_weight*dy_ov_dt;
    double crz_p = charge_weight*particles.momentum(2, ipart)/gf;


    // variable declaration
    double xpn, ypn;
    double delta, delta2;
    // arrays used for the Esirkepov projection method
    double  Sx0[5], Sx1[5], Sy0[5], Sy1[5], DSx[5], DSy[5], tmpJx[5];

    for (unsigned int i=0; i<5; i++) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
	// local array to accumulate Jx
	// Jx_p[i][j] = Jx_p[i-1][j] - crx_p * Wx[i-1][j];
	tmpJx[i] = 0.;
    }
    Sx0[0] = 0.;
    Sx0[4] = 0.;
    Sy0[0] = 0.;
    Sy0[4] = 0.;

    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------

    // locate the particle on the primal grid at former time-step & calculate coeff. S0
    xpn = particles.position_old(0, ipart) * dx_inv_;
    int ipo = round(xpn);
    delta  = xpn - (double)ipo;
    delta2 = delta*delta;
    Sx0[1] = 0.5 * (delta2-delta+0.25);
    Sx0[2] = 0.75-delta2;
    Sx0[3] = 0.5 * (delta2+delta+0.25);

    ypn = particles.position_old(1, ipart) * dy_inv_;
    int jpo = round(ypn);
    delta  = ypn - (double)jpo;
    delta2 = delta*delta;
    Sy0[1] = 0.5 * (delta2-delta+0.25);
    Sy0[2] = 0.75-delta2;
    Sy0[3] = 0.5 * (delta2+delta+0.25);


    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position(0, ipart) * dx_inv_;
    int ip = round(xpn);
    int ip_m_ipo = ip-ipo;
    delta  = xpn - (double)ip;
    delta2 = delta*delta;
    Sx1[ip_m_ipo+1] = 0.5 * (delta2-delta+0.25);
    Sx1[ip_m_ipo+2] = 0.75-delta2;
    Sx1[ip_m_ipo+3] = 0.5 * (delta2+delta+0.25);

    ypn = particles.position(1, ipart) * dy_inv_;
    int jp = round(ypn);
    int jp_m_jpo = jp-jpo;
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

    // ---------------------------
    // Calculate the total current
    // ---------------------------
    ipo -= i_domain_begin + bin;
    jpo -= j_domain_begin;
    // i =0
    {
	iloc = (ipo-2)*b_dim1;
	jloc = iloc+jpo-2; 
	tmp2 = 0.5*Sx1[0];
	tmp3 =     Sx1[0];
	Jz[jloc]  += crz_p * one_third * ( Sy1[0]*tmp3 );
	rho[jloc] += charge_weight * Sx1[0]*Sy1[0];	
	tmp = 0;
	tmpY = Sx0[0] + 0.5*DSx[0];
	for (unsigned int j=1 ; j<5 ; j++) {
	    jloc = iloc+j+jpo-2; 
	    tmp -= cry_p * DSy[j-1] * tmpY;
	    Jy[jloc]  += tmp;
	    Jz[jloc]  += crz_p * one_third * ( Sy0[j]*tmp2 + Sy1[j]*tmp3 );
	    rho[jloc] += charge_weight * Sx1[0]*Sy1[j];
	}

    }//i

    for (unsigned int i=1 ; i<5 ; i++) {
        iloc = (i+ipo-2)*b_dim1;
	jloc = iloc+jpo-2; 
	tmpJx[0] -= crx_p *  DSx[i-1] * (0.5*DSy[0]);
	Jx[jloc]  += tmpJx[0];
        tmp2 = 0.5*Sx1[i] + Sx0[i];
        tmp3 = 0.5*Sx0[i] + Sx1[i];
	Jz[jloc]  += crz_p * one_third * ( Sy1[0]*tmp3 );
	rho[jloc] += charge_weight * Sx1[i]*Sy1[0];	
	tmp = 0;
	tmpY = Sx0[i] + 0.5*DSx[i];
        for (unsigned int j=1 ; j<5 ; j++) {
            jloc = iloc+j+jpo-2; 
            tmpJx[j] -= crx_p * DSx[i-1] * (Sy0[j] + 0.5*DSy[j]);
	    Jx[jloc]  += tmpJx[j];
	    tmp -= cry_p * DSy[j-1] * tmpY;
            Jy[jloc]  += tmp;
            Jz[jloc]  += crz_p * one_third * ( Sy0[j]*tmp2 + Sy1[j]*tmp3 );
	    rho[jloc] += charge_weight * Sx1[i]*Sy1[j];
        }

    }//i

} // END Project local current densities (sort)


// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities (ionize)
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D2Order::operator() (Field* Jx, Field* Jy, Field* Jz, Particles &particles, int ipart, LocalFields Jion)
{
    ERROR("Projection of ionization current not yet defined for 2D 2nd order");

} // END Project global current densities (ionize)
