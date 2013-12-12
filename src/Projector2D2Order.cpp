#include "Projector2D2Order.h"
#include "ElectroMagn.h"
#include "Field2D.h"
#include "Particle.h"
#include "Tools.h"
#include "SmileiMPI_Cart2D.h"

#include <iostream>
#include <cmath>

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Projector2D2Order
// ---------------------------------------------------------------------------------------------------------------------
Projector2D2Order::Projector2D2Order (PicParams* params, SmileiMPI* smpi) : Projector2D(params, smpi)
{
	SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);

	dx_inv_   = 1.0/params->cell_length[0];
	dx_ov_dt  = params->cell_length[0] / params->timestep;
    dy_inv_   = 1.0/params->cell_length[1];
	dy_ov_dt  = params->cell_length[1] / params->timestep;
    
    one_third = 1.0/3.0;
    
	i_domain_begin = smpi2D->getCellStartingGlobalIndex(0);
    j_domain_begin = smpi2D->getCellStartingGlobalIndex(1);

	DEBUG("cell_length "<< params->cell_length[0]);

}



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Projector2D2Order
// ---------------------------------------------------------------------------------------------------------------------
Projector2D2Order::~Projector2D2Order()
{
}



// ---------------------------------------------------------------------------------------------------------------------
// 2nd order projection in 1d3v simulations
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D2Order::operator() (ElectroMagn* EMfields, Particle* part, double gf)
{
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    
    // static cast of the currents
    Field2D* Jx2D  = static_cast<Field2D*>(EMfields->Jx_);
	Field2D* Jy2D  = static_cast<Field2D*>(EMfields->Jy_);
	Field2D* Jz2D  = static_cast<Field2D*>(EMfields->Jz_);
    
    
    // (x,y,z) components of the current density for the macro-particle
	double charge_weight = (double)(part->charge())*part->weight();
    double crx_p = charge_weight*dx_ov_dt;
    double cry_p = charge_weight*dy_ov_dt;
    double crz_p = charge_weight*part->momentum(2)/gf;
    
    
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
    xpn = part->position_old(0) * dx_inv_;
    unsigned int ipo = round(xpn);
    delta  = xpn - (double)ipo;
    delta2 = delta*delta;
    Sx0[1] = 0.5 * (delta2-delta+0.25);
	Sx0[2] = 0.75-delta2;
	Sx0[3] = 0.5 * (delta2+delta+0.25);
    
    ypn = part->position_old(1) * dy_inv_;
    unsigned int jpo = round(ypn);
    delta  = ypn - (double)jpo;
    delta2 = delta*delta;
    Sy0[1] = 0.5 * (delta2-delta+0.25);
	Sy0[2] = 0.75-delta2;
	Sy0[3] = 0.5 * (delta2+delta+0.25);
    
    
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = part->position(0) * dx_inv_;
    unsigned int ip = round(xpn);
    unsigned int ip_m_ipo = ip-ipo;
    delta  = xpn - (double)ip;
    delta2 = delta*delta;
    Sx1[ip_m_ipo+1] = 0.5 * (delta2-delta+0.25);
	Sx1[ip_m_ipo+2] = 0.75-delta2;
	Sx1[ip_m_ipo+3] = 0.5 * (delta2+delta+0.25);
    
    ypn = part->position(1) * dy_inv_;
    unsigned int jp = round(ypn);
    unsigned int jp_m_jpo = jp-jpo;
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
//#pragma omp atomic
          (*Jx2D)(iloc,jloc) += Jx_p[i][j];
//#pragma omp atomic
          (*Jy2D)(iloc,jloc) += Jy_p[i][j];
//#pragma omp atomic
         (*Jz2D)(iloc,jloc) += crz_p * Wz[i][j];
        }
        
	}//i

    
} // END Projector2D2Order


void Projector2D2Order::operator() (Field* rho, Particle* part)
{
    
    //Static cast of the total charge density
    Field2D* rho2D  = static_cast<Field2D*>(rho);
	
	//Declaration of local variables
	double delta, delta2;
	double rho_p = part->charge()*part->weight();   // charge density of the macro-particle
    double Sx[3], Sy[3];             // projection coefficient arrays
    
	//Locate particle on the primal grid & calculate the projection coefficients
	double       xpn = part->position(0) * dx_inv_;  // normalized distance to the first node
	unsigned int ic  = round(xpn);                   // index of the central node
	delta  = xpn - (double)ic;                       // normalized distance to the nearest grid point
	delta2 = delta*delta;                            // square of the normalized distance to the nearest grid point
    Sx[0]  = 0.5 * (delta2-delta+0.25);
    Sx[1]  = 0.75-delta2;
    Sx[2]  = 0.5 * (delta2+delta+0.25);
    
    double       ypn = part->position(1) * dy_inv_;  // normalized distance to the first node
	unsigned int jc   = round(ypn);                  // index of the central node
	delta  = ypn - (double)jc;                       // normalized distance to the nearest grid point
	delta2 = delta*delta;                            // square of the normalized distance to the nearest grid point
    Sy[0]  = 0.5 * (delta2-delta+0.25);
    Sy[1]  = 0.75-delta2;
    Sy[2]  = 0.5 * (delta2+delta+0.25);
    
	//cout << "Pos = " << part->position(0) << " - i global = " << i << " - i local = " << i-index_domain_begin <<endl;
    
	unsigned int i = ic-i_domain_begin-1; // index of first point for projection in x
    unsigned int j = jc-j_domain_begin-1; // index of first point for projection in y
	
	// 2nd order projection for the total charge density
    for (unsigned int iloc=0 ; iloc<3 ; iloc++) {
        for (unsigned int jloc=0 ; jloc<3 ; jloc++) {
            (*rho2D)(i+iloc,j+jloc) += Sx[iloc]*Sy[jloc]*rho_p;
        }
    }
    
}//END TotalChargeDensityProjection

