#include "Projector2D2Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field2D.h"
#include "Particles.h"
#include "Tools.h"
#include "Patch.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Projector2D2Order
// ---------------------------------------------------------------------------------------------------------------------
Projector2D2Order::Projector2D2Order( Params &params, Patch *patch ) : Projector2D( params, patch )
{
    dx_inv_   = 1.0/params.cell_length[0];
    dx_ov_dt  = params.cell_length[0] / params.timestep;
    dy_inv_   = 1.0/params.cell_length[1];
    dy_ov_dt  = params.cell_length[1] / params.timestep;
    
    i_domain_begin = patch->getCellStartingGlobalIndex( 0 );
    j_domain_begin = patch->getCellStartingGlobalIndex( 1 );
    
    nprimy = params.n_space[1] + 2*params.oversize[1] + 1;
    
    pxr = !params.is_pxr;
    
    dt             = params.timestep;
    dts2           = params.timestep/2.;
    dts4           = params.timestep/4.;
    
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Projector2D2Order
// ---------------------------------------------------------------------------------------------------------------------
Projector2D2Order::~Projector2D2Order()
{
}


// ---------------------------------------------------------------------------------------------------------------------
//! Project current densities : main projector
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D2Order::currents( double *Jx, double *Jy, double *Jz, Particles &particles, unsigned int ipart, double invgf, int *iold, double *deltaold )
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
    double delta, delta2;
    // arrays used for the Esirkepov projection method
    double  Sx0[5], Sx1[5], Sy0[5], Sy1[5], DSx[5], DSy[5], tmpJx[5];
    
    for( unsigned int i=0; i<5; i++ ) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
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
    delta = deltaold[0*nparts];
    delta2 = delta*delta;
    Sx0[1] = 0.5 * ( delta2-delta+0.25 );
    Sx0[2] = 0.75-delta2;
    Sx0[3] = 0.5 * ( delta2+delta+0.25 );
    
    delta = deltaold[1*nparts];
    delta2 = delta*delta;
    Sy0[1] = 0.5 * ( delta2-delta+0.25 );
    Sy0[2] = 0.75-delta2;
    Sy0[3] = 0.5 * ( delta2+delta+0.25 );
    
    
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position( 0, ipart ) * dx_inv_;
    int ip = round( xpn );
    int ipo = iold[0*nparts];
    int ip_m_ipo = ip-ipo-i_domain_begin;
    delta  = xpn - ( double )ip;
    delta2 = delta*delta;
    Sx1[ip_m_ipo+1] = 0.5 * ( delta2-delta+0.25 );
    Sx1[ip_m_ipo+2] = 0.75-delta2;
    Sx1[ip_m_ipo+3] = 0.5 * ( delta2+delta+0.25 );
    
    ypn = particles.position( 1, ipart ) * dy_inv_;
    int jp = round( ypn );
    int jpo = iold[1*nparts];
    int jp_m_jpo = jp-jpo-j_domain_begin;
    delta  = ypn - ( double )jp;
    delta2 = delta*delta;
    // cerr << " ipart: " << ipart
    //      << " nparts: " << nparts
    //      << " iold[0]: " << iold[0*nparts]
    //      << " iold[1]: " << iold[1*nparts]
    //      << " ip_m_ipo: " << ip_m_ipo
    //      << " jp_m_jpo: " << jp_m_jpo
    //      << " ypn: " << ypn
    //      << endl;
    
    Sy1[jp_m_jpo+1] = 0.5 * ( delta2-delta+0.25 );
    Sy1[jp_m_jpo+2] = 0.75-delta2;
    Sy1[jp_m_jpo+3] = 0.5 * ( delta2+delta+0.25 );
    
    for( unsigned int i=0; i < 5; i++ ) {
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
    ipo -= 2; //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
    jpo -= 2;
    // i =0
    {
        iloc = ipo*nprimy+jpo;
        tmp2 = 0.5*Sx1[0];
        tmp3 =     Sx1[0];
        Jz[iloc]  += crz_p * ( Sy1[0]*tmp3 );
        tmp = 0;
        tmpY = Sx0[0] + 0.5*DSx[0];
        for( unsigned int j=1 ; j<5 ; j++ ) {
            tmp -= cry_p * DSy[j-1] * tmpY;
            Jy[iloc+j+pxr*ipo]  += tmp; //Because size of Jy in Y is nprimy+1.
            Jz[iloc+j]  += crz_p * ( Sy0[j]*tmp2 + Sy1[j]*tmp3 );
            // cerr << " iloc+j+ipo: " << iloc+j+ipo
            //      << " iloc: " << iloc
            //      << " ipo: " << ipo
            //      << endl;
        }
    }//i
    
    for( unsigned int i=1 ; i<5 ; i++ ) {
        iloc = ( i+ipo )*nprimy+jpo;
        // cerr << " iloc: " << iloc
        //      << endl;
        tmpJx[0] -= crx_p *  DSx[i-1] * ( 0.5*DSy[0] );
        Jx[iloc]  += tmpJx[0];
        tmp2 = 0.5*Sx1[i] + Sx0[i];
        tmp3 = 0.5*Sx0[i] + Sx1[i];
        Jz[iloc]  += crz_p * ( Sy1[0]*tmp3 );
        tmp = 0;
        tmpY = Sx0[i] + 0.5*DSx[i];
        for( unsigned int j=1 ; j<5 ; j++ ) {
            tmpJx[j] -= crx_p * DSx[i-1] * ( Sy0[j] + 0.5*DSy[j] );
            Jx[iloc+j]  += tmpJx[j];
            tmp -= cry_p * DSy[j-1] * tmpY;
            Jy[iloc+j+pxr*(i+ipo)]  += tmp; //Because size of Jy in Y is nprimy+1.
            Jz[iloc+j]  += crz_p * ( Sy0[j]*tmp2 + Sy1[j]*tmp3 );
        }
    }//i
    
} // END Project local current densities (sort)


// ---------------------------------------------------------------------------------------------------------------------
//!  Project current densities & charge : diagFields timstep
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D2Order::currentsAndDensity( double *Jx, double *Jy, double *Jz, double *rho, Particles &particles, unsigned int ipart, double invgf, int *iold, double *deltaold )
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
    double delta, delta2;
    // arrays used for the Esirkepov projection method
    double  Sx0[5], Sx1[5], Sy0[5], Sy1[5], DSx[5], DSy[5], tmpJx[5];
    
    for( unsigned int i=0; i<5; i++ ) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
        // local array to accumulate Jx
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
    delta = deltaold[0*nparts];
    delta2 = delta*delta;
    Sx0[1] = 0.5 * ( delta2-delta+0.25 );
    Sx0[2] = 0.75-delta2;
    Sx0[3] = 0.5 * ( delta2+delta+0.25 );
    
    delta = deltaold[1*nparts];
    delta2 = delta*delta;
    Sy0[1] = 0.5 * ( delta2-delta+0.25 );
    Sy0[2] = 0.75-delta2;
    Sy0[3] = 0.5 * ( delta2+delta+0.25 );
    
    
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position( 0, ipart ) * dx_inv_;
    int ip = round( xpn );
    int ipo = iold[0*nparts];
    int ip_m_ipo = ip-ipo-i_domain_begin;
    delta  = xpn - ( double )ip;
    delta2 = delta*delta;
    Sx1[ip_m_ipo+1] = 0.5 * ( delta2-delta+0.25 );
    Sx1[ip_m_ipo+2] = 0.75-delta2;
    Sx1[ip_m_ipo+3] = 0.5 * ( delta2+delta+0.25 );
    
    ypn = particles.position( 1, ipart ) * dy_inv_;
    int jp = round( ypn );
    int jpo = iold[1*nparts];
    int jp_m_jpo = jp-jpo-j_domain_begin;
    delta  = ypn - ( double )jp;
    delta2 = delta*delta;
    Sy1[jp_m_jpo+1] = 0.5 * ( delta2-delta+0.25 );
    Sy1[jp_m_jpo+2] = 0.75-delta2;
    Sy1[jp_m_jpo+3] = 0.5 * ( delta2+delta+0.25 );
    
    for( unsigned int i=0; i < 5; i++ ) {
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
    ipo -= 2; //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
    jpo -= 2;
    // case i =0
    {
        iloc = ipo*nprimy+jpo;
        tmp2 = 0.5*Sx1[0];
        tmp3 =     Sx1[0];
        Jz[iloc]  += crz_p * ( Sy1[0]*tmp3 );
        rho[iloc] += charge_weight * Sx1[0]*Sy1[0];
        tmp = 0;
        tmpY = Sx0[0] + 0.5*DSx[0];
        for( unsigned int j=1 ; j<5 ; j++ ) {
            tmp -= cry_p * DSy[j-1] * tmpY;
            Jy[iloc+j+pxr*ipo]  += tmp; //Because size of Jy in Y is nprimy+1.
            Jz[iloc+j]  += crz_p * ( Sy0[j]*tmp2 + Sy1[j]*tmp3 );
            rho[iloc+j] += charge_weight * Sx1[0]*Sy1[j];
        }
        
    }//end i=0 case
    
    // case i> 0
    for( unsigned int i=1 ; i<5 ; i++ ) {
        iloc = ( i+ipo )*nprimy+jpo;
        tmpJx[0] -= crx_p *  DSx[i-1] * ( 0.5*DSy[0] );
        Jx[iloc]  += tmpJx[0];
        tmp2 = 0.5*Sx1[i] + Sx0[i];
        tmp3 = 0.5*Sx0[i] + Sx1[i];
        Jz[iloc]  += crz_p * ( Sy1[0]*tmp3 );
        rho[iloc] += charge_weight * Sx1[i]*Sy1[0];
        tmp = 0;
        tmpY = Sx0[i] + 0.5*DSx[i];
        for( unsigned int j=1 ; j<5 ; j++ ) {
            tmpJx[j] -= crx_p * DSx[i-1] * ( Sy0[j] + 0.5*DSy[j] );
            Jx[iloc+j]  += tmpJx[j];
            tmp -= cry_p * DSy[j-1] * tmpY;
                Jy[iloc+j+pxr*(i+ipo)]  += tmp; //Because size of Jy in Y is nprimy+1.
            Jz[iloc+j]  += crz_p * ( Sy0[j]*tmp2 + Sy1[j]*tmp3 );
            rho[iloc+j] += charge_weight * Sx1[i]*Sy1[j];
        }
        
    }//i
} // END Project local current densities (sort)


// ---------------------------------------------------------------------------------------------------------------------
//! Project charge : frozen & diagFields timstep
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D2Order::basic( double *rhoj, Particles &particles, unsigned int ipart, unsigned int type )
{
    //Warning : this function is used for frozen species only. It is assumed that position = position_old !!!
    
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    
    int iloc, ny( nprimy );
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
            ny++;
        } else {
            charge_weight *= particles.momentum( 2, ipart );
        }
    }
    
    // variable declaration
    double xpn, ypn;
    double delta, delta2;
    double Sx1[5], Sy1[5]; // arrays used for the Esirkepov projection method
    
    // Initialize all current-related arrays to zero
    for( unsigned int i=0; i<5; i++ ) {
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
    Sx1[1] = 0.5 * ( delta2-delta+0.25 );
    Sx1[2] = 0.75-delta2;
    Sx1[3] = 0.5 * ( delta2+delta+0.25 );
    
    ypn = particles.position( 1, ipart ) * dy_inv_;
    int jp = round( ypn + 0.5*( type==2 ) );
    delta  = ypn - ( double )jp;
    delta2 = delta*delta;
    Sy1[1] = 0.5 * ( delta2-delta+0.25 );
    Sy1[2] = 0.75-delta2;
    Sy1[3] = 0.5 * ( delta2+delta+0.25 );
    
    // ---------------------------
    // Calculate the total current
    // ---------------------------
    ip -= i_domain_begin + 2;
    jp -= j_domain_begin + 2;
    
    for( unsigned int i=0 ; i<5 ; i++ ) {
        iloc = ( i+ip )*ny+jp;
        for( unsigned int j=0 ; j<5 ; j++ ) {
            rhoj[iloc+j] += charge_weight * Sx1[i]*Sy1[j];
        }
    }//i
} // END Project local current densities (sort)


// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities : ionization
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D2Order::ionizationCurrents( Field *Jx, Field *Jy, Field *Jz, Particles &particles, int ipart, LocalFields Jion )
{
    Field2D *Jx2D  = static_cast<Field2D *>( Jx );
    Field2D *Jy2D  = static_cast<Field2D *>( Jy );
    Field2D *Jz2D  = static_cast<Field2D *>( Jz );
    
    
    //Declaration of local variables
    int ip, id, jp, jd;
    double xpn, xpmxip, xpmxip2, xpmxid, xpmxid2;
    double ypn, ypmyjp, ypmyjp2, ypmyjd, ypmyjd2;
    double Sxp[3], Sxd[3], Syp[3], Syd[3];
    
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
    
    // x-dual index
    id      = round( xpn+0.5 );              // x-index of the central node
    xpmxid  = xpn - ( double )id + 0.5;      // normalized distance to the nearest grid point
    xpmxid2 = xpmxid*xpmxid;                 // square of the normalized distance to the nearest grid point
    
    // y-primal index
    jp      = round( ypn );                  // y-index of the central node
    ypmyjp  = ypn - ( double )jp;            // normalized distance to the nearest grid point
    ypmyjp2 = ypmyjp*ypmyjp;                 // square of the normalized distance to the nearest grid point
    
    // y-dual index
    jd      = round( ypn+0.5 );              // y-index of the central node
    ypmyjd  = ypn - ( double )jd + 0.5;      // normalized distance to the nearest grid point
    ypmyjd2 = ypmyjd*ypmyjd;                 // square of the normalized distance to the nearest grid point
    
    Sxp[0] = 0.5 * ( xpmxip2-xpmxip+0.25 );
    Sxp[1] = ( 0.75-xpmxip2 );
    Sxp[2] = 0.5 * ( xpmxip2+xpmxip+0.25 );
    
    Sxd[0] = 0.5 * ( xpmxid2-xpmxid+0.25 );
    Sxd[1] = ( 0.75-xpmxid2 );
    Sxd[2] = 0.5 * ( xpmxid2+xpmxid+0.25 );
    
    Syp[0] = 0.5 * ( ypmyjp2-ypmyjp+0.25 );
    Syp[1] = ( 0.75-ypmyjp2 );
    Syp[2] = 0.5 * ( ypmyjp2+ypmyjp+0.25 );
    
    Syd[0] = 0.5 * ( ypmyjd2-ypmyjd+0.25 );
    Syd[1] = ( 0.75-ypmyjd2 );
    Syd[2] = 0.5 * ( ypmyjd2+ypmyjd+0.25 );
    
    ip  -= i_domain_begin;
    id  -= i_domain_begin;
    jp  -= j_domain_begin;
    jd  -= j_domain_begin;
    
    for( unsigned int i=0 ; i<3 ; i++ ) {
        int iploc=ip+i-1;
        int idloc=id+i-1;
        for( unsigned int j=0 ; j<3 ; j++ ) {
            int jploc=jp+j-1;
            int jdloc=jd+j-1;
            // Jx^(d,p)
            ( *Jx2D )( idloc, jploc ) += Jx_ion * Sxd[i]*Syp[j];
            // Jy^(p,d)
            ( *Jy2D )( iploc, jdloc ) += Jy_ion * Sxp[i]*Syd[j];
            // Jz^(p,p)
            ( *Jz2D )( iploc, jploc ) += Jz_ion * Sxp[i]*Syp[j];
        }
    }//i
    
    
} // END Project global current densities (ionize)


// ---------------------------------------------------------------------------------------------------------------------
//! Wrapper for projection
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D2Order::currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int icell, int ipart_ref )
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
                // cerr << ipart << endl;
                // cerr << ( *iold )[ipart] << endl;
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
void Projector2D2Order::susceptibility( ElectroMagn *EMfields, Particles &particles, double species_mass, SmileiMPI *smpi, int istart, int iend,  int ithread, int icell, int ipart_ref )

{
    double *Chi_envelope = &( *EMfields->Env_Chi_ )( 0 );
    
    std::vector<double> *Epart       = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Phipart     = &( smpi->dynamics_PHIpart[ithread] );
    std::vector<double> *GradPhipart = &( smpi->dynamics_GradPHIpart[ithread] );
    std::vector<double> *inv_gamma_ponderomotive = &( smpi->dynamics_inv_gamma_ponderomotive[ithread] );
    
    int iloc;
    
    double momentum[3];
    
    double gamma_ponderomotive, gamma0, gamma0_sq;
    double charge_over_mass_dts2, charge_sq_over_mass_sq_dts4, charge_sq_over_mass_sq;
    double pxsm, pysm, pzsm;
    double one_over_mass=1./species_mass;
    
    int nparts = particles.size();
    double *Ex       = &( ( *Epart )[0*nparts] );
    double *Ey       = &( ( *Epart )[1*nparts] );
    double *Ez       = &( ( *Epart )[2*nparts] );
    double *Phi      = &( ( *Phipart )[0*nparts] );
    double *GradPhix = &( ( *GradPhipart )[0*nparts] );
    double *GradPhiy = &( ( *GradPhipart )[1*nparts] );
    double *GradPhiz = &( ( *GradPhipart )[2*nparts] );
    
    for( int ipart=istart ; ipart<iend; ipart++ ) {//Loop on bin particles
    
        charge_over_mass_dts2       = ( double )( particles.charge( ipart ) )*dts2*one_over_mass;
        // ! ponderomotive force is proportional to charge squared and the field is divided by 4 instead of 2
        charge_sq_over_mass_sq_dts4 = ( double )( particles.charge( ipart ) )*( double )( particles.charge( ipart ) )*dts4*one_over_mass*one_over_mass;
        // (charge over mass)^2
        charge_sq_over_mass_sq      = ( double )( particles.charge( ipart ) )*( double )( particles.charge( ipart ) )*one_over_mass*one_over_mass;
        
        for( int i = 0 ; i<3 ; i++ ) {
            momentum[i] = particles.momentum( i, ipart );
        }
        
        // compute initial ponderomotive gamma
        gamma0_sq = 1. + momentum[0]*momentum[0]+ momentum[1]*momentum[1] + momentum[2]*momentum[2] + *( Phi+ipart )*charge_sq_over_mass_sq ;
        gamma0    = sqrt( gamma0_sq ) ;
        
        // ( electric field + ponderomotive force for ponderomotive gamma advance ) scalar multiplied by momentum
        pxsm = ( gamma0 * charge_over_mass_dts2*( *( Ex+ipart ) ) - charge_sq_over_mass_sq_dts4*( *( GradPhix+ipart ) ) ) * momentum[0] / gamma0_sq;
        pysm = ( gamma0 * charge_over_mass_dts2*( *( Ey+ipart ) ) - charge_sq_over_mass_sq_dts4*( *( GradPhiy+ipart ) ) ) * momentum[1] / gamma0_sq;
        pzsm = ( gamma0 * charge_over_mass_dts2*( *( Ez+ipart ) ) - charge_sq_over_mass_sq_dts4*( *( GradPhiz+ipart ) ) ) * momentum[2] / gamma0_sq;
        
        // update of gamma ponderomotive
        gamma_ponderomotive = gamma0 + ( pxsm+pysm+pzsm )*0.5 ;
        // buffer inverse of ponderomotive gamma to use it in ponderomotive momentum pusher
        ( *inv_gamma_ponderomotive )[ipart] = 1./gamma_ponderomotive;
        
        // susceptibility for the macro-particle
        double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*( double )( particles.charge( ipart ) )*particles.weight( ipart )*one_over_mass/gamma_ponderomotive;
        
        // variable declaration
        double xpn, ypn;
        double delta, delta2;
        double Sx1[5], Sy1[5]; // arrays used for the Esirkepov projection method
        
        // Initialize all current-related arrays to zero
        for( unsigned int i=0; i<5; i++ ) {
            Sx1[i] = 0.;
            Sy1[i] = 0.;
        }
        
        // --------------------------------------------------------
        // Locate particles & Calculate Esirkepov coef. S, DS and W
        // --------------------------------------------------------
        
        // locate the particle on the primal grid at current time-step & calculate coeff. S1
        xpn = particles.position( 0, ipart ) * dx_inv_;
        int ip = round( xpn );
        delta  = xpn - ( double )ip;
        delta2 = delta*delta;
        Sx1[1] = 0.5 * ( delta2-delta+0.25 );
        Sx1[2] = 0.75-delta2;
        Sx1[3] = 0.5 * ( delta2+delta+0.25 );
        
        ypn = particles.position( 1, ipart ) * dy_inv_;
        int jp = round( ypn );
        delta  = ypn - ( double )jp;
        delta2 = delta*delta;
        Sy1[1] = 0.5 * ( delta2-delta+0.25 );
        Sy1[2] = 0.75-delta2;
        Sy1[3] = 0.5 * ( delta2+delta+0.25 );
        
        
        // ---------------------------
        // Calculate the total susceptibility
        // ---------------------------
        ip -= i_domain_begin + 2;
        jp -= j_domain_begin + 2;
        
        for( unsigned int i=0 ; i<5 ; i++ ) {
            iloc = ( i+ip )*nprimy+jp;
            for( unsigned int j=0 ; j<5 ; j++ ) {
                Chi_envelope[iloc+j] += charge_weight * Sx1[i]*Sy1[j];
            }
        }//i
        
        
    }
    
}
