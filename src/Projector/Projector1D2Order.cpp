#include "Projector1D2Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field1D.h"
#include "Particles.h"
#include "Tools.h"
#include "Patch.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Projector1D2Order
// ---------------------------------------------------------------------------------------------------------------------
Projector1D2Order::Projector1D2Order( Params &params, Patch *patch ) : Projector1D( params, patch )
{
    dx_inv_  = 1.0/params.cell_length[0];
    dx_ov_dt = params.cell_length[0] / params.timestep;
    
    index_domain_begin = patch->getCellStartingGlobalIndex( 0 );
    
    dt             = params.timestep;
    dts2           = params.timestep/2.;
    dts4           = params.timestep/4.;
    
}


Projector1D2Order::~Projector1D2Order()
{
}


// ---------------------------------------------------------------------------------------------------------------------
//! Project current densities : main projector
// ---------------------------------------------------------------------------------------------------------------------
void Projector1D2Order::currents( double *Jx, double *Jy, double *Jz, Particles &particles, unsigned int ipart, double invgf, int *iold, double *delta )
{
    // Declare local variables
    int ipo, ip;
    int ip_m_ipo;
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    double xjn, xj_m_xipo, xj_m_xipo2, xj_m_xip, xj_m_xip2;
    double crx_p = charge_weight*dx_ov_dt;                // current density for particle moving in the x-direction
    double cry_p = charge_weight*particles.momentum( 1, ipart )*invgf;  // current density in the y-direction of the macroparticle
    double crz_p = charge_weight*particles.momentum( 2, ipart )*invgf;  // current density allow the y-direction of the macroparticle
    double S0[5], S1[5], Wl[5], Wt[5], Jx_p[5];            // arrays used for the Esirkepov projection method
    
    // Initialize variables
    for( unsigned int i=0; i<5; i++ ) {
        S0[i]=0.;
        S1[i]=0.;
        Wl[i]=0.;
        Wt[i]=0.;
        Jx_p[i]=0.;
    }//i
    
    
    // Locate particle old position on the primal grid
    xj_m_xipo  = *delta;                              // normalized distance to the nearest grid point already stored
    xj_m_xipo2 = xj_m_xipo*xj_m_xipo;                 // square of the normalized distance to the nearest grid point
    
    // Locate particle new position on the primal grid
    xjn       = particles.position( 0, ipart ) * dx_inv_;
    ip        = round( xjn );                         // index of the central node
    xj_m_xip  = xjn - ( double )ip;                   // normalized distance to the nearest grid point
    xj_m_xip2 = xj_m_xip*xj_m_xip;                    // square of the normalized distance to the nearest grid point
    
    
    // coefficients 2nd order interpolation on 3 nodes
    S0[1] = 0.5 * ( xj_m_xipo2-xj_m_xipo+0.25 );
    S0[2] = ( 0.75-xj_m_xipo2 );
    S0[3] = 0.5 * ( xj_m_xipo2+xj_m_xipo+0.25 );
    
    // coefficients 2nd order interpolation on 3 nodes
    ipo        = *iold;                          // index of the central node
    ip_m_ipo = ip-ipo-index_domain_begin;
    S1[ip_m_ipo+1] = 0.5 * ( xj_m_xip2-xj_m_xip+0.25 );
    S1[ip_m_ipo+2] = ( 0.75-xj_m_xip2 );
    S1[ip_m_ipo+3] = 0.5 * ( xj_m_xip2+xj_m_xip+0.25 );
    
    // coefficients used in the Esirkepov method
    for( unsigned int i=0; i<5; i++ ) {
        Wl[i] = S0[i] - S1[i];           // for longitudinal current (x)
        Wt[i] = 0.5 * ( S0[i] + S1[i] ); // for transverse currents (y,z)
    }//i
    
    // local current created by the particle
    // calculate using the charge conservation equation
    for( unsigned int i=1; i<5; i++ ) {
        Jx_p[i] = Jx_p[i-1] + crx_p * Wl[i-1];
    }
    
    
    // 2nd order projection for the total currents & charge density
    ipo -= 2; // At the 2nd order, oversize = 2.
    for( unsigned int i=0; i<5; i++ ) {
        Jx[i + ipo ]  += Jx_p[i];
        Jy[i + ipo ]  += cry_p * Wt[i];
        Jz[i + ipo ]  += crz_p * Wt[i];
    }//i
}


// ---------------------------------------------------------------------------------------------------------------------
//!  Project current densities & charge : diagFields timstep
// ---------------------------------------------------------------------------------------------------------------------
void Projector1D2Order::currentsAndDensity( double *Jx, double *Jy, double *Jz, double *rho, Particles &particles, unsigned int ipart, double invgf, int *iold, double *delta )
{
    // Declare local variables
    int ipo, ip;
    int ip_m_ipo;
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    double xjn, xj_m_xipo, xj_m_xipo2, xj_m_xip, xj_m_xip2;
    double crx_p = charge_weight*dx_ov_dt;                // current density for particle moving in the x-direction
    double cry_p = charge_weight*particles.momentum( 1, ipart )*invgf;  // current density in the y-direction of the macroparticle
    double crz_p = charge_weight*particles.momentum( 2, ipart )*invgf;  // current density allow the y-direction of the macroparticle
    double S0[5], S1[5], Wl[5], Wt[5], Jx_p[5];            // arrays used for the Esirkepov projection method
    
    // Initialize variables
    for( unsigned int i=0; i<5; i++ ) {
        S0[i]=0.;
        S1[i]=0.;
        Wl[i]=0.;
        Wt[i]=0.;
        Jx_p[i]=0.;
    }//i
    
    
    // Locate particle old position on the primal grid
    xj_m_xipo  = *delta;                   // normalized distance to the nearest grid point
    xj_m_xipo2 = xj_m_xipo*xj_m_xipo;                 // square of the normalized distance to the nearest grid point
    
    // Locate particle new position on the primal grid
    xjn       = particles.position( 0, ipart ) * dx_inv_;
    ip        = round( xjn );                         // index of the central node
    xj_m_xip  = xjn - ( double )ip;                   // normalized distance to the nearest grid point
    xj_m_xip2 = xj_m_xip*xj_m_xip;                    // square of the normalized distance to the nearest grid point
    
    
    // coefficients 2nd order interpolation on 3 nodes
    S0[1] = 0.5 * ( xj_m_xipo2-xj_m_xipo+0.25 );
    S0[2] = ( 0.75-xj_m_xipo2 );
    S0[3] = 0.5 * ( xj_m_xipo2+xj_m_xipo+0.25 );
    
    // coefficients 2nd order interpolation on 3 nodes
    ipo = *iold;
    ip_m_ipo = ip-ipo-index_domain_begin;
    S1[ip_m_ipo+1] = 0.5 * ( xj_m_xip2-xj_m_xip+0.25 );
    S1[ip_m_ipo+2] = ( 0.75-xj_m_xip2 );
    S1[ip_m_ipo+3] = 0.5 * ( xj_m_xip2+xj_m_xip+0.25 );
    
    // coefficients used in the Esirkepov method
    for( unsigned int i=0; i<5; i++ ) {
        Wl[i] = S0[i] - S1[i];           // for longitudinal current (x)
        Wt[i] = 0.5 * ( S0[i] + S1[i] ); // for transverse currents (y,z)
    }//i
    
    // local current created by the particle
    // calculate using the charge conservation equation
    for( unsigned int i=1; i<5; i++ ) {
        Jx_p[i] = Jx_p[i-1] + crx_p * Wl[i-1];
    }
    
    
    // 2nd order projection for the total currents & charge density
    ipo -= 2;// At the 2nd order, oversize = 2.
    for( unsigned int i=0; i<5; i++ ) {
        Jx[i + ipo]  += Jx_p[i];
        Jy[i + ipo]  += cry_p * Wt[i];
        Jz[i + ipo]  += crz_p * Wt[i];
        rho[i + ipo] += charge_weight * S1[i];
    }//i
    
    
} // END Project local current densities (sort)


// ---------------------------------------------------------------------------------------------------------------------
//! Project charge : frozen & diagFields timstep
// ---------------------------------------------------------------------------------------------------------------------
void Projector1D2Order::basic( double *rhoj, Particles &particles, unsigned int ipart, unsigned int type )
{

    //Warning : this function is used for frozen species or initialization only and doesn't use the standard scheme.
    //rho type = 0
    //Jx type = 1
    //Jy type = 2
    //Jz type = 3
    
    // The variable bin received is  number of bin * cluster width.
    // Declare local variables
    int ip;
    double xjn, xj_m_xip, xj_m_xip2;
    double S1[5];            // arrays used for the Esirkepov projection method
    
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    if( type > 0 ) {
        charge_weight *= 1./sqrt( 1.0 + particles.momentum( 0, ipart )*particles.momentum( 0, ipart )
                                  + particles.momentum( 1, ipart )*particles.momentum( 1, ipart )
                                  + particles.momentum( 2, ipart )*particles.momentum( 2, ipart ) );
                                  
        if( type == 1 ) {
            charge_weight *= particles.momentum( 0, ipart );
        } else if( type == 2 ) {
            charge_weight *= particles.momentum( 1, ipart );
        } else {
            charge_weight *= particles.momentum( 2, ipart );
        }
    }
    
    // Initialize variables
    for( unsigned int i=0; i<5; i++ ) {
        S1[i]=0.;
    }//i
    
    // Locate particle new position on the primal grid
    xjn       = particles.position( 0, ipart ) * dx_inv_;
    ip        = round( xjn + 0.5 * ( type==1 ) );                       // index of the central node
    xj_m_xip  = xjn - ( double )ip;                   // normalized distance to the nearest grid point
    xj_m_xip2 = xj_m_xip*xj_m_xip;                    // square of the normalized distance to the nearest grid point
    
    // coefficients 2nd order interpolation on 3 nodes
    //ip_m_ipo = ip-ipo;
    S1[1] = 0.5 * ( xj_m_xip2-xj_m_xip+0.25 );
    S1[2] = ( 0.75-xj_m_xip2 );
    S1[3] = 0.5 * ( xj_m_xip2+xj_m_xip+0.25 );
    
    ip -= index_domain_begin + 2;
    
    // 2nd order projection for charge density
    // At the 2nd order, oversize = 2.
    for( unsigned int i=0; i<5; i++ ) {
        rhoj[i + ip ] += charge_weight * S1[i];
    }//i
    
}

// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities : ionization
// ---------------------------------------------------------------------------------------------------------------------
void Projector1D2Order::ionizationCurrents( Field *Jx, Field *Jy, Field *Jz, Particles &particles, int ipart, LocalFields Jion )
{
    Field1D *Jx1D  = static_cast<Field1D *>( Jx );
    Field1D *Jy1D  = static_cast<Field1D *>( Jy );
    Field1D *Jz1D  = static_cast<Field1D *>( Jz );
    
    
    //Declaration of local variables
    int i, im1, ip1;
    double xjn, xjmxi, xjmxi2;
    double cim1, ci, cip1;
    
    // weighted currents
    double weight = inv_cell_volume * particles.weight( ipart );
    double Jx_ion = Jion.x * weight;
    double Jy_ion = Jion.y * weight;
    double Jz_ion = Jion.z * weight;
    
    //Locate particle on the grid
    xjn    = particles.position( 0, ipart ) * dx_inv_; // normalized distance to the first node
    
    
    // Compute Jx_ion on the dual grid
    // -------------------------------
    
    i      = round( xjn+0.5 );             // index of the central node
    xjmxi  = xjn - ( double )i + 0.5;      // normalized distance to the nearest grid point
    xjmxi2 = xjmxi*xjmxi;                  // square of the normalized distance to the nearest grid point
    
    i  -= index_domain_begin;
    im1 = i-1;
    ip1 = i+1;
    
    cim1 = 0.5 * ( xjmxi2-xjmxi+0.25 );
    ci   = ( 0.75-xjmxi2 );
    cip1 = 0.5 * ( xjmxi2+xjmxi+0.25 );
    
    // Jx
    ( *Jx1D )( im1 )  += cim1 * Jx_ion;
    ( *Jx1D )( i )  += ci   * Jx_ion;
    ( *Jx1D )( ip1 )  += cip1 * Jx_ion;
    
    
    // Compute Jy_ion & Jz_ion on the primal grid
    // ------------------------------------------
    
    i      = round( xjn );                 // index of the central node
    xjmxi  = xjn - ( double )i;            // normalized distance to the nearest grid point
    xjmxi2 = xjmxi*xjmxi;                  // square of the normalized distance to the nearest grid point
    
    i  -= index_domain_begin;
    im1 = i-1;
    ip1 = i+1;
    
    cim1 = 0.5 * ( xjmxi2-xjmxi+0.25 );
    ci   = ( 0.75-xjmxi2 );
    cip1 = 0.5 * ( xjmxi2+xjmxi+0.25 );
    
    // Jy
    ( *Jy1D )( im1 )  += cim1 * Jy_ion;
    ( *Jy1D )( i )  += ci   * Jy_ion;
    ( *Jy1D )( ip1 )  += cip1 * Jy_ion;
    
    // Jz
    ( *Jz1D )( im1 )  += cim1 * Jz_ion;
    ( *Jz1D )( i )  += ci   * Jz_ion;
    ( *Jz1D )( ip1 )  += cip1 * Jz_ion;
    
} // END Project global current densities (ionize)

void Projector1D2Order::currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int icell, int ipart_ref )
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
                currents( Jx_, Jy_, Jz_, particles,  ipart, ( *invgf )[ipart], &( *iold )[ipart], &( *delta )[ipart] );
            }
        } else {
            for( int ipart=istart ; ipart<iend; ipart++ ) {
                currentsAndDensity( Jx_, Jy_, Jz_, rho_, particles,  ipart, ( *invgf )[ipart], &( *iold )[ipart], &( *delta )[ipart] );
            }
        }
        // Otherwise, the projection may apply to the species-specific arrays
    } else {
        double *b_Jxs  = EMfields->Jx_s [ispec] ? &( *EMfields->Jx_s [ispec] )( 0 ) : &( *EMfields->Jx_ )( 0 ) ;
        double *b_Jys  = EMfields->Jy_s [ispec] ? &( *EMfields->Jy_s [ispec] )( 0 ) : &( *EMfields->Jy_ )( 0 ) ;
        double *b_Jzs  = EMfields->Jz_s [ispec] ? &( *EMfields->Jz_s [ispec] )( 0 ) : &( *EMfields->Jz_ )( 0 ) ;
        double *b_rhos = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
        for( int ipart=istart ; ipart<iend; ipart++ ) {
            currentsAndDensity( b_Jxs, b_Jys, b_Jzs, b_rhos, particles,  ipart, ( *invgf )[ipart], &( *iold )[ipart], &( *delta )[ipart] );
        }
    }
}

// Projector for susceptibility used as source term in envelope equation
void Projector1D2Order::susceptibility( ElectroMagn *EMfields, Particles &particles, double species_mass, SmileiMPI *smpi, int istart, int iend,  int ithread, int icell, int ipart_ref )

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
        double xpn;
        double delta, delta2;
        double Sx1[5]; // arrays used for the Esirkepov projection method
        
        // Initialize all current-related arrays to zero
        for( unsigned int i=0; i<5; i++ ) {
            Sx1[i] = 0.;
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
        
        // ---------------------------
        // Calculate the total susceptibility
        // ---------------------------
        ip -= index_domain_begin + 2;
        
        for( unsigned int i=0 ; i<5 ; i++ ) {
            iloc = ( i+ip );
            Chi_envelope[iloc] += charge_weight * Sx1[i];
        }//i
        
        
    }
    
}
