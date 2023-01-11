
#include "LaserEnvelope.h"

#include "Params.h"
#include "Patch.h"
#include "cField1D.h"
#include "Field1D.h"
#include "ElectroMagn.h"
#include "Profile.h"
#include "ElectroMagnFactory.h"
#include "EnvelopeBC.h"
#include "EnvelopeBC_Factory.h"
#include <complex>
#include "SimWindow.h"


using namespace std;

LaserEnvelope1D::LaserEnvelope1D( Params &params, Patch *patch )
    : LaserEnvelope( params, patch )
{
    std::vector<unsigned int>  dimPrim( params.nDim_field );
    // Dimension of the primal and dual grids
    for( size_t i=0 ; i<params.nDim_field ; i++ ) {
        // Standard scheme
        dimPrim[i] = params.patch_size_[i]+1;
        // + Ghost domain
        dimPrim[i] += 2*params.oversize[i];
    }
    
    
    A_  = new cField1D( dimPrim, "A" );
    A0_ = new cField1D( dimPrim, "Aold" );
    
    Phi_         = new Field1D( dimPrim, "Phi" );
    Phi_m        = new Field1D( dimPrim, "Phi_m" );
    
    GradPhix_    = new Field1D( dimPrim, "GradPhix" );
    GradPhix_m   = new Field1D( dimPrim, "GradPhix_m" );
    
    GradPhiy_    = new Field1D( dimPrim, "GradPhiy" );
    GradPhiy_m   = new Field1D( dimPrim, "GradPhiy_m" );
    
    GradPhiz_    = new Field1D( dimPrim, "GradPhiz" );
    GradPhiz_m   = new Field1D( dimPrim, "GradPhiz_m" );
    
}


LaserEnvelope1D::LaserEnvelope1D( LaserEnvelope *envelope, Patch *patch, Params &params, unsigned int n_moved )
    : LaserEnvelope( envelope, patch, params, n_moved )
{
    A_           = new cField1D( envelope->A_->dims_, "A" );
    A0_          = new cField1D( envelope->A0_->dims_, "Aold" );
    
    Phi_         = new Field1D( envelope->Phi_->dims_ );
    Phi_m        = new Field1D( envelope->Phi_m->dims_ );
    
    GradPhix_    = new Field1D( envelope->GradPhix_->dims_ );
    GradPhix_m   = new Field1D( envelope->GradPhix_m->dims_ );
    
    GradPhiy_    = new Field1D( envelope->GradPhiy_->dims_ );
    GradPhiy_m   = new Field1D( envelope->GradPhiy_m->dims_ );
    
    GradPhiz_    = new Field1D( envelope->GradPhiz_->dims_ );
    GradPhiz_m   = new Field1D( envelope->GradPhiz_m->dims_ );
    
}


void LaserEnvelope1D::initEnvelope( Patch *patch, ElectroMagn *EMfields )
{
    cField1D *A1D          = static_cast<cField1D *>( A_ );
    cField1D *A01D         = static_cast<cField1D *>( A0_ );
    Field1D *Env_Aabs1D    = static_cast<Field1D *>( EMfields->Env_A_abs_ );
    Field1D *Env_Eabs1D    = static_cast<Field1D *>( EMfields->Env_E_abs_ );
    //Field1D *Env_Exabs1D   = static_cast<Field1D *>( EMfields->Env_Ex_abs_ );
    
    Field1D *Phi1D         = static_cast<Field1D *>( Phi_ );
    Field1D *Phi_m1D       = static_cast<Field1D *>( Phi_m );
    
    Field1D *GradPhix1D    = static_cast<Field1D *>( GradPhix_ );
    Field1D *GradPhix_m1D  = static_cast<Field1D *>( GradPhix_m );
    
    
    vector<double> position( 1, 0 );
    double t;
    double t_previous_timestep;
    
    
    // position[0]: x coordinate
    // t: time coordinate --> x/c for the envelope initialization
    
    position[0]           = cell_length[0]*( ( double )( patch->getCellStartingGlobalIndex( 0 ) )+( A1D->isDual( 0 )?-0.5:0. ) );
    t                     = position[0];          // x-ct     , t=0
    t_previous_timestep   = position[0]+timestep; // x-c(t-dt), t=0
    
    
    // UNSIGNED INT LEADS TO PB IN PERIODIC BCs
    for( unsigned int i=0 ; i<A_->dims_[0] ; i++ ) { // x loop
    
        // init envelope through Python function
        ( *A1D )( i )      += profile_->complexValueAt( position, t );
        ( *A01D )( i )     += profile_->complexValueAt( position, t_previous_timestep );
        
        // |A|
        ( *Env_Aabs1D )( i )= std::abs( ( *A1D )( i ) );
        // |E envelope| = |-(dA/dt-ik0cA)|
        ( *Env_Eabs1D )( i )= std::abs( ( ( *A1D )( i )-( *A01D )( i ) )/timestep - i1*omega*( *A1D )( i ) );
        // |Ex envelope| = 0 in 1D
        ( *Env_Eabs1D )( i )= 0;
        // compute ponderomotive potential at timestep n
        ( *Phi1D )( i )     = ellipticity_factor*std::abs( ( *A1D )( i ) ) * std::abs( ( *A1D )( i ) ) * 0.5;
        // compute ponderomotive potential at timestep n-1
        ( *Phi_m1D )( i )   = ellipticity_factor*std::abs( ( *A01D )( i ) ) * std::abs( ( *A01D )( i ) ) * 0.5;
        // interpolate in time
        ( *Phi_m1D )( i )   = 0.5*( ( *Phi_m1D )( i )+( *Phi1D )( i ) );
        
        position[0]          += cell_length[0];
        t                     = position[0];
        t_previous_timestep   = position[0]+timestep;
    } // end x loop
    
    // Compute gradient of ponderomotive potential
    for( unsigned int i=1 ; i<A_->dims_[0]-1 ; i++ ) { // x loop
        // gradient in x direction
        ( *GradPhix1D )( i ) = ( ( *Phi1D )( i+1 )-( *Phi1D )( i-1 ) ) * one_ov_2dx;
        ( *GradPhix_m1D )( i ) = ( ( *Phi_m1D )( i+1 )-( *Phi_m1D )( i-1 ) ) * one_ov_2dx;
    } // end x loop
    
}


LaserEnvelope1D::~LaserEnvelope1D()
{
}

void LaserEnvelope1D::updateEnvelope( Patch *patch )
{
    //// solves envelope equation in lab frame (see doc):
    // full_laplacian(A)+2ik0*(dA/dz+(1/c)*dA/dt)-d^2A/dt^2*(1/c^2)=Chi*A
    // where Chi is the plasma susceptibility [= sum(q^2*rho/mass/gamma_ponderomotive) for all species]
    // gamma_ponderomotive=sqrt(1+p^2+|A|^2/2) in normalized units
    
    // the following explicit finite difference scheme is obtained through centered finite difference derivatives
    // e.g. (dA/dx) @ time n and indices i = (A^n    _{i+1} - A^n    _{i-1}) /2/dx
    //      (dA/dt) @ time n and indices i = (A^{n+1}_{i  } - A^{n-1}_{i  }) /2/dt
    // A0 is A^{n-1}
    //      (d^2A/dx^2) @ time n and indices i = (A^{n}_{i+1,j,k}-2*A^{n}_{i}+A^{n}_{i-1})/dx^2
    
    cField1D *A1D          = static_cast<cField1D *>( A_ );               // the envelope at timestep n
    cField1D *A01D         = static_cast<cField1D *>( A0_ );              // the envelope at timestep n-1
    Field1D *Env_Chi1D     = static_cast<Field1D *>( patch->EMfields->Env_Chi_ ); // source term of envelope equation
    
    
    // temporary variable for updated envelope
    cField1D *A1Dnew;
    A1Dnew  = new cField1D( A_->dims_ );
    
    //// explicit solver
    for( unsigned int i=1 ; i <A_->dims_[0]-1; i++ ) { // x loop
        ( *A1Dnew )( i ) -= ( *Env_Chi1D )( i )*( *A1D )( i ); // subtract here source term Chi*A from plasma
        // A1Dnew = laplacian - source term
        ( *A1Dnew )( i ) += ( ( *A1D )( i-1 )-2.*( *A1D )( i )+( *A1D )( i+1 ) )*one_ov_dx_sq; // x part
        
        // A1Dnew = A1Dnew+2ik0*dA/dx
        ( *A1Dnew )( i ) += i1_2k0_over_2dx*( ( *A1D )( i+1 )-( *A1D )( i-1 ) );
        // A1Dnew = A1Dnew*dt^2
        ( *A1Dnew )( i )  = ( *A1Dnew )( i )*dt_sq;
        // A1Dnew = A1Dnew + 2/c^2 A1D - (1+ik0cdt)A01D/c^2
        ( *A1Dnew )( i ) += 2.*( *A1D )( i )-one_plus_ik0dt*( *A01D )( i );
        // A1Dnew = A1Dnew * (1+ik0dct)/(1+k0^2c^1Dt^2)
        ( *A1Dnew )( i )  = ( *A1Dnew )( i )*one_plus_ik0dt_ov_one_plus_k0sq_dtsq;
    } // end x loop
    
    for( unsigned int i=1 ; i <A_->dims_[0]-1; i++ ) { // x loop
    
        // final back-substitution
        ( *A01D )( i )       = ( *A1D )( i );
        ( *A1D )( i )        = ( *A1Dnew )( i );
        
    } // end x loop
    
    delete A1Dnew;
} // end LaserEnvelope1D::updateEnvelope

void LaserEnvelope1D::updateEnvelopeReducedDispersion( Patch *patch )
{
    //// solves envelope equation in lab frame (see doc):
    // full_laplacian(A)+2ik0*(dA/dz+(1/c)*dA/dt)-d^2A/dt^2*(1/c^2)=Chi*A
    // where Chi is the plasma susceptibility [= sum(q^2*rho/mass/gamma_ponderomotive) for all species]
    // gamma_ponderomotive=sqrt(1+p^2+|A|^2/2) in normalized units
    
    // the following explicit finite difference scheme is obtained through centered finite difference derivatives
    // e.g. (dA/dx) @ time n and indices i = (A^n    _{i+1} - A^n    _{i-1}) /2/dx
    //      (dA/dt) @ time n and indices i = (A^{n+1}_{i  } - A^{n-1}_{i  }) /2/dt
    // A0 is A^{n-1}
    //      (d^2A/dx^2) @ time n and indices i = (A^{n}_{i+1,j,k}-2*A^{n}_{i}+A^{n}_{i-1})/dx^2

    // An optimized form for the derivatives along x has been proposed in D. Terzani, P. Londrillo, JCP 2019
    // to reduce the numerical dispersion for the envelope solver.
    // The derivatives along x of the reduced dispersion scheme are defined as follows:
    // delta= [1-(dt/dx)^2]/3,
    // (dA/dx)_opt = (1+delta)*(dA/dx) - delta*(A_{i+2,j,k}-A_{i-2,j,k})/4/dx
    // (d^2A/dx^2)_opt = (1+delta)*(d^2A/dx^2) - delta*(A_{i+2,j,k}-2*A_{i,j,k}+A_{i-2,j,k})/(4dx^2)
    
    
    cField1D *A1D          = static_cast<cField1D *>( A_ );               // the envelope at timestep n
    cField1D *A01D         = static_cast<cField1D *>( A0_ );              // the envelope at timestep n-1
    Field1D *Env_Chi1D     = static_cast<Field1D *>( patch->EMfields->Env_Chi_ ); // source term of envelope equation
    
    
    // temporary variable for updated envelope
    cField1D *A1Dnew;
    A1Dnew  = new cField1D( A_->dims_ );
    
    //// explicit solver
    for( unsigned int i=2 ; i <A_->dims_[0]-2; i++ ) { // x loop
        ( *A1Dnew )( i ) -= ( *Env_Chi1D )( i )*( *A1D )( i ); // subtract here source term Chi*A from plasma
        // A1Dnew = laplacian - source term
        ( *A1Dnew )( i ) += (1.+delta)*( ( *A1D )( i-1 ) -2.*( *A1D )( i ) +( *A1D )( i+1 )   )*one_ov_dx_sq; // x part with optimized derivative
        ( *A1Dnew )( i ) -= delta*     ( ( *A1D )( i-2 ) -2.*( *A1D )( i ) +( *A1D )( i+2 )   )*0.25*one_ov_dx_sq;

        // A1Dnew = A1Dnew+2ik0*dA/dx, where dA/dx uses the optimized form
        ( *A1Dnew )( i ) += i1_2k0_over_2dx*(1.+delta)*( ( *A1D )( i+1 )-( *A1D )( i-1 ) );
        ( *A1Dnew )( i ) -= i1_2k0_over_2dx*delta*0.5 *( ( *A1D )( i+2 )-( *A1D )( i-2 ) );

        // A1Dnew = A1Dnew*dt^2
        ( *A1Dnew )( i )  = ( *A1Dnew )( i )*dt_sq;
        // A1Dnew = A1Dnew + 2/c^2 A1D - (1+ik0cdt)A01D/c^2
        ( *A1Dnew )( i ) += 2.*( *A1D )( i )-one_plus_ik0dt*( *A01D )( i );
        // A1Dnew = A1Dnew * (1+ik0dct)/(1+k0^2c^1Dt^2)
        ( *A1Dnew )( i )  = ( *A1Dnew )( i )*one_plus_ik0dt_ov_one_plus_k0sq_dtsq;
    } // end x loop
    
    for( unsigned int i=2 ; i <A_->dims_[0]-2; i++ ) { // x loop
        // final back-substitution
        ( *A01D )( i )       = ( *A1D )( i );
        ( *A1D )( i )        = ( *A1Dnew )( i );      
    } // end x loop
    
    delete A1Dnew;
} // end LaserEnvelope1D::updateEnvelopeReducedDispersion

void LaserEnvelope1D::computePhiEnvAEnvE( ElectroMagn *EMfields )
{

    // computes Phi=|A|^2/2 (the ponderomotive potential), new values immediately after the envelope update
    cField1D *A1D          = static_cast<cField1D *>( A_ );                   // the envelope at timestep n
    cField1D *A01D         = static_cast<cField1D *>( A0_ );                  // the envelope at timestep n-1
    Field1D *Phi1D         = static_cast<Field1D *>( Phi_ );                  //Phi=|A|^2/2 is the ponderomotive potential
    Field1D *Env_Aabs1D    = static_cast<Field1D *>( EMfields->Env_A_abs_ );  // field for diagnostic and ionization
    Field1D *Env_Eabs1D    = static_cast<Field1D *>( EMfields->Env_E_abs_ );  // field for diagnostic and ionization
    //Field1D *Env_Exabs1D   = static_cast<Field1D *>( EMfields->Env_Ex_abs_ ); // field for diagnostic and ionization
    
    // Compute ponderomotive potential Phi=|A|^2/2, at timesteps n+1, including ghost cells
    for( unsigned int i=0 ; i <A_->dims_[0]-1; i++ ) { // x loop
        ( *Phi1D )( i )      = ellipticity_factor*std::abs( ( *A1D )( i ) ) * std::abs( ( *A1D )( i ) ) * 0.5;
        ( *Env_Aabs1D )( i ) = std::abs( ( *A1D )( i ) );
        // |E envelope| = |-(dA/dt-ik0cA)|, forward finite difference for the time derivative
        ( *Env_Eabs1D )( i ) = std::abs( ( ( *A1D )( i )-( *A01D )( i ) )/timestep - i1*omega*( *A1D )( i ) );
    } // end x loop
    
} // end LaserEnvelope1D::computePhiEnvAEnvE


void LaserEnvelope1D::computeGradientPhi( ElectroMagn * )
{

    // computes gradient of Phi=|A|^2/2 (the ponderomotive potential), new values immediately after the envelope update
    Field1D *GradPhix1D    = static_cast<Field1D *>( GradPhix_ );
    Field1D *Phi1D         = static_cast<Field1D *>( Phi_ );      //Phi=|A|^2/2 is the ponderomotive potential
    
    
    // Compute gradients of Phi, at timesteps n
    for( unsigned int i=1 ; i <A_->dims_[0]-1; i++ ) { // x loop
    
        // gradient in x direction
        ( *GradPhix1D )( i ) = ( ( *Phi1D )( i+1 )-( *Phi1D )( i-1 ) ) * one_ov_2dx;
        
    } // end x loop
    
} // end LaserEnvelope1D::computeGradientPhi


void LaserEnvelope1D::savePhiAndGradPhi()
{
    // Static cast of the fields
    Field1D *Phi1D         = static_cast<Field1D *>( Phi_ );
    Field1D *Phi_m1D       = static_cast<Field1D *>( Phi_m );
    
    Field1D *GradPhix1D    = static_cast<Field1D *>( GradPhix_ );
    Field1D *GradPhix_m1D  = static_cast<Field1D *>( GradPhix_m );
    
    for( unsigned int i=0 ; i <A_->dims_[0]-1; i++ ) { // x loop
    
        // ponderomotive potential Phi=ellipticity_factor*|A|^2/2
        ( *Phi_m1D )( i )       = ( *Phi1D )( i );
        // gradient of ponderomotive potential
        ( *GradPhix_m1D )( i )  = ( *GradPhix1D )( i );
        
    } // end x loop
    
    
}//END savePhiAndGradPhi


void LaserEnvelope1D::centerPhiAndGradPhi()
{
    // Static cast of the fields
    Field1D *Phi1D         = static_cast<Field1D *>( Phi_ );
    Field1D *Phi_m1D       = static_cast<Field1D *>( Phi_m );
    
    Field1D *GradPhix1D    = static_cast<Field1D *>( GradPhix_ );
    Field1D *GradPhix_m1D  = static_cast<Field1D *>( GradPhix_m );
    
    // Phi_m and GradPhi_m quantities now contain values at timestep n
    
    for( unsigned int i=0 ; i <A_->dims_[0]-1; i++ ) { // x loop
    
        // ponderomotive potential Phi=|A|^2/2
        ( *Phi_m1D )( i )       = 0.5*( ( *Phi_m1D )( i )+( *Phi1D )( i ) );
        
        // gradient of ponderomotive potential
        ( *GradPhix_m1D )( i )  = 0.5*( ( *GradPhix_m1D )( i )+( *GradPhix1D )( i ) );
        
    } // end x loop
    
    // Phi_m and GradPhi_m quantities now contain values interpolated at timestep n+1/2
    // these are used for the ponderomotive position advance
    
    
}//END centerPhiAndGradPhi


