
#include "LaserEnvelope.h"

#include "Params.h"
#include "Patch.h"
#include "cField2D.h"
#include "Field2D.h"
#include "ElectroMagn.h"
#include "Profile.h"
#include "ElectroMagnFactory.h"
#include "EnvelopeBC.h"
#include "EnvelopeBC_Factory.h"
#include <complex>
#include "SimWindow.h"


using namespace std;

LaserEnvelopeAM::LaserEnvelopeAM( Params &params, Patch *patch, ElectroMagn *EMfields )
    : LaserEnvelope( params, patch, EMfields )
{
    std::vector<unsigned int>  dimPrim( params.nDim_field );
    // Dimension of the primal and dual grids
    for( size_t i=0 ; i<params.nDim_field ; i++ ) {
        // Standard scheme
        dimPrim[i] = params.n_space[i]+1;
        // + Ghost domain
        dimPrim[i] += 2*params.oversize[i];
    }
    
    
    A_  = new cField2D( dimPrim, "A" );
    A0_ = new cField2D( dimPrim, "Aold" );
    
    Phi_         = new Field2D( dimPrim, "Phi" );
    Phi_m        = new Field2D( dimPrim, "Phi_m" );
    
    GradPhix_    = new Field2D( dimPrim, "GradPhix" );
    GradPhix_m   = new Field2D( dimPrim, "GradPhix_m" );
    
    GradPhir_    = new Field2D( dimPrim, "GradPhir" );
    GradPhir_m   = new Field2D( dimPrim, "GradPhir_m" );

    // in cylindrical symmetry the gradient along theta is zero
    
    
}


LaserEnvelopeAM::LaserEnvelopeAM( LaserEnvelope *envelope, Patch *patch, ElectroMagn *EMfields, Params &params, unsigned int n_moved )
    : LaserEnvelope( envelope, patch, EMfields, params, n_moved )
{
    A_           = new cField2D( envelope->A_->dims_, "A" );
    A0_          = new cField2D( envelope->A0_->dims_, "Aold" );
    
    Phi_         = new Field2D( envelope->Phi_->dims_ );
    Phi_m        = new Field2D( envelope->Phi_m->dims_ );
    
    GradPhix_    = new Field2D( envelope->GradPhix_->dims_ );
    GradPhir_m   = new Field2D( envelope->GradPhix_m->dims_ );
    
    GradPhir_    = new Field2D( envelope->GradPhir_->dims_ );
    GradPhir_m   = new Field2D( envelope->GradPhir_m->dims_ );
 
    // in cylindrical symmetry the gradient along theta is zero
    
    
}


void LaserEnvelopeAM::initEnvelope( Patch *patch, ElectroMagn *EMfields )
{
    cField2D *A2Dcyl          = static_cast<cField2D *>( A_ );
    cField2D *A02Dcyl         = static_cast<cField2D *>( A0_ );
    Field2D *Env_Aabs2Dcyl    = static_cast<Field2D *>( EMfields->Env_A_abs_ );
    Field2D *Env_Eabs2Dcyl    = static_cast<Field2D *>( EMfields->Env_E_abs_ );
    
    Field2D *Phi2Dcyl         = static_cast<Field2D *>( Phi_ );
    Field2D *Phi_m2Dcyl       = static_cast<Field2D *>( Phi_m );
    
    Field2D *GradPhix2Dcyl    = static_cast<Field2D *>( GradPhix_ );
    Field2D *GradPhix_m2Dcyl  = static_cast<Field2D *>( GradPhix_m );
    
    Field2D *GradPhir2Dcyl    = static_cast<Field2D *>( GradPhir_ );
    Field2D *GradPhir_m2Dcyl  = static_cast<Field2D *>( GradPhir_m );
    
    
    vector<double> position( 2, 0 );
    double t;
    double t_previous_timestep;
    
    complex<double>     i1 = std::complex<double>( 0., 1 );
    
    //! 1/(2dx), where dx is the spatial step dx for 2D3V cylindrical simulations
    double one_ov_2dx=1./2./cell_length[0];
    //! 1/(2dr), where dr is the spatial step dr for 2D3V cylindrical simulations
    double one_ov_2dr=1./2./cell_length[1];
    
    // position[0]: x coordinate
    // position[1]: r coordinate
    // t: time coordinate --> x/c for the envelope initialization
    
    position[0]           = cell_length[0]*( ( double )( patch->getCellStartingGlobalIndex( 0 ) )+( A2Dcyl->isDual( 0 )?-0.5:0. ) );
    t                     = position[0];          // x-ct     , t=0
    t_previous_timestep   = position[0]+timestep; // x-c(t-dt), t=0
    double pos1 = cell_length[1]*( ( double )( patch->getCellStartingGlobalIndex( 1 ) )+( A2Dcyl->isDual( 1 )?-0.5:0. ) );
    
    // UNSIGNED INT LEADS TO PB IN PERIODIC BCs
    for( unsigned int i=0 ; i<A_->dims_[0] ; i++ ) { // x loop
        position[1] = pos1;
        for( unsigned int j=0 ; j<A_->dims_[1] ; j++ ) { // r loop
            // init envelope through Python function
            ( *A2Dcyl )( i, j )      += profile_->complexValueAt( position, t );
            ( *A02Dcyl )( i, j )     += profile_->complexValueAt( position, t_previous_timestep );
            
            // |A|
            ( *Env_Eabs2Dcyl )( i, j )= std::abs( ( *A2Dcyl )( i, j ) );
            // |E envelope| = |-(dA/dt-ik0cA)|
            ( *Env_Eabs2Dcyl )( i, j )= std::abs( ( ( *A2Dcyl )( i, j )-( *A02Dcyl )( i, j ) )/timestep - i1*( *A2Dcyl )( i, j ) );
            // compute ponderomotive potential at timestep n
            ( *Phi2Dcyl )( i, j )     = std::abs( ( *A2Dcyl )( i, j ) ) * std::abs( ( *A2Dcyl )( i, j ) ) * 0.5;
            // compute ponderomotive potential at timestep n-1
            ( *Phi_m2Dcyl )( i, j )   = std::abs( ( *A02Dcyl )( i, j ) ) * std::abs( ( *A02Dcyl )( i, j ) ) * 0.5;
            // interpolate in time
            ( *Phi_m2Dcyl )( i, j )   = 0.5*( ( *Phi_m2Dcyl )( i, j )+( *Phi2Dcyl )( i, j ) );
            
            position[1] += cell_length[1];
        } // end y loop
        position[0]          += cell_length[0];
        t                     = position[0];
        t_previous_timestep   = position[0]+timestep;
    } // end x loop
    
    // Compute gradient of ponderomotive potential
    for( unsigned int i=1 ; i<A_->dims_[0]-1 ; i++ ) { // x loop
        for( unsigned int j=1 ; j<A_->dims_[1]-1 ; j++ ) { // r loop
            // gradient in x direction
            ( *GradPhix2Dcyl )( i, j ) = ( ( *Phi2Dcyl )( i+1, j )-( *Phi2D )( i-1, j ) ) * one_ov_2dx;
            ( *GradPhix_m2Dcyl )( i, j ) = ( ( *Phi_m2D )( i+1, j )-( *Phi_m2D )( i-1, j ) ) * one_ov_2dx;
            // gradient in thet direction
            ( *GradPhir2Dcyl )( i, j ) = ( ( *Phi2D )( i, j+1 )-( *Phi2D )( i, j-1 ) ) * one_ov_2dr;
            ( *GradPhir_m2Dcyl )( i, j ) = ( ( *Phi_m2D )( i, j+1 )-( *Phi_m2D )( i, j-1 ) ) * one_ov_2dr;

            // in cylindrical symmetry the gradient along theta is zero

        } // end r loop
    } // end x loop
    
}


LaserEnvelopeAM::~LaserEnvelopeAM()
{
}

void LaserEnvelopeAM::compute( ElectroMagn *EMfields )
{
    //// solves envelope equation in lab frame (see doc):
    // full_laplacian(A)+2ik0*(dA/dz+(1/c)*dA/dt)-d^2A/dt^2*(1/c^2)=Chi*A
    // where Chi is the plasma susceptibility [= sum(q^2*rho/mass/gamma_ponderomotive) for all species]
    // gamma_ponderomotive=sqrt(1+p^2+|A|^2/2) in normalized units
    
    // For an envelope moving from right to left, replace the imaginary unit i with its opposite (-i)
    // if using an envelope moving to the left, change the sign of the phase in the envelope initialization
    
    // the following explicit finite difference scheme is obtained through centered finite difference derivatives
    // e.g. (dA/dx) @ time n and indices ijk = (A^n    _{i+1,j,k} - A^n    _{i-1,j,k}) /2/dx
    //      (dA/dt) @ time n and indices ijk = (A^{n+1}_{i  ,j,k} - A^{n-1}_{i  ,j,k}) /2/dt
    // A0 is A^{n-1}
    //      (d^2A/dx^2) @ time n and indices ijk = (A^{n}_{i+1,j,k}-2*A^{n}_{i,j,k}+A^{n}_{i-1,j,k})/dx^2
    
    
    
    //// auxiliary quantities
    
    //! 1/dt^2, where dt is the temporal step
    double           dt_sq = timestep*timestep;
    // imaginary unit
    complex<double>     i1 = std::complex<double>( 0., 1 );
    
    //! 1/dx^2, 1/dy^2, 1/dz^2, where dx,dy,dz are the spatial step dx for 2D3V cartesian simulations
    double one_ov_dx_sq    = 1./cell_length[0]/cell_length[0];
    double one_ov_dy_sq    = 1./cell_length[1]/cell_length[1];
    
    
    cField2D *A2Dcyl       = static_cast<cField2D *>( A_ );               // the envelope at timestep n
    cField2D *A02Dcyl      = static_cast<cField2D *>( A0_ );              // the envelope at timestep n-1
    Field2D *Env_Chi2Dcyl  = static_cast<Field2D *>( EMfields->Env_Chi_ ); // source term of envelope equation
    Field2D *Env_Aabs2Dcyl = static_cast<Field2D *>( EMfields->Env_A_abs_ ); // field for diagnostic
    Field2D *Env_Eabs2Dcyl = static_cast<Field2D *>( EMfields->Env_E_abs_ ); // field for diagnostic
    
    
    //! 1/(2dx), where dx is the spatial step dx for 2D3V cartesian simulations
    double one_ov_2dt      = 1./2./timestep;
    
    // temporary variable for updated envelope
    cField2D *A2Dcylnew;
    A2Dcylnew  = new cField2D( A_->dims_ );
    
    //// explicit solver
    for( unsigned int i=1 ; i <A_->dims_[0]-1; i++ ) { // x loop
        for( unsigned int j=1 ; j < A_->dims_[1]-1 ; j++ ) { // r loop
            ( *A2Dcylnew )( i, j ) -= ( *Env_Chi2Dcyl )( i, j )*( *A2Dcyl )( i, j ); // subtract here source term Chi*A from plasma
            // A2Dcylnew = laplacian - source term
            ( *A2Dcylnew )( i, j ) += ( ( *A2Dcyl )( i-1, j )-2.*( *A2Dcyl )( i, j )+( *A2Dcyl )( i+1, j ) )*one_ov_dx_sq; // x part
            ( *A2Dcylnew )( i, j ) += ( ( *A2Dcyl )( i, j-1 )-2.*( *A2Dcyl )( i, j )+( *A2Dcyl )( i, j+1 ) )*one_ov_dr_sq; // r part
            
            // A2Dcylnew = A2Dcylnew+2ik0*dA/dx
            ( *A2Dcylnew )( i, j ) += i1_2k0_over_2dx*( ( *A2Dcyl )( i+1, j )-( *A2Dcyl )( i-1, j ) );
            // A2Dcylnew = A2Dcylnew*dt^2
            ( *A2Dcylnew )( i, j )  = ( *A2Dcylnew )( i, j )*dt_sq;
            // A2Dcylnew = A2Dcylnew + 2/c^2 A2Dcyl - (1+ik0cdt)A02Dcyl/c^2
            ( *A2Dcylnew )( i, j ) += 2.*( *A2Dcyl )( i, j )-one_plus_ik0dt*( *A02Dcyl )( i, j );
            // A2Dcylnew = A2Dcylnew * (1+ik0dct)/(1+k0^2c^2dt^2)
            ( *A2Dcylnew )( i, j )  = ( *A2Dcylnew )( i, j )*one_plus_ik0dt_ov_one_plus_k0sq_dtsq;
        } // end y loop
    } // end x loop
    
    for( unsigned int i=1 ; i <A_->dims_[0]-1; i++ ) { // x loop
        for( unsigned int j=1 ; j < A_->dims_[1]-1 ; j++ ) { // y loop
        
            // final back-substitution
            // |E envelope| = |-(dA/dt-ik0cA)|
            ( *Env_Eabs2Dcyl )( i, j ) = std::abs( ( ( *A2Dcylnew )( i, j )-( *A02Dcyl )( i, j ) )*one_ov_2dt - i1*( *A2Dcyl )( i, j ) );
            ( *A02Dcyl )( i, j )       = ( *A2Dcyl )( i, j );
            ( *A2Dcyl )( i, j )        = ( *A2Dcylnew )( i, j );
            ( *Env_Aabs2D )( i, j ) = std::abs( ( *A2D )( i, j ) );
            
        } // end y loop
    } // end x loop
    
    delete A2Dcylnew;
} // end LaserEnvelopeAM::compute


void LaserEnvelopeAM::compute_Phi( ElectroMagn *EMfields )
{

    // computes Phi=|A|^2/2 (the ponderomotive potential), new values immediately after the envelope update
    
    cField2D *A2Dcyl          = static_cast<cField2D *>( A_ );       // the envelope at timestep n
    
    Field2D *Phi2Dcyl         = static_cast<Field2D *>( Phi_ );      //Phi=|A|^2/2 is the ponderomotive potential
    
    
    
    // Compute ponderomotive potential Phi=|A|^2/2, at timesteps n+1, including ghost cells
    for( unsigned int i=1 ; i <A_->dims_[0]-1; i++ ) { // x loop
        for( unsigned int j=1 ; j < A_->dims_[1]-1; j++ ) { // r loop
            ( *Phi2Dcyl )( i, j )       = std::abs( ( *A2Dcyl )( i, j ) ) * std::abs( ( *A2Dcyl )( i, j ) ) * 0.5;
        } // end y loop
    } // end x loop
    
} // end LaserEnvelopeAM::compute_Phi


void LaserEnvelopeAM::compute_gradient_Phi( ElectroMagn *EMfields )
{

    // computes gradient of Phi=|A|^2/2 (the ponderomotive potential), new values immediately after the envelope update
    Field2D *GradPhix2Dcyl    = static_cast<Field2D *>( GradPhix_ );
    Field2D *GradPhir2Dcyl    = static_cast<Field2D *>( GradPhir_ );
    Field2D *Phi2Dcyl         = static_cast<Field2D *>( Phi_ );      //Phi=|A|^2/2 is the ponderomotive potential
    
    
    //! 1/(2dx), where dx is the spatial step dx for 2D3V cylindrical simulations
    double one_ov_2dx=1./2./cell_length[0];
    //! 1/(2dr), where dy is the spatial step dy for 2D3V cylindrical simulations
    double one_ov_2dr=1./2./cell_length[1];
    
    
    // Compute gradients of Phi, at timesteps n
    for( unsigned int i=1 ; i <A_->dims_[0]-1; i++ ) { // x loop
        for( unsigned int j=1 ; j < A_->dims_[1]-1 ; j++ ) { // y loop
        
            // gradient in x direction
            ( *GradPhix2Dcyl )( i, j ) = ( ( *Phi2Dcyl )( i+1, j )-( *Phi2Dcyl )( i-1, j ) ) * one_ov_2dx;
            // gradient in y direction
            ( *GradPhir2Dcyl )( i, j ) = ( ( *Phi2Dcyl )( i, j+1 )-( *Phi2Dcyl )( i, j-1 ) ) * one_ov_2dr;
            
        } // end r loop
    } // end x loop
    
} // end LaserEnvelopeAM::compute_gradient_Phi


void LaserEnvelopeAM::savePhi_and_GradPhi()
{
    // Static cast of the fields
    Field2D *Phi2Dcyl         = static_cast<Field2D *>( Phi_ );
    Field2D *Phi_m2Dcyl       = static_cast<Field2D *>( Phi_m );
    
    Field2D *GradPhix2Dcyl    = static_cast<Field2D *>( GradPhix_ );
    Field2D *GradPhix_m2Dcyl  = static_cast<Field2D *>( GradPhix_m );
    
    Field2D *GradPhiy2Dcyl    = static_cast<Field2D *>( GradPhir_ );
    Field2D *GradPhiy_m2Dcyl  = static_cast<Field2D *>( GradPhir_m );
    
    for( unsigned int i=0 ; i <A_->dims_[0]-1; i++ ) { // x loop
        for( unsigned int j=0 ; j < A_->dims_[1]-1 ; j++ ) { // r loop
        
            // ponderomotive potential Phi=|A|^2/2
            ( *Phi_m2Dcyl )( i, j )       = ( *Phi2Dcyl )( i, j );
            
            // gradient of ponderomotive potential
            ( *GradPhix_m2Dcyl )( i, j )  = ( *GradPhix2Dcyl )( i, j );
            ( *GradPhiy_m2Dcyl )( i, j )  = ( *GradPhir2Dcyl )( i, j );
            
        } // end r loop
    } // end x loop
    
    
}//END savePhi_and_GradPhi


void LaserEnvelopeAM::centerPhi_and_GradPhi()
{
    // Static cast of the fields
    Field2D *Phi2Dcyl         = static_cast<Field2D *>( Phi_ );
    Field2D *Phi_m2Dcyl       = static_cast<Field2D *>( Phi_m );
    
    Field2D *GradPhix2Dcyl    = static_cast<Field2D *>( GradPhix_ );
    Field2D *GradPhix_m2Dcyl  = static_cast<Field2D *>( GradPhix_m );
    
    Field2D *GradPhiy2Dcyl    = static_cast<Field2D *>( GradPhiy_ );
    Field2D *GradPhiy_m2Dcyl  = static_cast<Field2D *>( GradPhiy_m );
    
    // Phi_m and GradPhi_m quantities now contain values at timestep n
    
    for( unsigned int i=0 ; i <A_->dims_[0]-1; i++ ) { // x loop
        for( unsigned int j=0 ; j < A_->dims_[1]-1 ; j++ ) { // r loop
        
            // ponderomotive potential Phi=|A|^2/2
            ( *Phi_m2Dcyl )( i, j )       = 0.5*( ( *Phi_m2Dcyl )( i, j )+( *Phi2Dcyl )( i, j ) );
            
            // gradient of ponderomotive potential
            ( *GradPhix_m2Dcyl )( i, j )  = 0.5*( ( *GradPhix_m2Dcyl )( i, j )+( *GradPhix2Dcyl )( i, j ) );
            ( *GradPhiy_m2Dcyl )( i, j )  = 0.5*( ( *GradPhiy_m2Dcyl )( i, j )+( *GradPhiy2Dcyl )( i, j ) );
            
        } // end r loop
    } // end x loop
    
    // Phi_m and GradPhi_m quantities now contain values interpolated at timestep n+1/2
    // these are used for the ponderomotive position advance
    
    
}//END centerPhi_and_GradPhi


