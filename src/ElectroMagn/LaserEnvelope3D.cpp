
#include "LaserEnvelope.h"

#include "Params.h"
#include "Patch.h"
#include "cField3D.h"
#include "Field3D.h"
#include "ElectroMagn.h"
#include "Profile.h"
#include "ElectroMagnFactory.h"
#include "EnvelopeBC.h"
#include "EnvelopeBC_Factory.h"
#include <complex>
#include "SimWindow.h"


using namespace std;


LaserEnvelope3D::LaserEnvelope3D( Params &params, Patch *patch, ElectroMagn *EMfields )
    : LaserEnvelope( params, patch, EMfields )
{

    one_ov_dy_sq    = 1./cell_length[1]/cell_length[1];
    one_ov_2dy      = 1./2./cell_length[1];
    one_ov_dz_sq    = 1./cell_length[2]/cell_length[2];
    one_ov_2dz      = 1./2./cell_length[2];

    std::vector<unsigned int>  dimPrim( params.nDim_field );
    // Dimension of the primal and dual grids
    for( size_t i=0 ; i<params.nDim_field ; i++ ) {
        // Standard scheme
        dimPrim[i] = params.n_space[i]+1;
        // + Ghost domain
        dimPrim[i] += 2*params.oversize[i];
    }
    
    
    A_  = new cField3D( dimPrim, "A" );
    A0_ = new cField3D( dimPrim, "Aold" );
    
    Phi_         = new Field3D( dimPrim, "Phi" );
    Phi_m        = new Field3D( dimPrim, "Phi_m" );
    
    GradPhix_    = new Field3D( dimPrim, "GradPhix" );
    GradPhix_m   = new Field3D( dimPrim, "GradPhix_m" );
    
    GradPhiy_    = new Field3D( dimPrim, "GradPhiy" );
    GradPhiy_m   = new Field3D( dimPrim, "GradPhiy_m" );
    
    GradPhiz_    = new Field3D( dimPrim, "GradPhiz" );
    GradPhiz_m   = new Field3D( dimPrim, "GradPhiz_m" );
    
}


LaserEnvelope3D::LaserEnvelope3D( LaserEnvelope *envelope, Patch *patch, ElectroMagn *EMfields, Params &params, unsigned int n_moved )
    : LaserEnvelope( envelope, patch, EMfields, params, n_moved )
{
    A_           = new cField3D( envelope->A_->dims_, "A" );
    A0_          = new cField3D( envelope->A0_->dims_, "Aold" );
    
    Phi_         = new Field3D( envelope->Phi_->dims_ );
    Phi_m        = new Field3D( envelope->Phi_m->dims_ );
    
    GradPhix_    = new Field3D( envelope->GradPhix_->dims_ );
    GradPhix_m   = new Field3D( envelope->GradPhix_m->dims_ );
    
    GradPhiy_    = new Field3D( envelope->GradPhiy_->dims_ );
    GradPhiy_m   = new Field3D( envelope->GradPhiy_m->dims_ );
    
    GradPhiz_    = new Field3D( envelope->GradPhiz_->dims_ );
    GradPhiz_m   = new Field3D( envelope->GradPhiz_m->dims_ );
    
}


void LaserEnvelope3D::initEnvelope( Patch *patch, ElectroMagn *EMfields )
{
    cField3D *A3D          = static_cast<cField3D *>( A_ );
    cField3D *A03D         = static_cast<cField3D *>( A0_ );
    Field3D *Env_Aabs3D    = static_cast<Field3D *>( EMfields->Env_A_abs_ );
    Field3D *Env_Eabs3D    = static_cast<Field3D *>( EMfields->Env_E_abs_ );
    Field3D *Env_Exabs3D   = static_cast<Field3D *>( EMfields->Env_Ex_abs_ );
    
    Field3D *Phi3D         = static_cast<Field3D *>( Phi_ );
    Field3D *Phi_m3D       = static_cast<Field3D *>( Phi_m );
    
    Field3D *GradPhix3D    = static_cast<Field3D *>( GradPhix_ );
    Field3D *GradPhix_m3D  = static_cast<Field3D *>( GradPhix_m );
    
    Field3D *GradPhiy3D    = static_cast<Field3D *>( GradPhiy_ );
    Field3D *GradPhiy_m3D  = static_cast<Field3D *>( GradPhiy_m );
    
    Field3D *GradPhiz3D    = static_cast<Field3D *>( GradPhiz_ );
    Field3D *GradPhiz_m3D  = static_cast<Field3D *>( GradPhiz_m );
    
    
    complex<double>     i1 = std::complex<double>( 0., 1 );
    
    // dx is the spatial step dx for 3D3V cartesian simulations
    double dx=cell_length[0];
    // dy is the spatial step dy for 3D3V cartesian simulations
    double dy=cell_length[1];
    // dz is the spatial step dz for 3D3V cartesian simulations
    double dz=cell_length[2];
    
    
    // position[0]: x coordinate
    // position[1]: y coordinate
    // position[2]: z coordinate
    // t: time coordinate --> x/c for the envelope initialization
    
    vector<double> position( 3 );
    position[0]      = dx*( ( double )( patch->getCellStartingGlobalIndex( 0 ) )+( A3D->isDual( 0 )?-0.5:0. ) );
    double pos1 = dy*( ( double )( patch->getCellStartingGlobalIndex( 1 ) )+( A3D->isDual( 1 )?-0.5:0. ) );
    double pos2 = dz*( ( double )( patch->getCellStartingGlobalIndex( 2 ) )+( A3D->isDual( 2 )?-0.5:0. ) );
    int N0 = ( int )A3D->dims()[0];
    int N1 = ( int )A3D->dims()[1];
    int N2 = ( int )A3D->dims()[2];
    
    // Create the x,y,z,t maps where profiles will be evaluated
    vector<Field *> xyz( 3 );
    Field *t;
    Field *t_previous_timestep;
    vector<unsigned int> n_space_to_create( 3 );
    n_space_to_create[0] = N0;
    n_space_to_create[1] = N1;
    n_space_to_create[2] = N2;
    
    for( unsigned int idim=0 ; idim<3 ; idim++ ) {
        xyz[idim] = new Field3D( n_space_to_create );
    }
    t                   = new Field3D( n_space_to_create );
    t_previous_timestep = new Field3D( n_space_to_create );
    
    for( int i=0 ; i<N0 ; i++ ) {
        position[1] = pos1;
        for( int j=0 ; j<N1 ; j++ ) {
            position[2] = pos2;
            for( int k=0 ; k<N2 ; k++ ) {
                //(*field3D)(i,j,k) += profile->valueAt(pos);
                for( unsigned int idim=0 ; idim<3 ; idim++ ) {
                    ( *xyz[idim] )( i, j, k ) = position[idim];
                }
                ( *t )( i, j, k )                   = position[0];
                ( *t_previous_timestep )( i, j, k ) = ( *t )( i, j, k )+timestep;
                position[2] += dz;
            }
            position[1] += dy;
        }
        position[0] += dx;
    }
    
    profile_->complexValuesAt( xyz, t, *A3D );
    profile_->complexValuesAt( xyz, t_previous_timestep, *A03D );
    
    for( unsigned int idim=0 ; idim<3 ; idim++ ) {
        delete xyz[idim];
    }
    
    delete t;
    delete t_previous_timestep;
    
    // UNSIGNED INT LEADS TO PB IN PERIODIC BCs
    for( unsigned int i=0 ; i<A_->dims_[0] ; i++ ) { // x loop
        //position[1] = pos1;
        for( unsigned int j=0 ; j<A_->dims_[1] ; j++ ) { // y loop
            //position[2] = pos2;
            for( unsigned int k=0 ; k<A_->dims_[2] ; k++ ) { // z loop
                // |A|
                ( *Env_Aabs3D )( i, j, k )= std::abs( ( *A3D )( i, j, k ) );
                // |E envelope| = |-(dA/dt-ik0cA)|
                ( *Env_Eabs3D )( i, j, k )= std::abs( ( ( *A3D )( i, j, k )-( *A03D )( i, j, k ) )/timestep - i1*( *A3D )( i, j, k ) );
                // compute ponderomotive potential at timestep n
                ( *Phi3D )( i, j, k )     = ellipticity_factor*std::abs( ( *A3D )( i, j, k ) ) * std::abs( ( *A3D )( i, j, k ) ) * 0.5;
                // compute ponderomotive potential at timestep n-1
                ( *Phi_m3D )( i, j, k )   = ellipticity_factor*std::abs( ( *A03D )( i, j, k ) ) * std::abs( ( *A03D )( i, j, k ) ) * 0.5;
                // interpolate in time
                ( *Phi_m3D )( i, j, k )   = 0.5*( ( *Phi_m3D )( i, j, k )+( *Phi3D )( i, j, k ) );
            }  // end z loop
        } // end y loop
    } // end x loop
    
    // Compute gradient of ponderomotive potential
    for( unsigned int i=1 ; i<A_->dims_[0]-1 ; i++ ) { // x loop
        for( unsigned int j=1 ; j<A_->dims_[1]-1 ; j++ ) { // y loop
            for( unsigned int k=1 ; k<A_->dims_[2]-1 ; k++ ) { // z loop
                // gradient in x direction
                ( *GradPhix3D )( i, j, k ) = ( ( *Phi3D )( i+1, j, k )-( *Phi3D )( i-1, j, k ) ) * one_ov_2dx;
                ( *GradPhix_m3D )( i, j, k ) = ( ( *Phi_m3D )( i+1, j, k )-( *Phi_m3D )( i-1, j, k ) ) * one_ov_2dx;
                // gradient in y direction
                ( *GradPhiy3D )( i, j, k ) = ( ( *Phi3D )( i, j+1, k )-( *Phi3D )( i, j-1, k ) ) * one_ov_2dy;
                ( *GradPhiy_m3D )( i, j, k ) = ( ( *Phi_m3D )( i, j+1, k )-( *Phi_m3D )( i, j-1, k ) ) * one_ov_2dy;
                // gradient in z direction
                ( *GradPhiz3D )( i, j, k ) = ( ( *Phi3D )( i, j, k+1 )-( *Phi3D )( i, j, k-1 ) ) * one_ov_2dz;
                ( *GradPhiz_m3D )( i, j, k ) = ( ( *Phi_m3D )( i, j, k+1 )-( *Phi_m3D )( i, j, k-1 ) ) * one_ov_2dz;
            }  // end z loop
        } // end y loop
    } // end x loop
    
}


LaserEnvelope3D::~LaserEnvelope3D()
{
}

void LaserEnvelope3D::updateEnvelope( ElectroMagn *EMfields )
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
    
    cField3D *A3D          = static_cast<cField3D *>( A_ );               // the envelope at timestep n
    cField3D *A03D         = static_cast<cField3D *>( A0_ );              // the envelope at timestep n-1
    Field3D *Env_Chi3D     = static_cast<Field3D *>( EMfields->Env_Chi_ ); // source term of envelope equation
    
    // temporary variable for updated envelope
    cField3D *A3Dnew;
    A3Dnew  = new cField3D( A_->dims_ );
    
    //// explicit solver
    for( unsigned int i=1 ; i <A_->dims_[0]-1; i++ ) { // x loop
        for( unsigned int j=1 ; j < A_->dims_[1]-1 ; j++ ) { // y loop
            for( unsigned int k=1 ; k < A_->dims_[2]-1; k++ ) { // z loop
                ( *A3Dnew )( i, j, k ) -= ( *Env_Chi3D )( i, j, k )*( *A3D )( i, j, k ); // subtract here source term Chi*A from plasma
                // A3Dnew = laplacian - source term
                ( *A3Dnew )( i, j, k ) += ( ( *A3D )( i-1, j, k )-2.*( *A3D )( i, j, k )+( *A3D )( i+1, j, k ) )*one_ov_dx_sq; // x part
                ( *A3Dnew )( i, j, k ) += ( ( *A3D )( i, j-1, k )-2.*( *A3D )( i, j, k )+( *A3D )( i, j+1, k ) )*one_ov_dy_sq; // y part
                ( *A3Dnew )( i, j, k ) += ( ( *A3D )( i, j, k-1 )-2.*( *A3D )( i, j, k )+( *A3D )( i, j, k+1 ) )*one_ov_dz_sq; // z part
                // A3Dnew = A3Dnew+2ik0*dA/dx
                ( *A3Dnew )( i, j, k ) += i1_2k0_over_2dx*( ( *A3D )( i+1, j, k )-( *A3D )( i-1, j, k ) );
                // A3Dnew = A3Dnew*dt^2
                ( *A3Dnew )( i, j, k )  = ( *A3Dnew )( i, j, k )*dt_sq;
                // A3Dnew = A3Dnew + 2/c^2 A3D - (1+ik0cdt)A03D/c^2
                ( *A3Dnew )( i, j, k ) += 2.*( *A3D )( i, j, k )-one_plus_ik0dt*( *A03D )( i, j, k );
                // A3Dnew = A3Dnew * (1+ik0dct)/(1+k0^2c^2dt^2)
                ( *A3Dnew )( i, j, k )  = ( *A3Dnew )( i, j, k )*one_plus_ik0dt_ov_one_plus_k0sq_dtsq;
            } // end z loop
        } // end y loop
    } // end x loop
    
    for( unsigned int i=1 ; i <A_->dims_[0]-1; i++ ) { // x loop
        for( unsigned int j=1 ; j < A_->dims_[1]-1 ; j++ ) { // y loop
            for( unsigned int k=1 ; k < A_->dims_[2]-1; k++ ) { // z loop
                // final back-substitution
                ( *A03D )( i, j, k )       = ( *A3D )( i, j, k );
                ( *A3D )( i, j, k )        = ( *A3Dnew )( i, j, k );
            } // end z loop
        } // end y loop
    } // end x loop
    
    delete A3Dnew;
} // end LaserEnvelope3D::updateEnvelope

void LaserEnvelope3D::updateEnvelopeReducedDispersion( ElectroMagn *EMfields )
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
    
    // An optimized form for the derivatives along x has been proposed in D. Terzani, P. Londrillo, JCP 2019
    // to reduce the numerical dispersion for the envelope solver.
    // The derivatives along x of the reduced dispersion scheme are defined as follows:
    // delta= [1-(dt/dx)^2]/3,
    // (dA/dx)_opt = (1+delta)*(dA/dx) - delta*(A_{i+2,j,k}-A_{i-2,j,k})/4/dx
    // (d^2A/dx^2)_opt = (1+delta)*(d^2A/dx^2) - delta*(A_{i+2,j,k}-2*A_{i,j,k}+A_{i-2,j,k})/(4dx^2)
    
    cField3D *A3D          = static_cast<cField3D *>( A_ );               // the envelope at timestep n
    cField3D *A03D         = static_cast<cField3D *>( A0_ );              // the envelope at timestep n-1
    Field3D *Env_Chi3D     = static_cast<Field3D *>( EMfields->Env_Chi_ ); // source term of envelope equation
  
    
    // temporary variable for updated envelope
    cField3D *A3Dnew;
    A3Dnew  = new cField3D( A_->dims_ );
    
    //// explicit solver
    for( unsigned int i=2 ; i <A_->dims_[0]-2; i++ ) { // x loop
        for( unsigned int j=1 ; j < A_->dims_[1]-1 ; j++ ) { // y loop
            for( unsigned int k=1 ; k < A_->dims_[2]-1; k++ ) { // z loop
                ( *A3Dnew )( i, j, k ) -= ( *Env_Chi3D )( i, j, k )*( *A3D )( i, j, k ); // subtract here source term Chi*A from plasma
                // A3Dnew = laplacian - source term
                ( *A3Dnew )( i, j, k ) += (1.+delta)*    ( ( *A3D )( i-1, j, k )-2.*( *A3D )( i, j, k )+( *A3D )( i+1, j, k ) )*one_ov_dx_sq; // x part with optimized derivative
                ( *A3Dnew )( i, j, k ) -= delta*         ( ( *A3D )( i-2, j, k )-2.*( *A3D )( i, j, k )+( *A3D )( i+2, j, k ) )*one_ov_dx_sq*0.25;
                ( *A3Dnew )( i, j, k ) +=                ( ( *A3D )( i, j-1, k )-2.*( *A3D )( i, j, k )+( *A3D )( i, j+1, k ) )*one_ov_dy_sq; // y part
                ( *A3Dnew )( i, j, k ) +=                ( ( *A3D )( i, j, k-1 )-2.*( *A3D )( i, j, k )+( *A3D )( i, j, k+1 ) )*one_ov_dz_sq; // z part
                // A3Dnew = A3Dnew+2ik0*dA/dx, where dA/dx uses the optimized form
                ( *A3Dnew )( i, j, k ) += i1_2k0_over_2dx*(1.+delta)*( ( *A3D )( i+1, j, k )-( *A3D )( i-1, j, k ) );
                ( *A3Dnew )( i, j, k ) -= i1_2k0_over_2dx*delta*0.5 *( ( *A3D )( i+2, j, k )-( *A3D )( i-2, j, k ) );
                // A3Dnew = A3Dnew*dt^2
                ( *A3Dnew )( i, j, k )  = ( *A3Dnew )( i, j, k )*dt_sq;
                // A3Dnew = A3Dnew + 2/c^2 A3D - (1+ik0cdt)A03D/c^2
                ( *A3Dnew )( i, j, k ) += 2.*( *A3D )( i, j, k )-one_plus_ik0dt*( *A03D )( i, j, k );
                // A3Dnew = A3Dnew * (1+ik0dct)/(1+k0^2c^2dt^2)
                ( *A3Dnew )( i, j, k )  = ( *A3Dnew )( i, j, k )*one_plus_ik0dt_ov_one_plus_k0sq_dtsq;
            } // end z loop
        } // end y loop
    } // end x loop
    
    for( unsigned int i=2 ; i <A_->dims_[0]-2; i++ ) { // x loop
        for( unsigned int j=1 ; j < A_->dims_[1]-1 ; j++ ) { // y loop
            for( unsigned int k=1 ; k < A_->dims_[2]-1; k++ ) { // z loop
                // final back-substitution
                ( *A03D )( i, j, k )       = ( *A3D )( i, j, k );
                ( *A3D )( i, j, k )        = ( *A3Dnew )( i, j, k );
            } // end z loop
        } // end y loop
    } // end x loop
    
    delete A3Dnew;
} // end LaserEnvelope3D::updateEnvelopeReducedDispersion


void LaserEnvelope3D::computePhiEnvAEnvE( ElectroMagn *EMfields )
{

    // computes Phi=|A|^2/2 (the ponderomotive potential), new values immediately after the envelope update
    
    cField3D *A3D          = static_cast<cField3D *>( A_ );                   // the envelope at timestep n

    cField3D *A03D         = static_cast<cField3D *>( A0_ );                  // the envelope at timestep n-1
    
    Field3D *Phi3D         = static_cast<Field3D *>( Phi_ );                  //Phi=|A|^2/2 is the ponderomotive potential

    Field3D *Env_Aabs3D    = static_cast<Field3D *>( EMfields->Env_A_abs_ );  // field for diagnostic and ionization
     
    Field3D *Env_Eabs3D    = static_cast<Field3D *>( EMfields->Env_E_abs_ );  // field for diagnostic and ionization
 
    Field3D *Env_Exabs3D   = static_cast<Field3D *>( EMfields->Env_Ex_abs_ ); // field for diagnostic and ionization
    
    
    // Compute ponderomotive potential Phi=|A|^2/2, at timesteps n+1, including ghost cells
    for( unsigned int i=1 ; i <A_->dims_[0]-1; i++ ) { // x loop
        for( unsigned int j=1 ; j < A_->dims_[1]-1; j++ ) { // y loop
            for( unsigned int k=1 ; k < A_->dims_[2]-1; k++ ) { // z loop
                ( *Phi3D )( i, j, k )      = ellipticity_factor*std::abs( ( *A3D )( i, j, k ) ) * std::abs( ( *A3D )( i, j, k ) ) * 0.5;
                ( *Env_Aabs3D )( i, j, k ) = std::abs( ( *A3D )( i, j, k ) );
                // |E envelope| = |-(dA/dt-ik0cA)|, forward finite differences for the time derivative
                ( *Env_Eabs3D )( i, j, k ) = std::abs( ( ( *A3D )( i, j, k )-( *A03D )( i, j, k ) )/timestep - i1*( *A3D )( i, j, k ) );
                // |Ex envelope| = |-(dA/dy|, central finite difference for the space derivative
                ( *Env_Exabs3D )( i, j,k ) = std::abs( ( ( *A3D )( i, j+1,k)-( *A3D  )( i, j-1,k) )*one_ov_2dy );
            } // end z loop
        } // end y loop
    } // end x loop
    
} // end LaserEnvelope3D::computePhiEnvAEnvE


void LaserEnvelope3D::computeGradientPhi( ElectroMagn *EMfields )
{

    // computes gradient of Phi=|A|^2/2 (the ponderomotive potential), new values immediately after the envelope update
    Field3D *GradPhix3D    = static_cast<Field3D *>( GradPhix_ );
    Field3D *GradPhiy3D    = static_cast<Field3D *>( GradPhiy_ );
    Field3D *GradPhiz3D    = static_cast<Field3D *>( GradPhiz_ );
    Field3D *Phi3D         = static_cast<Field3D *>( Phi_ );      //Phi=|A|^2/2 is the ponderomotive potential
    
    
    // Compute gradients of Phi, at timesteps n
    for( unsigned int i=1 ; i <A_->dims_[0]-1; i++ ) { // x loop
        for( unsigned int j=1 ; j < A_->dims_[1]-1 ; j++ ) { // y loop
            for( unsigned int k=1 ; k < A_->dims_[2]-1; k++ ) { // z loop
            
                // gradient in x direction
                ( *GradPhix3D )( i, j, k ) = ( ( *Phi3D )( i+1, j, k )-( *Phi3D )( i-1, j, k ) ) * one_ov_2dx;
                // gradient in y direction
                ( *GradPhiy3D )( i, j, k ) = ( ( *Phi3D )( i, j+1, k )-( *Phi3D )( i, j-1, k ) ) * one_ov_2dy;
                // gradient in z direction
                ( *GradPhiz3D )( i, j, k ) = ( ( *Phi3D )( i, j, k+1 )-( *Phi3D )( i, j, k-1 ) ) * one_ov_2dz;
            } // end z loop
        } // end y loop
    } // end x loop
    
} // end LaserEnvelope3D::computeGradientPhi


void LaserEnvelope3D::savePhiAndGradPhi()
{
    // Static cast of the fields
    Field3D *Phi3D         = static_cast<Field3D *>( Phi_ );
    Field3D *Phi_m3D       = static_cast<Field3D *>( Phi_m );
    
    Field3D *GradPhix3D    = static_cast<Field3D *>( GradPhix_ );
    Field3D *GradPhix_m3D = static_cast<Field3D *>( GradPhix_m );
    
    Field3D *GradPhiy3D    = static_cast<Field3D *>( GradPhiy_ );
    Field3D *GradPhiy_m3D  = static_cast<Field3D *>( GradPhiy_m );
    
    Field3D *GradPhiz3D    = static_cast<Field3D *>( GradPhiz_ );
    Field3D *GradPhiz_m3D  = static_cast<Field3D *>( GradPhiz_m );
    
    
    for( unsigned int i=0 ; i <A_->dims_[0]-1; i++ ) { // x loop
        for( unsigned int j=0 ; j < A_->dims_[1]-1 ; j++ ) { // y loop
            for( unsigned int k=0 ; k < A_->dims_[2]-1; k++ ) { // z loop
            
                // ponderomotive potential Phi=|A|^2/2
                ( *Phi_m3D )( i, j, k )       = ( *Phi3D )( i, j, k );
                
                // gradient of ponderomotive potential
                ( *GradPhix_m3D )( i, j, k )  = ( *GradPhix3D )( i, j, k );
                ( *GradPhiy_m3D )( i, j, k )  = ( *GradPhiy3D )( i, j, k );
                ( *GradPhiz_m3D )( i, j, k )  = ( *GradPhiz3D )( i, j, k );
                
            } // end z loop
        } // end y loop
    } // end x loop
    
    
}//END savePhiAndGradPhi


void LaserEnvelope3D::centerPhiAndGradPhi()
{
    // Static cast of the fields
    Field3D *Phi3D         = static_cast<Field3D *>( Phi_ );
    Field3D *Phi_m3D       = static_cast<Field3D *>( Phi_m );
    
    Field3D *GradPhix3D    = static_cast<Field3D *>( GradPhix_ );
    Field3D *GradPhix_m3D  = static_cast<Field3D *>( GradPhix_m );
    
    Field3D *GradPhiy3D    = static_cast<Field3D *>( GradPhiy_ );
    Field3D *GradPhiy_m3D  = static_cast<Field3D *>( GradPhiy_m );
    
    Field3D *GradPhiz3D    = static_cast<Field3D *>( GradPhiz_ );
    Field3D *GradPhiz_m3D  = static_cast<Field3D *>( GradPhiz_m );
    
    // Phi_m and GradPhi_m quantities now contain values at timestep n
    
    for( unsigned int i=0 ; i <A_->dims_[0]-1; i++ ) { // x loop
        for( unsigned int j=0 ; j < A_->dims_[1]-1 ; j++ ) { // y loop
            for( unsigned int k=0 ; k < A_->dims_[2]-1; k++ ) { // z loop
            
                // ponderomotive potential Phi=|A|^2/2
                ( *Phi_m3D )( i, j, k )       = 0.5*( ( *Phi_m3D )( i, j, k )+( *Phi3D )( i, j, k ) );
                
                // gradient of ponderomotive potential
                ( *GradPhix_m3D )( i, j, k )  = 0.5*( ( *GradPhix_m3D )( i, j, k )+( *GradPhix3D )( i, j, k ) );
                ( *GradPhiy_m3D )( i, j, k )  = 0.5*( ( *GradPhiy_m3D )( i, j, k )+( *GradPhiy3D )( i, j, k ) );
                ( *GradPhiz_m3D )( i, j, k )  = 0.5*( ( *GradPhiz_m3D )( i, j, k )+( *GradPhiz3D )( i, j, k ) );
                
            } // end z loop
        } // end y loop
    } // end x loop
    
    // Phi_m and GradPhi_m quantities now contain values interpolated at timestep n+1/2
    // these are used for the ponderomotive position advance
    
    
}//END centerPhiAndGradPhi


