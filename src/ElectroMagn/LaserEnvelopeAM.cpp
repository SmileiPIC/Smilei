
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

    one_ov_dl_sq    = 1./cell_length[0]/cell_length[0];
    one_ov_2dl      = 1./2./cell_length[0];
    one_ov_dr_sq    = 1./cell_length[1]/cell_length[1];
    one_ov_2dr      = 1./2./cell_length[1];
    dr              = cell_length[1];

    std::vector<unsigned int>  dimPrim( params.nDim_field );
    // Dimension of the primal and dual grids
    for( size_t i=0 ; i<params.nDim_field ; i++ ) {
        // Standard scheme
        dimPrim[i] = params.n_space[i]+1;
        // + Ghost domain
        dimPrim[i] += 2*params.oversize[i];
    }
    
    
    A_  = new cField2D( dimPrim, "A_mode_0" );
    A0_ = new cField2D( dimPrim, "Aold_mode_0" );
    
    Phi_         = new Field2D( dimPrim, "Phi_mode_0" );
    Phi_m        = new Field2D( dimPrim, "Phi_m_mode_0" );
    
    GradPhil_    = new Field2D( dimPrim, "GradPhil_mode_0" );
    GradPhil_m   = new Field2D( dimPrim, "GradPhil_m_mode_0" );
    
    GradPhir_    = new Field2D( dimPrim, "GradPhir_mode_0" );
    GradPhir_m   = new Field2D( dimPrim, "GradPhir_m_mode_0" );

    // in cylindrical symmetry the gradient along theta is zero
    
    
}


LaserEnvelopeAM::LaserEnvelopeAM( LaserEnvelope *envelope, Patch *patch, ElectroMagn *EMfields, Params &params, unsigned int n_moved )
    : LaserEnvelope( envelope, patch, EMfields, params, n_moved )
{
    A_           = new cField2D( envelope->A_->dims_, "A_mode_0" );
    A0_          = new cField2D( envelope->A0_->dims_, "Aold_mode_0" );
    
    Phi_         = new Field2D( envelope->Phi_->dims_ );
    Phi_m        = new Field2D( envelope->Phi_m->dims_ );
    
    GradPhil_    = new Field2D( envelope->GradPhil_->dims_ );
    GradPhil_m   = new Field2D( envelope->GradPhil_m->dims_ );
    
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
    Field2D *Env_Exabs2Dcyl   = static_cast<Field2D *>( EMfields->Env_Ex_abs_ );
    
    Field2D *Phi2Dcyl         = static_cast<Field2D *>( Phi_ );
    Field2D *Phi_m2Dcyl       = static_cast<Field2D *>( Phi_m );
    
    Field2D *GradPhil2Dcyl    = static_cast<Field2D *>( GradPhil_ );
    Field2D *GradPhil_m2Dcyl  = static_cast<Field2D *>( GradPhil_m );
    
    Field2D *GradPhir2Dcyl    = static_cast<Field2D *>( GradPhir_ );
    Field2D *GradPhir_m2Dcyl  = static_cast<Field2D *>( GradPhir_m );

    bool isYmin = ( static_cast<ElectroMagnAM *>( EMfields ) )->isYmin;
    int  j_glob = ( static_cast<ElectroMagnAM *>( EMfields ) )->j_glob_;
    
    
    vector<double> position( 2, 0 );
    double t;
    double t_previous_timestep;    
    
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
            ( *Env_Aabs2Dcyl )( i, j )= std::abs( ( *A2Dcyl )( i, j ) );
            // |E envelope| = |-(dA/dt-ik0cA)|
            ( *Env_Eabs2Dcyl )( i, j )= std::abs( ( ( *A2Dcyl )( i, j )-( *A02Dcyl )( i, j ) )/timestep - i1*( *A2Dcyl )( i, j ) );
            // compute ponderomotive potential at timestep n
            ( *Phi2Dcyl )( i, j )     = ellipticity_factor*std::abs( ( *A2Dcyl )( i, j ) ) * std::abs( ( *A2Dcyl )( i, j ) ) * 0.5;
            // compute ponderomotive potential at timestep n-1
            ( *Phi_m2Dcyl )( i, j )   = ellipticity_factor*std::abs( ( *A02Dcyl )( i, j ) ) * std::abs( ( *A02Dcyl )( i, j ) ) * 0.5;
            // interpolate in time
            ( *Phi_m2Dcyl )( i, j )   = 0.5*( ( *Phi_m2Dcyl )( i, j )+( *Phi2Dcyl )( i, j ) );
            
            position[1] += cell_length[1];
        } // end y loop
        position[0]          += cell_length[0];
        t                     = position[0];
        t_previous_timestep   = position[0]+timestep;
    } // end x loop
    
    // Compute gradient of ponderomotive potential and |Ex|
    for( unsigned int i=1 ; i<A_->dims_[0]-1 ; i++ ) { // x loop
        for( unsigned int j=std::max(isYmin*3,1) ; j<A_->dims_[1]-1 ; j++ ) { // r loop
            // gradient in x direction
            ( *GradPhil2Dcyl )( i, j ) = ( ( *Phi2Dcyl )( i+1, j )-( *Phi2Dcyl )( i-1, j ) ) * one_ov_2dl;
            ( *GradPhil_m2Dcyl )( i, j ) = ( ( *Phi_m2Dcyl )( i+1, j )-( *Phi_m2Dcyl )( i-1, j ) ) * one_ov_2dl;
            // gradient in thet direction
            ( *GradPhir2Dcyl )( i, j ) = ( ( *Phi2Dcyl )( i, j+1 )-( *Phi2Dcyl )( i, j-1 ) ) * one_ov_2dr;
            ( *GradPhir_m2Dcyl )( i, j ) = ( ( *Phi_m2Dcyl )( i, j+1 )-( *Phi_m2Dcyl )( i, j-1 ) ) * one_ov_2dr;

            // in cylindrical symmetry the gradient along theta is zero

            // |Ex envelope| = |-dA/dr|, central finite difference for the space derivative  
            ( *Env_Exabs2Dcyl )( i, j ) =  std::abs( (( *A2Dcyl )( i, j+1)-( *A2Dcyl )( i, j-1) )*one_ov_2dr );

        } // end r loop
    } // end l loop
    
    // 
    if (isYmin){ // axis BC
        for( unsigned int i=1 ; i <A_->dims_[0]-1; i++ ) { // l loop
            unsigned int j = 2;  // j_p=2 corresponds to r=0    
            ( *Phi2Dcyl )      ( i, j )   = ellipticity_factor*std::abs( ( *A2Dcyl )( i, j ) ) * std::abs( ( *A2Dcyl )( i, j ) ) * 0.5;  
            ( *Env_Aabs2Dcyl ) ( i, j )   = std::abs( ( *A2Dcyl )( i, j ) );
            ( *Env_Exabs2Dcyl )( i, j )   = 0.;
            // Axis BC on |A|
            ( *Env_Aabs2Dcyl )( i, j-1 )  = ( *Env_Aabs2Dcyl )( i, j+1 );
            ( *Env_Aabs2Dcyl )( i, j-2 )  = ( *Env_Aabs2Dcyl )( i, j+2 );
            // Axis BC on |E|
            ( *Env_Eabs2Dcyl )( i, j-1 )  = ( *Env_Eabs2Dcyl )( i, j+1 );
            ( *Env_Eabs2Dcyl )( i, j-2 )  = ( *Env_Eabs2Dcyl )( i, j+2 );
            // Axis BC on |Ex|
            ( *Env_Exabs2Dcyl )( i, j-1 ) = ( *Env_Eabs2Dcyl )( i, j+1 );
            ( *Env_Exabs2Dcyl )( i, j-2 ) = ( *Env_Eabs2Dcyl )( i, j+2 );
            
                  
        } // end l loop
    }
}


LaserEnvelopeAM::~LaserEnvelopeAM()
{
}

void LaserEnvelopeAM::updateEnvelope( ElectroMagn *EMfields )
{
    //// solves envelope equation in lab frame (see doc):
    // full_laplacian(A)+2ik0*(dA/dl+(1/c)*dA/dt)-d^2A/dt^2*(1/c^2)=Chi*A
    // where Chi is the plasma susceptibility [= sum(q^2*rho/mass/gamma_ponderomotive) for all species]
    // gamma_ponderomotive=sqrt(1+p^2+|A|^2/2) in normalized units
    
    // For an envelope moving from right to left, replace the imaginary unit i with its opposite (-i)
    // if using an envelope moving to the left, change the sign of the phase in the envelope initialization
    
    // the following explicit finite difference scheme is obtained through centered finite difference derivatives
    // e.g. (dA/dl) @ time n and indices ij = (A^n    _{i+1,j} - A^n    _{i-1,j}) /2/dl
    //      (dA/dt) @ time n and indices ij = (A^{n+1}_{i  ,j} - A^{n-1}_{i  ,j}) /2/dt
    // A0 is A^{n-1}
    //      (d^2A/dl^2) @ time n and indices ij = (A^{n}_{i+1,j}-2*A^{n}_{i,j}+A^{n}_{i-1,j})/dl^2
    
   
    cField2D *A2Dcyl       = static_cast<cField2D *>( A_ );               // the envelope at timestep n
    cField2D *A02Dcyl      = static_cast<cField2D *>( A0_ );              // the envelope at timestep n-1
    Field2D *Env_Chi2Dcyl  = static_cast<Field2D *>( EMfields->Env_Chi_ ); // source term of envelope equation
  
    int  j_glob = ( static_cast<ElectroMagnAM *>( EMfields ) )->j_glob_;
    bool isYmin = ( static_cast<ElectroMagnAM *>( EMfields ) )->isYmin;
  

    // temporary variable for updated envelope
    cField2D *A2Dcylnew;
    A2Dcylnew  = new cField2D( A_->dims_ );
 
    //// explicit solver
    for( unsigned int i=1 ; i <A_->dims_[0]-1; i++ ) { // l loop
        for( unsigned int j=std::max(3*isYmin,1) ; j < A_->dims_[1]-1 ; j++ ) { // r loop
            ( *A2Dcylnew )( i, j ) -= ( *Env_Chi2Dcyl )( i, j )*( *A2Dcyl )( i, j ); // subtract here source term Chi*A from plasma
            // A2Dcylnew = laplacian - source term
            ( *A2Dcylnew )( i, j ) += ( ( *A2Dcyl )( i-1, j )-2.*( *A2Dcyl )( i, j )+( *A2Dcyl )( i+1, j ) )*one_ov_dl_sq; // l part
            ( *A2Dcylnew )( i, j ) += ( ( *A2Dcyl )( i, j-1 )-2.*( *A2Dcyl )( i, j )+( *A2Dcyl )( i, j+1 ) )*one_ov_dr_sq; // r part
            ( *A2Dcylnew )( i, j ) += ( ( *A2Dcyl )( i, j+1 )-( *A2Dcyl )( i, j-1 ) ) * one_ov_2dr / ( ( double )( j_glob+j )*dr ); // r part         

            // A2Dcylnew = A2Dcylnew+2ik0*dA/dl
            ( *A2Dcylnew )( i, j ) += i1_2k0_over_2dl*( ( *A2Dcyl )( i+1, j )-( *A2Dcyl )( i-1, j ) );
            // A2Dcylnew = A2Dcylnew*dt^2
            ( *A2Dcylnew )( i, j )  = ( *A2Dcylnew )( i, j )*dt_sq;
            // A2Dcylnew = A2Dcylnew + 2/c^2 A2Dcyl - (1+ik0cdt)A02Dcyl/c^2
            ( *A2Dcylnew )( i, j ) += 2.*( *A2Dcyl )( i, j )-one_plus_ik0dt*( *A02Dcyl )( i, j );
            // A2Dcylnew = A2Dcylnew * (1+ik0dct)/(1+k0^2c^2dt^2)
            ( *A2Dcylnew )( i, j )  = ( *A2Dcylnew )( i, j )*one_plus_ik0dt_ov_one_plus_k0sq_dtsq;
        } // end r loop
    } // end l loop

    if (isYmin){ // axis BC
        for( unsigned int i=1 ; i <A_->dims_[0]-1; i++ ) { // l loop
           unsigned int j = 2; // j_p = 2 corresponds to r=0

           ( *A2Dcylnew )( i, j ) -= ( *Env_Chi2Dcyl )( i, j )*( *A2Dcyl )( i, j ); // subtract here source term Chi*A from plasma
           ( *A2Dcylnew )( i, j ) += ( ( *A2Dcyl )( i-1, j )-2.*( *A2Dcyl )( i, j )+( *A2Dcyl )( i+1, j ) )*one_ov_dl_sq; // l part
           ( *A2Dcylnew )( i, j ) += 4. * ( ( *A2Dcyl )( i, j+1 )-( *A2Dcyl )( i, j ) ) * one_ov_dr_sq; // r part

           // A2Dcylnew = A2Dcylnew+2ik0*dA/dl
           ( *A2Dcylnew )( i, j ) += i1_2k0_over_2dl*( ( *A2Dcyl )( i+1, j )-( *A2Dcyl )( i-1, j ) );
           // A2Dcylnew = A2Dcylnew*dt^2
           ( *A2Dcylnew )( i, j )  = ( *A2Dcylnew )( i, j )*dt_sq;
           // A2Dcylnew = A2Dcylnew + 2/c^2 A2Dcyl - (1+ik0cdt)A02Dcyl/c^2
           ( *A2Dcylnew )( i, j ) += 2.*( *A2Dcyl )( i, j )-one_plus_ik0dt*( *A02Dcyl )( i, j );
           // A2Dcylnew = A2Dcylnew * (1+ik0dct)/(1+k0^2c^2dt^2)
           ( *A2Dcylnew )( i, j )  = ( *A2Dcylnew )( i, j )*one_plus_ik0dt_ov_one_plus_k0sq_dtsq;
           
        } 
    }    

    for( unsigned int i=1 ; i <A_->dims_[0]-1; i++ ) { // x loop
        for( unsigned int j=isYmin*2 ; j < A_->dims_[1]-1 ; j++ ) { // r loop
            // final back-substitution
            ( *A02Dcyl )( i, j )       = ( *A2Dcyl )( i, j );
            ( *A2Dcyl )( i, j )        = ( *A2Dcylnew )( i, j );
        } // end r loop
    } // end l loop


    delete A2Dcylnew;
} // end LaserEnvelopeAM::updateEnvelope

void LaserEnvelopeAM::updateEnvelopeReducedDispersion( ElectroMagn *EMfields )
{
    //// solves envelope equation in lab frame (see doc):
    // full_laplacian(A)+2ik0*(dA/dz+(1/c)*dA/dt)-d^2A/dt^2*(1/c^2)=Chi*A
    // where Chi is the plasma susceptibility [= sum(q^2*rho/mass/gamma_ponderomotive) for all species]
    // gamma_ponderomotive=sqrt(1+p^2+|A|^2/2) in normalized units
    
    // For an envelope moving from right to left, replace the imaginary unit i with its opposite (-i)
    // if using an envelope moving to the left, change the sign of the phase in the envelope initialization
    
    // the following explicit finite difference scheme is obtained through centered finite difference derivatives
    // e.g. (dA/dl) @ time n and indices ij = (A^n    _{i+1,j} - A^n    _{i-1,j}) /2/dl
    //      (dA/dt) @ time n and indices ij = (A^{n+1}_{i  ,j} - A^{n-1}_{i  ,j}) /2/dt
    // A0 is A^{n-1}
    //      (d^2A/dl^2) @ time n and indices ij = (A^{n}_{i+1,j,k}-2*A^{n}_{i,j}+A^{n}_{i-1,j})/dl^2
    
    // An optimized form for the derivatives along x has been proposed in D. Terzani, P. Londrillo, JCP 2019
    // to reduce the numerical dispersion for the envelope solver.
    // The derivatives along l of the reduced dispersion scheme are defined as follows:
    // delta= [1-(dt/dl)^2]/3,
    // (dA/dl)_opt = (1+delta)*(dA/dl) - delta*(A_{i+2,j,k}-A_{i-2,j,k})/4/dl
    // (d^2A/dl^2)_opt = (1+delta)*(d^2A/dl^2) - delta*(A_{i+2,j,k}-2*A_{i,j,k}+A_{i-2,j,k})/(4dl^2)
  
   
    cField2D *A2Dcyl       = static_cast<cField2D *>( A_ );                  // the envelope at timestep n
    cField2D *A02Dcyl      = static_cast<cField2D *>( A0_ );                 // the envelope at timestep n-1
    Field2D *Env_Chi2Dcyl  = static_cast<Field2D *>( EMfields->Env_Chi_ );   // source term of envelope equation
  
    int  j_glob = ( static_cast<ElectroMagnAM *>( EMfields ) )->j_glob_;
    bool isYmin = ( static_cast<ElectroMagnAM *>( EMfields ) )->isYmin;


    // temporary variable for updated envelope
    cField2D *A2Dcylnew;
    A2Dcylnew  = new cField2D( A_->dims_ );
 
    //// explicit solver
    for( unsigned int i=2 ; i <A_->dims_[0]-2; i++ ) { // l loop
        for( unsigned int j=std::max(3*isYmin,1) ; j < A_->dims_[1]-1 ; j++ ) { // r loop
            ( *A2Dcylnew )( i, j ) -= ( *Env_Chi2Dcyl )( i, j )*( *A2Dcyl )( i, j ); // subtract here source term Chi*A from plasma
            // A2Dcylnew = laplacian - source term
            ( *A2Dcylnew )( i, j ) += (1.+delta)*( ( *A2Dcyl )( i-1, j )-2.*( *A2Dcyl )( i, j )+( *A2Dcyl )( i+1, j ) )*one_ov_dl_sq; // l part with optimized derivative
            ( *A2Dcylnew )( i, j ) -= delta*     ( ( *A2Dcyl )( i-2, j )-2.*( *A2Dcyl )( i, j )+( *A2Dcyl )( i+2, j ) )*one_ov_dl_sq*0.25;
            ( *A2Dcylnew )( i, j ) +=            ( ( *A2Dcyl )( i, j-1 )-2.*( *A2Dcyl )( i, j )+( *A2Dcyl )( i, j+1 ) )*one_ov_dr_sq; // r part
            ( *A2Dcylnew )( i, j ) +=            ( ( *A2Dcyl )( i, j+1 )-   ( *A2Dcyl )( i, j-1 ) ) * one_ov_2dr / ( ( double )( j_glob+j )*dr ); // r part         

            // A2Dcylnew = A2Dcylnew+2ik0*dA/dl, where dA/dl uses the optimized form
            ( *A2Dcylnew )( i, j ) += i1_2k0_over_2dl*(1.+delta)*( ( *A2Dcyl )( i+1, j )-( *A2Dcyl )( i-1, j ) );
            ( *A2Dcylnew )( i, j ) -= i1_2k0_over_2dl*delta*0.5 *( ( *A2Dcyl )( i+2, j )-( *A2Dcyl )( i-2, j ) );
            // A2Dcylnew = A2Dcylnew*dt^2
            ( *A2Dcylnew )( i, j )  = ( *A2Dcylnew )( i, j )*dt_sq;
            // A2Dcylnew = A2Dcylnew + 2/c^2 A2Dcyl - (1+ik0cdt)A02Dcyl/c^2
            ( *A2Dcylnew )( i, j ) += 2.*( *A2Dcyl )( i, j )-one_plus_ik0dt*( *A02Dcyl )( i, j );
            // A2Dcylnew = A2Dcylnew * (1+ik0dct)/(1+k0^2c^2dt^2)
            ( *A2Dcylnew )( i, j )  = ( *A2Dcylnew )( i, j )*one_plus_ik0dt_ov_one_plus_k0sq_dtsq;
        } // end r loop
    } // end l loop

    if (isYmin){ // axis BC
        for( unsigned int i=2 ; i <A_->dims_[0]-2; i++ ) { // l loop
           unsigned int j = 2; // j_p = 2 corresponds to r=0

           ( *A2Dcylnew )( i, j ) -= ( *Env_Chi2Dcyl )( i, j )*( *A2Dcyl )( i, j ); // subtract here source term Chi*A from plasma
           ( *A2Dcylnew )( i, j ) += (1.+delta)*( ( *A2Dcyl )( i-1, j )-2.*( *A2Dcyl )( i, j )+( *A2Dcyl )( i+1, j ) )*one_ov_dl_sq; // l part with optimized derivative
           ( *A2Dcylnew )( i, j ) -= delta*     ( ( *A2Dcyl )( i-2, j )-2.*( *A2Dcyl )( i, j )+( *A2Dcyl )( i+2, j ) )*one_ov_dl_sq*0.25;
           ( *A2Dcylnew )( i, j ) += 4. *       ( ( *A2Dcyl )( i, j+1 )   -( *A2Dcyl )( i, j ) ) * one_ov_dr_sq; // r part

           // A2Dcylnew = A2Dcylnew+2ik0*dA/dl, where dA/dl uses the optimized form
           ( *A2Dcylnew )( i, j ) += i1_2k0_over_2dl*(1.+delta)*( ( *A2Dcyl )( i+1, j )-( *A2Dcyl )( i-1, j ) );
           ( *A2Dcylnew )( i, j ) -= i1_2k0_over_2dl*delta*0.5 *( ( *A2Dcyl )( i+2, j )-( *A2Dcyl )( i-2, j ) );
           // A2Dcylnew = A2Dcylnew*dt^2
           ( *A2Dcylnew )( i, j )  = ( *A2Dcylnew )( i, j )*dt_sq;
           // A2Dcylnew = A2Dcylnew + 2/c^2 A2Dcyl - (1+ik0cdt)A02Dcyl/c^2
           ( *A2Dcylnew )( i, j ) += 2.*( *A2Dcyl )( i, j )-one_plus_ik0dt*( *A02Dcyl )( i, j );
           // A2Dcylnew = A2Dcylnew * (1+ik0dct)/(1+k0^2c^2dt^2)
           ( *A2Dcylnew )( i, j )  = ( *A2Dcylnew )( i, j )*one_plus_ik0dt_ov_one_plus_k0sq_dtsq;
           
        } 
    }    

    for( unsigned int i=2 ; i <A_->dims_[0]-2; i++ ) { // x loop
        for( unsigned int j=isYmin*2 ; j < A_->dims_[1]-1 ; j++ ) { // r loop
            // final back-substitution
            ( *A02Dcyl )( i, j )       = ( *A2Dcyl )( i, j );
            ( *A2Dcyl )( i, j )        = ( *A2Dcylnew )( i, j );            
        } // end r loop
    } // end l loop


    delete A2Dcylnew;
} // end LaserEnvelopeAM::updateEnvelopeReducedDispersion


void LaserEnvelopeAM::computePhiEnvAEnvE( ElectroMagn *EMfields )
{

    // computes Phi=|A|^2/2 (the ponderomotive potential), new values immediately after the envelope update
    
    cField2D *A2Dcyl          = static_cast<cField2D *>( A_ );               // the envelope at timestep n

    cField2D *A02Dcyl      = static_cast<cField2D *>( A0_ );                 // the envelope at timestep n-1
    
    Field2D *Phi2Dcyl         = static_cast<Field2D *>( Phi_ );              //Phi=|A|^2/2 is the ponderomotive potential
    
    Field2D *Env_Aabs2Dcyl = static_cast<Field2D *>( EMfields->Env_A_abs_ ); // field for diagnostic and ionization

    Field2D *Env_Eabs2Dcyl = static_cast<Field2D *>( EMfields->Env_E_abs_ ); // field for diagnostic and ionization

    Field2D *Env_Exabs2Dcyl = static_cast<Field2D *>( EMfields->Env_Ex_abs_ ); // field for diagnostic and ionization

    bool isYmin = ( static_cast<ElectroMagnAM *>( EMfields ) )->isYmin;

    int  j_glob = ( static_cast<ElectroMagnAM *>( EMfields ) )->j_glob_;
    
    
    // Compute ponderomotive potential Phi=|A|^2/2, at timesteps n+1, including ghost cells
    for( unsigned int i=0 ; i <A_->dims_[0]-1; i++ ) { // x loop
        for( unsigned int j=std::max(3*isYmin,1) ; j < A_->dims_[1]-1; j++ ) { // r loop
            ( *Phi2Dcyl )( i, j )       = ellipticity_factor*std::abs( ( *A2Dcyl )( i, j ) ) * std::abs( ( *A2Dcyl )( i, j ) ) * 0.5;
            ( *Env_Aabs2Dcyl )( i, j )  = std::abs( ( *A2Dcyl )( i, j ) );
            // |E envelope| = |-(dA/dt-ik0cA)|, forward finite difference for the time derivative
            ( *Env_Eabs2Dcyl )( i, j )  = std::abs( ( ( *A2Dcyl )( i, j )-( *A02Dcyl )( i, j ) )/timestep - i1*( *A2Dcyl )( i, j ) );
            // |Ex envelope| = |-dA/dr|, central finite difference for the space derivative
            ( *Env_Exabs2Dcyl )( i, j ) =  std::abs( (( *A2Dcyl )( i, j+1)-( *A2Dcyl )( i, j-1) )*one_ov_2dr );
        } // end r loop
    } // end l loop

    // 
    if (isYmin){ // axis BC
        for( unsigned int i=1 ; i <A_->dims_[0]-1; i++ ) { // l loop
            unsigned int j = 2;  // j_p=2 corresponds to r=0    
            ( *Phi2Dcyl )      ( i, j )   = std::abs( ( *A2Dcyl )( i, j ) ) * std::abs( ( *A2Dcyl )( i, j ) ) * 0.5;  
            ( *Env_Aabs2Dcyl ) ( i, j )   = std::abs( ( *A2Dcyl )( i, j ) );
            ( *Env_Eabs2Dcyl ) ( i, j )   = std::abs( ( ( *A2Dcyl )( i, j )-( *A02Dcyl )( i, j ) )/timestep - i1*( *A2Dcyl )( i, j ) );
            ( *Env_Exabs2Dcyl )( i, j )   = 0.;
            // Axis BC on |A|
            ( *Env_Aabs2Dcyl )( i, j-1 )  = ( *Env_Aabs2Dcyl )( i, j+1 );
            ( *Env_Aabs2Dcyl )( i, j-2 )  = ( *Env_Aabs2Dcyl )( i, j+2 );
            // Axis BC on |E|
            ( *Env_Eabs2Dcyl )( i, j-1 )  = ( *Env_Eabs2Dcyl )( i, j+1 );
            ( *Env_Eabs2Dcyl )( i, j-2 )  = ( *Env_Eabs2Dcyl )( i, j+2 );
            // Axis BC on |Ex|
            ( *Env_Exabs2Dcyl)( i, j-1 )  = ( *Env_Exabs2Dcyl)( i, j+1 );
            ( *Env_Exabs2Dcyl)( i, j-2 )  = ( *Env_Exabs2Dcyl)( i, j+2 );
            
                  
        } // end l loop
    }
    
} // end LaserEnvelopeAM::computePhiEnvAEnvE


void LaserEnvelopeAM::computeGradientPhi( ElectroMagn *EMfields )
{

    // computes gradient of Phi=|A|^2/2 (the ponderomotive potential), new values immediately after the envelope update
    Field2D *GradPhil2Dcyl    = static_cast<Field2D *>( GradPhil_ );
    Field2D *GradPhir2Dcyl    = static_cast<Field2D *>( GradPhir_ );
    Field2D *Phi2Dcyl         = static_cast<Field2D *>( Phi_ );      //Phi=|A|^2/2 is the ponderomotive potential
    bool isYmin = ( static_cast<ElectroMagnAM *>( EMfields ) )->isYmin;
    
    
    // Compute gradients of Phi, at timesteps n
    for( unsigned int i=1 ; i <A_->dims_[0]-1; i++ ) { // l loop
        for( unsigned int j=std::max(3*isYmin,1) ; j < A_->dims_[1]-1 ; j++ ) { // r loop
        
            // gradient in l direction
            ( *GradPhil2Dcyl )( i, j ) = ( ( *Phi2Dcyl )( i+1, j )-( *Phi2Dcyl )( i-1, j ) ) * one_ov_2dl;
            // gradient in r direction
            ( *GradPhir2Dcyl )( i, j ) = ( ( *Phi2Dcyl )( i, j+1 )-( *Phi2Dcyl )( i, j-1 ) ) * one_ov_2dr;
            
        } // end r loop
    } // end l loop

    
    if (isYmin){ // axis BC
        for( unsigned int i=1 ; i <A_->dims_[0]-1; i++ ) { // l loop
            unsigned int j = 2;  // j_p=2 corresponds to r=0      
      
            // gradient in x direction
            ( *GradPhil2Dcyl )( i, j ) = ( ( *Phi2Dcyl )( i+1, j )-( *Phi2Dcyl )( i-1, j ) ) * one_ov_2dl;
            // gradient in r direction, identically zero on r = 0
            ( *GradPhir2Dcyl )( i, j ) = 0. ; // ( ( *Phi2Dcyl )( i, j+1 )-( *Phi2Dcyl )( i, j-1 ) ) * one_ov_2dr;

            // Axis BC on gradient in x direction
            ( *GradPhil2Dcyl )( i, j-1 ) = ( *GradPhil2Dcyl )( i, j+1 );
            ( *GradPhil2Dcyl )( i, j-2 ) = ( *GradPhil2Dcyl )( i, j+2 );
            // Axis BC on gradient in r direction
            ( *GradPhir2Dcyl )( i, j-1 ) = ( *GradPhir2Dcyl )( i, j+1 );
            ( *GradPhir2Dcyl )( i, j-2 ) = ( *GradPhir2Dcyl )( i, j+2 );
                            
        } // end l loop
    }

} // end LaserEnvelopeAM::computeGradientPhi


void LaserEnvelopeAM::savePhiAndGradPhi()
{
    // Static cast of the fields
    Field2D *Phi2Dcyl         = static_cast<Field2D *>( Phi_ );
    Field2D *Phi_m2Dcyl       = static_cast<Field2D *>( Phi_m );
    
    Field2D *GradPhil2Dcyl    = static_cast<Field2D *>( GradPhil_ );
    Field2D *GradPhil_m2Dcyl  = static_cast<Field2D *>( GradPhil_m );
    
    Field2D *GradPhir2Dcyl    = static_cast<Field2D *>( GradPhir_ );
    Field2D *GradPhir_m2Dcyl  = static_cast<Field2D *>( GradPhir_m );
    
    for( unsigned int i=0 ; i <A_->dims_[0]-1; i++ ) { // l loop
        for( unsigned int j=0 ; j < A_->dims_[1]-1 ; j++ ) { // r loop
        
            // ponderomotive potential Phi=|A|^2/2
            ( *Phi_m2Dcyl )( i, j )       = ( *Phi2Dcyl )( i, j );
            
            // gradient of ponderomotive potential
            ( *GradPhil_m2Dcyl )( i, j )  = ( *GradPhil2Dcyl )( i, j );
            ( *GradPhir_m2Dcyl )( i, j )  = ( *GradPhir2Dcyl )( i, j );
            
        } // end r loop
    } // end l loop
    
    
}//END savePhiAndGradPhi


void LaserEnvelopeAM::centerPhiAndGradPhi()
{
    // Static cast of the fields
    Field2D *Phi2Dcyl         = static_cast<Field2D *>( Phi_ );
    Field2D *Phi_m2Dcyl       = static_cast<Field2D *>( Phi_m );
    
    Field2D *GradPhil2Dcyl    = static_cast<Field2D *>( GradPhil_ );
    Field2D *GradPhil_m2Dcyl  = static_cast<Field2D *>( GradPhil_m );
    
    Field2D *GradPhir2Dcyl    = static_cast<Field2D *>( GradPhir_ );
    Field2D *GradPhir_m2Dcyl  = static_cast<Field2D *>( GradPhir_m );
    
    // Phi_m and GradPhi_m quantities now contain values at timestep n
    
    for( unsigned int i=0 ; i <A_->dims_[0]-1; i++ ) { // l loop
        for( unsigned int j=0 ; j < A_->dims_[1]-1 ; j++ ) { // r loop
        
            // ponderomotive potential Phi=|A|^2/2
            ( *Phi_m2Dcyl )( i, j )       = 0.5*( ( *Phi_m2Dcyl )( i, j )+( *Phi2Dcyl )( i, j ) );
            
            // gradient of ponderomotive potential
            ( *GradPhil_m2Dcyl )( i, j )  = 0.5*( ( *GradPhil_m2Dcyl )( i, j )+( *GradPhil2Dcyl )( i, j ) );
            ( *GradPhir_m2Dcyl )( i, j )  = 0.5*( ( *GradPhir_m2Dcyl )( i, j )+( *GradPhir2Dcyl )( i, j ) );
            
        } // end r loop
    } // end l loop
    
    // Phi_m and GradPhi_m quantities now contain values interpolated at timestep n+1/2
    // these are used for the ponderomotive position advance
    
    
}//END centerPhiAndGradPhi


