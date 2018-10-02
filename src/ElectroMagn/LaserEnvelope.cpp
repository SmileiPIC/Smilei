
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

LaserEnvelope::LaserEnvelope( Params& params, Patch* patch, ElectroMagn* EMfields ) :
    profile_   ( NULL               ),
    cell_length( params.cell_length ),
    timestep   ( params.timestep    )
    
{
    
    PyObject * profile = NULL;
    if (!PyTools::extract_pyProfile("envelope_profile",profile,"LaserEnvelope"))
        MESSAGE("No envelope profile set !");
    profile_ = new Profile(profile, params.nDim_field+1, "envelope");
    
    // params.Laser_Envelope_model = true;
    
    ostringstream name("");
    name << "Laser Envelope " << endl;
    ostringstream info("");
    
    // Read laser envelope parameters
    std:: string envelope_solver  = "explicit"; // default value
    
    PyTools::extract("envelope_solver",envelope_solver,"LaserEnvelope");
    

    
    // double omega_value(0);
    // PyTools::extract("omega",omega_value,"LaserEnvelope");

    //double a0_laser;
    //PyTools::extract("a0",a0_laser,"LaserEnvelope");

    //double waist_laser;
    //PyTools::extract("waist",waist_laser,"LaserEnvelope");
  
    
    info << "\t Laser Envelope parameters: "<< endl;
    // envelope solver
    info << "\t\tEnvelope solver    : " << envelope_solver << endl;

    
    // a0
    //info << "\t\t\ta0                 : " << a0_laser << endl;
    // waist
    //info << "\t\t\twaist              : " << waist_laser << endl;
    
    // Display info
    if( patch->isMaster() ) {
        MESSAGE( info.str() );
    }
    
    EnvBoundCond = EnvelopeBC_Factory::create(params, patch);

    if( patch->isMaster() ) {
    for ( unsigned int i=0 ; i<params.grid_length.size() ; i++ ){
        MESSAGE(1,"\tdimension " << i);
        MESSAGE(1,"\t- Envelope boundary conditions: " << "(" << params.Env_BCs[i][0] << ", " << params.Env_BCs[i][1] << ")");
        }
    }
}

// Cloning constructor
LaserEnvelope::LaserEnvelope( LaserEnvelope *envelope, Patch* patch, ElectroMagn* EMfields, Params& params, unsigned int n_moved ) :
cell_length    ( envelope->cell_length ), timestep( envelope->timestep)
{
    if (n_moved ==0)
        profile_ = new Profile( envelope->profile_ );
    EnvBoundCond = EnvelopeBC_Factory::create(params, patch);  
}


LaserEnvelope::~LaserEnvelope()
{
    // Pb wih clones, know problem
    //if (profile_ != NULL) {
    //    delete profile_;
    //    profile_ = NULL;
    //}

    if( A0_ )         delete A0_;
    if( A_  )         delete A_;

    if( Phi_)         delete Phi_;
    if( Phiold_)      delete Phiold_;
  
    if( GradPhix_)    delete GradPhix_;
    if( GradPhixold_) delete GradPhixold_;

    if( GradPhiy_)    delete GradPhiy_;
    if( GradPhiyold_) delete GradPhiyold_;

    if( GradPhiz_)    delete GradPhiz_;
    if( GradPhizold_) delete GradPhizold_;
    
    int nBC = EnvBoundCond.size();
    for ( int i=0 ; i<nBC ;i++ )
        if (EnvBoundCond[i]!=NULL) delete EnvBoundCond[i];
}


LaserEnvelope3D::LaserEnvelope3D( Params& params, Patch* patch, ElectroMagn* EMfields )
    : LaserEnvelope(params, patch, EMfields )
{
    std::vector<unsigned int>  dimPrim( params.nDim_field );
    // Dimension of the primal and dual grids
    for (size_t i=0 ; i<params.nDim_field ; i++) {
        // Standard scheme
        dimPrim[i] = params.n_space[i]+1;
        // + Ghost domain
        dimPrim[i] += 2*params.oversize[i];
    }


    A_  = new cField3D( dimPrim, "A" );
    A0_ = new cField3D( dimPrim, "Aold" );

    Phi_         = new Field3D( dimPrim, "Phi" );
    Phiold_      = new Field3D( dimPrim, "Phiold" );

    GradPhix_    = new Field3D( dimPrim, "GradPhix" );
    GradPhixold_ = new Field3D( dimPrim, "GradPhixold" );

    GradPhiy_    = new Field3D( dimPrim, "GradPhiy" );
    GradPhiyold_ = new Field3D( dimPrim, "GradPhiyold" );

    GradPhiz_    = new Field3D( dimPrim, "GradPhiz" );
    GradPhizold_ = new Field3D( dimPrim, "GradPhizold" );

}


LaserEnvelope3D::LaserEnvelope3D( LaserEnvelope *envelope, Patch* patch,ElectroMagn* EMfields, Params& params, unsigned int n_moved )
    : LaserEnvelope(envelope,patch,EMfields,params, n_moved)
{
    A_           = new cField3D( envelope->A_->dims_ , "A"    );
    A0_          = new cField3D( envelope->A0_->dims_, "Aold" );

    Phi_         = new Field3D( envelope->Phi_->dims_ );
    Phiold_      = new Field3D( envelope->Phiold_->dims_ );

    GradPhix_    = new Field3D( envelope->GradPhix_->dims_ );
    GradPhixold_ = new Field3D( envelope->GradPhixold_->dims_ );

    GradPhiy_    = new Field3D( envelope->GradPhiy_->dims_ );
    GradPhiyold_ = new Field3D( envelope->GradPhiyold_->dims_ );

    GradPhiz_    = new Field3D( envelope->GradPhiz_->dims_ );
    GradPhizold_ = new Field3D( envelope->GradPhizold_->dims_ );

}


void LaserEnvelope3D::initEnvelope( Patch* patch,ElectroMagn* EMfields )
{
    cField3D* A3D          = static_cast<cField3D*>(A_);
    cField3D* A03D         = static_cast<cField3D*>(A0_);
    Field3D* Env_Aabs3D      = static_cast<Field3D*>(EMfields->Env_A_abs_);
    Field3D* Env_Eabs3D      = static_cast<Field3D*>(EMfields->Env_E_abs_);

    Field3D* Phi3D         = static_cast<Field3D*>(Phi_);
    Field3D* Phiold3D      = static_cast<Field3D*>(Phiold_); 

    Field3D* GradPhix3D    = static_cast<Field3D*>(GradPhix_);
    Field3D* GradPhixold3D = static_cast<Field3D*>(GradPhixold_); 

    Field3D* GradPhiy3D    = static_cast<Field3D*>(GradPhiy_);
    Field3D* GradPhiyold3D = static_cast<Field3D*>(GradPhiyold_); 
  
    Field3D* GradPhiz3D    = static_cast<Field3D*>(GradPhiz_);
    Field3D* GradPhizold3D = static_cast<Field3D*>(GradPhizold_); 

    vector<double> position(3,0);
    double t;
    double t_previous_timestep;

    complex<double>     i1 = std::complex<double>(0., 1);

    //! 1/(2dx), where dx is the spatial step dx for 3D3V cartesian simulations
    double one_ov_2dx=1./2./cell_length[0];
    //! 1/(2dy), where dy is the spatial step dy for 3D3V cartesian simulations
    double one_ov_2dy=1./2./cell_length[1];
    //! 1/(2dz), where dz is the spatial step dz for 3D3V cartesian simulations
    double one_ov_2dz=1./2./cell_length[2];

    // position_time[0]: x coordinate
    // position_time[1]: y coordinate
    // position_time[2]: z coordinate
    // t: time coordinate --> x/c for the envelope initialization
    
    position[0]           = cell_length[0]*((double)(patch->getCellStartingGlobalIndex(0))+(A3D->isDual(0)?-0.5:0.));
    t                     = position[0];          // x-ct     , t=0
    t_previous_timestep   = position[0]+timestep; // x-c(t-dt), t=0
    double pos1 = cell_length[1]*((double)(patch->getCellStartingGlobalIndex(1))+(A3D->isDual(1)?-0.5:0.));
    double pos2 = cell_length[2]*((double)(patch->getCellStartingGlobalIndex(2))+(A3D->isDual(2)?-0.5:0.));
    // UNSIGNED INT LEADS TO PB IN PERIODIC BCs
    for (unsigned int i=0 ; i<A_->dims_[0] ; i++) { // x loop
        position[1] = pos1;
        for (unsigned int j=0 ; j<A_->dims_[1] ; j++) { // y loop
            position[2] = pos2;
            for (unsigned int k=0 ; k<A_->dims_[2] ; k++) { // z loop
                (*A3D)(i,j,k)      += profile_->complexValueAt(position,t);
                (*A03D)(i,j,k)     += profile_->complexValueAt(position,t_previous_timestep);

                (*Env_Aabs3D)(i,j,k)= std::abs ((*A3D)(i,j,k));
                // |E envelope| = |-(dA/dt-ik0cA)|
                (*Env_Eabs3D)(i,j,k)= std::abs ( ((*A3D)(i,j,k)-(*A03D)(i,j,k))/timestep - i1*(*A3D)(i,j,k)    );

                (*Phi3D)(i,j,k)     = std::abs((*A3D) (i,j,k)) * std::abs((*A3D) (i,j,k)) * 0.5;
                (*Phiold3D)(i,j,k)  = std::abs((*A03D)(i,j,k)) * std::abs((*A03D)(i,j,k)) * 0.5;

                position[2] += cell_length[2];
            }  // end z loop
            position[1] += cell_length[1];
        } // end y loop
        position[0]          += cell_length[0];
        t                     = position[0];
        t_previous_timestep   = position[0]+timestep;
    } // end x loop

    // Compute gradients
    for (unsigned int i=1 ; i<A_->dims_[0]-1 ; i++) { // x loop
        for (unsigned int j=1 ; j<A_->dims_[1]-1 ; j++) { // y loop
            for (unsigned int k=1 ; k<A_->dims_[2]-1 ; k++) { // z loop
                // gradient in x direction
                (*GradPhix3D)   (i,j,k) = ( (*Phi3D)   (i+1,j  ,k  )-(*Phi3D)   (i-1,j  ,k  ) ) * one_ov_2dx;
                (*GradPhixold3D)(i,j,k) = ( (*Phiold3D)(i+1,j  ,k  )-(*Phiold3D)(i-1,j  ,k  ) ) * one_ov_2dx;
                // gradient in y direction
                (*GradPhiy3D)   (i,j,k) = ( (*Phi3D)   (i  ,j+1,k  )-(*Phi3D)   (i  ,j-1,k  ) ) * one_ov_2dy;
                (*GradPhiyold3D)(i,j,k) = ( (*Phiold3D)(i  ,j+1,k  )-(*Phiold3D)(i  ,j-1,k  ) ) * one_ov_2dy;
                // gradient in z direction
                (*GradPhiz3D)   (i,j,k) = ( (*Phi3D)   (i  ,j  ,k+1)-(*Phi3D)   (i  ,j  ,k-1) ) * one_ov_2dz;
                (*GradPhizold3D)(i,j,k) = ( (*Phiold3D)(i  ,j  ,k+1)-(*Phiold3D)(i  ,j  ,k-1) ) * one_ov_2dz; 
            }  // end z loop
        } // end y loop
    } // end x loop

}


LaserEnvelope3D::~LaserEnvelope3D()
{
}

void LaserEnvelope3D::compute(ElectroMagn* EMfields)
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
    //! laser wavenumber, i.e. omega0/c
    double              k0 = 1.;
    //! laser wavenumber times the temporal step, i.e. omega0/c * dt
    double           k0_dt = 1.*timestep;
    //! 1/dt^2, where dt is the temporal step
    double           dt_sq = timestep*timestep; 
    // imaginary unit
    complex<double>     i1 = std::complex<double>(0., 1);
  
    //! 1/dx^2, 1/dy^2, 1/dz^2, where dx,dy,dz are the spatial step dx for 3D3V cartesian simulations
    double one_ov_dx_sq    = 1./cell_length[0]/cell_length[0];
    double one_ov_dy_sq    = 1./cell_length[1]/cell_length[1];
    double one_ov_dz_sq    = 1./cell_length[2]/cell_length[2];
  
    cField3D* A3D          = static_cast<cField3D*>(A_);                  // the envelope at timestep n
    cField3D* A03D         = static_cast<cField3D*>(A0_);                 // the envelope at timestep n-1
    Field3D* Env_Chi3D     = static_cast<Field3D*>(EMfields->Env_Chi_);   // source term of envelope equation
    Field3D* Env_Aabs3D    = static_cast<Field3D*>(EMfields->Env_A_abs_); // field for diagnostic
    Field3D* Env_Eabs3D    = static_cast<Field3D*>(EMfields->Env_E_abs_); // field for diagnostic
    

    //! 1/(2dx), where dx is the spatial step dx for 3D3V cartesian simulations
    double one_ov_2dx      = 1./2./cell_length[0];

    //! 1/(2dx), where dx is the spatial step dx for 3D3V cartesian simulations
    double one_ov_2dt      = 1./2./timestep;

    // temporary variable for updated envelope
    cField3D* A3Dnew;
    A3Dnew  = new cField3D( A_->dims_  );

    //// explicit solver 
    for (unsigned int i=1 ; i <A_->dims_[0]-1; i++){ // x loop
        for (unsigned int j=1 ; j < A_->dims_[1]-1 ; j++){ // y loop
            for (unsigned int k=1 ; k < A_->dims_[2]-1; k++){ // z loop
                (*A3Dnew)(i,j,k) -= (*Env_Chi3D)(i,j,k)*(*A3D)(i,j,k); // subtract here source term Chi*A from plasma
                // A3Dnew = laplacian - source term
                (*A3Dnew)(i,j,k) += ((*A3D)(i-1,j  ,k  )-2.*(*A3D)(i,j,k)+(*A3D)(i+1,j  ,k  ))*one_ov_dx_sq; // x part
                (*A3Dnew)(i,j,k) += ((*A3D)(i  ,j-1,k  )-2.*(*A3D)(i,j,k)+(*A3D)(i  ,j+1,k  ))*one_ov_dy_sq; // y part
                (*A3Dnew)(i,j,k) += ((*A3D)(i  ,j  ,k-1)-2.*(*A3D)(i,j,k)+(*A3D)(i  ,j  ,k+1))*one_ov_dz_sq; // z part
                // A3Dnew = A3Dnew+2ik0*dA/dx
                (*A3Dnew)(i,j,k) += 2.*i1*k0*((*A3D)(i+1,j,k)-(*A3D)(i-1,j,k))*one_ov_2dx;
                // A3Dnew = A3Dnew*dt^2
                (*A3Dnew)(i,j,k)  = (*A3Dnew)(i,j,k)*dt_sq;
                // A3Dnew = A3Dnew + 2/c^2 A3D - (1+ik0cdt)A03D/c^2
                (*A3Dnew)(i,j,k) += 2.*(*A3D)(i,j,k)-(1.+i1*k0_dt)*(*A03D)(i,j,k);
                // A3Dnew = A3Dnew * (1+ik0dct)/(1+k0^2c^2dt^2)
                (*A3Dnew)(i,j,k)  = (*A3Dnew)(i,j,k)*(1.+i1*k0_dt)/(1.+k0_dt*k0_dt);
            } // end z loop
        } // end y loop
    } // end x loop

    for (unsigned int i=1 ; i <A_->dims_[0]-1; i++){ // x loop
        for (unsigned int j=1 ; j < A_->dims_[1]-1 ; j++){ // y loop
            for (unsigned int k=1 ; k < A_->dims_[2]-1; k++){ // z loop
             // final back-substitution
             // |E envelope| = |-(dA/dt-ik0cA)|
             (*Env_Eabs3D)(i,j,k) = std::abs ( ((*A3Dnew)(i,j,k)-(*A03D)(i,j,k))*one_ov_2dt - i1*(*A3D)(i,j,k) ); 
             (*A03D)(i,j,k)       = (*A3D)(i,j,k);
             (*A3D)(i,j,k)        = (*A3Dnew)(i,j,k); 
             (*Env_Aabs3D)(i,j,k) = std::abs((*A3D)(i,j,k));
            } // end z loop
        } // end y loop
    } // end x loop

    delete A3Dnew;
} // end LaserEnvelope3D::compute


void LaserEnvelope::boundaryConditions(int itime, double time_dual, Patch* patch, Params &params, SimWindow* simWindow)
{
    // Compute Envelope Bcs
    if ( ! (simWindow && simWindow->isMoving(time_dual)) ) {
        if (EnvBoundCond[0]!=NULL) { // <=> if !periodic
            EnvBoundCond[0]->apply(this, time_dual, patch);
            EnvBoundCond[1]->apply(this, time_dual, patch);
        }
    }
    if (EnvBoundCond.size()>2) {
        if (EnvBoundCond[2]!=NULL) {// <=> if !periodic
            EnvBoundCond[2]->apply(this, time_dual, patch);
            EnvBoundCond[3]->apply(this, time_dual, patch);
        }
    }
    if (EnvBoundCond.size()>4) {
        if (EnvBoundCond[4]!=NULL) {// <=> if !periodic
            EnvBoundCond[4]->apply(this, time_dual, patch);
            EnvBoundCond[5]->apply(this, time_dual, patch);
        }
    }

} // end LaserEnvelope::boundaryConditions


void LaserEnvelope3D::compute_Phi(ElectroMagn* EMfields){

    // computes Phi=|A|^2/2 (the ponderomotive potential) old and present values

    cField3D* A3D          = static_cast<cField3D*>(A_);          // the envelope at timestep n
    cField3D* A03D         = static_cast<cField3D*>(A0_);         // the envelope at timestep n-1

    Field3D* Phi3D         = static_cast<Field3D*>(Phi_);         //Phi=|A|^2/2 is the ponderomotive potential
    Field3D* Phiold3D      = static_cast<Field3D*>(Phiold_); 


    // Compute ponderomotive potential Phi=|A|^2/2, at timesteps n and n-1, including ghost cells and 
    // and old value of Phi
    for (unsigned int i=1 ; i <A_->dims_[0]-1; i++){ // x loop
        for (unsigned int j=1 ; j < A_->dims_[1]-1; j++){ // y loop
            for (unsigned int k=1 ; k < A_->dims_[2]-1; k++){ // z loop
                (*Phi3D)   (i,j,k)       = std::abs((*A3D) (i,j,k)) * std::abs((*A3D) (i,j,k)) * 0.5;
                (*Phiold3D)(i,j,k)       = std::abs((*A03D)(i,j,k)) * std::abs((*A03D)(i,j,k)) * 0.5;
            } // end z loop
        } // end y loop
    } // end x loop

} // end LaserEnvelope3D::compute_Phi


void LaserEnvelope3D::compute_gradient_Phi(ElectroMagn* EMfields){

    // computes gradient of Phi=|A|^2/2 (the ponderomotive potential), old and present values

    Field3D* GradPhix3D    = static_cast<Field3D*>(GradPhix_);
    Field3D* GradPhixold3D = static_cast<Field3D*>(GradPhixold_); 

    Field3D* GradPhiy3D    = static_cast<Field3D*>(GradPhiy_);
    Field3D* GradPhiyold3D = static_cast<Field3D*>(GradPhiyold_); 

    Field3D* GradPhiz3D    = static_cast<Field3D*>(GradPhiz_);
    Field3D* GradPhizold3D = static_cast<Field3D*>(GradPhizold_); 

    Field3D* Phi3D         = static_cast<Field3D*>(Phi_);         //Phi=|A|^2/2 is the ponderomotive potential


    //! 1/(2dx), where dx is the spatial step dx for 3D3V cartesian simulations
    double one_ov_2dx=1./2./cell_length[0];
    //! 1/(2dy), where dy is the spatial step dy for 3D3V cartesian simulations
    double one_ov_2dy=1./2./cell_length[1];
    //! 1/(2dz), where dz is the spatial step dz for 3D3V cartesian simulations
    double one_ov_2dz=1./2./cell_length[2];

    // Compute gradients of Phi, at timesteps n 
    for (unsigned int i=1 ; i <A_->dims_[0]-1; i++){ // x loop
        for (unsigned int j=1 ; j < A_->dims_[1]-1 ; j++){ // y loop
            for (unsigned int k=1 ; k < A_->dims_[2]-1; k++){ // z loop

             // old gradient values
             (*GradPhixold3D)(i,j,k)  = (*GradPhix3D)(i,j,k);
             (*GradPhiyold3D)(i,j,k)  = (*GradPhiy3D)(i,j,k);
             (*GradPhizold3D)(i,j,k)  = (*GradPhiy3D)(i,j,k);

             // gradient in x direction
             (*GradPhix3D)   (i,j,k) = ( (*Phi3D)   (i+1,j  ,k  )-(*Phi3D)   (i-1,j  ,k  ) ) * one_ov_2dx;
             // gradient in y direction
             (*GradPhiy3D)   (i,j,k) = ( (*Phi3D)   (i  ,j+1,k  )-(*Phi3D)   (i  ,j-1,k  ) ) * one_ov_2dy;
             // gradient in z direction
             (*GradPhiz3D)   (i,j,k) = ( (*Phi3D)   (i  ,j  ,k+1)-(*Phi3D)   (i  ,j  ,k-1) ) * one_ov_2dz;
            } // end z loop
        } // end y loop
    } // end x loop

} // end LaserEnvelope3D::compute_gradient_Phi
