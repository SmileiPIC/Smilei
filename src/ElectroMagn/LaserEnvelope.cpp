
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
cell_length    ( params.cell_length) ,timestep( params.timestep)
{
    
    PyObject * profile;
    if (!PyTools::extract_pyProfile("envelope_profile",profile,"LaserEnvelope"))
        MESSAGE("No envelope profile set !");
    profile_ = new Profile(profile, params.nDim_field, "envelope");

    int ienvlaser = 0;
    ostringstream name("");
    name << "Laser Envelope " << endl;
    ostringstream info("");

    bool  envelope_solver_read;
    // Read laser envelope parameters
    envelope_solver_read          = false; // default value
    std:: string envelope_solver  = "explicit"; // default value
    envelope_solver_read          = PyTools::extract("envelope_solver",envelope_solver,"LaserEnvelope");

    bool omega;
    double omega_value(0);

    omega      = PyTools::extract("omega",omega_value,"LaserEnvelope");

    info << "\t Laser Envelope parameters: "<< endl;
    // omega
    info << "\t\t\tomega              : " << omega_value << endl;
    // envelope solver
    info << "\t\t\tenvelope solver    : " << envelope_solver << endl;
  
    // Display info
    if( patch->isMaster() ) {
        MESSAGE( info.str() );
      }
      
    EnvBoundCond = EnvelopeBC_Factory::create(params, patch);
}


LaserEnvelope::LaserEnvelope( LaserEnvelope *envelope, Patch* patch, ElectroMagn* EMfields ) :
cell_length    ( envelope->cell_length ), timestep( envelope->timestep)
{
    profile_ = envelope->profile_;
}


LaserEnvelope::~LaserEnvelope()
{
    delete profile_;

    delete A0_;
    delete A_;
  
    int nBC = EnvBoundCond.size();
    for ( int i=0 ; i<nBC ;i++ )
        if (EnvBoundCond[i]!=NULL) delete EnvBoundCond[i];
}


LaserEnvelope3D::LaserEnvelope3D( Params& params, Patch* patch, ElectroMagn* EMfields )
    : LaserEnvelope(params, patch, EMfields)
{
    std::vector<unsigned int>  dimPrim( params.nDim_field );
    // Dimension of the primal and dual grids
    for (size_t i=0 ; i<params.nDim_field ; i++) {
        // Standard scheme
        dimPrim[i] = params.n_space[i]+1;
        // + Ghost domain
        dimPrim[i] += 2*params.oversize[i];
    }


    A_  = new cField3D( dimPrim );
    A0_ = new cField3D( dimPrim );

    initEnvelope( patch,EMfields );
}


LaserEnvelope3D::LaserEnvelope3D( LaserEnvelope *envelope, Patch* patch,ElectroMagn* EMfields )
    : LaserEnvelope(envelope,patch,EMfields)
{
    A_  = new cField3D( envelope->A_->dims_  );
    A0_ = new cField3D( envelope->A0_->dims_ );

    initEnvelope( patch,EMfields );
}


void LaserEnvelope3D::initEnvelope( Patch* patch,ElectroMagn* EMfields )
{
    cField3D* A3D = static_cast<cField3D*>(A_);
    cField3D* A03D = static_cast<cField3D*>(A0_);
    Field3D* Env_Ar3D = static_cast<Field3D*>(EMfields->Env_Ar_);
    Field3D* Env_Ai3D = static_cast<Field3D*>(EMfields->Env_Ai_);
    vector<double> pos(3,0);
    vector<double> pos_previous_timestep(3,0);
    pos[0]      = cell_length[0]*((double)(patch->getCellStartingGlobalIndex(0))+(A3D->isDual(0)?-0.5:0.));
    pos_previous_timestep[0]      = pos[0]-timestep;
    double pos1 = cell_length[1]*((double)(patch->getCellStartingGlobalIndex(1))+(A3D->isDual(1)?-0.5:0.));
    double pos2 = cell_length[2]*((double)(patch->getCellStartingGlobalIndex(2))+(A3D->isDual(2)?-0.5:0.));
    // UNSIGNED INT LEADS TO PB IN PERIODIC BCs
    for (unsigned int i=0 ; i<A_->dims_[0] ; i++) {
        pos[1] = pos1;
        pos_previous_timestep[1] = pos1;
        for (unsigned int j=0 ; j<A_->dims_[1] ; j++) {
            pos[2] = pos2;
            pos_previous_timestep[2] = pos2;
            for (unsigned int k=0 ; k<A_->dims_[2] ; k++) {
                (*A3D)(i,j,k) += profile_->valueAt(pos);
                (*A03D)(i,j,k) += profile_->valueAt(pos_previous_timestep);
                (*Env_Ar3D)(i,j,k)=std::real((*A3D)(i,j,k));
                pos[2] += cell_length[2];
                pos_previous_timestep[2] += cell_length[2];
            }
            pos[1] += cell_length[1];
            pos_previous_timestep[1] += cell_length[1];
        }
        pos[0] += cell_length[0];
        pos_previous_timestep[0] += cell_length[0];
    }
}


LaserEnvelope3D::~LaserEnvelope3D()
{
    delete A_;
    delete A0_;
}

void LaserEnvelope3D::compute(ElectroMagn* EMfields)
{
    //// solves envelope equation in lab frame:
    //full_laplacian(a)+2ik0*(da/dz+(1/c)*da/dt)-d^2a/dt^2*(1/c^2)=kp^2 n/n0 a/gamma
    // where kp^2 n/n0 a/gamma is gathered in charge deposition
    
    //// auxiliary quantities
    //! laser wavenumber, i.e. omega0/c
    double k0=1.;
    //! laser wavenumber times the temporal step, i.e. omega0/c * dt
    double k0_dt=1.*timestep;
    //! 1/dt^2, where dt is the temporal step
    double dt_sq = timestep*timestep; 
    // imaginary unit
    complex<double> i1=std::complex<double>(0., 1);
  
    //! 1/dx^2, 1/dy^2, 1/dz^2, where dx,dy,dz are the spatial step dx for 3D3V cartesian simulations
    double one_ov_dx_sq=1./cell_length[0]/cell_length[0];
    double one_ov_dy_sq=1./cell_length[1]/cell_length[1];
    double one_ov_dz_sq=1./cell_length[2]/cell_length[2];
    
    //! 1/dx, where dx is the spatial step dx for 3D3V cartesian simulations
    double one_ov_2dx=1./2./cell_length[0];
  
    //->rho_e- ???;
    cField3D* A3D = static_cast<cField3D*>(A_);   // the envelope at timestep n
    cField3D* A03D = static_cast<cField3D*>(A0_); // the envelope at timestep n-1
    Field3D* Env_Ar3D = static_cast<Field3D*>(EMfields->Env_Ar_);
    Field3D* Env_Ai3D = static_cast<Field3D*>(EMfields->Env_Ai_);
    complex<double> A3Dnew;
    
    // find e_idx in all species
    int e_idx = 0;
    Field3D* rho_e = static_cast<Field3D*>(EMfields->rho_s[e_idx]);

    //// explicit solver 
    for (unsigned int i=1 ; i <A_->dims_[0]-1; i++){
        for (unsigned int j=1 ; j < A_->dims_[1]-1 ; j++){
            for (unsigned int k=1 ; k < A_->dims_[2]-1; k++){
                //(*A3D)(i,j,k) = (*A3D)(i,j,k);
                A3Dnew = 0.; // subtract here source term from plasma
                // A3Dnew = laplacian - source term
                A3Dnew += ((*A3D)(i-1,j,k)-2.*(*A3D)(i,j,k)+(*A3D)(i+1,j,k))*one_ov_dx_sq; // x part
                A3Dnew += ((*A3D)(i,j-1,k)-2.*(*A3D)(i,j,k)+(*A3D)(i,j+1,k))*one_ov_dy_sq; // y part
                A3Dnew += ((*A3D)(i,j,k-1)-2.*(*A3D)(i,j,k)+(*A3D)(i,j,k+1))*one_ov_dz_sq; // z part
                // A3Dnew = A3Dnew+2ik0*dA/dx
                A3Dnew += 2.*i1*k0*((*A3D)(i+1,j,k)-(*A3D)(i-1,j,k))*one_ov_2dx;
                // A3Dnew = A3Dnew*dt^2
                A3Dnew  = A3Dnew*dt_sq;
                // A3Dnew = A3Dnew + 2/c^2 A3D + (1+ik0cdt)A03D/c^2
                A3Dnew += 2.*(*A3D)(i,j,k)-(1.+k0_dt)*(*A03D)(i+1,j,k);
                // A3Dnew = A3Dnew * (1+ik0dct)/(1+k0^2c^2dt^2)
                A3Dnew  = A3Dnew* (1.+i1*k0_dt)/(1.+k0_dt*k0_dt);
                // final substitution
                (*A3D)(i,j,k)  = A3Dnew;
                (*A03D)(i,j,k) = (*A3D)(i,j,k);
                //(*A3D)(i,j,k) = (*A03D)(i,j,k);
                (*Env_Ar3D)(i,j,k) =std::real((*A3D)(i,j,k));
            }
        }
    }
}


void LaserEnvelope::boundaryConditions(int itime, double time_dual, Patch* patch, Params &params, SimWindow* simWindow)
{
    // Compute EM Bcs
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

}



