
#include "LaserEnvelope.h"

#include "Params.h"
#include "Patch.h"
#include "cField3D.h"
#include "Field3D.h"
#include "ElectroMagn.h"
#include "Profile.h"

using namespace std;

LaserEnvelope::LaserEnvelope( Params& params, Patch* patch ) :
cell_length    ( params.cell_length)
{
    //PyObject * profile;
    //if (!PyTools::extract_pyProfile("envelope_profile",profile,"LaserEnvelope"))
    //    MESSAGE("No envelope profile set !");
    //profile_ = new Profile(profile, params.nDim_field, "envelope");

    int ienvlaser = 0;
    ostringstream name("");
    name << "Laser Envelope " << endl;
    ostringstream info("");

    bool  envelope_solver_read;
    // Read laser envelope parameters
    envelope_solver_read          = false; // default value
    std:: string envelope_solver  = "explicit"; // default value
    envelope_solver_read          = PyTools::extract("envelope_solver",envelope_solver,"LaserEnvelope");

    PyObject *time_profile=nullptr;
    vector<PyObject*>  space_profile;
    bool time, space, omega;
    double omega_value(0);

    omega      = PyTools::extract("omega",omega_value,"LaserEnvelope",ienvlaser);
    time       = PyTools::extract_pyProfile("time_envelope" , time_profile , "LaserEnvelope",ienvlaser);
    space      = PyTools::extract2EnvelopeProfiles ("space_envelope", ienvlaser, space_profile     );

    info << "\t Laser Envelope parameters: "<< endl;
    // omega
    info << "\t\t\tomega              : " << omega_value << endl;
    // envelope solver
    info << "\t\t\tenvelope solver    : " << envelope_solver << endl;

    // Display info
    if( patch->isMaster() ) {
        MESSAGE( info.str() );
      }

}


LaserEnvelope::LaserEnvelope( LaserEnvelope *envelope, Patch* patch ) :
cell_length    ( envelope->cell_length )
{
    //profile_ = envelope->profile_;
}


LaserEnvelope::~LaserEnvelope()
{
    //delete profile_;

    delete A0_;
    delete A_;
}


LaserEnvelope3D::LaserEnvelope3D( Params& params, Patch* patch )
    : LaserEnvelope(params, patch)
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

    initEnvelope( patch );
}


LaserEnvelope3D::LaserEnvelope3D( LaserEnvelope *envelope, Patch* patch )
    : LaserEnvelope(envelope, patch)
{
    A_  = new cField3D( envelope->A_->dims_  );
    A0_ = new cField3D( envelope->A0_->dims_ );

    initEnvelope( patch );
}


void LaserEnvelope3D::initEnvelope( Patch* patch )
{
    cField3D* A3D = static_cast<cField3D*>(A_);
    vector<double> pos(3,0);
    pos[0]      = cell_length[0]*((double)(patch->getCellStartingGlobalIndex(0))+(A3D->isDual(0)?-0.5:0.));
    double pos1 = cell_length[1]*((double)(patch->getCellStartingGlobalIndex(1))+(A3D->isDual(1)?-0.5:0.));
    double pos2 = cell_length[2]*((double)(patch->getCellStartingGlobalIndex(2))+(A3D->isDual(2)?-0.5:0.));
    // UNSIGNED INT LEADS TO PB IN PERIODIC BCs
    for (int i=0 ; i<A_->dims_[0] ; i++) {
        pos[1] = pos1;
        for (int j=0 ; j<A_->dims_[1] ; j++) {
            pos[2] = pos2;
            for (int k=0 ; k<A_->dims_[2] ; k++) {
                //(*A3D)(i,j,k) += profile_->valueAt(pos);
                pos[2] += cell_length[2];
            }
            pos[1] += cell_length[1];
        }
        pos[0] += cell_length[0];
    }
}


LaserEnvelope3D::~LaserEnvelope3D()
{
    delete A_;
    delete A0_;
}

void LaserEnvelope3D::compute(ElectroMagn* EMfields)
{
    //->rho_e- ???;
    cField3D* A3D = static_cast<cField3D*>(A_);
    cField3D* A03D = static_cast<cField3D*>(A0_);

    // find e_idx in all species
    int e_idx = 0;
    Field3D* rho_e = static_cast<Field3D*>(EMfields->rho_s[e_idx]);

    for (unsigned int i=0 ; i <A_->dims_[0]; i++)
        for (unsigned int j=0 ; j < A_->dims_[1] ; j++)
            for (unsigned int k=0 ; k < A_->dims_[2]; k++)
                (*A3D)(i,j,k) = (*A03D)(i,j,k);

}
