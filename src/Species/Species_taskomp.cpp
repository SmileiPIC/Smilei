
#include "Species_taskomp.h"

#include <cmath>
#include <ctime>
#include <cstdlib>

#include <iostream>

#include <omp.h>

// IDRIS
#include <cstring>
// IDRIS
#include "PusherFactory.h"
#include "IonizationFactory.h"
#include "RadiationFactory.h"
#include "MultiphotonBreitWheelerFactory.h"
#include "MergingFactory.h"
#include "PartBoundCond.h"
#include "PartWall.h"
#include "BoundaryConditionType.h"

#include "ElectroMagn.h"
#include "Interpolator.h"
#include "InterpolatorFactory.h"
#include "ProjectorFactory.h"
#include "Profile.h"
#include "ElectroMagnAM.h"
#include "Projector.h"
#include "ProjectorFactory.h"
#include "ParticleCreator.h"

#include "SimWindow.h"
#include "Patch.h"

// #include "Field.h"
#include "Field1D.h"
#include "Field2D.h"
#include "Field3D.h"
#include "Tools.h"

#include "DiagnosticTrack.h"

// necessary for the static_cast
#include "ElectroMagn2D.h"
#include "Projector2D2Order.h"
#include "ElectroMagn3D.h"
#include "Projector3D2Order.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Species_taskomp
// input: simulation parameters & Species index
// ---------------------------------------------------------------------------------------------------------------------
Species_taskomp::Species_taskomp( Params &params, Patch *patch )
    : Species( params, patch )
{

    // Init tags for the task dependencies of the particle operations
    bin_has_pushed            = new int[particles->first_index.size()];
    bin_has_done_particles_BC = new int[particles->first_index.size()];
    bin_has_projected         = new int[particles->first_index.size()];

    //! buffers for currents and charge
    b_Jx.resize(particles->first_index.size());
    b_Jy.resize(particles->first_index.size());
    b_Jz.resize(particles->first_index.size());
    b_rho.resize(particles->first_index.size());


    for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) {
        // allocate current-buffers, then put to zero their content
        b_Jx[ibin]  = new double[size_proj_buffer_Jx ];
        b_Jy[ibin]  = new double[size_proj_buffer_Jy ];
        b_Jz[ibin]  = new double[size_proj_buffer_Jz ];
        b_rho[ibin] = new double[size_proj_buffer_rho];
    }

}//END Species creator



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species
// ---------------------------------------------------------------------------------------------------------------------
Species_taskomp::~Species_taskomp()
{
    // delete Push;
    // delete Interp;
    // delete Proj;
    // 
    // 
    // if( partBoundCond ) {
    //     delete partBoundCond;
    // }
    // if( particles_per_cell_profile_ ) {
    //     delete particles_per_cell_profile_;
    // }
    // if( charge_profile_ ) {
    //     delete charge_profile_;
    // }
    // if( density_profile_ ) {
    //     delete density_profile_;
    // }
    // for( unsigned int i=0; i<velocity_profile_.size(); i++ ) {
    //     delete velocity_profile_[i];
    // }
    // for( unsigned int i=0; i<temperature_profile_.size(); i++ ) {
    //     delete temperature_profile_[i];
    // }
    // if( ionization_rate_!=Py_None ) {
    //     Py_DECREF( ionization_rate_ );
    // }

    if (bin_has_pushed != NULL){
        delete bin_has_pushed;
        delete bin_has_done_particles_BC;
        delete bin_has_projected;
    }

    for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) {
        // delete buffers
        delete[] b_Jx[ibin];
        delete[] b_Jy[ibin];
        delete[] b_Jz[ibin];
        delete[] b_rho[ibin];
    }

}

void Species_taskomp::dynamics( double time_dual, unsigned int ispec,
                        ElectroMagn *EMfields,
                        Params &params, bool diag_flag,
                        PartWalls *partWalls,
                        Patch *patch, SmileiMPI *smpi,
                        RadiationTables &RadiationTables,
                        MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables,
                        vector<Diagnostic *> &localDiags )
{
}
