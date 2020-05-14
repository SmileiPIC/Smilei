#include "SpeciesVAdaptiveMixedSort.h"

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
#include "PartBoundCond.h"
#include "PartWall.h"
#include "BoundaryConditionType.h"

#include "ElectroMagn.h"
#include "Interpolator.h"
#include "InterpolatorFactory.h"
#include "Profile.h"

#include "Projector.h"
#include "ProjectorFactory.h"

#include "SimWindow.h"
#include "Patch.h"

// #include "Field.h"
#include "Field1D.h"
#include "Field2D.h"
#include "Field3D.h"
#include "Tools.h"

#include "DiagnosticTrack.h"

#include "SpeciesMetrics.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Species
// input: simulation parameters & Species index
// ---------------------------------------------------------------------------------------------------------------------
SpeciesVAdaptiveMixedSort::SpeciesVAdaptiveMixedSort( Params &params, Patch *patch ) :
    SpeciesV( params, patch )
{
    initCluster( params );
    npack_ = 0 ;
    packsize_ = 0;

}//END SpeciesVAdaptiveMixedSort creator

// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species
// ---------------------------------------------------------------------------------------------------------------------
SpeciesVAdaptiveMixedSort::~SpeciesVAdaptiveMixedSort()
{
}


void SpeciesVAdaptiveMixedSort::resizeCluster( Params &params )
{

    if( vectorized_operators ) {

        // We recompute the number of cells
        unsigned int ncells = ( params.n_space[0]+1 );
        for( unsigned int i=1; i < params.nDim_field; i++ ) {
            ncells *= ( params.n_space[i]+1 );
        }

        // We keep the current number of particles
        // int npart = particles->last_index[particles->last_index.size()-1];
        // int size = params.n_space[0]/clrw;

        particles->last_index.resize( ncells, 0 );
        particles->first_index.resize( ncells, 0 );
        //count.resize(ncells,0);

        particles->first_index[0] = 0;
        for( unsigned int ic=1; ic < ncells; ic++ ) {
            particles->first_index[ic] = particles->first_index[ic-1] + count[ic-1];
            particles->last_index[ic-1]= particles->first_index[ic];
        }
        //New total number of particles is stored as last element of particles->last_index
        particles->last_index[ncells-1] = particles->last_index[ncells-2] + count.back() ;

    } else {

        Species::resizeCluster( params );

    }
}// end resizeCluster



// -----------------------------------------------------------------------------
//! Compute part_cell_keys at patch creation.
//! This operation is normally done in the pusher to avoid additional particles pass.
// -----------------------------------------------------------------------------
void SpeciesVAdaptiveMixedSort::computeParticleCellKeys( Params &params )
{

    unsigned int ip, nparts;
    int IX;
    double X;

    //Number of particles before exchange
    nparts = particles->size();

    // Cell_keys is resized at the current number of particles
    particles->cell_keys.resize( nparts );

    // Reinitialize count to 0
    for( unsigned int ic=0; ic < count.size() ; ic++ ) {
        count[ic] = 0 ;
    }

    #pragma omp simd
    for( ip=0; ip < nparts ; ip++ ) {
        // Counts the # of particles in each cell (or sub_cell) and store it in sparticles->last_index.
        for( unsigned int ipos=0; ipos < nDim_particle ; ipos++ ) {
            X = particles->position( ipos, ip )-min_loc_vec[ipos];
            IX = round( X * dx_inv_[ipos] );
            particles->cell_keys[ip] = particles->cell_keys[ip] * this->length_[ipos] + IX;
        }
    }

    // Reduction of the number of particles per cell in count
    for( ip=0; ip < nparts ; ip++ ) {
        count[particles->cell_keys[ip]] ++ ;
    }

}

void SpeciesVAdaptiveMixedSort::importParticles( Params &params, Patch *patch, Particles &source_particles, vector<Diagnostic *> &localDiags )
{

    if( vectorized_operators ) {
        importParticles( params, patch, source_particles, localDiags );
    } else {
        Species::importParticles( params, patch, source_particles, localDiags );
    }
}

void SpeciesVAdaptiveMixedSort::sortParticles( Params &params, Patch * patch )
{
    if( vectorized_operators ) {
        SpeciesV::sortParticles( params , patch);
    } else {
        Species::sortParticles( params, patch );
    }
}

// -----------------------------------------------------------------------------
//! This function configures the type of species according
//! to the default mode regardless the particle distribution
//! params object containing global Parameters
//! patch object containing the current patch data and properties
// -----------------------------------------------------------------------------
void SpeciesVAdaptiveMixedSort::defaultConfigure( Params &params, Patch *patch )
{
    // Setup the vectorization state
    this->vectorized_operators = ( params.adaptive_default_mode == "on" );


    // Destroy and reconfigure operators
    this->reconfigure_operators( params, patch );

    // If we switch from non-vectorized to vectozied,
    // we have to reactivate the cell-sorting algorithm

    // We resize the bins
    resizeCluster( params );

    // We perform the sorting
    this->sortParticles( params , patch);

    // Reconfigure species to be imported
    this->reconfigure_particle_importation();

}

// -----------------------------------------------------------------------------
//! This function configures the type of species according
//! to the vectorization mode
//! params object containing global Parameters
//! patch object containing the current patch data and properties
// -----------------------------------------------------------------------------
void SpeciesVAdaptiveMixedSort::configuration( Params &params, Patch *patch )
{
    //float ratio_number_of_vecto_cells;
    float vecto_time = 0.;
    float scalar_time = 0.;

    // We first compute cell_keys: the number of particles per cell
    this->computeParticleCellKeys( params );

    // Species with particles
    if( particles->size() > 0 ) {

        // --------------------------------------------------------------------
        // Metrics 2 - based on the evaluation of the computational time
        SpeciesMetrics::get_computation_time( this->count,
                                              vecto_time,
                                              scalar_time );

        if( vecto_time <= scalar_time ) {
            this->vectorized_operators = true;
        } else if( vecto_time > scalar_time ) {
            this->vectorized_operators = false;
        }
    }
    // Default mode where there is no particles
    else {
        this->vectorized_operators = ( params.adaptive_default_mode == "on" );
    }

    // --------------------------------------------------------------------

#ifdef  __DEBUG
    std::cerr << "  > Species " << this->name_ << " configuration (" << this->vectorized_operators
              << ") default: " << params.adaptive_default_mode
              << " in patch (" << patch->Pcoordinates[0] << "," <<  patch->Pcoordinates[1] << "," <<  patch->Pcoordinates[2] << ")"
              << " of MPI process " << patch->MPI_me_
              << " (vecto time: " << vecto_time
              << ", scalar time: " << scalar_time
              << ", particle number: " << particles->size()
              << ")" << '\n';
#endif

    // Destroy and reconfigure operators
    this->reconfigure_operators( params, patch );

    // If we switch from non-vectorized to vectozied,
    // we have to reactivate the cell-sorting algorithm

    // We resize the bins
    resizeCluster( params );

    // We perform the sorting
    this->sortParticles( params , patch);

    // Reconfigure species to be imported
    this->reconfigure_particle_importation();

}

// -----------------------------------------------------------------------------
//! This function reconfigures the type of species according
//! to the vectorization mode
//! params object containing global Parameters
//! patch object containing the current patch data and properties
// -----------------------------------------------------------------------------
void SpeciesVAdaptiveMixedSort::reconfiguration( Params &params, Patch *patch )
{

    //unsigned int ncell;
    bool reasign_operators = false;
    //float ratio_number_of_vecto_cells;
    float vecto_time = 0;
    float scalar_time = 0;

    //split cell into smaller sub_cells for refined sorting
    //ncell = (params.n_space[0]+1);
    //for ( unsigned int i=1; i < params.nDim_field; i++) ncell *= (params.n_space[i]+1);

    // We first compute cell_keys: the number of particles per cell
    // if the current mode is without vectorization
    if( !this->vectorized_operators ) {
        this->computeParticleCellKeys( params );
    }

    // --------------------------------------------------------------------
    // Metrics 1 - based on the ratio of vectorized cells
    // Compute the number of cells that contain more than 8 particles
    //ratio_number_of_vecto_cells = SpeciesMetrics::get_ratio_number_of_vecto_cells(count,8);
    
    // Test metrics, if necessary we reasign operators
    //if ( (ratio_number_of_vecto_cells > 0.5 && this->vectorized_operators == false)
    //  || (ratio_number_of_vecto_cells < 0.5 && this->vectorized_operators == true))
    //{
    //    reasign_operators = true;
    //}
    // --------------------------------------------------------------------
    
    // --------------------------------------------------------------------
    // Metrics 2 - based on the evaluation of the computational time
    SpeciesMetrics::get_computation_time( count,
                                          vecto_time,
                                          scalar_time );

    if( ( vecto_time < scalar_time && this->vectorized_operators == false )
            || ( vecto_time > scalar_time && this->vectorized_operators == true ) ) {
        reasign_operators = true;
    }
    // --------------------------------------------------------------------

    // Operator reasignment if required by the metrics
    if( reasign_operators ) {
    
        // The type of operator is changed
        this->vectorized_operators = !this->vectorized_operators;

#ifdef  __DEBUG
        std::cerr << "  > Species " << this->name_ << " reconfiguration (" << this->vectorized_operators
                  << ") in patch (" << patch->Pcoordinates[0] << "," <<  patch->Pcoordinates[1] << "," <<  patch->Pcoordinates[2] << ")"
                  << " of MPI process " << patch->MPI_me_
                  << " (vecto time: " << vecto_time
                  << ", scalar time: " << scalar_time
                  << ", particle number: " << particles->size()
                  << ")" << '\n';
#endif

        // Destroy and reconfigure operators
        this->reconfigure_operators( params, patch );

        // If we switch from non-vectorized to vectozied,
        // we have to reactivate the cell-sorting algorithm

        // We resize the bins
        resizeCluster( params );

        // We perform the sorting
        this->sortParticles( params, patch );

        // Reconfigure species to be imported
        this->reconfigure_particle_importation();
    }
}

// -----------------------------------------------------------------------------
//! This function reconfigures the operators
// -----------------------------------------------------------------------------
void SpeciesVAdaptiveMixedSort::reconfigure_operators( Params &params, Patch *patch )
{
    // Destroy current operators
    delete Interp;
    delete Push;
    if( Push_ponderomotive_position ) {
        delete Push_ponderomotive_position;
    }
    delete Proj;
    
    // Reassign the correct Interpolator
    Interp = InterpolatorFactory::create( params, patch, this->vectorized_operators );
    // Reassign the correct Pusher to Push
    Push = PusherFactory::create( params, this );
    // Reassign the correct Ponderomotive Pusher if used
    if( Push_ponderomotive_position ) {
        Push_ponderomotive_position = PusherFactory::create_ponderomotive_position_updater( params, this );
    }
    // Reassign the correct Projector
    Proj = ProjectorFactory::create( params, patch, this->vectorized_operators );
}

// -----------------------------------------------------------------------------
//! This function reconfigures the operators
// -----------------------------------------------------------------------------
void SpeciesVAdaptiveMixedSort::reconfigure_particle_importation()
{
    // Local species for importation
    if( this->Ionize ) {
        this->electron_species->vectorized_operators = this->vectorized_operators;
    }
    if( this->Radiate ) {
        this->photon_species_->vectorized_operators = this->vectorized_operators;
    }
    if( this->Multiphoton_Breit_Wheeler_process ) {
        for( int k=0; k<2; k++ ) {
            this->mBW_pair_species[k]->vectorized_operators = this->vectorized_operators;
        }
    }
}
