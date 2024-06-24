#include "PyTools.h"

#include <string>
#include <sstream>

#include "ParticleData.h"
#include "PeekAtSpecies.h"
#include "DiagnosticTrack.h"
#include "VectorPatch.h"
#include "Params.h"

using namespace std;

DiagnosticTrack::DiagnosticTrack( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, unsigned int iDiagTrackParticles, unsigned int idiag, OpenPMDparams &oPMD ) :
    DiagnosticParticleList( params, smpi, vecPatches, "DiagTrackParticles", "TrackParticlesDisordered_", iDiagTrackParticles, oPMD ),
    IDs_done( params.restart )
{
    write_id_ = true;
    
    // Inform each patch about this diag
    for( unsigned int ipatch=0; ipatch<vecPatches.size(); ipatch++ ) {
        vecPatches( ipatch )->vecSpecies[species_index_]->tracking_diagnostic = idiag;
    }
    
    // Obtain the approximate number of particles in the species
    if( params.print_expected_disk_usage ) {
        PeekAtSpecies peek( params, species_index_ );
        npart_total = peek.totalNumberofParticles();
    } else {
        npart_total = 0;
    }
}

DiagnosticTrack::~DiagnosticTrack()
{
    closeFile();
}

void DiagnosticTrack::openFile( Params &, SmileiMPI *smpi )
{
    // Create HDF5 file
    file_ = new H5Write( filename, &smpi->world() );
    file_->attr( "name", diag_name_ );
    
    // Attributes for openPMD
    openPMD_->writeRootAttributes( *file_, "no_meshes", "particles/" );
    
    data_group_ = new H5Write( file_, "data" );
    
    file_->flush();
}

void DiagnosticTrack::closeFile()
{
    if( file_ ) {
        delete data_group_;
        delete file_;
        file_ = NULL;
    }
}

void DiagnosticTrack::init( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches )
{
    // Set the IDs of the particles
    if( ! IDs_done ) {
        latest_Id = smpi->getRank() * 4294967296; // 2^32
        
        for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
            setIDs( vecPatches( ipatch ) );
        }
        
        IDs_done = true;
    }
    
    // create the file
    openFile( params, smpi );
}


H5Space * DiagnosticTrack::prepareH5( SimWindow *simWindow, SmileiMPI *smpi, int itime, uint32_t nParticles_local, uint64_t nParticles_global, uint64_t offset )
{
    // Make a new group for this iteration
    ostringstream t( "" );
    t << setfill( '0' ) << setw( 10 ) << itime;
    // Create groups for openPMD compatibility
    H5Write iteration_group = data_group_->group( t.str() );
    H5Write particles_group = iteration_group.group( "particles" );
    H5Write * species_group = new H5Write( &particles_group, species_name_ );
    loc_id_ = species_group;
    loc_charge_ = species_group;
    loc_weight_ = species_group;
    loc_chi_ = species_group;
    if( write_any_position_ ) {
        fill( loc_position_.begin(), loc_position_.end(), new H5Write( species_group, "position" ) );
        openPMD_->writeRecordAttributes( *loc_position_[0], SMILEI_UNIT_POSITION );
    }
    if( write_any_momentum_ ) {
        fill( loc_momentum_.begin(), loc_momentum_.end(), new H5Write( species_group, "momentum" ) );
        openPMD_->writeRecordAttributes( *loc_momentum_[0], SMILEI_UNIT_MOMENTUM );
    }
    if( write_any_E_ ) {
        fill( loc_E_.begin(), loc_E_.end(), new H5Write( species_group, "E" ) );
        openPMD_->writeRecordAttributes( *loc_E_[0], SMILEI_UNIT_EFIELD );
    }
    if( write_any_B_ ) {
        fill( loc_B_.begin(), loc_B_.end(), new H5Write( species_group, "B" ) );
        openPMD_->writeRecordAttributes( *loc_B_[0], SMILEI_UNIT_BFIELD );
    }
    if( write_any_W_ ) {
        fill( loc_W_.begin(), loc_W_.end(), new H5Write( species_group, "W" ) );
        openPMD_->writeRecordAttributes( *loc_W_[0], SMILEI_UNIT_ENERGY );
    }
    
    // PositionOffset (for OpenPMD)
    string xyz = "xyz";
    H5Write positionoffset_group = species_group->group( "positionOffset" );
    openPMD_->writeRecordAttributes( positionoffset_group, SMILEI_UNIT_POSITION );
    vector<uint64_t> np = {nParticles_global};
    for( unsigned int idim=0; idim<nDim_particle; idim++ ) {
        H5Write xyz_group = positionoffset_group.group( xyz.substr( idim, 1 ) );
        openPMD_->writeComponentAttributes( xyz_group, SMILEI_UNIT_POSITION );
        xyz_group.attr( "value", 0. );
        xyz_group.attr( "shape", np, H5T_NATIVE_UINT64 );
    }
    
    // Attributes for openPMD
    openPMD_->writeBasePathAttributes( iteration_group, itime );
    openPMD_->writeParticlesAttributes( particles_group );
    openPMD_->writeSpeciesAttributes( *species_group );
    
    // Write x_moved
    iteration_group.attr( "x_moved", simWindow ? simWindow->getXmoved() : 0. );

    // Create the "latest_IDs" dataset
    // Create file space and select one element for each proc
    iteration_group.vect( "latest_IDs", latest_Id, smpi->getSize(), H5T_NATIVE_UINT64, smpi->getRank(), 1 );
    
    // Filespace and chunks
    hsize_t chunk = 0;
    if( nParticles_global>0 ) {
        // Set the chunk size
        unsigned int maximum_chunk_size = 100000000;
        unsigned int number_of_chunks = nParticles_global/maximum_chunk_size;
        if( nParticles_global%maximum_chunk_size != 0 ) {
            number_of_chunks++;
        }
        if( number_of_chunks <= 1 ) {
            chunk = 0;
        } else {
            unsigned int chunk_size = nParticles_global/number_of_chunks;
            if( nParticles_global%number_of_chunks != 0 ) {
                chunk_size++;
            }
            chunk = chunk_size;
        }
    }
    return new H5Space( nParticles_global, offset, nParticles_local, chunk );
}

void DiagnosticTrack::deleteH5()
{
    delete loc_position_[0];
    delete loc_momentum_[0];
    delete loc_E_[0];
    delete loc_B_[0];
    delete loc_W_[0];
    delete loc_id_;
}

void DiagnosticTrack::modifyFiltered( VectorPatch &vecPatches, unsigned int ipatch )
{
    Particles *p = getParticles( vecPatches( ipatch ) );
    for( auto ipart: patch_selection[ipatch] ) {
        // If particle not tracked before ( the 7 first bytes (ID<2^56) == 0 ), then set its ID
        if(( p->id( ipart ) & 72057594037927935) == 0 ) {
            p->id( ipart ) += ++latest_Id;
        }
    }
}

void DiagnosticTrack::setIDs( Patch *patch )
{
    // If filter, IDs are set on-the-fly
    if( has_filter ) {
        return;
    }
    unsigned int s = patch->vecSpecies[species_index_]->getNbrOfParticles();
    for( unsigned int iPart=0; iPart<s; iPart++ ) {
        patch->vecSpecies[species_index_]->particles->id( iPart ) = ++latest_Id;
    }
#if defined( SMILEI_ACCELERATOR_GPU )
    patch->vecSpecies[species_index_]->particles->initializeIDsOnDevice();
#endif
}


void DiagnosticTrack::setIDs( Particles &particles )
{
    // If filter, IDs are set on-the-fly
    if( has_filter ) {
        return;
    }
    unsigned int s = particles.numberOfParticles();
    #pragma omp critical
    {
        for( unsigned int iPart=0; iPart<s; iPart++ ) {
            particles.id( iPart ) = ++latest_Id;
        }
    }
}


void DiagnosticTrack::write_scalar_uint64( H5Write * location, string name, uint64_t &buffer, H5Space *file_space, H5Space *mem_space, unsigned int unit_type )
{
    H5Write a = location->array( name, buffer, H5T_NATIVE_UINT64, file_space, mem_space );
    openPMD_->writeRecordAttributes( a, unit_type );
    openPMD_->writeComponentAttributes( a, unit_type );
}
void DiagnosticTrack::write_scalar_short( H5Write * location, string name, short &buffer, H5Space *file_space, H5Space *mem_space, unsigned int unit_type )
{
    H5Write a = location->array( name, buffer, H5T_NATIVE_SHORT, file_space, mem_space );
    openPMD_->writeRecordAttributes( a, unit_type );
    openPMD_->writeComponentAttributes( a, unit_type );
}
void DiagnosticTrack::write_scalar_double( H5Write * location, string name, double &buffer, H5Space *file_space, H5Space *mem_space, unsigned int unit_type )
{
    H5Write a = location->array( name, buffer, H5T_NATIVE_DOUBLE, file_space, mem_space );
    openPMD_->writeRecordAttributes( a, unit_type );
    openPMD_->writeComponentAttributes( a, unit_type );
}

void DiagnosticTrack::write_component_uint64( H5Write * location, string name, uint64_t &buffer, H5Space *file_space, H5Space *mem_space, unsigned int unit_type )
{
    H5Write a = location->array( name, buffer, H5T_NATIVE_UINT64, file_space, mem_space );
    openPMD_->writeComponentAttributes( a, unit_type );
}
void DiagnosticTrack::write_component_short( H5Write * location, string name, short &buffer, H5Space *file_space, H5Space *mem_space, unsigned int unit_type )
{
    H5Write a = location->array( name, buffer, H5T_NATIVE_SHORT, file_space, mem_space );
    openPMD_->writeComponentAttributes( a, unit_type );
}
void DiagnosticTrack::write_component_double( H5Write * location, string name, double &buffer, H5Space *file_space, H5Space *mem_space, unsigned int unit_type )
{
    H5Write a = location->array( name, buffer, H5T_NATIVE_DOUBLE, file_space, mem_space );
    openPMD_->writeComponentAttributes( a, unit_type );
}



// SUPPOSED TO BE EXECUTED ONLY BY MASTER MPI
uint64_t DiagnosticTrack::getDiskFootPrint( int istart, int istop, Patch * )
{
    uint64_t footprint = 0;
    
    // Calculate the number of dumps between istart and istop
    uint64_t ndumps = timeSelection->howManyTimesBefore( istop ) - timeSelection->howManyTimesBefore( istart );
    
    // Calculate the number of written parameters
    int nparams = 6 + nDim_particle;
    
    // Add necessary global headers approximately
    footprint += 2500;
    
    // Add necessary timestep headers approximately
    footprint += ndumps * 11250;
    
    // Add size of each parameter
    footprint += ndumps * ( uint64_t )( nparams * npart_total * 8 );
    
    return footprint;
}

Particles * DiagnosticTrack::getParticles( Patch * patch )
{
    return patch->vecSpecies[species_index_]->particles;
}
