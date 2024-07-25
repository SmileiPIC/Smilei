#include "PyTools.h"

#include <string>
#include <sstream>

#include "ParticleData.h"
#include "PeekAtSpecies.h"
#include "DiagnosticNewParticles.h"
#include "VectorPatch.h"
#include "Params.h"

using namespace std;

DiagnosticNewParticles::DiagnosticNewParticles( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, unsigned int iDiagnosticNewParticles, unsigned int, OpenPMDparams &oPMD ) :
    DiagnosticParticleList( params, smpi, vecPatches, "DiagNewParticles", "NewParticles_", iDiagnosticNewParticles, oPMD )
{
    write_id_ = vecPatches.species( 0, species_index_ )->particles->tracked;
    
    // Inform each patch about this diag
    for( unsigned int ipatch=0; ipatch<vecPatches.size(); ipatch++ ) {
        Species * s = vecPatches.species( ipatch, species_index_ );
        s->birth_records_ = new BirthRecords( *s->particles );
        // Find out the other species that may create this one by ionization
        for( auto s: vecPatches( ipatch )->vecSpecies ) {
            if( s->Ionize && s->electron_species_index == species_index_ ) {
                s->Ionize->save_ion_charge_ = true;
            }
        }
    }
}

DiagnosticNewParticles::~DiagnosticNewParticles()
{
    closeFile();
}

void DiagnosticNewParticles::openFile( Params &params, SmileiMPI *smpi )
{
    // Create HDF5 file
    file_ = new H5Write( filename, &smpi->world() );
    file_->attr( "name", diag_name_ );
    
    // Groups for openPMD
    H5Write data_group( file_, "data" );
    H5Write iteration_group( &data_group, "0" );
    H5Write particles_group( &iteration_group, "particles" );
    H5Write species_group( &particles_group, species_name_ );
    
    // Attributes for openPMD
    openPMD_->writeRootAttributes( *file_, "no_meshes", "particles/" );
    openPMD_->writeBasePathAttributes( iteration_group, 0 );
    openPMD_->writeParticlesAttributes( particles_group );
    openPMD_->writeSpeciesAttributes( species_group );
    
    // PositionOffset (for OpenPMD)
    string xyz = "xyz";
    H5Write positionoffset_group = species_group.group( "positionOffset" );
    openPMD_->writeRecordAttributes( positionoffset_group, SMILEI_UNIT_POSITION );
    vector<uint64_t> np = {0};
    for( unsigned int idim=0; idim<nDim_particle; idim++ ) {
        H5Write xyz_group = positionoffset_group.group( xyz.substr( idim, 1 ) );
        openPMD_->writeComponentAttributes( xyz_group, SMILEI_UNIT_POSITION );
        xyz_group.attr( "value", 0. );
        xyz_group.attr( "shape", np, H5T_NATIVE_UINT64 );
    }
    
    // Make empty datasets
    H5Space file_space( 0, 0, 0, 1000, true );
    if( write_id_ ) {
        loc_id_ = newDataset( species_group, "id", H5T_NATIVE_UINT64, file_space, SMILEI_UNIT_NONE );
        openPMD_->writeRecordAttributes( *loc_id_, SMILEI_UNIT_NONE );
    }
    if( write_charge_ ) {
        loc_charge_ = newDataset( species_group, "charge", H5T_NATIVE_SHORT, file_space, SMILEI_UNIT_CHARGE );
        openPMD_->writeRecordAttributes( *loc_charge_, SMILEI_UNIT_CHARGE );
    }
    if( write_any_position_ ) {
        H5Write position_group( &species_group, "position" );
        openPMD_->writeRecordAttributes( position_group, SMILEI_UNIT_POSITION );
        for( unsigned int idim=0; idim<nDim_particle; idim++ ) {
            if( write_position_[idim] ) {
                loc_position_[idim] = newDataset( position_group, xyz.substr( idim, 1 ).c_str(), H5T_NATIVE_DOUBLE, file_space, SMILEI_UNIT_POSITION );
            }
        }
    }
    if( write_any_momentum_ ) {
        H5Write momentum_group( &species_group, "momentum" );
        openPMD_->writeRecordAttributes( momentum_group, SMILEI_UNIT_MOMENTUM );
        for( unsigned int idim=0; idim<3; idim++ ) {
            if( write_momentum_[idim] ) {
                loc_momentum_[idim] = newDataset( momentum_group, xyz.substr( idim, 1 ).c_str(), H5T_NATIVE_DOUBLE, file_space, SMILEI_UNIT_MOMENTUM );
            }
        }
    }
    if( write_weight_ ) {
        loc_weight_ = newDataset( species_group, "weight", H5T_NATIVE_DOUBLE, file_space, SMILEI_UNIT_WEIGHT );
        openPMD_->writeRecordAttributes( *loc_weight_, SMILEI_UNIT_WEIGHT );
    }
    if( write_chi_ ) {
        loc_chi_ = newDataset( species_group, "chi", H5T_NATIVE_DOUBLE, file_space, SMILEI_UNIT_NONE );
        openPMD_->writeRecordAttributes( *loc_chi_, SMILEI_UNIT_NONE );
    }
    if( write_any_E_ ) {
        H5Write E_group( &species_group, "E" );
        openPMD_->writeRecordAttributes( E_group, SMILEI_UNIT_EFIELD );
        for( unsigned int idim=0; idim<3; idim++ ) {
            if( write_E_[idim] ) {
                loc_E_[idim] = newDataset( E_group, xyz.substr( idim, 1 ).c_str(), H5T_NATIVE_DOUBLE, file_space, SMILEI_UNIT_EFIELD );
            }
        }
    }
    if( write_any_B_ ) {
        H5Write B_group( &species_group, "B" );
        openPMD_->writeRecordAttributes( B_group, SMILEI_UNIT_BFIELD );
        for( unsigned int idim=0; idim<3; idim++ ) {
            if( write_B_[idim] ) {
                loc_B_[idim] = newDataset( B_group, xyz.substr( idim, 1 ).c_str(), H5T_NATIVE_DOUBLE, file_space, SMILEI_UNIT_BFIELD );
            }
        }
    }
    if( write_any_W_ ) {
        H5Write W_group( &species_group, "W" );
        openPMD_->writeRecordAttributes( W_group, SMILEI_UNIT_ENERGY );
        for( unsigned int idim=0; idim<3; idim++ ) {
            if( write_W_[idim] ) {
                loc_W_[idim] = newDataset( W_group, xyz.substr( idim, 1 ).c_str(), H5T_NATIVE_DOUBLE, file_space, SMILEI_UNIT_ENERGY );
            }
        }
    }
    loc_birth_time_ = newDataset( species_group, "birth_time", H5T_NATIVE_DOUBLE, file_space, SMILEI_UNIT_TIME );
    openPMD_->writeRecordAttributes( *loc_birth_time_, SMILEI_UNIT_TIME );
    
    // Make a dataset supposed to contain the number of particles written at every output iteration
    hsize_t ntimes = timeSelection->howManyTimesBefore( params.n_time ) + 1;
    H5Space fs( {0, 2}, {0, 0}, {0, 2}, {ntimes, 2}, {true, false} );
    iteration_npart_ = new H5Write( file_, "iteration_npart", H5T_NATIVE_INT64, &fs );
    
    file_->flush();
}

void DiagnosticNewParticles::closeFile()
{
    if( file_ ) {
        for( auto d : loc_position_ ) delete d;
        for( auto d : loc_momentum_ ) delete d;
        delete loc_id_;
        delete loc_charge_;
        delete loc_weight_;
        delete loc_chi_;
        for( auto d : loc_E_ ) delete d;
        for( auto d : loc_B_ ) delete d;
        for( auto d : loc_W_ ) delete d;
        delete loc_birth_time_;
        delete iteration_npart_;
        delete file_;
        file_ = nullptr;
    }
}

void DiagnosticNewParticles::init( Params &params, SmileiMPI *smpi, VectorPatch & )
{
    // create the file
    openFile( params, smpi );
}


H5Space * DiagnosticNewParticles::prepareH5( SimWindow *, SmileiMPI *smpi, int itime, uint32_t nParticles_local, uint64_t nParticles_global, uint64_t offset )
{
    // Resize datasets
    hsize_t new_size = nParticles_written + nParticles_global;
    loc_birth_time_->extend( new_size );
    if( write_id_ ) {
        loc_id_->extend( new_size );
    }
    if( write_charge_ ) {
        loc_charge_->extend( new_size );
    }
    if( write_any_position_ ) {
        for( unsigned int idim=0; idim<nDim_particle; idim++ ) {
            if( write_position_[idim] ) {
                loc_position_[idim]->extend( new_size );
            }
        }
    }
    if( write_any_momentum_ ) {
        for( unsigned int idim=0; idim<3; idim++ ) {
            if( write_momentum_[idim] ) {
                loc_momentum_[idim]->extend( new_size );
            }
        }
    }
    if( write_weight_ ) {
        loc_weight_->extend( new_size );
    }
    if( write_chi_ ) {
        loc_chi_->extend( new_size );
    }
    if( write_any_E_ ) {
        for( unsigned int idim=0; idim<3; idim++ ) {
            if( write_E_[idim] ) {
                loc_E_[idim]->extend( new_size );
            }
        }
    }
    if( write_any_B_ ) {
        for( unsigned int idim=0; idim<3; idim++ ) {
            if( write_B_[idim] ) {
                loc_B_[idim]->extend( new_size );
            }
        }
    }
    if( write_any_W_ ) {
        for( unsigned int idim=0; idim<3; idim++ ) {
            if( write_W_[idim] ) {
                loc_W_[idim]->extend( new_size );
            }
        }
    }
    // update iteration_npart_
    iteration_npart_->extend( {nTimes_written+1, 2} );
    if( smpi->isMaster() ) {
        H5Space filespace( {nTimes_written+1, 2}, {nTimes_written, 0}, {1, 2} );
        H5Space memspace( 2 );
        uint64_t i_n[2] = { (uint64_t) itime, new_size };
        iteration_npart_->write( i_n[0], H5T_NATIVE_UINT64, &filespace, &memspace, true );
    }
    nTimes_written += 1;
    
    // Filespace
    hsize_t full_offset = nParticles_written + offset;
    nParticles_written = new_size;
    return new H5Space( new_size, full_offset, nParticles_local );
}

void DiagnosticNewParticles::writeOther( VectorPatch &vecPatches, size_t iprop, H5Space *file_space, H5Space *mem_space )
{
    fill_buffer( vecPatches, iprop, data_double );
    #pragma omp master
    write_scalar_double( loc_birth_time_, "birth_time", data_double[0], file_space, mem_space, SMILEI_UNIT_NONE );
    
    // Delete the records
    #pragma omp for schedule(static)
    for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
        vecPatches.species( ipatch, species_index_ )->birth_records_->clear();
    }
}


void DiagnosticNewParticles::write_scalar_uint64( H5Write * location, string /*name*/, uint64_t &buffer, H5Space *file_space, H5Space *mem_space, unsigned int /*unit_type*/ )
{
    location->write( buffer, H5T_NATIVE_UINT64, file_space, mem_space );
}
void DiagnosticNewParticles::write_scalar_short( H5Write * location, string /*name*/, short &buffer, H5Space *file_space, H5Space *mem_space, unsigned int /*unit_type*/ )
{
    location->write( buffer, H5T_NATIVE_SHORT, file_space, mem_space );
}
void DiagnosticNewParticles::write_scalar_double( H5Write * location, string /*name*/, double &buffer, H5Space *file_space, H5Space *mem_space, unsigned int /*unit_type*/ )
{
    location->write( buffer, H5T_NATIVE_DOUBLE, file_space, mem_space );
}

void DiagnosticNewParticles::write_component_uint64( H5Write * location, string /*name*/, uint64_t &buffer, H5Space *file_space, H5Space *mem_space, unsigned int /*unit_type*/ )
{
    location->write( buffer, H5T_NATIVE_UINT64, file_space, mem_space );
}
void DiagnosticNewParticles::write_component_short( H5Write * location, string /*name*/, short &buffer, H5Space *file_space, H5Space *mem_space, unsigned int /*unit_type*/ )
{
    location->write( buffer, H5T_NATIVE_SHORT, file_space, mem_space );
}
void DiagnosticNewParticles::write_component_double( H5Write * location, string /*name*/, double &buffer, H5Space *file_space, H5Space *mem_space, unsigned int /*unit_type*/ )
{
    location->write( buffer, H5T_NATIVE_DOUBLE, file_space, mem_space );
}


// SUPPOSED TO BE EXECUTED ONLY BY MASTER MPI
uint64_t DiagnosticNewParticles::getDiskFootPrint( int, int, Patch * )
{
    return 0;
}

Particles * DiagnosticNewParticles::getParticles( Patch * patch )
{
    return &patch->vecSpecies[species_index_]->birth_records_->p_;
}
