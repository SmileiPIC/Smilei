#include "PyTools.h"

#include <string>
#include <sstream>

#include "ParticleData.h"
#include "PeekAtSpecies.h"
#include "DiagnosticParticleList.h"
#include "VectorPatch.h"
#include "Params.h"

using namespace std;

DiagnosticParticleList::DiagnosticParticleList( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, string diag_type, string file_prefix, unsigned int idiag_of_this_type, OpenPMDparams &oPMD ) :
    Diagnostic( &oPMD, diag_type, idiag_of_this_type ),
    nDim_particle( params.nDim_particle )
{
    
    // Extract the species
    PyTools::extract( "species", species_name_, diag_type, idiag_of_this_type );
    vector<string> species_names = {species_name_};
    vector<unsigned int> species_ids = Params::FindSpecies( vecPatches( 0 )->vecSpecies, species_names );
    if( species_ids.size() > 1 ) {
        ERROR( diag_type << " #" << idiag_of_this_type << " corresponds to more than 1 species" );
    }
    if( species_ids.size() < 1 ) {
        ERROR( diag_type << " #" << idiag_of_this_type << " does not correspond to any existing species" );
    }
    species_index_ = species_ids[0];
    
    ostringstream name( "" );
    name << diag_type << " with species '" << species_name_ << "'";
    
    // Get parameter "every" which describes an iteration selection
    timeSelection = new TimeSelection( PyTools::extract_py( "every", diag_type, idiag_of_this_type ), name.str() );
    
    // Get parameter "flush_every" which decides the file flushing time selection
    flush_timeSelection = new TimeSelection( PyTools::extract_py( "flush_every", diag_type, idiag_of_this_type ), name.str() );
    
    // Get parameter "filter" which gives a python function to select particles
    filter = PyTools::extract_py( "filter", diag_type, idiag_of_this_type );
    has_filter = ( filter != Py_None );
    if( has_filter ) {
#ifdef SMILEI_USE_NUMPY
        // Test the filter with temporary, "fake" particles
        name << " filter:";
        bool *dummy = NULL;
        ParticleData test( nDim_particle, filter, name.str(), dummy );
#else
        ERROR( name.str() << " with a filter requires the numpy package" );
#endif
    }
    
    // Get the parameter "attributes": a list of attribute name that must be written
    vector<string> attributes( 0 );
    if( !PyTools::extractV( "attributes", attributes, diag_type, idiag_of_this_type ) ) {
        ERROR( diag_type << " #" << idiag_of_this_type << ": argument `attribute` must be a list of strings" );
    }
    if( attributes.size() == 0 ) {
        ERROR( diag_type << " #" << idiag_of_this_type << ": argument `attribute` must have at least one element" );
    }
    ostringstream attr_list( "" );
    InterpolatedFields * interpolated_fields = vecPatches( 0 )->vecSpecies[species_index_]->particles->interpolated_fields_;
    for( unsigned int i=0; i<attributes.size(); i++ ) {
        if( attributes[i] == "id" ) {
            write_id_ = true;
        } else if( attributes[i] == "x" ) {
            write_position_[0] = true;
        } else if( attributes[i] == "y" ) {
            if( nDim_particle>1 ) {
                write_position_[1] = true;
            } else {
                continue;
            }
        } else if( attributes[i] == "z" ) {
            if( nDim_particle>2 ) {
                write_position_[2] = true;
            } else {
                continue;
            }
        } else if( attributes[i] == "px" ) {
            write_momentum_[0] = true;
        } else if( attributes[i] == "py" ) {
            write_momentum_[1] = true;
        } else if( attributes[i] == "pz" ) {
            write_momentum_[2] = true;
        } else if( attributes[i] == "charge" || attributes[i] == "q" ) {
            write_charge_      = true;
        } else if( attributes[i] == "weight" || attributes[i] == "w" ) {
            write_weight_      = true;
        } else if( attributes[i] == "chi" ) {
            write_chi_         = true;
        } else if( attributes[i] == "Ex" ) {
            write_E_[0]        = true;
            interpolate_ = interpolate_ || !( interpolated_fields && interpolated_fields->mode_[0] > 0 );
        } else if( attributes[i] == "Ey" ) {
            write_E_[1]        = true;
            interpolate_ = interpolate_ || !( interpolated_fields && interpolated_fields->mode_[1] > 0 );
        } else if( attributes[i] == "Ez" ) {
            write_E_[2]        = true;
            interpolate_ = interpolate_ || !( interpolated_fields && interpolated_fields->mode_[2] > 0 );
        } else if( attributes[i] == "Bx" ) {
            write_B_[0]        = true;
            interpolate_ = interpolate_ || !( interpolated_fields && interpolated_fields->mode_[3] > 0 );
        } else if( attributes[i] == "By" ) {
            write_B_[1]        = true;
            interpolate_ = interpolate_ || !( interpolated_fields && interpolated_fields->mode_[4] > 0 );
        } else if( attributes[i] == "Bz" ) {
            write_B_[2]        = true;
            interpolate_ = interpolate_ || !( interpolated_fields && interpolated_fields->mode_[5] > 0 );
        } else if( attributes[i] == "Wx" ) {
            write_W_[0]        = true;
            if( ! interpolated_fields || interpolated_fields->mode_[6] == 0 ) {
                ERROR( diag_type << " #" << idiag_of_this_type << ": attribute Wx requires `keep_interpolated_fields` to include 'Wx'" );
            }
        } else if( attributes[i] == "Wy" ) {
            write_W_[1]        = true;
            if( ! interpolated_fields || interpolated_fields->mode_[7] == 0 ) {
                ERROR( diag_type << " #" << idiag_of_this_type << ": attribute Wy requires `keep_interpolated_fields` to include 'Wy'" );
            }
        } else if( attributes[i] == "Wz" ) {
            write_W_[2]        = true;
            if( ! interpolated_fields || interpolated_fields->mode_[8] == 0 ) {
                ERROR( diag_type << " #" << idiag_of_this_type << ": attribute Wz requires `keep_interpolated_fields` to include 'Wz'" );
            }
        } else {
            ERROR( diag_type << " #" << idiag_of_this_type << ": attribute `" << attributes[i] << "` unknown" );
        }
        attr_list << "," << attributes[i];
    }
    write_any_position_ = write_position_[0] || write_position_[1] || write_position_[2];
    write_any_momentum_ = write_momentum_[0] || write_momentum_[1] || write_momentum_[2];
    write_any_E_ = write_E_[0] || write_E_[1] || write_E_[2];
    write_any_B_ = write_B_[0] || write_B_[1] || write_B_[2];
    write_any_W_ = write_W_[0] || write_W_[1] || write_W_[2];
    if( write_chi_ && ! vecPatches( 0 )->vecSpecies[species_index_]->particles->has_quantum_parameter ) {
        ERROR( diag_type << " #" << idiag_of_this_type << ": attribute `chi` not available for this species" );
    }
    
    // Create the filename
    ostringstream hdf_filename( "" );
    hdf_filename << file_prefix << species_name_  << ".h5" ;
    filename = hdf_filename.str();
    
    // Print some info
    if( smpi->isMaster() ) {
        MESSAGE( 1, "Created " << diag_type << " #" << idiag_of_this_type << ": species " << species_name_ );
        MESSAGE( 2, attr_list.str() );
    }
}

DiagnosticParticleList::~DiagnosticParticleList()
{
    delete timeSelection;
    delete flush_timeSelection;
    Py_DECREF( filter );
}

bool DiagnosticParticleList::prepare( int itime )
{
    return timeSelection->theTimeIsNow( itime );
}


void DiagnosticParticleList::run( SmileiMPI *smpi, VectorPatch &vecPatches, int itime, SimWindow *simWindow, Timers & )
{
    uint64_t nParticles_global = 0;
    string xyz = "xyz";
    
    H5Space *file_space=NULL, *mem_space=NULL;
    #pragma omp master
    {
        // Obtain the particle partition of all the patches in this MPI
        nParticles_local = 0;
        patch_start.resize( vecPatches.size() );
        
        if( has_filter ) {
#ifdef SMILEI_USE_NUMPY
            patch_selection.resize( vecPatches.size() );
            PyArrayObject *ret;
            ParticleData particleData( 0 );
            for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
                patch_selection[ipatch].resize( 0 );
                Particles *p = getParticles( vecPatches( ipatch ) );
                unsigned int npart = p->numberOfParticles();
                if( npart > 0 ) {
                    // Expose particle data as numpy arrays
                    particleData.resize( npart );
                    particleData.set( p );
                    // run the filter function
                    ret = ( PyArrayObject * )PyObject_CallFunctionObjArgs( filter, particleData.get(), NULL );
                    PyTools::checkPyError();
                    particleData.clear();
                    if( ret == NULL ) {
                        ERROR( "A DiagTrackParticles filter has not provided a correct result" );
                    }
                    // Loop the return value and store the selected particle index
                    bool *arr = ( bool * ) PyArray_GETPTR1( ret, 0 );
                    for( unsigned int i=0; i<npart; i++ ) {
                        if( arr[i] ) {
                            patch_selection[ipatch].push_back( i );
                        }
                    }
                    // Apply changes to filtered particles
                    modifyFiltered( vecPatches, ipatch );
                    Py_DECREF( ret );
                }
                patch_start[ipatch] = nParticles_local;
                nParticles_local += patch_selection[ipatch].size();
            }
#endif
        } else {
            for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
                patch_start[ipatch] = nParticles_local;
                nParticles_local += getParticles( vecPatches( ipatch ) )->numberOfParticles();
            }
        }
        
        // Specify the memory dataspace (the size of the local buffer)
        mem_space = new H5Space( (hsize_t)nParticles_local );
        
        // Get the number of offset for this MPI rank
        uint64_t np_local = nParticles_local, offset;
        MPI_Scan( &np_local, &offset, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD );
        nParticles_global = offset;
        offset -= np_local;
        MPI_Bcast( &nParticles_global, 1, MPI_UNSIGNED_LONG_LONG, smpi->getSize()-1, MPI_COMM_WORLD );
        
        // Prepare all HDF5 groups and datasets
        file_space = prepareH5( simWindow, smpi, itime, nParticles_local, nParticles_global, offset );
    }
    
    // Id
    if( write_id_ ) {
        #pragma omp master
        data_uint64.resize( nParticles_local );
        fill_buffer( vecPatches, 0, data_uint64 );
        #pragma omp master
        {
            write_scalar_uint64( loc_id_, "id", data_uint64[0], file_space, mem_space, SMILEI_UNIT_NONE );
            data_uint64.resize( 0 );
        }
    }
    
    // Charge
    if( write_charge_ ) {
        #pragma omp master
        data_short.resize( nParticles_local );
        fill_buffer( vecPatches, 0, data_short );
        #pragma omp master
        {
            write_scalar_short( loc_charge_, "charge", data_short[0], file_space, mem_space, SMILEI_UNIT_CHARGE );
            data_short.resize( 0 );
        }
    }
    
    #pragma omp master
    data_double.resize( nParticles_local );
    
    // Position
    if( write_any_position_ ) {
        for( unsigned int idim=0; idim<nDim_particle; idim++ ) {
            if( write_position_[idim] ) {
                fill_buffer( vecPatches, idim, data_double );
                #pragma omp master
                write_component_double( loc_position_[idim], xyz.substr( idim, 1 ).c_str(), data_double[0], file_space, mem_space, SMILEI_UNIT_POSITION );
            }
        }
    }
    
    // Momentum
    if( write_any_momentum_ ) {
        for( unsigned int idim=0; idim<3; idim++ ) {
            if( write_momentum_[idim] ) {
                fill_buffer( vecPatches, nDim_particle+idim, data_double );
                #pragma omp master
                {
                    // Multiply by the mass to obtain an actual momentum (except for photons (mass = 0))
                    if( vecPatches( 0 )->vecSpecies[species_index_]->mass_ != 1. &&
                        vecPatches( 0 )->vecSpecies[species_index_]->mass_ > 0) {
                        for( unsigned int ip=0; ip<nParticles_local; ip++ ) {
                            data_double[ip] *= vecPatches( 0 )->vecSpecies[species_index_]->mass_;
                        }
                    }
                    write_component_double( loc_momentum_[idim], xyz.substr( idim, 1 ).c_str(), data_double[0], file_space, mem_space, SMILEI_UNIT_MOMENTUM );
                }
            }
        }
    }
    
    // Properties included in Particles::double_prop are ordered a certain way that must be followed with care.
    // iprop is the index of the property currently being treated
    // positions (x, y, z) are indexed from 0 to nDim_particles
    // momenta (px, py, pz) are indexed from nDim_particles to nDim_particles+3
    // etc
    size_t iprop = nDim_particle+3;
    
    // Weight
    if( write_weight_ ) {
        fill_buffer( vecPatches, iprop, data_double );
        #pragma omp master
        write_scalar_double( loc_weight_, "weight", data_double[0], file_space, mem_space, SMILEI_UNIT_WEIGHT );
    }
    
    iprop++;
    // If position_old exist, skip its components
    if( ! getParticles( vecPatches( 0 ) )->Position_old.empty() ) {
        iprop += nDim_particle;
    }
    
    // Chi - quantum parameter
    if( write_chi_ ) {
        fill_buffer( vecPatches, iprop, data_double );
        #pragma omp master
        write_scalar_double( loc_chi_, "chi", data_double[0], file_space, mem_space, SMILEI_UNIT_NONE );
    }
    
    #pragma omp barrier
    if( getParticles( vecPatches( 0 ) )->has_quantum_parameter ) {
        iprop++;
    }
    if( getParticles( vecPatches( 0 ) )->has_Monte_Carlo_process ) {
        iprop++;
    }
    
    InterpolatedFields * interpolated_fields = getParticles( vecPatches( 0 ) )->interpolated_fields_;
    
    // If field interpolation necessary
    if( interpolate_ ) {
        
        #pragma omp master
        data_double.resize( nParticles_local*6 );
        
        // Do the interpolation
        unsigned int nPatches = vecPatches.size();
        #pragma omp barrier
        
        if( has_filter ) {
            #pragma omp for schedule(static)
            for( unsigned int ipatch=0 ; ipatch<nPatches ; ipatch++ ) {
                vecPatches.species( ipatch, species_index_ )->Interp->fieldsSelection(
                    vecPatches.emfields( ipatch ),
                    *getParticles( vecPatches( ipatch ) ),
                    &data_double[patch_start[ipatch]],
                    ( int ) nParticles_local,
                    &patch_selection[ipatch]
                );
            }
        } else {
            #pragma omp for schedule(static)
            for( unsigned int ipatch=0 ; ipatch<nPatches ; ipatch++ ) {
                vecPatches.species( ipatch, species_index_ )->Interp->fieldsSelection(
                    vecPatches.emfields( ipatch ),
                    *getParticles( vecPatches( ipatch ) ),
                    &data_double[patch_start[ipatch]],
                    ( int ) nParticles_local,
                    NULL
                );
            }
        }
        #pragma omp barrier
        
        // Write out the fields
        #pragma omp master
        {
            if( write_any_E_ ) {
                for( unsigned int idim=0; idim<3; idim++ ) {
                    if( write_E_[idim] ) {
                        write_component_double( loc_E_[idim], xyz.substr( idim, 1 ).c_str(), data_double[idim*nParticles_local], file_space, mem_space, SMILEI_UNIT_EFIELD );
                    }
                }
            }
            
            if( write_any_B_ ) {
                for( unsigned int idim=0; idim<3; idim++ ) {
                    if( write_B_[idim] ) {
                        write_component_double( loc_B_[idim], xyz.substr( idim, 1 ).c_str(), data_double[( 3+idim )*nParticles_local], file_space, mem_space, SMILEI_UNIT_BFIELD );
                    }
                }
            }
        }
        
        if( interpolated_fields ) {
            iprop += count( interpolated_fields->mode_.begin(), interpolated_fields->mode_.end(), 1 );
        }
        
    // If field interpolation not necessary, maybe fields have been stored using `keep_interpolated_fields`
    } else if( write_any_E_ || write_any_B_ ) {
        
        for( unsigned int idim=0; idim<3; idim++ ) {
            if( write_E_[idim] ) {
                fill_buffer( vecPatches, iprop, data_double );
                #pragma omp master
                write_component_double( loc_E_[idim], xyz.substr( idim, 1 ).c_str(), data_double[0], file_space, mem_space, SMILEI_UNIT_EFIELD );
            }
            if( interpolated_fields->mode_[idim] > 0 ) {
                iprop++;
            }
        }
        
        for( unsigned int idim=0; idim<3; idim++ ) {
            if( write_B_[idim] ) {
                fill_buffer( vecPatches, iprop, data_double );
                #pragma omp master
                write_component_double( loc_B_[idim], xyz.substr( idim, 1 ).c_str(), data_double[0], file_space, mem_space, SMILEI_UNIT_BFIELD );
            }
            if( interpolated_fields->mode_[3+idim] > 0 ) {
                iprop++;
            }
        }
    
    // If no fields required, just count the offset they create
    } else if( interpolated_fields ) {
        iprop += count( interpolated_fields->mode_.begin(), interpolated_fields->mode_.end(), 1 );
    }
    
    if( write_any_W_ ) {
        for( unsigned int idim=0; idim<3; idim++ ) {
            if( write_W_[idim] ) {
                fill_buffer( vecPatches, iprop, data_double );
                #pragma omp master
                write_component_double( loc_W_[idim], xyz.substr( idim, 1 ).c_str(), data_double[0], file_space, mem_space, SMILEI_UNIT_ENERGY );
            }
            if( interpolated_fields->mode_[6+idim] > 0 ) {
                iprop++;
            }
        }
    } else if( interpolated_fields ) {
        iprop += count( interpolated_fields->mode_.begin()+6, interpolated_fields->mode_.end(), 2 );
    }
    
    writeOther( vecPatches, iprop, file_space, mem_space );
    
    #pragma omp master
    {
        data_double.resize( 0 );
        
        // Close and flush
        patch_selection.resize( 0 );
        
        delete file_space;
        delete mem_space;
        deleteH5();
        
        if( flush_timeSelection->theTimeIsNow( itime ) ) {
            file_->flush();
        }
    }
    #pragma omp barrier
}

