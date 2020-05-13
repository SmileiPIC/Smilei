// -----------------------------------------------------------------------------
//
//! \file ParticleCreator.cpp
//
//! \brief Class with functions to create particles
//
// -----------------------------------------------------------------------------

#include "ParticleCreator.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
//! constructor for ParticleCreator
// ---------------------------------------------------------------------------------------------------------------------
ParticleCreator::ParticleCreator()
{
    position_initialization_on_species_ = false;
    initialized_in_species_ = true;
    time_profile_ = NULL;
}

// ---------------------------------------------------------------------------------------------------------------------
//! Destructor for ParticleCreator
// ---------------------------------------------------------------------------------------------------------------------
ParticleCreator::~ParticleCreator() {
    // for( unsigned int i=0; i<velocity_profile_.size(); i++ ) {
    //     delete velocity_profile_[i];
    // }
    // for( unsigned int i=0; i<temperature_profile_.size(); i++ ) {
    //     delete temperature_profile_[i];
    // }
    // if( density_profile_ ) {
    //     delete density_profile_;
    // }
    // if( time_profile_ ) {
    //     delete time_profile_;
    // }
    // if( particles_per_cell_profile_ ) {
    //     delete particles_per_cell_profile_;
    // }
}

// ---------------------------------------------------------------------------------------------------------------------
//! \brief Associate this particle creator object to the specified particle injector
//! \param particle_injector : particle injector to associate with
//! \param particles : particles object to inject particles
//! \param species : associated species
// ---------------------------------------------------------------------------------------------------------------------
void ParticleCreator::associate( ParticleInjector * particle_injector, Particles * particles, Species * species)
{
    species_ = species;
    
    particles_ = particles;
    
    // If we do not use the particles object associated to the species
    if (&particles_ != &species->particles)
    {
        initialized_in_species_ = false;
    }
    
    position_initialization_ = particle_injector->position_initialization_;
    position_initialization_on_species_ = particle_injector->position_initialization_on_injector_;
    momentum_initialization_ = particle_injector->momentum_initialization_;
    velocity_profile_.resize(particle_injector->velocity_profile_.size());
    for (unsigned int i = 0 ; i < velocity_profile_.size() ; i++) {
        velocity_profile_[i] = particle_injector->velocity_profile_[i];
    }
    temperature_profile_.resize(particle_injector->temperature_profile_.size());
    for (unsigned int i = 0 ; i < temperature_profile_.size() ; i++) {
        temperature_profile_[i] = particle_injector->temperature_profile_[i];
    }
    density_profile_ = particle_injector->density_profile_;
    density_profile_type_ = particle_injector->density_profile_type_;
    time_profile_ = particle_injector->time_profile_;
    particles_per_cell_profile_ = particle_injector->particles_per_cell_profile_;
}

// ---------------------------------------------------------------------------------------------------------------------
//! \brief Associate this particle creator object to the specified particle injector
//! \param species :species to associate with
// ---------------------------------------------------------------------------------------------------------------------
void ParticleCreator::associate( Species * species)
{
    species_ = species;
    
    particles_ = species->particles;
    
    position_initialization_ = species->position_initialization_;
    position_initialization_on_species_ = species->position_initialization_on_species_;
    position_initialization_on_species_type_ = species->position_initialization_on_species_type_;
    momentum_initialization_ = species->momentum_initialization_;
    velocity_profile_.resize(species->velocity_profile_.size());
    for (unsigned int i = 0 ; i < velocity_profile_.size() ; i++) {
        velocity_profile_[i] = species->velocity_profile_[i];
    }
    temperature_profile_.resize(species->temperature_profile_.size());
    for (unsigned int i = 0 ; i < temperature_profile_.size() ; i++) {
        temperature_profile_[i] = species->temperature_profile_[i];
    }
    density_profile_ = species->density_profile_;
    density_profile_type_ = species->density_profile_type_;
    particles_per_cell_profile_ = species->particles_per_cell_profile_;
}

// ---------------------------------------------------------------------------------------------------------------------
//! \brief Creation of the particle properties in the given particle vector `particles`
//! \param particles : vector of particles
//! \param species : species object necessary for some properties
//! \param n_space_to_create
//! \param params : general parameters
// ---------------------------------------------------------------------------------------------------------------------
int ParticleCreator::create( std::vector<unsigned int> n_space_to_create,
                             Params &params,
                             Patch *patch,
                             int new_cell_idx,
                             unsigned int itime)
{
    
    // n_space_to_create_generalized = n_space_to_create, + copy of 2nd direction data among 3rd direction
    // same for local species_::cell_length[2]
    std::vector<unsigned int> n_space_to_create_generalized( n_space_to_create );
    unsigned int nPart, i, j, k;
    unsigned int npart_effective = 0 ;
    double *momentum[3], *position[species_->nDim_particle], *weight_arr;
    std::vector<int> my_particles_indices;
    std::vector<Field *> xyz( species_->nDim_field );
    
    // Create particles_ in a space starting at cell_position
    std::vector<double> cell_position( 3, 0 );
    std::vector<double> cell_index( 3, 0 );
    for( unsigned int idim=0 ; idim<species_->nDim_field ; idim++ ) {
        //if (params.cell_length[idim]!=0) { // Useless, nDim_field defined for (params.cell_length[idim>=nDim_field]==0)
        cell_position[idim] = patch->getDomainLocalMin( idim );
        cell_index   [idim] = ( double ) patch->getCellStartingGlobalIndex( idim );
        xyz[idim] = new Field3D( n_space_to_create_generalized );
        //}
    }
    // Create the x,y,z maps where profiles will be evaluated
    std::vector<double> ijk( 3 );
    for( ijk[0]=0; ijk[0]<n_space_to_create_generalized[0]; ijk[0]++ ) {
        for( ijk[1]=0; ijk[1]<n_space_to_create_generalized[1]; ijk[1]++ ) {
            for( ijk[2]=0; ijk[2]<n_space_to_create_generalized[2]; ijk[2]++ ) {
                for( unsigned int idim=0 ; idim<species_->nDim_field ; idim++ ) {
                    ( *xyz[idim] )( ijk[0], ijk[1], ijk[2] ) = cell_position[idim] + ( ijk[idim]+0.5 )*species_->cell_length[idim];
                }
                ( *xyz[0] )( ijk[0], ijk[1], ijk[2] ) += new_cell_idx*species_->cell_length[0];
            }
        }
    }
    
    // ---------------------------------------------------------
    // Calculate density and number of particles_ for the species_
    // ---------------------------------------------------------

    // field containing the charge distribution (always 3d)
    Field3D charge( n_space_to_create_generalized );
    species_->max_charge_ = 0.;

    // field containing the number of particles in each cell
    Field3D n_part_in_cell( n_space_to_create_generalized );

    // field containing the density distribution (always 3d)
    Field3D density( n_space_to_create_generalized );

    // field containing the temperature distribution along all 3 momentum coordinates (always 3d * 3)
    Field3D temperature[3];
    // field containing the temperature distribution along all 3 momentum coordinates (always 3d * 3)
    Field3D velocity[3];
    
    if( species_->momentum_initialization_array_ != NULL ) {
        for( unsigned int idim = 0; idim < 3; idim++ ) {
            momentum[idim] = &( species_->momentum_initialization_array_[idim*species_->n_numpy_particles_] );
        }
    } else {
        //Initialize velocity and temperature profiles
        for( i=0; i<3; i++ ) {
            velocity[i].allocateDims( n_space_to_create_generalized );
            temperature[i].allocateDims( n_space_to_create_generalized );
        }
        // Evaluate profiles
        for( unsigned int m=0; m<3; m++ ) {
            if( temperature_profile_[m] ) {
                temperature_profile_[m]->valuesAt( xyz, temperature[m] );
            } else {
                temperature[m].put_to( 0.0000000001 ); // default value
            }

            if( velocity_profile_[m] ) {
                velocity_profile_[m]   ->valuesAt( xyz, velocity   [m] );
            } else {
                velocity[m].put_to( 0.0 ); //default value
            }
        }
    } // end if momentum_initialization_array_
    
    // Initialize charge profile
    if( species_->mass_ > 0 ) {
        species_->charge_profile_->valuesAt( xyz, charge );
    }
    if( species_->position_initialization_array_ != NULL ) {
        for( unsigned int idim = 0; idim < species_->nDim_particle; idim++ ) {
            position[idim] = &( species_->position_initialization_array_[idim*species_->n_numpy_particles_] );
        }
        weight_arr =         &( species_->position_initialization_array_[species_->nDim_particle*species_->n_numpy_particles_] );
        //Idea to speed up selection, provides xmin, xmax of the bunch and check if there is an intersection with the patch instead of going through all particles for all patches.
        for( unsigned int ip = 0; ip < species_->n_numpy_particles_; ip++ ) {
            //If the particle belongs to this patch
            if (params.geometry!="AMcylindrical") {
                if( position[0][ip] >= patch->getDomainLocalMin( 0 ) && position[0][ip] < patch->getDomainLocalMax( 0 )
                    && ( species_->nDim_particle < 2  || ( position[1][ip] >= patch->getDomainLocalMin( 1 ) && position[1][ip] < patch->getDomainLocalMax( 1 ) ) )
                    && ( species_->nDim_particle < 3  || ( position[2][ip] >= patch->getDomainLocalMin( 2 ) && position[2][ip] < patch->getDomainLocalMax( 2 ) ) ) ) {
                    my_particles_indices.push_back( ip ); //This vector stores particles initially sittinig in the current patch.
                }
            }
            else {
                double distance =  sqrt( position[1][ip]*position[1][ip]+position[2][ip]*position[2][ip] );
                if( position[0][ip] >= patch->getDomainLocalMin( 0 ) && position[0][ip] < patch->getDomainLocalMax( 0 )
                    && ( distance >= patch->getDomainLocalMin( 1 ) && distance < patch->getDomainLocalMax( 1 ) ) ) {
                    my_particles_indices.push_back( ip ); //This vector stores particles initially sittinig in the current patch.
                }
            }
        }
        npart_effective = my_particles_indices.size();
    } else {
        //Initialize density and ppc profiles
        density_profile_->valuesAt( xyz, density );
        particles_per_cell_profile_->valuesAt( xyz, n_part_in_cell );
        // Take into account the time profile
        double time_amplitude;
        if (time_profile_) {
            time_amplitude = time_profile_->valueAt(itime*params.timestep);
        } else {
            time_amplitude = 1;
        }
        weight_arr = NULL;
        //Now compute number of particles per cell
        double remainder, nppc;
        for( i=0; i<n_space_to_create_generalized[0]; i++ ) {
            for( j=0; j<n_space_to_create_generalized[1]; j++ ) {
                for( k=0; k<n_space_to_create_generalized[2]; k++ ) {

                    // Obtain the number of particles per cell
                    nppc = n_part_in_cell( i, j, k );
                    
                    n_part_in_cell( i, j, k ) = floor( nppc );
                    
                    // If not a round number, then we need to decide how to round
                    double intpart;
                    if( modf( nppc, &intpart ) > 0 ) {
                        remainder = pow( nppc - floor( nppc ), -species_->inv_nDim_particles);
                        if( fmod( cell_index[0]+( double )i, remainder ) < 1.
                                && fmod( cell_index[1]+( double )j, remainder ) < 1.
                                && fmod( cell_index[2]+( double )k, remainder ) < 1. ) {
                            n_part_in_cell( i, j, k )++;
                        }
                    }

                    // assign charge its correct value in the cell
                    if( species_->mass_ > 0 ) {
                        if( charge( i, j, k )>species_->max_charge_ ) {
                            species_->max_charge_=charge( i, j, k );
                        }
                    }

                    // If zero or less, zero particles
                    if( n_part_in_cell( i, j, k )<=0. || density( i, j, k )==0. ) {
                        n_part_in_cell( i, j, k ) = 0.;
                        density( i, j, k ) = 0.;
                        continue;
                    }

                    // assign density its correct value in the cell
                    if( density_profile_type_=="charge" ) {
                        if( charge( i, j, k )==0. ) {
                            ERROR( "Encountered non-zero charge density and zero charge at the same location" );
                        }
                        density( i, j, k ) /= charge( i, j, k );
                    }
                    
                    density( i, j, k ) = abs( density( i, j, k ) );
                    
                    // Time amplitude
                    density( i, j, k ) *= time_amplitude;
                    
                    // multiply by the cell volume
                    density( i, j, k ) *= params.cell_volume;

                    // increment the effective number of particle by n_part_in_cell(i,j,k)
                    // for each cell with as non-zero density
                    npart_effective += ( unsigned int ) n_part_in_cell( i, j, k );

                }//i
            }//j
        }//k end the loop on all cells
    }
    
    unsigned int n_existing_particles = particles_->size();
    
    //if (!n_existing_particles)
    particles_->initialize( n_existing_particles+npart_effective, species_->nDim_particle );

    // Initialization of the particles properties
    // ------------------------------------------
    unsigned int iPart=n_existing_particles;
    double *indexes=new double[species_->nDim_particle];
    double *temp=new double[3];
    double *vel=new double[3];
    
    if( species_->position_initialization_array_ == NULL ) {
        for( i=0; i<n_space_to_create_generalized[0]; i++ ) {
            if(( !n_existing_particles )&&( i%species_->clrw == 0 )&&( initialized_in_species_ )) {
                species_->particles->first_index[(new_cell_idx+i)/species_->clrw] = iPart;
            }
            for( j=0; j<n_space_to_create_generalized[1]; j++ ) {
                for( k=0; k<n_space_to_create_generalized[2]; k++ ) {
                    // initialize particles in meshes where the density is non-zero
                    if( density( i, j, k )>0 ) {

                        vel[0]  = velocity[0]( i, j, k );
                        vel[1]  = velocity[1]( i, j, k );
                        vel[2]  = velocity[2]( i, j, k );
                        temp[0] = temperature[0]( i, j, k );
                        temp[1] = temperature[1]( i, j, k );
                        temp[2] = temperature[2]( i, j, k );
                        nPart = n_part_in_cell( i, j, k );
                        
                        //cerr << vel[0] << endl;
                        
                        //if (n_existing_particles) {
                        //    iPart = n_existing_particles;
                        //    iPart = particles->first_index[(new_cell_idx+i)/clrw];
                        //    particles->last_index[(new_cell_idx+i)/clrw] += nPart;
                        //    for ( int idx=(new_cell_idx+i)/clrw+1 ; idx<particles->last_index.size() ; idx++ ) {
                        //        particles->first_index[idx] += nPart;
                        //        particles->last_index[idx] += nPart;
                        //    }
                        //    particles->createParticles( nPart, iPart );
                        //
                        //}

                        indexes[0]=i*species_->cell_length[0]+cell_position[0] + new_cell_idx*species_->cell_length[0];;
                        if( species_->nDim_particle > 1 ) {
                            indexes[1]=j*species_->cell_length[1]+cell_position[1];
                            if( species_->nDim_particle > 2 ) {
                                indexes[2]=k*species_->cell_length[2]+cell_position[2];
                            }
                        }
                        if( !position_initialization_on_species_ ) {
                            ParticleCreator::createPosition( position_initialization_, particles_, species_, nPart, iPart, indexes, params );
                        }
                        ParticleCreator::createMomentum( momentum_initialization_, particles_, species_,  nPart, iPart, temp, vel );
                        
                        ParticleCreator::createWeight( position_initialization_, particles_, nPart, iPart, density( i, j, k ), params );

                        ParticleCreator::createCharge( particles_, species_, nPart, iPart, charge( i, j, k ) );

                        //if (n_existing_particles) {
                        //    // operate filter
                        //    for ( int ip = nPart-1 ; ip >= 0 ; ip-- ){
                        //        double lf=1;
                        //        for (int iDim=0;iDim<3;iDim++)
                        //            lf += particles->Momentum[iDim][iPart+ip]*particles->Momentum[iDim][iPart+ip];
                        //        lf =sqrt( lf );
                        //        double X = particles->Position[0][iPart+ip] - cell_length[0]+params.timestep*particles->Momentum[0][iPart+ip]/lf;
                        //        if ( patch->isXmin() ) {
                        //            if ( ( X < 0. ) ) {
                        //                nPart--; // ne sert à rien ici
                        //                particles->eraseParticle(iPart+ip);
                        //                particles->last_index[(new_cell_idx+i)/clrw]--;
                        //                for ( int idx=(new_cell_idx+i)/clrw+1 ; idx<particles->last_index.size() ; idx++ ) {
                        //                    particles->first_index[idx]--;
                        //                    particles->last_index[idx]--;
                        //                }
                        //            }
                        //            else
                        //                particles->Position[0][iPart+ip] = X;
                        //        }
                        //        else if ( patch->isXmax() ) {
                        //            X = particles->Position[0][iPart+ip] + cell_length[0]+params.timestep*particles->Momentum[0][iPart+ip]/lf;
                        //            if (  X > params.grid_length[0] ) {
                        //                //cout << "params.grid_length[0]+ cell_length[0]*vel[0] = " << params.grid_length[0]+ cell_length[0]*vel[0] << endl;
                        //                //cout << "params.grid_length[0]                        = " << params.grid_length[0] << endl;
                        //                nPart--; // ne sert à rien ici
                        //                particles->eraseParticle(iPart+ip);
                        //                particles->last_index[(new_cell_idx+i)/clrw]--;
                        //                for ( int idx=(new_cell_idx+i)/clrw+1 ; idx<particles->last_index.size() ; idx++ ) {
                        //                    particles->first_index[idx]--;
                        //                    particles->last_index[idx]--;
                        //                }
                        //            }
                        //            else
                        //                particles->Position[0][iPart+ip] = X;
                        //        }
                        //
                        //    }
                        //}

                        
                        iPart+=nPart;
                    }//END if density > 0
                }//k end the loop on all cells
            }//j
            if((!n_existing_particles)&&( i%species_->clrw == species_->clrw -1 ) &&(initialized_in_species_)) {
                 species_->particles->last_index[(new_cell_idx+i)/species_->clrw] = iPart;
            }

        }//i
    } else if( n_existing_particles == 0 ) {
        // Here position are created from a numpy array.
        // Do not recreate particles from numpy array again after initialization. Is this condition enough ?
        // Initializing particles from numpy array and based on a count sort to comply with initial sorting.
        int nbins = species_->particles->first_index.size();
        int indices[nbins];
        double one_ov_dbin = 1. / ( species_->cell_length[0] * species_->clrw ) ;

        for( int ibin=0; ibin < nbins ; ibin++ ) {
            indices[ibin] = 0 ;
        }

        ///Compute proper indices for particle susing a count sort
        for( unsigned int ipart = 0; ipart < npart_effective ; ipart++ ) {
            unsigned int ip = my_particles_indices[ipart];
            double x = position[0][ip]-species_->min_loc ;
            int ix = int( x * one_ov_dbin ) ;
            indices[ix] ++;
        }
        unsigned int tot=0;
        for( int ibin=0; ibin < nbins; ibin++ ) {
            unsigned int oc = indices[ibin];
            indices[ibin] = tot;
            tot += oc;
        }
        for( int ibin=0; ibin < nbins   ; ibin++ ) {
            species_->particles->first_index[ibin] = indices[ibin] ;
        }
        for( int ibin=0; ibin < nbins-1 ; ibin++ ) {
            species_->particles->last_index[ibin] = species_->particles->first_index[ibin+1] ;
        }
        species_->particles->last_index[nbins-1] = npart_effective ;

        //Now initialize particles at their proper indices
        for( unsigned int ipart = 0; ipart < npart_effective ; ipart++ ) {
            unsigned int ippy = my_particles_indices[ipart];//Indice of the particle in the python array.
            double x = position[0][ippy]-species_->min_loc ;
            unsigned int ibin = int( x * one_ov_dbin ) ;
            int ip = indices[ibin] ; //Indice of the position of the particle in the particles array.

            unsigned int int_ijk[3] = {0, 0, 0};
            if ( params.geometry != "AMcylindrical") {
                for( unsigned int idim=0; idim<species_->nDim_particle; idim++ ) {
                    particles_->position( idim, ip ) = position[idim][ippy];
                    int_ijk[idim] = ( unsigned int )( ( particles_->position( idim, ip ) - species_->min_loc_vec[idim] )/species_->cell_length[idim] );
                }
            }
            else {
                for( unsigned int idim=0; idim<species_->nDim_particle; idim++ ) {
                    particles_->position( idim, ip ) = position[idim][ippy];
                }
                int_ijk[0] = ( unsigned int )( ( particles_->position( 0, ip ) - species_->min_loc_vec[0] )/species_->cell_length[0] );
                int_ijk[1] = ( unsigned int )( ( sqrt( position[1][ippy]*position[1][ippy]+position[2][ippy]*position[2][ippy] )
                                               - species_->min_loc_vec[1] )/species_->cell_length[1] );
            }
            if( !species_->momentum_initialization_array_ ) {
                vel [0] = velocity   [0]( int_ijk[0], int_ijk[1], int_ijk[2] );
                vel [1] = velocity   [1]( int_ijk[0], int_ijk[1], int_ijk[2] );
                vel [2] = velocity   [2]( int_ijk[0], int_ijk[1], int_ijk[2] );
                temp[0] = temperature[0]( int_ijk[0], int_ijk[1], int_ijk[2] );
                temp[1] = temperature[1]( int_ijk[0], int_ijk[1], int_ijk[2] );
                temp[2] = temperature[2]( int_ijk[0], int_ijk[1], int_ijk[2] );
                ParticleCreator::createMomentum( momentum_initialization_, particles_, species_, 1, ip, temp, vel );
            } else {
                for( unsigned int idim=0; idim < 3; idim++ ) {
                    particles_->momentum( idim, ip ) = momentum[idim][ippy]/species_->mass_ ;
                }
            }

            particles_->weight( ip ) = weight_arr[ippy] ;
            ParticleCreator::createCharge( particles_, species_, 1, ip, charge( int_ijk[0], int_ijk[1], int_ijk[2] ) );
            indices[ibin]++;
        }
    }

    // Delete map xyz.
    for( unsigned int idim=0 ; idim<species_->nDim_field ; idim++ ) {
        delete xyz[idim];
    }

    delete [] indexes;
    delete [] temp;
    delete [] vel;
    
    if( particles_->tracked ) {
        particles_->resetIds();
    }
    return npart_effective;
    
} // end create

// ---------------------------------------------------------------------------------------------------------------------
//! Creation of the position for all particles (nPart)
// ---------------------------------------------------------------------------------------------------------------------
void ParticleCreator::createPosition( std::string position_initialization,
                                    Particles * particles,
                                    Species * species,
                                    unsigned int nPart,
                                    unsigned int iPart,
                                    double *indexes,
                                    Params &params )
{
    if( position_initialization == "regular" ) {

        if ( species->regular_number_array_.size()!=0){
            if ( species->regular_number_array_.size() != species->nDim_particle){
                ERROR( "The number of particles required per cell per dimension (regular_number) must be of length " << species->nDim_particle << " in this geometry." );
            }
            unsigned int npart_check=1;
            for( unsigned int idim=0; idim<species->nDim_particle; idim++ ) {
                npart_check *= species->regular_number_array_[idim];
            }
            if( nPart != npart_check) {
                ERROR( "The number of particles required per cell and per dimension is not coherent with the total number of particles per cell." );
            }
        }

        if( params.geometry != "AMcylindrical" ) {
            double inv_coeff_array[3];
            int    coeff_array[3];

            if ( species->regular_number_array_.size()==0){
                double coeff = pow( ( double )nPart, species->inv_nDim_particles );
                if( nPart != ( unsigned int ) pow( round( coeff ), ( double ) species->nDim_particle ) ) {
                    ERROR( "Impossible to put "<<nPart<<" particles regularly spaced in one cell. Use a square number, or `position_initialization = 'random'`" );
                }
                for( unsigned int idim=0; idim<species->nDim_particle; idim++ ) {
                    coeff_array[idim] = coeff;
                    inv_coeff_array[idim] = 1./coeff;
                }
            } else{
               for( unsigned int idim=0; idim<species->nDim_particle; idim++ ) {
                   coeff_array[idim] = species->regular_number_array_[idim];
                   inv_coeff_array[idim] = 1./(double)coeff_array[idim];
               }
            }

            for( unsigned int  p=iPart; p<iPart+nPart; p++ ) {
                int i = ( int )( p-iPart );
                for( unsigned int idim=0; idim<species->nDim_particle; idim++ ) {
                    particles->position( idim, p ) = indexes[idim] + species->cell_length[idim] * 0.975 * inv_coeff_array[idim] * ( 0.5 + i%coeff_array[idim] );
                    i /= coeff_array[idim]; // integer division
                }
            }
        } else { // AM geometry

            unsigned int Np_array[species->nDim_particle];
            double dx, dr, dtheta, theta_offset;

            if ( species->regular_number_array_.size()==0){
                //Trick to derive number of particles per dimension from total number of particles per cell
                int Np = nPart;
                int counter = 0;
                unsigned int prime = 2;
                for( unsigned int idim=0; idim<species->nDim_particle; idim++ ) {
                    Np_array[idim] = 1;
                }

                while( prime <= 23 && Np > 1 ) {
                    if( Np%prime == 0 ) {
                        Np = Np/prime;
                        Np_array[counter%species->nDim_particle] *= prime;
                        counter++;
                    } else {
                        prime++;
                    }
                }
                Np_array[counter%species->nDim_particle] *= Np; //At that point, if Np is not equal to 1, it means that nPart has a prime divisor greater than 23.
                std::sort( Np_array, Np_array + species->nDim_particle ); //sort so that the largest number of particles per dimension is used along theta.

            } else{

                for( unsigned int idim=0; idim<species->nDim_particle; idim++ ) {
                    Np_array[idim] = species->regular_number_array_[idim];
                }
            }

            dx = species->cell_length[0]/Np_array[0];
            dr = species->cell_length[1]/Np_array[1];
            dtheta = 2.*M_PI   /Np_array[2];

            for( unsigned int ix = 0 ; ix < Np_array[0]; ix++ ) {
                double qx = indexes[0] + dx*( ix+0.5 );
                int nx = ix*( Np_array[2]*Np_array[1] );
                for( unsigned int ir = 0 ; ir < Np_array[1]; ir++ ) {
                    double qr = indexes[1] + dr*( ir+0.5 );
                    int nr = ir*( Np_array[2] );
                    theta_offset = Rand::uniform()*2.*M_PI;
                    for( unsigned int itheta = 0 ; itheta < Np_array[2]; itheta++ ) {
                        int p = nx+nr+itheta+iPart;
                        double theta = theta_offset + itheta*dtheta;
                        particles->position( 0, p ) = qx ;
                        particles->position( 1, p ) = qr*cos( theta );
                        particles->position( 2, p ) = qr*sin( theta );
                    }
                }
            }
        }

    } else if( position_initialization == "random" ) {
        if( params.geometry=="AMcylindrical" ) {
            double particles_r, particles_theta;
            for( unsigned int p= iPart; p<iPart+nPart; p++ ) {
                particles->position( 0, p )=indexes[0]+Rand::uniform()*species->cell_length[0];
                particles_r=sqrt( indexes[1]*indexes[1]+ 2.*Rand::uniform()*( indexes[1]+species->cell_length[1]*0.5 )*species->cell_length[1] );
                particles_theta=Rand::uniform()*2.*M_PI;
                particles->position( 2, p )=particles_r*sin( particles_theta );
                particles->position( 1, p )= particles_r*cos( particles_theta );
            }
        } else {
            for( unsigned int p= iPart; p<iPart+nPart; p++ ) {
                for( unsigned int i=0; i<species->nDim_particle ; i++ ) {
                    particles->position( i, p )=indexes[i]+Rand::uniform()*species->cell_length[i];
                }
            }
        }
    } else if( position_initialization == "centered" ) {

        for( unsigned int p=iPart; p<iPart+nPart; p++ )
            for( unsigned int i=0; i<species->nDim_particle ; i++ ) {
                particles->position( i, p )=indexes[i]+0.5*species->cell_length[i];
            }
    }
}

// ---------------------------------------------------------------------------------------------------------------------
//! Creation of the particle momentum
//! For all (np) particles in a mesh initialize their momentum
//!   - at zero (init_momentum_type = cold)
//!   - using random distribution (init_momentum_type = maxwell-juettner)
// ---------------------------------------------------------------------------------------------------------------------
void ParticleCreator::createMomentum( std::string momentum_initialization,
                                    Particles * particles,
                                    Species * species,
                                    unsigned int nPart,
                                    unsigned int iPart,
                                    double * temp,
                                    double * vel )
{
    // -------------------------------------------------------------------------
    // Particles
    // -------------------------------------------------------------------------
    if( species->mass_ > 0 ) {

        // Cold distribution
        if( momentum_initialization == "cold" ) {

            for( unsigned int p=iPart; p<iPart+nPart; p++ ) {
                particles->momentum( 0, p ) = 0.0;
                particles->momentum( 1, p ) = 0.0;
                particles->momentum( 2, p ) = 0.0;
            }

            // Maxwell-Juttner distribution
        } else if( momentum_initialization == "maxwell-juettner" ) {

            // Sample the energies in the MJ distribution
            std::vector<double> energies = maxwellJuttner( species, nPart, temp[0]/species->mass_ );

            // Sample angles randomly and calculate the momentum
            for( unsigned int p=iPart; p<iPart+nPart; p++ ) {
                double phi   = acos( -Rand::uniform2() );
                double theta = 2.0*M_PI*Rand::uniform();
                double psm = sqrt( pow( 1.0+energies[p-iPart], 2 )-1.0 );

                particles->momentum( 0, p ) = psm*cos( theta )*sin( phi );
                particles->momentum( 1, p ) = psm*sin( theta )*sin( phi );
                particles->momentum( 2, p ) = psm*cos( phi );
            }

            // Trick to have non-isotropic distribution (not good)
            double t1 = sqrt( temp[1]/temp[0] ), t2 = sqrt( temp[2]/temp[0] );
            if( t1!=1. || t2 !=1. ) {
                for( unsigned int p= iPart; p<iPart+nPart; p++ ) {
                    particles->momentum( 1, p ) *= t1;
                    particles->momentum( 2, p ) *= t2;
                }
            }

            // Rectangular distribution
        } else if( momentum_initialization == "rectangular" ) {

            double t0 = sqrt( temp[0]/species->mass_ ), t1 = sqrt( temp[1]/species->mass_ ), t2 = sqrt( temp[2]/species->mass_ );
            for( unsigned int p= iPart; p<iPart+nPart; p++ ) {
                particles->momentum( 0, p ) = Rand::uniform2() * t0;
                particles->momentum( 1, p ) = Rand::uniform2() * t1;
                particles->momentum( 2, p ) = Rand::uniform2() * t2;
            }
        }

        // Adding the mean velocity (using relativistic composition)
        // Also relies on the method proposed in Zenitani, Phys. Plasmas 22, 042116 (2015)
        // to ensure the correct properties of a boosted distribution function
        // -------------------------------------------------------------------------------
        double vx, vy, vz, v2, g, gm1, Lxx, Lyy, Lzz, Lxy, Lxz, Lyz, px, py, pz;
        double gamma, inverse_gamma;
        // mean-velocity
        vx  = -vel[0];
        vy  = -vel[1];
        vz  = -vel[2];
        v2  = vx*vx + vy*vy + vz*vz;
        if( v2>0. ) {
            
            if( v2>=1. ) {
                ERROR("The mean velocity should not be higher than the speed of light");
            }
            
            g   = 1.0/sqrt( 1.0-v2 );
            gm1 = g - 1.0;

            // compute the different component of the Matrix block of the Lorentz transformation
            Lxx = 1.0 + gm1 * vx*vx/v2;
            Lyy = 1.0 + gm1 * vy*vy/v2;
            Lzz = 1.0 + gm1 * vz*vz/v2;
            Lxy = gm1 * vx*vy/v2;
            Lxz = gm1 * vx*vz/v2;
            Lyz = gm1 * vy*vz/v2;

            // Volume transformation method (here is the correction by Zenitani)
            double Volume_Acc;
            double CheckVelocity;

            // Lorentz transformation of the momentum
            for( unsigned int p=iPart; p<iPart+nPart; p++ ) {
                gamma = sqrt( 1.0 + particles->momentum( 0, p )*particles->momentum( 0, p )
                           + particles->momentum( 1, p )*particles->momentum( 1, p )
                           + particles->momentum( 2, p )*particles->momentum( 2, p ) );
                inverse_gamma = 1./gamma;

                CheckVelocity = ( vx*particles->momentum( 0, p )
                              + vy*particles->momentum( 1, p )
                              + vz*particles->momentum( 2, p ) ) * inverse_gamma;
                Volume_Acc = Rand::uniform();
                if( CheckVelocity > Volume_Acc ) {

                    double Phi, Theta, vfl, vflx, vfly, vflz, vpx, vpy, vpz ;
                    Phi = atan2( sqrt( vx*vx +vy*vy ), vz );
                    Theta = atan2( vy, vx );

                    vpx = particles->momentum( 0, p )*inverse_gamma ;
                    vpy = particles->momentum( 1, p )*inverse_gamma ;
                    vpz = particles->momentum( 2, p )*inverse_gamma ;
                    vfl = vpx*cos( Theta )*sin( Phi ) +vpy*sin( Theta )*sin( Phi ) + vpz*cos( Phi ) ;
                    vflx = vfl*cos( Theta )*sin( Phi ) ;
                    vfly = vfl*sin( Theta )*sin( Phi ) ;
                    vflz = vfl*cos( Phi ) ;
                    vpx -= 2.*vflx ;
                    vpy -= 2.*vfly ;
                    vpz -= 2.*vflz ;
                    inverse_gamma = sqrt( 1.0 - vpx*vpx - vpy*vpy - vpz*vpz );
                    gamma = 1./inverse_gamma;
                    particles->momentum( 0, p ) = vpx*gamma ;
                    particles->momentum( 1, p ) = vpy*gamma ;
                    particles->momentum( 2, p ) = vpz*gamma ;

                }//here ends the corrections by Zenitani

                px = -gamma*g*vx + Lxx * particles->momentum( 0, p ) + Lxy * particles->momentum( 1, p ) + Lxz * particles->momentum( 2, p );
                py = -gamma*g*vy + Lxy * particles->momentum( 0, p ) + Lyy * particles->momentum( 1, p ) + Lyz * particles->momentum( 2, p );
                pz = -gamma*g*vz + Lxz * particles->momentum( 0, p ) + Lyz * particles->momentum( 1, p ) + Lzz * particles->momentum( 2, p );

                particles->momentum( 0, p ) = px;
                particles->momentum( 1, p ) = py;
                particles->momentum( 2, p ) = pz;
            }

        }//ENDif vel != 0

    }
    // -------------------------------------------------------------------------
    // Photons
    // -------------------------------------------------------------------------
    else if( species->mass_ == 0 ) {
        // Cold distribution
        if( momentum_initialization == "cold" ) {

            //double gamma =sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
            for( unsigned int p=iPart; p<iPart+nPart; p++ ) {
                particles->momentum( 0, p ) = vel[0];
                particles->momentum( 1, p ) = vel[1];
                particles->momentum( 2, p ) = vel[2];
            }

            // Rectangular distribution
        } else if( momentum_initialization == "rectangular" ) {

            //double gamma =sqrt(temp[0]*temp[0] + temp[1]*temp[1] + temp[2]*temp[2]);
            for( unsigned int p= iPart; p<iPart+nPart; p++ ) {
                particles->momentum( 0, p ) = Rand::uniform2()*temp[0];
                particles->momentum( 1, p ) = Rand::uniform2()*temp[1];
                particles->momentum( 2, p ) = Rand::uniform2()*temp[2];
            }

        }
    }
}


// ---------------------------------------------------------------------------------------------------------------------
//! For all (nPart) particles in a mesh initialize its numerical weight (equivalent to a number density)
// ---------------------------------------------------------------------------------------------------------------------
void ParticleCreator::createWeight( std::string position_initialization,
                                    Particles * particles, unsigned int nPart, unsigned int iPart, double n_real_particles,
                                    Params &params )
{
    double w = n_real_particles / nPart;
    for( unsigned int p= iPart; p<iPart+nPart; p++ ) {
        particles->weight( p ) = w ;
    }
}


// ---------------------------------------------------------------------------------------------------------------------
//! For all (nPart) particles in a mesh initialize its numerical weight (equivalent to a number density)
// ---------------------------------------------------------------------------------------------------------------------
void ParticleCreator::regulateWeightwithPositionAM( Particles * particles, std::string position_initialization_on_species_type_, double dr )
{
    unsigned int nParts = particles->Weight.size();

    if ( position_initialization_on_species_type_ == "regular" ){
        //Particles in regular have a weight proportional to their position along r.
        for (unsigned int ipart=0; ipart < nParts ; ipart++){
            double radius = sqrt(particles->position(1,ipart)*particles->position(1,ipart) + particles->position(2,ipart)*particles->position(2,ipart));
            particles->weight(ipart) *= radius;
        }
    } else {
        //Particles in AM have a weight proportional to their intial cell radius
        double dr_inv = 1./dr;
        for (unsigned int ipart=0; ipart < nParts ; ipart++){
            double cell_radius = dr * (floor ( sqrt(particles->position(1,ipart)*particles->position(1,ipart) + particles->position(2,ipart)*particles->position(2,ipart)) * dr_inv) + 0.5);
            particles->weight(ipart) *= cell_radius;
        }
    }
}


// ---------------------------------------------------------------------------------------------------------------------
// For all (np) particles in a mesh initialize its charge state
// ---------------------------------------------------------------------------------------------------------------------
void ParticleCreator::createCharge( Particles * particles, Species * species,
                                    unsigned int nPart, unsigned int iPart, double q )
{
    short Z = ( short )q;
    double r = q-( double )Z;

    // if charge is integer, then all particles have the same charge
    if( r == 0. ) {
        for( unsigned int p = iPart; p<iPart+nPart; p++ ) {
            particles->charge( p ) = Z;
        }
        // if charge is not integer, then particles can have two different charges
    } else {
        int tot = 0, Nm, Np;
        double rr=r/( 1.-r ), diff;
        Np = ( int )round( r*( double )nPart );
        Nm = ( int )nPart - Np;
        for( unsigned int p = iPart; p<iPart+nPart; p++ ) {
            if( Np > rr*Nm ) {
                particles->charge( p ) = Z+1;
                Np--;
            } else {
                particles->charge( p ) = Z;
                Nm--;
            }
            tot += particles->charge( p );
        }
        diff = q - ( ( double )tot )/( ( double )nPart ); // missing charge
        if( diff != 0. ) {
            WARNING( "Could not match exactly charge="<<q<<" for species "<< species->name_
                  << " (difference of "<<diff<<"). Try to add particles." );
        }
    }
}

// ---------------------------------------------------------------------------------------------------------------------
//! Provides a Maxwell-Juttner distribution of energies
// ---------------------------------------------------------------------------------------------------------------------
std::vector<double> ParticleCreator::maxwellJuttner( Species * species, unsigned int npoints, double temperature )
{
    if( temperature==0. ) {
        ERROR( "The species " << species->species_number_ << " is initializing its momentum with the following temperature : " << temperature );
    }
    std::vector<double> energies( npoints );

    // Classical case: Maxwell-Bolztmann
    if( temperature < 0.1 ) {
        double U, lnlnU, invF, I, remainder;
        const double invdU_F = 999./( 2.+19. );
        unsigned int index;
        // For each particle
        for( unsigned int i=0; i<npoints; i++ ) {
            // Pick a random number
            U = Rand::uniform();
            // Calculate the inverse of F
            lnlnU = log( -log( U ) );
            if( lnlnU>2. ) {
                invF = 3.*sqrt( M_PI )/4. * pow( U, 2./3. );
            } else if( lnlnU<-19. ) {
                invF = 1.;
            } else {
                I = ( lnlnU + 19. )*invdU_F;
                index = ( unsigned int )I;
                remainder = I - ( double )index;
                invF = exp( lnInvF[index] + remainder*( lnInvF[index+1]-lnInvF[index] ) );
            }
            // Store that value of the energy
            energies[i] = temperature * invF;
        }

        // Relativistic case: Maxwell-Juttner
    } else {
        double U, lnU, invH, I, remainder, gamma;
        const double invdU_H = 999./( 12.+30. );
        unsigned int index;
        // Calculate the constant H(1/T)
        double invT = 1./temperature;
        double H0 = -invT + log( 1. + invT + 0.5*invT*invT );
        // For each particle
        for( unsigned int i=0; i<npoints; i++ ) {
            do {
                // Pick a random number
                U = Rand::uniform();
                // Calculate the inverse of H at the point log(1.-U) + H0
                lnU = log( -log( 1.-U ) - H0 );
                if( lnU<-26. ) {
                    invH = pow( -6.*U, 1./3. );
                } else if( lnU>12. ) {
                    invH = -U + 11.35 * pow( -U, 0.06 );
                } else {
                    I = ( lnU + 30. )*invdU_H;
                    index = ( unsigned int )I;
                    remainder = I - ( double )index;
                    invH = exp( lnInvH[index] + remainder*( lnInvH[index+1]-lnInvH[index] ) );
                }
                // Make a first guess for the value of gamma
                gamma = temperature * invH;
                // We use the rejection method, so we pick another random number
                U = Rand::uniform();
                // And we are done only if U < beta, otherwise we try again
            } while( U >= sqrt( 1.-1./( gamma*gamma ) ) );
            // Store that value of the energy
            energies[i] = gamma - 1.;
        }
    }

    return energies;
}

// Array used in the Maxwell-Juttner sampling
const double ParticleCreator::lnInvF[1000] = {
    3.028113261189336214e+00, 3.027071073732334305e+00, 3.026027773526156039e+00, 3.024983358705674252e+00, 3.023937826984889554e+00, 3.022891176588271556e+00, 3.021843405329281751e+00,
    3.020794508115738797e+00, 3.019744485967266634e+00, 3.018693334699208197e+00, 3.017641052204572016e+00, 3.016587635695088032e+00, 3.015533082181309776e+00, 3.014477390819529479e+00,
    3.013420557495547936e+00, 3.012362581105449966e+00, 3.011303458598878713e+00, 3.010243186517333580e+00, 3.009181763646249674e+00, 3.008119187197372923e+00, 3.007055454521499804e+00,
    3.005990562743902750e+00, 3.004924509173909630e+00, 3.003857292634036114e+00, 3.002788908474254725e+00, 3.001719355680854129e+00, 3.000648630268972106e+00, 2.999576731991148826e+00,
    2.998503655125519529e+00, 2.997429399666421634e+00, 2.996353962705502472e+00, 2.995277339443297215e+00, 2.994199530000347664e+00, 2.993120530477835217e+00, 2.992040337402567474e+00,
    2.990958949033188929e+00, 2.989876362981995772e+00, 2.988792575310888822e+00, 2.987707585049156567e+00, 2.986621388358625229e+00, 2.985533982440079281e+00, 2.984445365362788039e+00,
    2.983355533727858777e+00, 2.982264485037745771e+00, 2.981172216817944864e+00, 2.980078725335808532e+00, 2.978984008843210685e+00, 2.977888064324651030e+00, 2.976790888518913381e+00,
    2.975692479479559172e+00, 2.974592833341854536e+00, 2.973491948129783236e+00, 2.972389820849905107e+00, 2.971286448441586625e+00, 2.970181828061679408e+00, 2.969075957208892724e+00,
    2.967968832275599933e+00, 2.966860451481935002e+00, 2.965750810995817499e+00, 2.964639908489009379e+00, 2.963527740903820096e+00, 2.962414304820172539e+00, 2.961299597945019624e+00,
    2.960183617061270755e+00, 2.959066359649807243e+00, 2.957947822143851102e+00, 2.956828002037547254e+00, 2.955706895985176885e+00, 2.954584501187757173e+00, 2.953460814660211931e+00,
    2.952335833265368414e+00, 2.951209554214652364e+00, 2.950081974153428543e+00, 2.948953090668825716e+00, 2.947822899669829244e+00, 2.946691399036428738e+00, 2.945558585590287937e+00,
    2.944424455781003314e+00, 2.943289006595699586e+00, 2.942152235293856943e+00, 2.941014138411835344e+00, 2.939874712888751240e+00, 2.938733955997462566e+00, 2.937591863890403943e+00,
    2.936448433824479842e+00, 2.935303662619723308e+00, 2.934157546896903224e+00, 2.933010083853386796e+00, 2.931861270010330145e+00, 2.930711101906631644e+00, 2.929559576625540007e+00,
    2.928406691054723510e+00, 2.927252441657371751e+00, 2.926096825262792578e+00, 2.924939838353341592e+00, 2.923781478255301547e+00, 2.922621741115138061e+00, 2.921460623900341336e+00,
    2.920298122948058239e+00, 2.919134235281933165e+00, 2.917968957511259731e+00, 2.916802286044811510e+00, 2.915634217678482187e+00, 2.914464748961894003e+00, 2.913293876608659794e+00,
    2.912121597106724913e+00, 2.910947907143853985e+00, 2.909772803068770841e+00, 2.908596281477778600e+00, 2.907418339173144961e+00, 2.906238972494791906e+00, 2.905058177761874472e+00,
    2.903875951883352791e+00, 2.902692291057715757e+00, 2.901507191842126243e+00, 2.900320650726790461e+00, 2.899132664073211352e+00, 2.897943228293437645e+00, 2.896752339956896627e+00,
    2.895559995385308838e+00, 2.894366190812813322e+00, 2.893170923001124883e+00, 2.891974188049351024e+00, 2.890775982289901069e+00, 2.889576302190630663e+00, 2.888375144068416223e+00,
    2.887172504172867971e+00, 2.885968378755603858e+00, 2.884762764158861792e+00, 2.883555656589881444e+00, 2.882347052565732426e+00, 2.881136948008295562e+00, 2.879925339228762926e+00,
    2.878712222535656728e+00, 2.877497594136535053e+00, 2.876281450075311774e+00, 2.875063786676820943e+00, 2.873844599923041532e+00, 2.872623886146674632e+00, 2.871401641221379641e+00,
    2.870177861549739085e+00, 2.868952543061442650e+00, 2.867725681847173025e+00, 2.866497273976656768e+00, 2.865267315541965676e+00, 2.864035802530226604e+00, 2.862802730991182099e+00,
    2.861568096884958390e+00, 2.860331896301832710e+00, 2.859094125094009620e+00, 2.857854779315351035e+00, 2.856613854776802519e+00, 2.855371347595562881e+00, 2.854127253560217792e+00,
    2.852881568520355682e+00, 2.851634288537852058e+00, 2.850385409199351461e+00, 2.849134926562906678e+00, 2.847882836492185543e+00, 2.846629134633422264e+00, 2.845373816843268511e+00,
    2.844116878950580407e+00, 2.842858316662027374e+00, 2.841598125822740961e+00, 2.840336302002615554e+00, 2.839072841028552396e+00, 2.837807738507422961e+00, 2.836540990284869057e+00,
    2.835272591834488765e+00, 2.834002538926763126e+00, 2.832730827126246798e+00, 2.831457452075073711e+00, 2.830182409297550272e+00, 2.828905694487738653e+00, 2.827627303066788222e+00,
    2.826347230686399925e+00, 2.825065472803526045e+00, 2.823782024967665283e+00, 2.822496882576258859e+00, 2.821210041215960196e+00, 2.819921496214795376e+00, 2.818631243084864124e+00,
    2.817339277197156377e+00, 2.816045593982217543e+00, 2.814750188733192271e+00, 2.813453056898449489e+00, 2.812154193791040147e+00, 2.810853594670671640e+00, 2.809551254881757387e+00,
    2.808247169700560875e+00, 2.806941334410369748e+00, 2.805633744179458322e+00, 2.804324394269614995e+00, 2.803013279897377252e+00, 2.801700396185725861e+00, 2.800385738398021740e+00,
    2.799069301561947221e+00, 2.797751080840848559e+00, 2.796431071354171127e+00, 2.795109268148242343e+00, 2.793785666321451533e+00, 2.792460260881505008e+00, 2.791133046852197985e+00,
    2.789804019239056299e+00, 2.788473173006982719e+00, 2.787140503101419142e+00, 2.785806004486786946e+00, 2.784469672054690204e+00, 2.783131500686589543e+00, 2.781791485274847542e+00,
    2.780449620648572484e+00, 2.779105901651405475e+00, 2.777760323073211968e+00, 2.776412879715631110e+00, 2.775063566300681739e+00, 2.773712377596409873e+00, 2.772359308309426229e+00,
    2.771004353150072763e+00, 2.769647506758373456e+00, 2.768288763781067807e+00, 2.766928118871522457e+00, 2.765565566608112924e+00, 2.764201101585431974e+00, 2.762834718338050610e+00,
    2.761466411416351630e+00, 2.760096175316944844e+00, 2.758724004532232765e+00, 2.757349893504113858e+00, 2.755973836691637757e+00, 2.754595828493303422e+00, 2.753215863309641964e+00,
    2.751833935474973902e+00, 2.750450039363840027e+00, 2.749064169258609969e+00, 2.747676319470597317e+00, 2.746286484264893080e+00, 2.744894657857348097e+00, 2.743500834488827422e+00,
    2.742105008327913929e+00, 2.740707173548328157e+00, 2.739307324283253742e+00, 2.737905454637300284e+00, 2.736501558705558335e+00, 2.735095630545556489e+00, 2.733687664185743493e+00,
    2.732277653640771131e+00, 2.730865592871666081e+00, 2.729451475853627240e+00, 2.728035296494928374e+00, 2.726617048691977185e+00, 2.725196726325556096e+00, 2.723774323245006901e+00,
    2.722349833245849116e+00, 2.720923250132553761e+00, 2.719494567664118456e+00, 2.718063779572909677e+00, 2.716630879545971489e+00, 2.715195861279537048e+00, 2.713758718405862691e+00,
    2.712319444555017167e+00, 2.710878033300977208e+00, 2.709434478216621756e+00, 2.707988772820614454e+00, 2.706540910616433759e+00, 2.705090885087406694e+00, 2.703638689660375238e+00,
    2.702184317739106945e+00, 2.700727762725688930e+00, 2.699269017956015926e+00, 2.697808076748571260e+00, 2.696344932399573402e+00, 2.694879578164103062e+00, 2.693412007259300633e+00,
    2.691942212892081354e+00, 2.690470188220285053e+00, 2.688995926374097234e+00, 2.687519420457483044e+00, 2.686040663537480278e+00, 2.684559648638043861e+00, 2.683076368776276865e+00,
    2.681590816908736130e+00, 2.680102985979682106e+00, 2.678612868886856013e+00, 2.677120458510173773e+00, 2.675625747673846089e+00, 2.674128729184321429e+00, 2.672629395805679042e+00,
    2.671127740279562790e+00, 2.669623755300685186e+00, 2.668117433529207272e+00, 2.666608767603057206e+00, 2.665097750114018726e+00, 2.663584373613864020e+00, 2.662068630634735200e+00,
    2.660550513662186756e+00, 2.659030015139278724e+00, 2.657507127491555377e+00, 2.655981843095368333e+00, 2.654454154286059353e+00, 2.652924053374173585e+00, 2.651391532625470671e+00,
    2.649856584268291737e+00, 2.648319200494160253e+00, 2.646779373459224871e+00, 2.645237095280490003e+00, 2.643692358034541279e+00, 2.642145153752421649e+00, 2.640595474442429591e+00,
    2.639043312061721824e+00, 2.637488658533055919e+00, 2.635931505730037205e+00, 2.634371845500786513e+00, 2.632809669643025430e+00, 2.631244969918466570e+00, 2.629677738044042368e+00,
    2.628107965700841486e+00, 2.626535644519062629e+00, 2.624960766101144038e+00, 2.623383321996961115e+00, 2.621803303713478694e+00, 2.620220702721562489e+00, 2.618635510446848169e+00,
    2.617047718273066703e+00, 2.615457317536296067e+00, 2.613864299536005298e+00, 2.612268655521309491e+00, 2.610670376693941641e+00, 2.609069454222949336e+00, 2.607465879225293826e+00,
    2.605859642771145790e+00, 2.604250735887037482e+00, 2.602639149556862819e+00, 2.601024874714857660e+00, 2.599407902248199953e+00, 2.597788222998665297e+00, 2.596165827762979106e+00,
    2.594540707285014403e+00, 2.592912852269290358e+00, 2.591282253365011723e+00, 2.589648901175708229e+00, 2.588012786255381670e+00, 2.586373899110266272e+00, 2.584732230196375102e+00,
    2.583087769921877275e+00, 2.581440508639612386e+00, 2.579790436659435304e+00, 2.578137544233525702e+00, 2.576481821567448982e+00, 2.574823258815181592e+00, 2.573161846073629189e+00,
    2.571497573394337266e+00, 2.569830430771919971e+00, 2.568160408149524176e+00, 2.566487495417556275e+00, 2.564811682411042959e+00, 2.563132958911872095e+00, 2.561451314647463118e+00,
    2.559766739291108983e+00, 2.558079222458634838e+00, 2.556388753712818929e+00, 2.554695322557521742e+00, 2.552998918443355247e+00, 2.551299530761480749e+00, 2.549597148847723815e+00,
    2.547891761978005931e+00, 2.546183359372441224e+00, 2.544471930189601050e+00, 2.542757463534301543e+00, 2.541039948445674046e+00, 2.539319373907536814e+00, 2.537595728841977483e+00,
    2.535869002110008985e+00, 2.534139182512810340e+00, 2.532406258786655151e+00, 2.530670219610361205e+00, 2.528931053597593070e+00, 2.527188749298794335e+00, 2.525443295203204652e+00,
    2.523694679733073709e+00, 2.521942891247395568e+00, 2.520187918042916753e+00, 2.518429748346313612e+00, 2.516668370323443149e+00, 2.514903772070240073e+00, 2.513135941616603031e+00,
    2.511364866925872352e+00, 2.509590535894033803e+00, 2.507812936346065502e+00, 2.506032056042026390e+00, 2.504247882668526604e+00, 2.502460403845353287e+00, 2.500669607121238425e+00,
    2.498875479972036739e+00, 2.497078009805216325e+00, 2.495277183953487299e+00, 2.493472989678181762e+00, 2.491665414167052894e+00, 2.489854444534345568e+00, 2.488040067821501555e+00,
    2.486222270992040517e+00, 2.484401040935984906e+00, 2.482576364467900643e+00, 2.480748228323957250e+00, 2.478916619164317492e+00, 2.477081523571096788e+00, 2.475242928047906243e+00,
    2.473400819019639929e+00, 2.471555182831365993e+00, 2.469706005747788424e+00, 2.467853273952659077e+00, 2.465996973548479687e+00, 2.464137090554903597e+00, 2.462273610910279853e+00,
    2.460406520468072511e+00, 2.458535804997896701e+00, 2.456661450184707274e+00, 2.454783441627811591e+00, 2.452901764841395327e+00, 2.451016405251206454e+00, 2.449127348196260101e+00,
    2.447234578927442339e+00, 2.445338082606219654e+00, 2.443437844304995998e+00, 2.441533849005145029e+00, 2.439626081597783713e+00, 2.437714526881108679e+00, 2.435799169561309707e+00,
    2.433879994251672230e+00, 2.431956985469836852e+00, 2.430030127640440352e+00, 2.428099405091363572e+00, 2.426164802054011638e+00, 2.424226302663136323e+00, 2.422283890954779473e+00,
    2.420337550866872078e+00, 2.418387266236796673e+00, 2.416433020802418064e+00, 2.414474798200034300e+00, 2.412512581962929836e+00, 2.410546355522440454e+00, 2.408576102205855829e+00,
    2.406601805235187630e+00, 2.404623447727224139e+00, 2.402641012692436018e+00, 2.400654483033482389e+00, 2.398663841545188191e+00, 2.396669070912619048e+00, 2.394670153711276228e+00,
    2.392667072405136430e+00, 2.390659809346837417e+00, 2.388648346774972619e+00, 2.386632666815606374e+00, 2.384612751478727422e+00, 2.382588582659103782e+00, 2.380560142134247936e+00,
    2.378527411564189453e+00, 2.376490372489827863e+00, 2.374449006332351342e+00, 2.372403294391666417e+00, 2.370353217846057792e+00, 2.368298757750526118e+00, 2.366239895036310603e+00,
    2.364176610508826215e+00, 2.362108884847826662e+00, 2.360036698605156857e+00, 2.357960032204726719e+00, 2.355878865940187250e+00, 2.353793179975101069e+00, 2.351702954340258334e+00,
    2.349608168934095964e+00, 2.347508803520371501e+00, 2.345404837727459668e+00, 2.343296251047237266e+00, 2.341183022833486227e+00, 2.339065132300874428e+00, 2.336942558523775748e+00,
    2.334815280435037277e+00, 2.332683276824467633e+00, 2.330546526337776925e+00, 2.328405007475015331e+00, 2.326258698589899421e+00, 2.324107577887416287e+00, 2.321951623423519795e+00,
    2.319790813103117078e+00, 2.317625124679225657e+00, 2.315454535750915088e+00, 2.313279023762432995e+00, 2.311098566001612564e+00, 2.308913139598376407e+00, 2.306722721523330577e+00,
    2.304527288586462053e+00, 2.302326817435369044e+00, 2.300121284553910517e+00, 2.297910666260743362e+00, 2.295694938707419031e+00, 2.293474077877455386e+00, 2.291248059584263252e+00,
    2.289016859469788390e+00, 2.286780453002789315e+00, 2.284538815477031637e+00, 2.282291922010251994e+00, 2.280039747541927397e+00, 2.277782266831583691e+00, 2.275519454457507251e+00,
    2.273251284814709283e+00, 2.270977732113204084e+00, 2.268698770376391671e+00, 2.266414373438985663e+00, 2.264124514945681899e+00, 2.261829168348788333e+00, 2.259528306906793294e+00,
    2.257221903682209874e+00, 2.254909931539875068e+00, 2.252592363145054843e+00, 2.250269170961193499e+00, 2.247940327248446390e+00, 2.245605804061153510e+00, 2.243265573246389533e+00,
    2.240919606441381440e+00, 2.238567875071938662e+00, 2.236210350349980391e+00, 2.233847003271789866e+00, 2.231477804615432881e+00, 2.229102724939111990e+00, 2.226721734578658296e+00,
    2.224334803645326541e+00, 2.221941902023694571e+00, 2.219542999369362501e+00, 2.217138065106579958e+00, 2.214727068426054046e+00, 2.212309978282390066e+00, 2.209886763392083342e+00,
    2.207457392230697035e+00, 2.205021833030698541e+00, 2.202580053778893543e+00, 2.200132022214010163e+00, 2.197677705824045091e+00, 2.195217071843857504e+00, 2.192750087252414382e+00,
    2.190276718770312492e+00, 2.187796932857083210e+00, 2.185310695708424511e+00, 2.182817973253639909e+00, 2.180318731152714573e+00, 2.177812934793672994e+00, 2.175300549289609364e+00,
    2.172781539476034585e+00, 2.170255869907644186e+00, 2.167723504855801675e+00, 2.165184408305282471e+00, 2.162638543951286962e+00, 2.160085875196527283e+00, 2.157526365148061398e+00,
    2.154959976614116091e+00, 2.152386672101062715e+00, 2.149806413810088301e+00, 2.147219163634004779e+00, 2.144624883153926298e+00, 2.142023533636032706e+00, 2.139415076028050144e+00,
    2.136799470955930147e+00, 2.134176678720404841e+00, 2.131546659293431567e+00, 2.128909372314629067e+00, 2.126264777087764291e+00, 2.123612832576956322e+00, 2.120953497403178734e+00,
    2.118286729840338722e+00, 2.115612487811585396e+00, 2.112930728885486609e+00, 2.110241410272010398e+00, 2.107544488818754225e+00, 2.104839921006786074e+00, 2.102127662946750242e+00,
    2.099407670374596080e+00, 2.096679898647592299e+00, 2.093944302739995766e+00, 2.091200837238794019e+00, 2.088449456339469545e+00, 2.085690113841469628e+00, 2.082922763143882694e+00,
    2.080147357240880179e+00, 2.077363848717170391e+00, 2.074572189743352446e+00, 2.071772332071244005e+00, 2.068964227029113978e+00, 2.066147825516855718e+00, 2.063323078001130018e+00,
    2.060489934510390864e+00, 2.057648344629890325e+00, 2.054798257496523561e+00, 2.051939621793771984e+00, 2.049072385746379066e+00, 2.046196497115118529e+00, 2.043311903191355583e+00,
    2.040418550791640140e+00, 2.037516386252159695e+00, 2.034605355423152684e+00, 2.031685403663211709e+00, 2.028756475833488615e+00, 2.025818516291934213e+00, 2.022871468887311508e+00,
    2.019915276953162753e+00, 2.016949883301822677e+00, 2.013975230218087109e+00, 2.010991259453104973e+00, 2.007997912217875047e+00, 2.004995129176902591e+00, 2.001982850441623718e+00,
    1.998961015563778254e+00, 1.995929563528624051e+00, 1.992888432748232796e+00, 1.989837561054456971e+00, 1.986776885691970529e+00, 1.983706343311098630e+00, 1.980625869960619401e+00,
    1.977535401080458666e+00, 1.974434871494165078e+00, 1.971324215401476510e+00, 1.968203366370548713e+00, 1.965072257330279015e+00, 1.961930820562358013e+00, 1.958778987693311269e+00,
    1.955616689686361820e+00, 1.952443856833217639e+00, 1.949260418745664580e+00, 1.946066304347138454e+00, 1.942861441864086824e+00, 1.939645758817238219e+00, 1.936419182012731222e+00,
    1.933181637533145425e+00, 1.929933050728310340e+00, 1.926673346206103199e+00, 1.923402447822974537e+00, 1.920120278674454450e+00, 1.916826761085404396e+00, 1.913521816600212144e+00,
    1.910205365972802882e+00, 1.906877329156460910e+00, 1.903537625293552082e+00, 1.900186172705079723e+00, 1.896822888880065117e+00, 1.893447690464737931e+00, 1.890060493251651819e+00,
    1.886661212168522672e+00, 1.883249761266976074e+00, 1.879826053711074696e+00, 1.876390001765686044e+00, 1.872941516784688565e+00, 1.869480509198935936e+00, 1.866006888504112782e+00,
    1.862520563248323491e+00, 1.859021441019548915e+00, 1.855509428432867036e+00, 1.851984431117487118e+00, 1.848446353703593781e+00, 1.844895099808936845e+00, 1.841330572025270618e+00,
    1.837752671904536061e+00, 1.834161299944834456e+00, 1.830556355576170802e+00, 1.826937737146010710e+00, 1.823305341904525756e+00, 1.819659065989680524e+00, 1.815998804412052525e+00,
    1.812324451039377227e+00, 1.808635898580905010e+00, 1.804933038571444825e+00, 1.801215761355188238e+00, 1.797483956069260147e+00, 1.793737510626989495e+00, 1.789976311700933076e+00,
    1.786200244705582918e+00, 1.782409193779852963e+00, 1.778603041769189819e+00, 1.774781670207480122e+00, 1.770944959298595300e+00, 1.767092787897660999e+00, 1.763225033492008320e+00,
    1.759341572181807090e+00, 1.755442278660384936e+00, 1.751527026194209302e+00, 1.747595686602526843e+00, 1.743648130236697957e+00, 1.739684225959130748e+00, 1.735703841121919000e+00,
    1.731706841545051834e+00, 1.727693091494336919e+00, 1.723662453658859839e+00, 1.719614789128142363e+00, 1.715549957368860401e+00, 1.711467816201181602e+00, 1.707368221774694383e+00,
    1.703251028543929291e+00, 1.699116089243453587e+00, 1.694963254862551727e+00, 1.690792374619444205e+00, 1.686603295935078739e+00, 1.682395864406460051e+00, 1.678169923779530448e+00,
    1.673925315921539703e+00, 1.669661880792991049e+00, 1.665379456419029180e+00, 1.661077878860387935e+00, 1.656756982183783089e+00, 1.652416598431809991e+00, 1.648056557592287863e+00,
    1.643676687567082295e+00, 1.639276814140357530e+00, 1.634856760946272747e+00, 1.630416349436097923e+00, 1.625955398844746824e+00, 1.621473726156702044e+00, 1.616971146071332965e+00,
    1.612447470967592000e+00, 1.607902510868064239e+00, 1.603336073402373163e+00, 1.598747963769924896e+00, 1.594137984701968103e+00, 1.589505936422963783e+00, 1.584851616611246961e+00,
    1.580174820358985821e+00, 1.575475340131367474e+00, 1.570752965725072725e+00, 1.566007484225963253e+00, 1.561238679965980536e+00, 1.556446334479266724e+00, 1.551630226457446637e+00,
    1.546790131704089966e+00, 1.541925823088317049e+00, 1.537037070497540237e+00, 1.532123640789306629e+00, 1.527185297742251313e+00, 1.522221802006102687e+00, 1.517232911050763189e+00,
    1.512218379114408595e+00, 1.507177957150604541e+00, 1.502111392774418741e+00, 1.497018430207492923e+00, 1.491898810222073379e+00, 1.486752270083952387e+00, 1.481578543494325384e+00,
    1.476377360530513050e+00, 1.471148447585529651e+00, 1.465891527306504516e+00, 1.460606318531867620e+00, 1.455292536227339673e+00, 1.449949891420651227e+00, 1.444578091134982145e+00,
    1.439176838321101659e+00, 1.433745831788158176e+00, 1.428284766133101957e+00, 1.422793331668710026e+00, 1.417271214350181552e+00, 1.411718095700257525e+00, 1.406133652732862727e+00,
    1.400517557875193386e+00, 1.394869478888252967e+00, 1.389189078785785014e+00, 1.383476015751558874e+00, 1.377729943054986528e+00, 1.371950508965018134e+00, 1.366137356662287861e+00,
    1.360290124149454716e+00, 1.354408444159718172e+00, 1.348491944063437087e+00, 1.342540245772842367e+00, 1.336552965644766111e+00, 1.330529714381354589e+00, 1.324470096928725082e+00,
    1.318373712373507534e+00, 1.312240153837215484e+00, 1.306069008368417661e+00, 1.299859856832624949e+00, 1.293612273799870760e+00, 1.287325827429903979e+00, 1.281000079354956078e+00,
    1.274634584560003780e+00, 1.268228891260492430e+00, 1.261782540777428796e+00, 1.255295067409807119e+00, 1.248765998304282920e+00, 1.242194853322038606e+00, 1.235581144902765605e+00,
    1.228924377925692424e+00, 1.222224049567603776e+00, 1.215479649157727993e+00, 1.208690658029483833e+00, 1.201856549368946370e+00, 1.194976788059981532e+00, 1.188050830525968715e+00,
    1.181078124568002430e+00, 1.174058109199509481e+00, 1.166990214477174659e+00, 1.159873861328088340e+00, 1.152708461373012083e+00, 1.145493416745687165e+00, 1.138228119908041958e+00,
    1.130911953461235075e+00, 1.123544289952409958e+00, 1.116124491677038755e+00, 1.108651910476764924e+00, 1.101125887532611536e+00, 1.093545753153434275e+00, 1.085910826559507436e+00,
    1.078220415661100162e+00, 1.070473816831928549e+00, 1.062670314677327976e+00, 1.054809181797020434e+00, 1.046889678542338942e+00, 1.038911052767741428e+00, 1.030872539576482394e+00,
    1.022773361060272501e+00, 1.014612726032773971e+00, 1.006389829756761278e+00, 9.981038536647748316e-01, 9.897539650731033145e-01, 9.813393168888950857e-01, 9.728590473102259883e-01,
    9.643122795189308372e-01, 9.556981213659938579e-01, 9.470156650493019024e-01, 9.382639867835408376e-01, 9.294421464620334916e-01, 9.205491873102766842e-01, 9.115841355309679539e-01,
    9.025459999402748457e-01, 8.934337715950928516e-01, 8.842464234110668508e-01, 8.749829097710923875e-01, 8.656421661240399912e-01, 8.562231085734303138e-01, 8.467246334557674281e-01,
    8.371456169082474030e-01, 8.274849144255297384e-01, 8.177413604052780061e-01, 8.079137676821201985e-01, 7.980009270497354645e-01, 7.880016067706957505e-01, 7.779145520737251740e-01,
    7.677384846380108652e-01, 7.574721020641970171e-01, 7.471140773316662376e-01, 7.366630582417227346e-01, 7.261176668462573369e-01, 7.154764988614864540e-01, 7.047381230662990159e-01,
    6.939010806848042723e-01, 6.829638847525810741e-01, 6.719250194661655629e-01, 6.607829395152855501e-01, 6.495360693973274424e-01, 6.381828027135106884e-01, 6.267215014462346190e-01,
    6.151504952170321339e-01, 6.034680805245669077e-01, 5.916725199620531672e-01, 5.797620414135221667e-01, 5.677348372282775557e-01, 5.555890633728943762e-01, 5.433228385600959998e-01,
    5.309342433538073447e-01, 5.184213192496778255e-01, 5.057820677303297430e-01, 4.930144492945795487e-01, 4.801163824598490004e-01, 4.670857427369426995e-01, 4.539203615763838240e-01,
    4.406180252854359769e-01, 4.271764739149248591e-01, 4.135934001149544903e-01, 3.998664479585762876e-01, 3.859932117324395562e-01, 3.719712346934307967e-01, 3.577980077902795375e-01,
    3.434709683490433374e-01, 3.289874987214279622e-01, 3.143449248947801089e-01, 2.995405150626198609e-01, 2.845714781545091165e-01, 2.694349623240440139e-01, 2.541280533937015362e-01,
    2.386477732552492736e-01, 2.229910782243804401e-01, 2.071548573482217825e-01, 1.911359306642606559e-01, 1.749310474093025702e-01, 1.585368841769239923e-01, 1.419500430219155229e-01,
    1.251670495101367131e-01, 1.081843507121787035e-01, 9.099831313917583486e-02, 7.360522061907438796e-02, 5.600127211162326396e-02, 3.818257946031319849e-02, 2.014516507940277909e-02,
    1.884959574231336163e-03, -1.660220070715516874e-02, -3.532057619297919449e-02, -5.427452667931363661e-02, -7.346851399632661761e-02, -9.290710475123310774e-02, -1.125949731501099393e-01,
    -1.325369039007213978e-01, -1.527377951985385285e-01, -1.732026617982476979e-01, -1.939366381729568933e-01, -2.149449817634470628e-01, -2.362330763198251526e-01, -2.578064353379159734e-01,
    -2.796707055928397967e-01, -3.018316707721606651e-01, -3.242952552111035835e-01, -3.470675277322518570e-01, -3.701547055922924656e-01, -3.935631585382081687e-01, -4.172994129754923343e-01,
    -4.413701562508508536e-01, -4.657822410518934197e-01, -4.905426899263066187e-01, -5.156586999229669788e-01, -5.411376473574273094e-01, -5.669870927042111042e-01, -5.932147856182231616e-01,
    -6.198286700876785016e-01, -6.468368897207085189e-01, -6.742477931678908520e-01, -7.020699396827745353e-01, -7.303121048224235912e-01, -7.589832862898669985e-01, -7.880927099202489350e-01,
    -8.176498358123210908e-01, -8.476643646068054982e-01, -8.781462439129084085e-01, -9.091056748842766266e-01, -9.405531189452329688e-01, -9.724993046681473796e-01, -1.004955234802481279e+00,
    -1.037932193455777918e+00, -1.071441753426678245e+00, -1.105495783689714573e+00, -1.140106457031379916e+00, -1.175286257836624326e+00, -1.211047990024640253e+00, -1.247404785132468108e+00,
    -1.284370110544520438e+00, -1.321957777865867323e+00, -1.360181951436657322e+00, -1.399057156984695904e+00, -1.438598290412803626e+00, -1.478820626717161746e+00, -1.519739829032453526e+00,
    -1.561371957799146948e+00, -1.603733480048015858e+00, -1.646841278796369101e+00, -1.690712662550236223e+00, -1.735365374906370661e+00, -1.780817604247509900e+00, -1.827087993524071807e+00,
    -1.874195650115166645e+00, -1.922160155761540201e+00, -1.971001576562910973e+00, -2.020740473031910955e+00, -2.071397910196988157e+00, -2.122995467746291887e+00, -2.175555250204936009e+00,
    -2.229099897138042685e+00, -2.283652593372224260e+00, -2.339237079228581706e+00, -2.395877660760723504e+00, -2.453599219991908331e+00, -2.512427225146172116e+00, -2.572387740869024775e+00,
    -2.633507438434537828e+00, -2.695813605936466573e+00, -2.759334158462542153e+00, -2.824097648252419823e+00, -2.890133274841317768e+00, -2.957470895193073268e+00, -3.026141033828172677e+00,
    -3.096174892954419722e+00, -3.167604362609613755e+00, -3.240462030828213713e+00, -3.314781193846027207e+00, -3.390595866359104260e+00, -3.467940791855854687e+00, -3.546851453043248448e+00,
    -3.627364082390851063e+00, -3.709515672818399779e+00, -3.793343988555180957e+00, -3.878887576201507592e+00, -3.966185776024640397e+00, -4.055278733523468127e+00, -4.146207411298033385e+00,
    -4.239013601260967157e+00, -4.333739937229747596e+00, -4.430429907939251954e+00, -4.529127870514209953e+00, -4.629879064442645742e+00, -4.732729626089497543e+00
};

// Array used in the Maxwell-Juttner sampling
const double ParticleCreator::lnInvH[1000] = {
    -9.402756483482924921e+00, -9.388318921346392898e+00, -9.374883046130019437e+00, -9.360818128319227327e+00, -9.346352500952018971e+00, -9.332609535599788231e+00, -9.318636378401940590e+00,
        -9.304555775711108367e+00, -9.290525531420621874e+00, -9.276495588326515218e+00, -9.262659647248140615e+00, -9.248642923366197977e+00, -9.234537986282660427e+00, -9.220440348595868585e+00,
        -9.206526581779595375e+00, -9.192429091051817380e+00, -9.178493014106759773e+00, -9.164404479064138798e+00, -9.150395754817472138e+00, -9.136394001864923453e+00, -9.122462891956370612e+00,
        -9.108479530895117193e+00, -9.094472838919676505e+00, -9.080484396623118570e+00, -9.066302396191021629e+00, -9.052305266120191263e+00, -9.038372060472962488e+00, -9.024257614905998537e+00,
        -9.010261730228865673e+00, -8.996299574888430683e+00, -8.982360935717023764e+00, -8.968261999420034769e+00, -8.954257835762128082e+00, -8.940289487588648498e+00, -8.926164230689391133e+00,
        -8.912199784045913731e+00, -8.898161019156789919e+00, -8.884171549684843683e+00, -8.870206113458378283e+00, -8.856192689655291161e+00, -8.842146098452111858e+00, -8.828133587340268207e+00,
        -8.814141950394821734e+00, -8.800128482639740568e+00, -8.786091791355046254e+00, -8.772065927576454314e+00, -8.758060948224377640e+00, -8.744020933568334542e+00, -8.729979616492244077e+00,
        -8.715998025917098602e+00, -8.701976304743247681e+00, -8.688005470915637574e+00, -8.673989982215227101e+00, -8.659985003093263245e+00, -8.645943134466779867e+00, -8.631941877630428195e+00,
        -8.617932722102422005e+00, -8.603896655344126287e+00, -8.589888145629942073e+00, -8.575831329584334028e+00, -8.561871862507890896e+00, -8.547821679997451128e+00, -8.533818689215104669e+00,
        -8.519833016046151286e+00, -8.505800768836330406e+00, -8.491788588790452508e+00, -8.477763718789315561e+00, -8.463751871535670546e+00, -8.449750573726415581e+00, -8.435726296060591878e+00,
        -8.421710291650438052e+00, -8.407696461704228241e+00, -8.393684105181973720e+00, -8.379666566421843044e+00, -8.365649358002741565e+00, -8.351631625725728370e+00, -8.337625365285765255e+00,
        -8.323606711928743351e+00, -8.309588119961361485e+00, -8.295570530937240505e+00, -8.281566039082182584e+00, -8.267550775100456661e+00, -8.253535360283727584e+00, -8.239518589650394631e+00,
        -8.225500054829783636e+00, -8.211478531743365394e+00, -8.197475804036692182e+00, -8.183450797927536158e+00, -8.169443383719595886e+00, -8.155428558470900313e+00, -8.141413336443854121e+00,
        -8.127399047429936019e+00, -8.113385870188610127e+00, -8.099365847720060074e+00, -8.085356566734459349e+00, -8.071338108809586132e+00, -8.057323083786892326e+00, -8.043308533608312771e+00,
        -8.029283462705111063e+00, -8.015277322954988293e+00, -8.001260912686870341e+00, -7.987248308454050871e+00, -7.973235083485832320e+00, -7.959214299713621266e+00, -7.945203032260913290e+00,
        -7.931188194720145468e+00, -7.917173759481055839e+00, -7.903154847403775385e+00, -7.889135491874881723e+00, -7.875121967064339756e+00, -7.861110899903347438e+00, -7.847093968442999667e+00,
        -7.833077727373989774e+00, -7.819062115094970622e+00, -7.805048389543351561e+00, -7.791032607063780091e+00, -7.777020408512951732e+00, -7.762999809997588763e+00, -7.748986096367973531e+00,
        -7.734970856036877507e+00, -7.720951970669967857e+00, -7.706938603231540341e+00, -7.692924736361460347e+00, -7.678906461019358254e+00, -7.664893161079090689e+00, -7.650876041121127180e+00,
        -7.636861215452872997e+00, -7.622845466222778477e+00, -7.608828024984916283e+00, -7.594812640914675228e+00, -7.580797680436783814e+00, -7.566782149947995251e+00, -7.552765744609319931e+00,
        -7.538750580869197471e+00, -7.524734998913311657e+00, -7.510718694153761810e+00, -7.496702604314066321e+00, -7.482686172118370393e+00, -7.468669842403806491e+00, -7.454654512935313448e+00,
        -7.440638616231906255e+00, -7.426621810922931388e+00, -7.412607158172065169e+00, -7.398589434939852438e+00, -7.384574035610326881e+00, -7.370557463576223434e+00, -7.356540657627733459e+00,
        -7.342525153318804065e+00, -7.328508954197332059e+00, -7.314492188771704484e+00, -7.300475443269295539e+00, -7.286460013886568277e+00, -7.272443057492030682e+00, -7.258426815444318336e+00,
        -7.244409972475176929e+00, -7.230393545186772464e+00, -7.216377394851891225e+00, -7.202360689208528122e+00, -7.188344086261769128e+00, -7.174327161166150546e+00, -7.160310137083985893e+00,
        -7.146293051212126102e+00, -7.132276902050797673e+00, -7.118259888890483111e+00, -7.104242846535062661e+00, -7.090226085231720710e+00, -7.076209195279213660e+00, -7.062192351651393807e+00,
        -7.048175543678467214e+00, -7.034158119357607930e+00, -7.020141117404175013e+00, -7.006123850412526721e+00, -6.992106579636755193e+00, -6.978089207294207341e+00, -6.964072209306565675e+00,
        -6.950054581777039608e+00, -6.936037293725791031e+00, -6.922019986767113942e+00, -6.908002377494772084e+00, -6.893984966736285358e+00, -6.879967314518558474e+00, -6.865949626461007149e+00,
        -6.851932031486784425e+00, -6.837914145871644145e+00, -6.823896411506032322e+00, -6.809878490100980564e+00, -6.795860633787632388e+00, -6.781842802936246528e+00, -6.767824770899844466e+00,
        -6.753806395345270275e+00, -6.739788309358444529e+00, -6.725770146734140198e+00, -6.711752028083568966e+00, -6.697733741195642132e+00, -6.683715263982263011e+00, -6.669697011456180213e+00,
        -6.655678448295303973e+00, -6.641659934976074986e+00, -6.627641256235382805e+00, -6.613622641520617407e+00, -6.599603827649048959e+00, -6.585584973605825176e+00, -6.571566171568333559e+00,
        -6.557547183105451261e+00, -6.543528196118535867e+00, -6.529509089637392627e+00, -6.515489894663369697e+00, -6.501470707639946056e+00, -6.487451315675093255e+00, -6.473431951649657456e+00,
        -6.459412536766182988e+00, -6.445392946609437068e+00, -6.431373347339028435e+00, -6.417353628364098839e+00, -6.403333843548605131e+00, -6.389314006174331873e+00, -6.375294070149028158e+00,
        -6.361274032230475051e+00, -6.347253896227756265e+00, -6.333233728056003820e+00, -6.319213446735220785e+00, -6.305193076729772805e+00, -6.291172602855333196e+00, -6.277152056531630109e+00,
        -6.263131434948522092e+00, -6.249110669730993273e+00, -6.235089834011095178e+00, -6.221068901366077597e+00, -6.207047871955993834e+00, -6.193026738061672809e+00, -6.179005530463898666e+00,
        -6.164984203796887385e+00, -6.150962772820938618e+00, -6.136941223435446346e+00, -6.122919590081422392e+00, -6.108897823041722575e+00, -6.094875982483746846e+00, -6.080853991215137810e+00,
        -6.066831914689460703e+00, -6.052809719325111359e+00, -6.038787398623465918e+00, -6.024764964099071030e+00, -6.010742424656190686e+00, -5.996719754664948887e+00, -5.982696958028578038e+00,
        -5.968674044105862997e+00, -5.954650997286028868e+00, -5.940627840583593944e+00, -5.926604542202445813e+00, -5.912581112044525078e+00, -5.898557539480430378e+00, -5.884533835062183194e+00,
        -5.870510003642364083e+00, -5.856486031202377873e+00, -5.842461916555514279e+00, -5.828437659701839024e+00, -5.814413251724160681e+00, -5.800388708424720541e+00, -5.786364015360013546e+00,
        -5.772339155289128776e+00, -5.758314155110069166e+00, -5.744288994507930290e+00, -5.730263684034876626e+00, -5.716238206091529328e+00, -5.702212571495541837e+00, -5.688186769701041534e+00,
        -5.674160802365624257e+00, -5.660134668080841536e+00, -5.646108355685827362e+00, -5.632081881824617220e+00, -5.618055220226613855e+00, -5.604028383548823378e+00, -5.590001365458268090e+00,
        -5.575974167980627172e+00, -5.561946776948706095e+00, -5.547919199106292609e+00, -5.533891426332606223e+00, -5.519863466651733219e+00, -5.505835303004738890e+00, -5.491806946614043561e+00,
        -5.477778380315056594e+00, -5.463749609942059493e+00, -5.449720632566440237e+00, -5.435691440821471154e+00, -5.421662039716173886e+00, -5.407632416902895756e+00, -5.393602573387089372e+00,
        -5.379572509719777074e+00, -5.365542216350482008e+00, -5.351511692335392034e+00, -5.337480935451307751e+00, -5.323449941825153076e+00, -5.309418707795056314e+00, -5.295387229476762769e+00,
        -5.281355504533746803e+00, -5.267323530078401816e+00, -5.253291301650716782e+00, -5.239258814888134275e+00, -5.225226066090911559e+00, -5.211193053448391233e+00, -5.197159771824590990e+00,
        -5.183126218151476117e+00, -5.169092385672625412e+00, -5.155058275382643274e+00, -5.141023879201224389e+00, -5.126989194548196238e+00, -5.112954218078122892e+00, -5.098918945828792459e+00,
        -5.084883371575669386e+00, -5.070847492041045790e+00, -5.056811303731988616e+00, -5.042774801179291444e+00, -5.028737980487945514e+00, -5.014700836583935839e+00, -5.000663365711495167e+00,
        -4.986625562378193877e+00, -4.972587422641337795e+00, -4.958548941040521463e+00, -4.944510113413161712e+00, -4.930470934645069470e+00, -4.916431399638748090e+00, -4.902391502425701653e+00,
        -4.888351239474742371e+00, -4.874310604571175709e+00, -4.860269592484651291e+00, -4.846228198304002532e+00, -4.832186416119078842e+00, -4.818144241376513648e+00, -4.804101666899612155e+00,
        -4.790058688038118184e+00, -4.776015299456037866e+00, -4.761971494014358264e+00, -4.747927267417144215e+00, -4.733882612267035661e+00, -4.719837522939142715e+00, -4.705791993126767991e+00,
        -4.691746016933713292e+00, -4.677699588149253351e+00, -4.663652699310194549e+00, -4.649605345185644723e+00, -4.635557518140910105e+00, -4.621509212121083898e+00, -4.607460419599018309e+00,
        -4.593411134326189860e+00, -4.579361349129068337e+00, -4.565311056852286775e+00, -4.551260249935657143e+00, -4.537208921127983352e+00, -4.523157063716806370e+00, -4.509104669361667206e+00,
        -4.495051730786109978e+00, -4.480998239915424008e+00, -4.466944189213031713e+00, -4.452889570797960772e+00, -4.438834376060603049e+00, -4.424778597381012979e+00, -4.410722225868971336e+00,
        -4.396665253509551619e+00, -4.382607671517379622e+00, -4.368549471322706879e+00, -4.354490643789567628e+00, -4.340431180308653047e+00, -4.326371071702013182e+00, -4.312310308837443706e+00,
        -4.298248882124207526e+00, -4.284186782349916456e+00, -4.270123999670713211e+00, -4.256060524440004933e+00, -4.241996346796963913e+00, -4.227931456807191068e+00, -4.213865843972508074e+00,
        -4.199799498259279673e+00, -4.185732409026496903e+00, -4.171664565769908073e+00, -4.157595957540261011e+00, -4.143526573541322477e+00, -4.129456402581863195e+00, -4.115385433401063331e+00,
        -4.101313654585375446e+00, -4.087241054608058199e+00, -4.073167621647699299e+00, -4.059093343771843720e+00, -4.045018208971966622e+00, -4.030942204774508930e+00, -4.016865318883647618e+00,
        -4.002787538662525790e+00, -3.988708851168091929e+00, -3.974629243506720311e+00, -3.960548702433660750e+00, -3.946467214545992253e+00, -3.932384766302702506e+00, -3.918301343838640260e+00,
        -3.904216933265594491e+00, -3.890131520343185656e+00, -3.876045090659561154e+00, -3.861957629639284129e+00, -3.847869122433445366e+00, -3.833779554071717666e+00, -3.819688909227958007e+00,
        -3.805597172458903277e+00, -3.791504328004591251e+00, -3.777410359985496235e+00, -3.763315252231650643e+00, -3.749218988318105961e+00, -3.735121551585395494e+00, -3.721022925173522733e+00,
        -3.706923091918343349e+00, -3.692822034423430377e+00, -3.678719735067865759e+00, -3.664616175927026820e+00, -3.650511338833094754e+00, -3.636405205371272764e+00, -3.622297756801203583e+00,
        -3.608188974152441997e+00, -3.594078838159660005e+00, -3.579967329266379839e+00, -3.565854427652642045e+00, -3.551740113183203373e+00, -3.537624365419402839e+00, -3.523507163639728823e+00,
        -3.509388486827117681e+00, -3.495268313597620136e+00, -3.481146622322163875e+00, -3.467023390990199516e+00, -3.452898597305658956e+00, -3.438772218629355226e+00, -3.424644231970350283e+00,
        -3.410514614028766545e+00, -3.396383341118600185e+00, -3.382250389237801524e+00, -3.368115734004724082e+00, -3.353979350675957427e+00, -3.339841214159893212e+00, -3.325701298958231078e+00,
        -3.311559579221938776e+00, -3.297416028694333701e+00, -3.283270620742043899e+00, -3.269123328322449229e+00, -3.254974123986548840e+00, -3.240822979894216527e+00, -3.226669867761286170e+00,
        -3.212514758894621458e+00, -3.198357624168206126e+00, -3.184198434015861778e+00, -3.170037158433746782e+00, -3.155873766953881443e+00, -3.141708228660941238e+00, -3.127540512170100762e+00,
        -3.113370585620815145e+00, -3.099198416671123280e+00, -3.085023972502690448e+00, -3.070847219786304816e+00, -3.056668124686660715e+00, -3.042486652874943243e+00, -3.028302769482326795e+00,
        -3.014116439119186808e+00, -2.999927625862666591e+00, -2.985736293235094152e+00, -2.971542404209721155e+00, -2.957345921201353089e+00, -2.943146806038983154e+00, -2.928945019976269659e+00,
        -2.914740523681580875e+00, -2.900533277212366645e+00, -2.886323240022723002e+00, -2.872110370940061053e+00, -2.857894628167629847e+00, -2.843675969261741709e+00, -2.829454351133278145e+00,
        -2.815229730025366806e+00, -2.801002061511853380e+00, -2.786771300482580038e+00, -2.772537401131430990e+00, -2.758300316947211339e+00, -2.744060000697740698e+00, -2.729816404423872633e+00,
        -2.715569479423635801e+00, -2.701319176241026909e+00, -2.687065444653787605e+00, -2.672808233660158184e+00, -2.658547491464986923e+00, -2.644283165469476504e+00, -2.630015202255621354e+00,
        -2.615743547571750316e+00, -2.601468146322303721e+00, -2.587188942549503601e+00, -2.572905879421945929e+00, -2.558618899217506293e+00, -2.544327943310193074e+00, -2.530032952155425452e+00,
        -2.515733865272170355e+00, -2.501430621229792539e+00, -2.487123157629747006e+00, -2.472811411091746070e+00, -2.458495317233948541e+00, -2.444174810658960695e+00, -2.429849824935295999e+00,
        -2.415520292579488970e+00, -2.401186145038158415e+00, -2.386847312671634747e+00, -2.372503724731980679e+00, -2.358155309347800266e+00, -2.343801993501943404e+00, -2.329443703013275524e+00,
        -2.315080362516850343e+00, -2.300711895443124710e+00, -2.286338223996969177e+00, -2.271959269136864634e+00, -2.257574950553179249e+00, -2.243185186646586615e+00, -2.228789894504906499e+00,
        -2.214388989881381331e+00, -2.199982387170757114e+00, -2.185569999385429174e+00, -2.171151738131902764e+00, -2.156727513586190081e+00, -2.142297234468256484e+00, -2.127860808017024041e+00,
        -2.113418139964494014e+00, -2.098969134508561929e+00, -2.084513694286517005e+00, -2.070051720347950930e+00, -2.055583112125752177e+00, -2.041107767408548490e+00, -2.026625582310815243e+00,
        -2.012136451243868418e+00, -1.997640266884977533e+00, -1.983136920147360760e+00, -1.968626300148024644e+00, -1.954108294176254423e+00, -1.939582787661000784e+00, -1.925049664137404415e+00,
        -1.910508805213203543e+00, -1.895960090533721276e+00, -1.881403397747692496e+00, -1.866838602470123343e+00, -1.852265578246936162e+00, -1.837684196517090696e+00, -1.823094326574594470e+00,
        -1.808495835530226525e+00, -1.793888588271750573e+00, -1.779272447423927428e+00, -1.764647273307479391e+00, -1.750012923897452444e+00, -1.735369254780609216e+00, -1.720716119112000175e+00,
        -1.706053367570993240e+00, -1.691380848315968466e+00, -1.676698406938626595e+00, -1.662005886417345257e+00, -1.647303127068936668e+00, -1.632589966500902712e+00, -1.617866239561131625e+00,
        -1.603131778288010123e+00, -1.588386411858628611e+00, -1.573629966536664426e+00, -1.558862265618827125e+00, -1.544083129380501962e+00, -1.529292375020237182e+00, -1.514489816603383998e+00,
        -1.499675265004333014e+00, -1.484848527848137811e+00, -1.470009409450440341e+00, -1.455157710756823208e+00, -1.440293229280589449e+00, -1.425415759039740937e+00, -1.410525090492460665e+00,
        -1.395621010471538304e+00, -1.380703302117652687e+00, -1.365771744811192123e+00, -1.350826114103020892e+00, -1.335866181643817541e+00, -1.320891715112248566e+00, -1.305902478141627965e+00,
        -1.290898230245432377e+00, -1.275878726741366043e+00, -1.260843718673973157e+00, -1.245792952735989001e+00, -1.230726171188135876e+00, -1.215643111777624297e+00, -1.200543507654948661e+00,
        -1.185427087289522730e+00, -1.170293574383549196e+00, -1.155142687784439515e+00, -1.139974141395854179e+00, -1.124787644086856497e+00, -1.109582899599844241e+00, -1.094359606456657197e+00,
        -1.079117457863106644e+00, -1.063856141612011585e+00, -1.048575339984358035e+00, -1.033274729649019674e+00, -1.017953981560683063e+00, -1.002612760856079754e+00, -9.872507267485509663e-01,
        -9.718675324209239408e-01, -9.564628249165176843e-01, -9.410362450286021696e-01, -9.255874271879587223e-01, -9.101159993487000222e-01, -8.946215828724155550e-01, -8.791037924103899392e-01,
        -8.635622357841905572e-01, -8.479965138643655864e-01, -8.324062204473884341e-01, -8.167909421308313656e-01, -8.011502581867100403e-01, -7.854837404331064254e-01, -7.697909531039377473e-01,
        -7.540714527170240267e-01, -7.383247879403610492e-01, -7.225504994566308570e-01, -7.067481198260057162e-01, -6.909171733472264654e-01, -6.750571759169931019e-01, -6.591676348876497338e-01,
        -6.432480489232650367e-01, -6.272979078540815712e-01, -6.113166925293762599e-01, -5.953038746687773219e-01, -5.792589167121133809e-01, -5.631812716677440100e-01, -5.470703829595495726e-01,
        -5.309256842725363912e-01, -5.147465993971619413e-01, -4.985325420724325274e-01, -4.822829158278582606e-01, -4.659971138243255151e-01, -4.496745186939796191e-01, -4.333145023792023820e-01,
        -4.169164259707845432e-01, -4.004796395453985025e-01, -3.840034820024536555e-01, -3.674872809005117480e-01, -3.509303522933207575e-01, -3.343320005656342242e-01, -3.176915182689690753e-01,
        -3.010081859574224028e-01, -2.842812720237357094e-01, -2.675100325357526176e-01, -2.506937110734680507e-01, -2.338315385668545687e-01, -2.169227331346592069e-01, -1.999664999244038510e-01,
        -1.829620309537795531e-01, -1.659085049536948575e-01, -1.488050872132367364e-01, -1.316509294267726449e-01, -1.144451695435026312e-01, -9.718693161971100891e-02, -7.987532567406357975e-02,
        -6.250944754623571908e-02, -4.508837875921357929e-02, -2.761118638561238861e-02, -1.007692291837952031e-02, 7.515373853751791157e-03, 2.516668096563849308e-02, 4.287799038955116687e-02,
        6.065030914300877790e-02, 7.848465939016274762e-02, 9.638207853628300015e-02, 1.143436193148676949e-01, 1.323703498668596101e-01, 1.504633538114835134e-01, 1.686237303081316863e-01,
        1.868525941087738618e-01, 2.051510756003103453e-01, 2.235203208362644411e-01, 2.419614915572445846e-01, 2.604757651994907275e-01, 2.790643348909225274e-01, 2.977284094339847642e-01,
        3.164692132746190767e-01, 3.352879864566974955e-01, 3.541859845611505797e-01, 3.731644786291101745e-01, 3.922247550683047868e-01, 4.113681155419415814e-01, 4.305958768393399749e-01,
        4.499093707274900988e-01, 4.693099437827778497e-01, 4.887989572020708939e-01, 5.083777865923346795e-01, 5.280478217380019101e-01, 5.478104663452385559e-01, 5.676671377623091486e-01,
        5.876192666752143579e-01, 6.076682967777744526e-01, 6.278156844153628402e-01, 6.480628982014593475e-01, 6.684114186062461993e-01, 6.888627375164465549e-01, 7.094183577656660855e-01,
        7.300797926344519961e-01, 7.508485653193812670e-01, 7.717262083704640174e-01, 7.927142630961924175e-01, 8.138142789356321849e-01, 8.350278127969287256e-01, 8.563564283617107753e-01,
        8.778016953548820611e-01, 8.993651887793497890e-01, 9.210484881153276904e-01, 9.428531764838573581e-01, 9.647808397743138364e-01, 9.868330657357101687e-01, 1.009011443031680555e+00,
        1.031317560259154220e+00, 1.053753004930747306e+00, 1.076319362421058745e+00, 1.099018214877129251e+00, 1.121851140093406807e+00, 1.144819710351726805e+00, 1.167925491226848189e+00,
        1.191170040358282112e+00, 1.214554906189223793e+00, 1.238081626673530833e+00, 1.261751727951857527e+00, 1.285566722998118916e+00, 1.309528110237661025e+00, 1.333637372138616994e+00,
        1.357895973778068477e+00, 1.382305361384801090e+00, 1.406866960860538418e+00, 1.431582176281727126e+00, 1.456452388384066099e+00, 1.481478953032097312e+00, 1.506663199676377740e+00,
        1.532006429800786274e+00, 1.557509915362758957e+00, 1.583174897229287392e+00, 1.609002583611719883e+00, 1.634994148502431299e+00, 1.661150730116619734e+00, 1.687473429342547293e+00,
        1.713963308203631852e+00, 1.740621388335936937e+00, 1.767448649484585133e+00, 1.794446028022785189e+00, 1.821614415497148087e+00, 1.848954657203021990e+00, 1.876467550793639250e+00,
        1.904153844926801176e+00, 1.932014237952910163e+00, 1.960049376648085095e+00, 1.988259854996070608e+00, 2.016646213022645018e+00, 2.045208935686086527e+00, 2.073948451827253514e+00,
        2.102865133182698631e+00, 2.131959293464111393e+00, 2.161231187507307805e+00, 2.190681010493759295e+00, 2.220308897247562196e+00, 2.250114921610536278e+00, 2.280099095897945016e+00,
        2.310261370437167727e+00, 2.340601633191360609e+00, 2.371119709469981629e+00, 2.401815361727784648e+00, 2.432688289453607489e+00, 2.463738129150097489e+00, 2.494964454405129484e+00,
        2.526366776055511743e+00, 2.557944542443212121e+00, 2.589697139764099543e+00, 2.621623892508883724e+00, 2.653724063995622107e+00, 2.685996856992998172e+00, 2.718441414433075298e+00,
        2.751056820212167686e+00, 2.783842100078050663e+00, 2.816796222601529287e+00, 2.849918100230090179e+00, 2.883206590421129256e+00, 2.916660496851998019e+00, 2.950278570703887215e+00,
        2.984059512016312254e+00, 3.018001971108888792e+00, 3.052104550066684130e+00, 3.086365804285480685e+00, 3.120784244073000924e+00, 3.155358336302044542e+00, 3.190086506111383446e+00,
        3.224967138650136622e+00, 3.259998580861284800e+00, 3.295179143299920455e+00, 3.330507101981744711e+00, 3.365980700257427571e+00, 3.401598150708227575e+00, 3.437357637058491289e+00,
        3.473257316100561987e+00, 3.509295319627706355e+00, 3.545469756370736469e+00, 3.581778713934076208e+00, 3.618220260727126991e+00, 3.654792447886851203e+00, 3.691493311187746151e+00,
        3.728320872935293906e+00, 3.765273143839341863e+00, 3.802348124863886181e+00, 3.839543809049948297e+00, 3.876858183308386607e+00, 3.914289230179682999e+00, 3.951834929557909604e+00,
        3.989493260376288930e+00, 4.027262202251907119e+00, 4.065139737087459793e+00, 4.103123850627900637e+00, 4.141212533970256793e+00, 4.179403785024969231e+00, 4.217695609927317513e+00,
        4.256086024397712109e+00, 4.294573055049793808e+00, 4.333154740645484004e+00, 4.371829133296307290e+00, 4.410594299610440139e+00, 4.449448321785244964e+00, 4.488389298644973380e+00,
        4.527415346623724801e+00, 4.566524600693744773e+00, 4.605715215239329474e+00, 4.644985364876748513e+00, 4.684333245220703823e+00, 4.723757073597968592e+00, 4.763255089708972712e+00,
        4.802825556238135896e+00, 4.842466759414006283e+00, 4.882177009520106914e+00, 4.921954641357675797e+00, 4.961798014661471079e+00, 5.001705514469860780e+00, 5.041675551450504500e+00,
        5.081706562182967346e+00, 5.121797009399650769e+00, 5.161945382186415188e+00, 5.202150196144439853e+00, 5.242409993514641542e+00, 5.282723343266237492e+00, 5.323088841150888761e+00,
        5.363505109723930353e+00, 5.403970798334178482e+00, 5.444484583083772478e+00, 5.485045166759582536e+00, 5.525651278737584882e+00, 5.566301674861654725e+00, 5.606995137298286913e+00,
        5.647730474368503728e+00, 5.688506520358434848e+00, 5.729322135309868003e+00, 5.770176204792100805e+00, 5.811067639656364747e+00, 5.851995375774075470e+00, 5.892958373760132318e+00,
        5.933955618682425381e+00, 5.974986119758676217e+00, 6.016048910041793540e+00, 6.057143046094681793e+00, 6.098267607655615485e+00, 6.139421697295129654e+00, 6.180604440065376792e+00,
        6.221814983142841982e+00, 6.263052495465290193e+00, 6.304316167363771761e+00, 6.345605210190464085e+00, 6.386918855943076956e+00, 6.428256356886620004e+00, 6.469616985173070489e+00,
        6.511000032459687503e+00, 6.552404809526525931e+00, 6.593830645893738129e+00, 6.635276889439174042e+00, 6.676742906016807311e+00, 6.718228079076448367e+00, 6.759731809285145054e+00,
        6.801253514150790380e+00, 6.842792627648170978e+00, 6.884348599847924355e+00, 6.925920896548682926e+00, 6.967508998912713913e+00, 7.009112403105341116e+00, 7.050730619938390120e+00,
        7.092363174517901214e+00, 7.134009605896315165e+00, 7.175669466729295287e+00, 7.217342322937444266e+00, 7.259027753372938818e+00, 7.300725349491324678e+00, 7.342434715028535841e+00,
        7.384155465683258868e+00, 7.425887228804723961e+00, 7.467629643086000080e+00, 7.509382358262848278e+00, 7.551145034818199875e+00, 7.592917343692232812e+00, 7.634698965998204301e+00,
        7.676489592743889645e+00, 7.718288924558750708e+00, 7.760096671426786941e+00, 7.801912552425056546e+00, 7.843736295467861552e+00, 7.885567637056555945e+00, 7.927406322034950215e+00,
        7.969252103350283001e+00, 8.011104741819659480e+00, 8.052964005902026656e+00, 8.094829671475483934e+00, 8.136701521619965405e+00, 8.178579346405202344e+00, 8.220462942683900209e+00,
        8.262352113890059968e+00, 8.304246669842363815e+00, 8.346146426552573772e+00, 8.388051206038806384e+00, 8.429960836143713721e+00, 8.471875150357364603e+00, 8.513793987644810102e+00,
        8.555717192278262218e+00, 8.597644613673752545e+00, 8.639576106232228270e+00, 8.681511529184973597e+00, 8.723450746443274895e+00, 8.765393626452272713e+00, 8.807340042048823037e+00,
        8.849289870323442386e+00, 8.891242992486070307e+00, 8.933199293735690105e+00, 8.975158663133653292e+00, 9.017120993480681079e+00, 9.059086181197388399e+00, 9.101054126208316220e+00,
        9.143024731829335394e+00, 9.184997904658382950e+00, 9.226973554469383387e+00, 9.268951594109404724e+00, 9.310931939398793489e+00, 9.352914509034373935e+00, 9.394899224495514289e+00,
        9.436886009953081356e+00, 9.478874792181123610e+00, 9.520865500471277443e+00, 9.562858066549781100e+00, 9.604852424497055452e+00, 9.646848510669737919e+00, 9.688846263625196187e+00,
        9.730845624048324538e+00, 9.772846534680663666e+00, 9.814848940251730269e+00, 9.856852787412515582e+00, 9.898858024671071121e+00, 9.940864602330142574e+00, 9.982872472426784327e+00,
        1.002488158867386936e+01, 1.006689190640353360e+01, 1.010890338251236287e+01, 1.015091597540838286e+01, 1.019292964495973486e+01, 1.023494435244503542e+01, 1.027696006050531352e+01,
        1.031897673309754104e+01, 1.036099433544965187e+01, 1.040301283401705135e+01, 1.044503219644050596e+01, 1.048705239150549318e+01, 1.052907338910280899e+01, 1.057109516019052187e+01,
        1.061311767675715778e+01, 1.065514091178612688e+01, 1.069716483922131189e+01, 1.073918943393381298e+01, 1.078121467168979386e+01, 1.082324052911941337e+01, 1.086526698368675881e+01,
        1.090729401366085582e+01, 1.094932159808758954e+01, 1.099134971676259376e+01, 1.103337835020506397e+01, 1.107540747963243177e+01, 1.111743708693592225e+01, 1.115946715465691952e+01,
        1.120149766596414409e+01, 1.124352860463161008e+01, 1.128555995501730358e+01, 1.132759170204265331e+01, 1.136962383117263720e+01, 1.141165632839659239e+01, 1.145368918020970028e+01,
        1.149572237359508264e+01, 1.153775589600653184e+01, 1.157978973535181844e+01, 1.162182387997659738e+01, 1.166385831864882938e+01, 1.170589304054380619e+01, 1.174792803522962892e+01,
        1.178996329265322451e+01, 1.183199880312684549e+01, 1.187403455731502611e+01, 1.191607054622200579e+01, 1.195810676117958415e+01, 1.200014319383540240e+01
    };
