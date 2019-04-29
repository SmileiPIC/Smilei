#include "SpeciesV.h"

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

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Species
// input: simulation parameters & Species index
// ---------------------------------------------------------------------------------------------------------------------
SpeciesV::SpeciesV( Params &params, Patch *patch ) :
    Species( params, patch )
{
    initCluster( params );
    npack_ = 0 ;
    packsize_ = 0;
    
    
}//END SpeciesV creator

// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species
// ---------------------------------------------------------------------------------------------------------------------
SpeciesV::~SpeciesV()
{
}


void SpeciesV::initCluster( Params &params )
{
    int ncells = 1;
    for( unsigned int iDim=0 ; iDim<nDim_particle ; iDim++ ) {
        ncells *= ( params.n_space[iDim]+1 );
    }
    last_index.resize( ncells, 0 );
    first_index.resize( ncells, 0 );
    count.resize( ncells, 0 );
    
    //Size in each dimension of the buffers on which each bin are projected
    //In 1D the particles of a given bin can be projected on 6 different nodes at the second order (oversize = 2)
    
    //Primal dimension of fields.
    f_dim0 =  params.n_space[0] + 2 * oversize[0] +1;
    f_dim1 =  params.n_space[1] + 2 * oversize[1] +1;
    f_dim2 =  params.n_space[2] + 2 * oversize[2] +1;
    
    b_dim.resize( params.nDim_field, 1 );
    if( nDim_particle == 1 ) {
        b_dim[0] = ( 1 + clrw ) + 2 * oversize[0];
        f_dim1 = 1;
        f_dim2 = 1;
    }
    if( nDim_particle == 2 ) {
        b_dim[0] = ( 1 + clrw ) + 2 * oversize[0]; // There is a primal number of bins.
        b_dim[1] =  f_dim1;
        f_dim2 = 1;
    }
    if( nDim_particle == 3 ) {
        b_dim[0] = ( 1 + clrw ) + 2 * oversize[0]; // There is a primal number of bins.
        b_dim[1] = f_dim1;
        b_dim[2] = f_dim2;
    }
    
    //Initialize specMPI
    MPIbuff.allocate( nDim_particle );
    
    //ener_tot = 0.;
    nrj_bc_lost = 0.;
    nrj_mw_lost = 0.;
    nrj_new_particles = 0.;
    
}//END initCluster


void SpeciesV::dynamics( double time_dual, unsigned int ispec,
                         ElectroMagn *EMfields, Params &params, bool diag_flag,
                         PartWalls *partWalls,
                         Patch *patch, SmileiMPI *smpi,
                         RadiationTables &RadiationTables,
                         MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables,
                         vector<Diagnostic *> &localDiags )
{
    int ithread;
#ifdef _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif
    
#ifdef  __DETAILED_TIMERS
    double timer;
#endif
    
    if( npack_==0 ) {
        npack_    = 1;
        packsize_ = ( f_dim1-2*oversize[1] );
        
        //if ( ( (long int)last_index.back() < (long int)60000 ) || (Radiate) || (Ionize) || (Multiphoton_Breit_Wheeler_process) )
        packsize_ *= ( f_dim0-2*oversize[0] );
        //else
        //    npack_ *= (f_dim0-2*oversize[0]);
        
        if( nDim_particle == 3 ) {
            packsize_ *= ( f_dim2-2*oversize[2] );
        }
    }
    
    unsigned int iPart;
    
    // Reset list of particles to exchange
    clearExchList();
    
    int tid( 0 );
    double ener_iPart( 0. );
    std::vector<double> nrj_lost_per_thd( 1, 0. );
    
    // -------------------------------
    // calculate the particle dynamics
    // -------------------------------
    if( time_dual>time_frozen ) { // moving particle
    
        smpi->dynamics_resize( ithread, nDim_field, last_index.back(), params.geometry=="AMcylindrical" );
        
        //Point to local thread dedicated buffers
        //Still needed for ionization
        vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
        
        //Prepare for sorting
        for( unsigned int i=0; i<count.size(); i++ ) {
            count[i] = 0;
        }
        
        for( unsigned int ipack = 0 ; ipack < npack_ ; ipack++ ) {
        
            int nparts_in_pack = last_index[( ipack+1 ) * packsize_-1 ];
            smpi->dynamics_resize( ithread, nDim_particle, nparts_in_pack );
            
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            
            // Interpolate the fields at the particle position
            for( unsigned int scell = 0 ; scell < packsize_ ; scell++ )
                Interp->fieldsWrapper( EMfields, *particles, smpi, &( first_index[ipack*packsize_+scell] ),
                                       &( last_index[ipack*packsize_+scell] ),
                                       ithread, first_index[ipack*packsize_] );
                                       
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[0] += MPI_Wtime() - timer;
#endif
            
            // Ionization
            if( Ionize ) {
#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif
                for( unsigned int scell = 0 ; scell < first_index.size() ; scell++ ) {
                    ( *Ionize )( particles, first_index[scell], last_index[scell], Epart, patch, Proj );
                }
#ifdef  __DETAILED_TIMERS
                patch->patch_timers[4] += MPI_Wtime() - timer;
#endif
            }
            
            // Radiation losses
            if( Radiate ) {
#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif
                for( unsigned int scell = 0 ; scell < first_index.size() ; scell++ ) {
                    // Radiation process
                    ( *Radiate )( *particles, this->photon_species, smpi,
                                  RadiationTables,
                                  first_index[scell], last_index[scell], ithread );
                                  
                    // Update scalar variable for diagnostics
                    nrj_radiation += Radiate->getRadiatedEnergy();
                    
                    // Update the quantum parameter chi
                    Radiate->computeParticlesChi( *particles,
                                                  smpi,
                                                  first_index[scell],
                                                  last_index[scell],
                                                  ithread );
                }
#ifdef  __DETAILED_TIMERS
                patch->patch_timers[5] += MPI_Wtime() - timer;
#endif
            }
            
            // Multiphoton Breit-Wheeler
            if( Multiphoton_Breit_Wheeler_process ) {
#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif
                for( unsigned int scell = 0 ; scell < first_index.size() ; scell++ ) {
                
                    // Pair generation process
                    ( *Multiphoton_Breit_Wheeler_process )( *particles,
                                                            smpi,
                                                            MultiphotonBreitWheelerTables,
                                                            first_index[scell], last_index[scell], ithread );
                                                            
                    // Update scalar variable for diagnostics
                    // We reuse nrj_radiation for the pairs
                    nrj_radiation += Multiphoton_Breit_Wheeler_process->getPairEnergy();
                    
                    // Update the photon quantum parameter chi of all photons
                    Multiphoton_Breit_Wheeler_process->compute_thread_chiph( *particles,
                            smpi,
                            first_index[scell],
                            last_index[scell],
                            ithread );
                            
                    // Suppression of the decayed photons into pairs
                    Multiphoton_Breit_Wheeler_process->decayed_photon_cleaning(
                        *particles, smpi, scell, first_index.size(), &first_index[0], &last_index[0], ithread );
                        
                }
#ifdef  __DETAILED_TIMERS
                patch->patch_timers[6] += MPI_Wtime() - timer;
#endif
            }
            
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            
            // Push the particles and the photons
            ( *Push )( *particles, smpi, first_index[ipack*packsize_],
                       last_index[ipack*packsize_+packsize_-1],
                       ithread, first_index[ipack*packsize_] );
                       
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[1] += MPI_Wtime() - timer;
            timer = MPI_Wtime();
#endif
            
            unsigned int length[3];
            length[0]=0;
            length[1]=params.n_space[1]+1;
            length[2]=params.n_space[2]+1;
            
            for( unsigned int scell = 0 ; scell < packsize_ ; scell++ ) {
                // Apply wall and boundary conditions
                if( mass>0 ) {
                    for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                        for( iPart=first_index[scell] ; ( int )iPart<last_index[scell]; iPart++ ) {
                            double dtgf = params.timestep * smpi->dynamics_invgf[ithread][iPart];
                            if( !( *partWalls )[iwall]->apply( *particles, iPart, this, dtgf, ener_iPart ) ) {
                                nrj_lost_per_thd[tid] += mass * ener_iPart;
                            }
                        }
                    }
                    
                    // Boundary Condition may be physical or due to domain decomposition
                    // apply returns 0 if iPart is not in the local domain anymore
                    
                    for( iPart=first_index[ipack*packsize_+scell] ; ( int )iPart<last_index[ipack*packsize_+scell]; iPart++ ) {
                        if( !partBoundCond->apply( *particles, iPart, this, ener_iPart ) ) {
                            addPartInExchList( iPart );
                            nrj_lost_per_thd[tid] += mass * ener_iPart;
                            particles->cell_keys[iPart] = -1;
                        } else {
                            //Compute cell_keys of remaining particles
                            for( unsigned int i = 0 ; i<nDim_particle; i++ ) {
                                particles->cell_keys[iPart] *= this->length_[i];
                                particles->cell_keys[iPart] += round( ( particles->position( i, iPart )-min_loc_vec[i] ) * dx_inv_[i] );
                            }
                            //First reduction of the count sort algorithm. Lost particles are not included.
                            count[particles->cell_keys[iPart]] ++;
                        }
                    }
                    
                } else if( mass==0 ) {
                    for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                        for( iPart=first_index[scell] ; ( int )iPart<last_index[scell]; iPart++ ) {
                            double dtgf = params.timestep * smpi->dynamics_invgf[ithread][iPart];
                            if( !( *partWalls )[iwall]->apply( *particles, iPart, this, dtgf, ener_iPart ) ) {
                                nrj_lost_per_thd[tid] += ener_iPart;
                            }
                        }
                    }
                    
                    // Boundary Condition may be physical or due to domain decomposition
                    // apply returns 0 if iPart is not in the local domain anymore
                    for( iPart=first_index[scell] ; ( int )iPart<last_index[scell]; iPart++ ) {
                        if( !partBoundCond->apply( *particles, iPart, this, ener_iPart ) ) {
                            addPartInExchList( iPart );
                            nrj_lost_per_thd[tid] += ener_iPart;
                            particles->cell_keys[iPart] = -1;
                        } else {
                            //Compute cell_keys of remaining particles
                            for( unsigned int i = 0 ; i<nDim_particle; i++ ) {
                                particles->cell_keys[iPart] *= length[i];
                                particles->cell_keys[iPart] += round( ( particles->position( i, iPart )-min_loc_vec[i] ) * dx_inv_[i] );
                            }
                            //First reduction of the count sort algorithm. Lost particles are not included.
                            count[particles->cell_keys[iPart]] ++;
                        }
                    }
                }
            }
            //START EXCHANGE PARTICLES OF THE CURRENT BIN ?
            
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[3] += MPI_Wtime() - timer;
#endif
            
            // Project currents if not a Test species and charges as well if a diag is needed.
            // Do not project if a photon
            if( ( !particles->is_test ) && ( mass > 0 ) )
#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif
                
            for( unsigned int scell = 0 ; scell < packsize_ ; scell++ )
                Proj->currentsAndDensityWrapper(
                    EMfields, *particles, smpi, first_index[ipack*packsize_+scell],
                    last_index[ipack*packsize_+scell],
                    ithread,
                    diag_flag, params.is_spectral,
                    ispec, ipack*packsize_+scell, first_index[ipack*packsize_]
                );
                
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[2] += MPI_Wtime() - timer;
#endif
            
            for( unsigned int ithd=0 ; ithd<nrj_lost_per_thd.size() ; ithd++ ) {
                nrj_bc_lost += nrj_lost_per_thd[tid];
            }
            
        }
        
    } else { // immobile particle (at the moment only project density)
        if( diag_flag &&( !particles->is_test ) ) {
            double *b_rho=nullptr;
            
            for( unsigned int scell = 0 ; scell < first_index.size() ; scell ++ ) { //Loop for projection on buffer_proj
            
                b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
                
                for( iPart=first_index[scell] ; ( int )iPart<last_index[scell]; iPart++ ) {
                    Proj->basic( b_rho, ( *particles ), iPart, 0 );
                } //End loop on particles
            }//End loop on bins
            
        }
    }//END if time vs. time_frozen
    
}//END dynamic


// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species
//   - increment the charge (projection)
//   - used at initialisation for Poisson (and diags if required, not for now dynamics )
// ---------------------------------------------------------------------------------------------------------------------
void SpeciesV::computeCharge( unsigned int ispec, ElectroMagn *EMfields )
{
    // -------------------------------
    // calculate the particle charge
    // -------------------------------
    if( ( !particles->is_test ) ) {
        double *b_rho=&( *EMfields->rho_ )( 0 );
        
        for( unsigned int iPart=first_index[0] ; ( int )iPart<last_index[last_index.size()-1]; iPart++ ) {
            Proj->basic( b_rho, ( *particles ), iPart, 0 );
        }
        
    }
    
}//END computeCharge


// ---------------------------------------------------------------------------------------------------------------------
// Sort particles
// ---------------------------------------------------------------------------------------------------------------------
void SpeciesV::sort_part( Params &params )
{
    unsigned int npart, ncell;
    int ip_dest, cell_target;
    unsigned int length[3];
    vector<int> buf_cell_keys[3][2];
    std::vector<unsigned int> cycle;
    unsigned int ip_src;
    
    length[0]=0;
    length[1]=params.n_space[1]+1;
    length[2]=params.n_space[2]+1;
    
    //Number of dual cells
    ncell = ( params.n_space[0]+1 );
    for( unsigned int i=1; i < params.nDim_field; i++ ) {
        ncell *= length[i];
    }
    
    //Number of particles before exchange
    npart = particles->size();
    
    //Loop over just arrived particles to compute their cell keys and contribution to count
    for( unsigned int idim=0; idim < nDim_particle ; idim++ ) {
        for( unsigned int ineighbor=0 ; ineighbor < 2 ; ineighbor++ ) {
            buf_cell_keys[idim][ineighbor].resize( MPIbuff.part_index_recv_sz[idim][ineighbor] );
            #pragma omp simd
            for( unsigned int ip=0; ip < MPIbuff.part_index_recv_sz[idim][ineighbor]; ip++ ) {
                for( unsigned int ipos=0; ipos < nDim_particle ; ipos++ ) {
                    double X = MPIbuff.partRecv[idim][ineighbor].position( ipos, ip )-min_loc_vec[ipos];
                    int IX = round( X * dx_inv_[ipos] );
                    buf_cell_keys[idim][ineighbor][ip] = buf_cell_keys[idim][ineighbor][ip] * length[ipos] + IX;
                }
            }
            //Can we vectorize this reduction ?
            for( unsigned int ip=0; ip < MPIbuff.part_index_recv_sz[idim][ineighbor]; ip++ ) {
                count[buf_cell_keys[idim][ineighbor][ip]] ++;
            }
        }
    }
    
    // second loop convert the count array in cumulative sum
    first_index[0]=0;
    for( unsigned int ic=1; ic < ncell; ic++ ) {
        first_index[ic] = first_index[ic-1] + count[ic-1];
        last_index[ic-1]= first_index[ic];
    }
    
    //New total number of particles is stored as last element of last_index
    last_index[ncell-1] = last_index[ncell-2] + count.back() ;
    
    //Now proceed to the cycle sort
    
    if( MPIbuff.partRecv[0][0].size() == 0 ) {
        MPIbuff.partRecv[0][0].initialize( 0, *particles );    //Is this correct ?
    }
    
    // Resize the particle vector
    if( ( unsigned int )last_index.back() > npart ) {
        particles->resize( last_index.back(), nDim_particle );
        particles->cell_keys.resize( last_index.back(), -1 ); // Merge this in particles.resize(..) ?
        for( unsigned int ipart = npart; ipart < ( unsigned int )last_index.back(); ipart ++ ) {
            addPartInExchList( ipart );
        }
    }
    
    //Copy all particles from MPI buffers back to the writable particles via cycle sort pass.
    for( unsigned int idim=0; idim < nDim_particle ; idim++ ) {
        for( unsigned int ineighbor=0 ; ineighbor < 2 ; ineighbor++ ) {
            for( unsigned int ip=0; ip < MPIbuff.part_index_recv_sz[idim][ineighbor]; ip++ ) {
                cycle.resize( 1 );
                cell_target = buf_cell_keys[idim][ineighbor][ip];
                ip_dest = first_index[cell_target];
                while( particles->cell_keys[ip_dest] == cell_target ) {
                    ip_dest++;
                }
                first_index[cell_target] = ip_dest + 1 ;
                cycle[0] = ip_dest;
                cell_target = particles->cell_keys[ip_dest];
                //As long as the particle is not erased, we can build up the cycle.
                while( cell_target != -1 ) {
                    ip_dest = first_index[cell_target];
                    while( particles->cell_keys[ip_dest] == cell_target ) {
                        ip_dest++;
                    }
                    first_index[cell_target] = ip_dest + 1 ;
                    cycle.push_back( ip_dest );
                    cell_target = particles->cell_keys[ip_dest];
                }
                //Last target_cell is -1, the particle must be erased:
                particles->translate_parts( cycle );
                //Eventually copy particle from the MPI buffer into the particle vector.
                MPIbuff.partRecv[idim][ineighbor].overwrite_part( ip, *particles, cycle[0] );
            }
        }
    }
    
    //Copy valid particles siting over last_index.back() back into the real particles array (happens when more particles are lost than received)
    for( unsigned int ip=( unsigned int )last_index.back(); ip < npart; ip++ ) {
        cell_target = particles->cell_keys[ip];
        if( cell_target == -1 ) {
            continue;
        }
        cycle.resize( 0 );
        cycle.push_back( ip );
        //As long as the particle is not erased, we can build up the cycle.
        while( cell_target != -1 ) {
            ip_dest = first_index[cell_target];
            while( particles->cell_keys[ip_dest] == cell_target ) {
                ip_dest++;
            }
            first_index[cell_target] = ip_dest + 1 ;
            cycle.push_back( ip_dest );
            cell_target = particles->cell_keys[ip_dest];
        }
        //Last target_cell is -1, the particle must be erased:
        particles->translate_parts( cycle );
    }
    
    // Resize the particle vector
    if( ( unsigned int )last_index.back() < npart ) {
        particles->resize( last_index.back(), nDim_particle );
        particles->cell_keys.resize( last_index.back() ); // Merge this in particles.resize(..) ?
    }
    
    
    //Loop over all cells
    for( int icell = 0 ; icell < ( int )ncell; icell++ ) {
        for( unsigned int ip=( unsigned int )first_index[icell]; ip < ( unsigned int )last_index[icell] ; ip++ ) {
            //update value of current cell 'icell' if necessary
            //if particle changes cell, build a cycle of exchange as long as possible. Treats all particles
            if( particles->cell_keys[ip] != icell ) {
                cycle.resize( 1 );
                cycle[0] = ip;
                ip_src = ip;
                //While the destination particle is not going out of the patch or back to the initial cell, keep building the cycle.
                while( particles->cell_keys[ip_src] != icell ) {
                    //Scan the next cell destination
                    ip_dest = first_index[particles->cell_keys[ip_src]];
                    while( particles->cell_keys[ip_dest] == particles->cell_keys[ip_src] ) {
                        ip_dest++;
                    }
                    //In the destination cell, if a particle is going out of this cell, add it to the cycle.
                    first_index[particles->cell_keys[ip_src]] = ip_dest + 1 ;
                    cycle.push_back( ip_dest );
                    ip_src = ip_dest; //Destination becomes source for the next iteration
                }
                //swap parts
                particles->swap_parts( cycle );
            }
        }
    } //end loop on cells
    // Restore first_index initial value
    first_index[0]=0;
    for( unsigned int ic=1; ic < ncell; ic++ ) {
        first_index[ic] = last_index[ic-1];
    }
    
}


void SpeciesV::compute_part_cell_keys( Params &params )
{
    //Compute part_cell_keys at patch creation. This operation is normally done in the pusher to avoid additional particles pass.
    
    unsigned int ip, npart;
    int IX;
    double X;
    
    npart = particles->size(); //Number of particles
    
    #pragma omp simd
    for( ip=0; ip < npart ; ip++ ) {
        // Counts the # of particles in each cell (or sub_cell) and store it in slast_index.
        for( unsigned int ipos=0; ipos < nDim_particle ; ipos++ ) {
            X = particles->position( ipos, ip )-min_loc_vec[ipos];
            IX = round( X * dx_inv_[ipos] );
            particles->cell_keys[ip] = particles->cell_keys[ip] * this->length_[ipos] + IX;
        }
    }
    for( ip=0; ip < npart ; ip++ ) {
        count[particles->cell_keys[ip]] ++ ;
    }
    
}

// -----------------------------------------------------------------------------
//! Compute cell_keys for the specified bin boundaries.
//! params object that contains the global parameters
//! istart first bin index
//! iend last bin index
// -----------------------------------------------------------------------------
void SpeciesV::compute_bin_cell_keys( Params &params, int istart, int iend )
{
    // Resize of cell_keys seems necessary here
    particles->cell_keys.resize( particles->size() );
    
    #pragma omp simd
    for( int ip=istart; ip < iend; ip++ ) {
        // Counts the # of particles in each cell (or sub_cell) and store it in slast_index.
        for( unsigned int ipos=0; ipos < nDim_particle ; ipos++ ) {
            particles->cell_keys[ip] *= this->length_[ipos];
            particles->cell_keys[ip] += round( ( particles->position( ipos, ip )-min_loc_vec[ipos] ) * dx_inv_[ipos] );
        }
    }
}

void SpeciesV::importParticles( Params &params, Patch *patch, Particles &source_particles, vector<Diagnostic *> &localDiags )
{

    unsigned int npart = source_particles.size(), scell, ii, nbin=first_index.size();
    
    // If this species is tracked, set the particle IDs
    if( particles->tracked ) {
        dynamic_cast<DiagnosticTrack *>( localDiags[tracking_diagnostic] )->setIDs( source_particles );
    }
    
    unsigned int length[3];
    length[0]=0;
    length[1]=params.n_space[1]+1;
    length[2]=params.n_space[2]+1;
    
    int IX;
    double X;
    
    // std::cerr << "SpeciesV::importParticles "
    //           << " for "<< this->name
    //           << " in patch (" << patch->Pcoordinates[0] << "," <<  patch->Pcoordinates[1] << "," <<  patch->Pcoordinates[2] << ") "
    //           << " mpi process " << patch->MPI_me_ << " - "
    //           << " mode: " << this->vectorized_operators << " - "
    //           << " nb bin: " << first_index.size() << " - "
    //           << " nbp: " << npart
    //           << std::endl;
    
    // Move particles
    for( unsigned int i=0; i<npart; i++ ) {
    
        // Compute the receiving bin index
        scell = 0;
        for( unsigned int ipos=0; ipos < nDim_particle ; ipos++ ) {
            X = source_particles.position( ipos, i )-min_loc_vec[ipos];
            IX = round( X * dx_inv_[ipos] );
            scell = scell * length[ipos] + IX;
        }
        
        // Copy particle to the correct bin
        source_particles.cp_particle( i, *particles, last_index[scell] );
        
        // Update the bin counts
        last_index[scell]++;
        for( ii=scell+1; ii<nbin; ii++ ) {
            first_index[ii]++;
            last_index[ii]++;
        }
        
        particles->cell_keys.insert( particles->cell_keys.begin() + first_index[scell] + count[scell], scell );
        count[scell] ++ ;
        
    }
    
    source_particles.clear();
}

// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species reacting to laser envelope
//   - interpolate the fields at the particle position
//   - deposit susceptibility
//   - calculate the new momentum
// ---------------------------------------------------------------------------------------------------------------------
void SpeciesV::ponderomotive_update_susceptibility_and_momentum( double time_dual, unsigned int ispec,
        ElectroMagn *EMfields,
        Params &params, bool diag_flag,
        Patch *patch, SmileiMPI *smpi,
        std::vector<Diagnostic *> &localDiags )
{

////////////////////////////// new vectorized
    int ithread;
#ifdef _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif
    
#ifdef  __DETAILED_TIMERS
    double timer;
#endif
    
    if( npack_==0 ) {
        npack_    = 1;
        packsize_ = ( f_dim1-2*oversize[1] );
        
        //if ( (long int)last_index.back() < (long int)60000 || (Radiate) || (Ionize) || (Multiphoton_Breit_Wheeler_process) )
        packsize_ *= ( f_dim0-2*oversize[0] );
        //else
        //    npack_ *= (f_dim0-2*oversize[0]);
        
        if( nDim_particle == 3 ) {
            packsize_ *= ( f_dim2-2*oversize[2] );
        }
    }
    
    
    // -------------------------------
    // calculate the particle dynamics
    // -------------------------------
    if( time_dual>time_frozen ) { // advance particle momentum
    
        for( unsigned int ipack = 0 ; ipack < npack_ ; ipack++ ) {
        
            // ipack start @ first_index [ ipack * packsize_ ]
            // ipack end   @ last_index [ ipack * packsize_ + packsize_ - 1 ]
            //int nparts_in_pack = last_index[ (ipack+1) * packsize_-1 ] - first_index [ ipack * packsize_ ];
            int nparts_in_pack = last_index[( ipack+1 ) * packsize_-1 ];
            smpi->dynamics_resize( ithread, nDim_particle, nparts_in_pack );
            
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            // Interpolate the fields at the particle position
            for( unsigned int scell = 0 ; scell < packsize_ ; scell++ ) {
                Interp->fieldsAndEnvelope( EMfields, *particles, smpi, &( first_index[ipack*packsize_+scell] ), &( last_index[ipack*packsize_+scell] ), ithread, first_index[ipack*packsize_] );
            }
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[7] += MPI_Wtime() - timer;
#endif
            
            // Project susceptibility, the source term of envelope equation
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            for( unsigned int scell = 0 ; scell < packsize_ ; scell++ ) {
                Proj->susceptibility( EMfields, *particles, mass, smpi, first_index[ipack*packsize_+scell], last_index[ipack*packsize_+scell], ithread, ipack*packsize_+scell, first_index[ipack*packsize_] );
            }
            
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[8] += MPI_Wtime() - timer;
#endif
            
            // Push the particles
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            ( *Push )( *particles, smpi, first_index[ipack*packsize_], last_index[ipack*packsize_+packsize_-1], ithread, first_index[ipack*packsize_] );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[9] += MPI_Wtime() - timer;
#endif
        }
        
    } else { // immobile particle (at the moment only project density)
    
    }//END if time vs. time_frozen
    
} // end ponderomotive_update_susceptibility_and_momentum

// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species reacting to laser envelope
//   - interpolate the fields at the particle position
//   - deposit susceptibility
// ---------------------------------------------------------------------------------------------------------------------
void SpeciesV::ponderomotive_project_susceptibility( double time_dual, unsigned int ispec,
        ElectroMagn *EMfields,
        Params &params, bool diag_flag,
        Patch *patch, SmileiMPI *smpi,
        std::vector<Diagnostic *> &localDiags )
{

////////////////////////////// new vectorized
    int ithread;
#ifdef _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif
    
#ifdef  __DETAILED_TIMERS
    double timer;
#endif
    
    if( npack_==0 ) {
        npack_    = 1;
        packsize_ = ( f_dim1-2*oversize[1] );
        
        //if ( (long int)last_index.back() < (long int)60000 || (Radiate) || (Ionize) || (Multiphoton_Breit_Wheeler_process) )
        packsize_ *= ( f_dim0-2*oversize[0] );
        //else
        //    npack_ *= (f_dim0-2*oversize[0]);
        
        if( nDim_particle == 3 ) {
            packsize_ *= ( f_dim2-2*oversize[2] );
        }
    }
    
    
    // -------------------------------
    // calculate the particle dynamics
    // -------------------------------
    if( time_dual>time_frozen ) { // advance particle momentum
    
        for( unsigned int ipack = 0 ; ipack < npack_ ; ipack++ ) {
        
            // ipack start @ first_index [ ipack * packsize_ ]
            // ipack end   @ last_index [ ipack * packsize_ + packsize_ - 1 ]
            //int nparts_in_pack = last_index[ (ipack+1) * packsize_-1 ] - first_index [ ipack * packsize_ ];
            int nparts_in_pack = last_index[( ipack+1 ) * packsize_-1 ];
            smpi->dynamics_resize( ithread, nDim_particle, nparts_in_pack );
            
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            // Interpolate the fields at the particle position
            for( unsigned int scell = 0 ; scell < packsize_ ; scell++ ) {
                Interp->fieldsAndEnvelope( EMfields, *particles, smpi, &( first_index[ipack*packsize_+scell] ), &( last_index[ipack*packsize_+scell] ), ithread, first_index[ipack*packsize_] );
            }
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[7] += MPI_Wtime() - timer;
#endif
            
            // Project susceptibility, the source term of envelope equation
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            for( unsigned int scell = 0 ; scell < packsize_ ; scell++ ) {
                Proj->susceptibility( EMfields, *particles, mass, smpi, first_index[ipack*packsize_+scell], last_index[ipack*packsize_+scell], ithread, ipack*packsize_+scell, first_index[ipack*packsize_] );
            }
            
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[8] += MPI_Wtime() - timer;
#endif
            
        }
        
    } else { // immobile particle (at the moment only project density)
    
    }//END if time vs. time_frozen
    
} // end ponderomotive_project_susceptibility

// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species reacting to laser envelope
//   - interpolate the ponderomotive potential and its gradient at the particle position, for present and previous timestep
//   - calculate the new particle position
//   - particles BC
//   - project charge and current density
// ---------------------------------------------------------------------------------------------------------------------
void SpeciesV::ponderomotive_update_position_and_currents( double time_dual, unsigned int ispec,
        ElectroMagn *EMfields,
        Params &params, bool diag_flag, PartWalls *partWalls,
        Patch *patch, SmileiMPI *smpi,
        std::vector<Diagnostic *> &localDiags )
{


    int ithread;
#ifdef _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif
    
#ifdef  __DETAILED_TIMERS
    double timer;
#endif
    
    unsigned int iPart;
    
    // Reset list of particles to exchange
    clearExchList();
    
    int tid( 0 );
    double ener_iPart( 0. );
    std::vector<double> nrj_lost_per_thd( 1, 0. );
    
    // -------------------------------
    // calculate the particle dynamics
    // -------------------------------
    if( time_dual>time_frozen ) { // moving particle
    
        //Prepare for sorting
        for( unsigned int i=0; i<count.size(); i++ ) {
            count[i] = 0;
        }
        
        for( unsigned int ipack = 0 ; ipack < npack_ ; ipack++ ) {
        
            //int nparts_in_pack = last_index[ (ipack+1) * packsize_-1 ] - first_index [ ipack * packsize_ ];
            int nparts_in_pack = last_index[( ipack+1 ) * packsize_-1 ];
            smpi->dynamics_resize( ithread, nDim_particle, nparts_in_pack );
            
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            // Interpolate the fields at the particle position
            for( unsigned int scell = 0 ; scell < packsize_ ; scell++ ) {
                Interp->timeCenteredEnvelope( EMfields, *particles, smpi, &( first_index[ipack*packsize_+scell] ), &( last_index[ipack*packsize_+scell] ), ithread, first_index[ipack*packsize_] );
            }
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[10] += MPI_Wtime() - timer;
#endif
            
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            // Push only the particle position
            ( *Push_ponderomotive_position )( *particles, smpi, first_index[ipack*packsize_], last_index[ipack*packsize_+packsize_-1], ithread, first_index[ipack*packsize_] );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[11] += MPI_Wtime() - timer;
            timer = MPI_Wtime();
#endif
            unsigned int length[3];
            length[0]=0;
            length[1]=params.n_space[1]+1;
            length[2]=params.n_space[2]+1;
            
            for( unsigned int scell = 0 ; scell < packsize_ ; scell++ ) {
                // Apply wall and boundary conditions
                if( mass>0 ) { // condition mass>0
                    for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                        for( iPart=first_index[scell] ; ( int )iPart<last_index[scell]; iPart++ ) {
                            double dtgf = params.timestep * smpi->dynamics_invgf[ithread][iPart];
                            if( !( *partWalls )[iwall]->apply( *particles, iPart, this, dtgf, ener_iPart ) ) {
                                nrj_lost_per_thd[tid] += mass * ener_iPart;
                            }
                        }
                    }
                    
                    // Boundary Condition may be physical or due to domain decomposition
                    // apply returns 0 if iPart is not in the local domain anymore
                    for( iPart=first_index[ipack*packsize_+scell] ; ( int )iPart<last_index[ipack*packsize_+scell]; iPart++ ) {
                        if( !partBoundCond->apply( *particles, iPart, this, ener_iPart ) ) {
                            addPartInExchList( iPart );
                            nrj_lost_per_thd[tid] += mass * ener_iPart;
                            particles->cell_keys[iPart] = -1;
                        } else {
                            //First reduction of the count sort algorithm. Lost particles are not included.
                            for( int i = 0 ; i<( int )nDim_particle; i++ ) {
                                particles->cell_keys[iPart] *= length[i];
                                particles->cell_keys[iPart] += round( ( particles->position( i, iPart )-min_loc_vec[i] ) * dx_inv_[i] );
                            }
                            count[particles->cell_keys[iPart]] ++; //First reduction of the count sort algorithm. Lost particles are not included.
                        }
                    }
                } else if( mass==0 ) { // condition mass=0
                    ERROR( "Particles with zero mass cannot interact with envelope" );
                }
            }
            //START EXCHANGE PARTICLES OF THE CURRENT BIN ?
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[3] += MPI_Wtime() - timer;
#endif
            
            // Project currents if not a Test species and charges as well if a diag is needed.
            // Do not project if a photon
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            if( ( !particles->is_test ) && ( mass > 0 ) )
                for( unsigned int scell = 0 ; scell < packsize_ ; scell++ ) {
                    Proj->currentsAndDensityWrapper( EMfields, *particles, smpi, first_index[ipack*packsize_+scell], last_index[ipack*packsize_+scell], ithread, diag_flag, params.is_spectral, ispec, ipack*packsize_+scell, first_index[ipack*packsize_] );
                }
                
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[12] += MPI_Wtime() - timer;
#endif
        }
        
        for( unsigned int ithd=0 ; ithd<nrj_lost_per_thd.size() ; ithd++ ) {
            nrj_bc_lost += nrj_lost_per_thd[tid];
        }
        
    } else { // immobile particle (at the moment only project density)
        if( diag_flag &&( !particles->is_test ) ) {
            double *b_rho=nullptr;
            for( unsigned int scell = 0 ; scell < first_index.size() ; scell ++ ) {
            
                if( nDim_field==2 ) {
                    b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
                }
                if( nDim_field==3 ) {
                    b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
                } else if( nDim_field==1 ) {
                    b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
                }
                for( iPart=first_index[scell] ; ( int )iPart<last_index[scell]; iPart++ ) {
                    Proj->basic( b_rho, ( *particles ), iPart, 0 );
                } //End loop on particles
            }
        }
    }//END if time vs. time_frozen
    
} // end ponderomotive_update_position_and_currents
