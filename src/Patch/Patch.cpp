
// CLOBAL COORDINATES:
//                           Patch_minGlobal                                                                      Patch_maxGlobal
//                      --------<===================================== gs ===================================>------------
//     GLOBAL INDICES:          0                                  .                                        nspace_global
//                           ix+oversize                                                                  ix+oversize
//                      ------------------------------------       .              ------------------------------------
//                      |   |   |     ...          |   |   |       .              |   |   |   |   ...    |   |   |   |
//                      |   |   |     ...          |   |   |       .              |   |   |   |   ...    |   |   |   |
//                      ------------------------------------       .              ------------------------------------
//                          Patch_minLocal    Patch_maxLocal       .             Patch_minLocal        Patch_maxLocal
//                                                 ----------------------------------------
//                                                 |   |   |       .              |   |   |
//                                                 |   |   |       .              |   |   |
//                                                 ----------------------------------------
// LOCAL COORDINATES:                             x(0) rlb        x(ix)             rub  x(nspace)
//                                                 ----<============= length =========>----
//     LOCAL INDICES:                              0   lb                            ub   nspace

#include "Patch.h"

#include <iostream>
#include <iomanip>

#include "DomainDecompositionFactory.h"
#include "Hilbert_functions.h"
#include "PatchesFactory.h"
#include "SpeciesFactory.h"
#include "ParticleInjectorFactory.h"
#include "Particles.h"
#include "ElectroMagnFactory.h"
#include "ElectroMagnBC_Factory.h"
#include "DiagnosticFactory.h"
#include "CollisionsFactory.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Patch constructor :
//   Called by PatchXD constructor which will finalize initialization
// ---------------------------------------------------------------------------------------------------------------------
Patch::Patch( Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch, unsigned int n_moved )
{

    hindex = ipatch;
    nDim_fields_ = params.nDim_field;

    initStep1( params );

#ifdef  __DETAILED_TIMERS
    // Initialize timers
    // 0 - Interpolation
    // 1 - Pusher
    // 2 - Projection
    // 3 - exchange init + cell_keys
    // 4 - ionization
    // 5 - radiation
    // 6 - Breit-Wheeler
    // 7 - Interp Fields_Env
    // 8 - Proj Susceptibility
    // 9 - Push Momentum
    // 10 - Interp Env_Old
    // 11 - Proj Currents
    // 12 - Push Pos
    // 13 - Sorting
    patch_timers.resize( 15, 0. );
#endif

} // END Patch::Patch




// Cloning patch constructor
Patch::Patch( Patch *patch, Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch, unsigned int n_moved, bool with_particles = true )
{

    hindex = ipatch;
    nDim_fields_ = patch->nDim_fields_;

    initStep1( params );

#ifdef  __DETAILED_TIMERS
    // Initialize timers
    patch_timers.resize( 15, 0. );
#endif

}

void Patch::initStep1( Params &params )
{
    // for nDim_fields = 1 : bug if Pcoordinates.size = 1 !!
    //Pcoordinates.resize(nDim_fields_);
    Pcoordinates.resize( 2 );

    nbNeighbors_ = 2;
    neighbor_.resize( nDim_fields_ );
    tmp_neighbor_.resize( nDim_fields_ );
    for( int iDim = 0 ; iDim < nDim_fields_ ; iDim++ ) {
        neighbor_[iDim].resize( 2, MPI_PROC_NULL );
        tmp_neighbor_[iDim].resize( 2, MPI_PROC_NULL );
    }
    MPI_neighbor_.resize( nDim_fields_ );
    tmp_MPI_neighbor_.resize( nDim_fields_ );
    for( int iDim = 0 ; iDim < nDim_fields_; iDim++ ) {
        MPI_neighbor_[iDim].resize( 2, MPI_PROC_NULL );
        tmp_MPI_neighbor_[iDim].resize( 2, MPI_PROC_NULL );
    }

    oversize.resize( nDim_fields_ );
    for( int iDim = 0 ; iDim < nDim_fields_; iDim++ ) {
        oversize[iDim] = params.oversize[iDim];
    }
    
    // Initialize the random number generator
    rand_ = new Random( params.random_seed );
    
    // Obtain the cell_volume
    cell_volume = params.cell_volume;
}


void Patch::initStep3( Params &params, SmileiMPI *smpi, unsigned int n_moved )
{
    // Compute MPI neighborood
    updateMPIenv( smpi );

    // Compute patch boundaries
    min_local.resize( params.nDim_field, 0. );
    max_local.resize( params.nDim_field, 0. );
    center   .resize( params.nDim_field, 0. );
    cell_starting_global_index.resize( params.nDim_field, 0 );
    radius = 0.;
    for( unsigned int i = 0 ; i<params.nDim_field ; i++ ) {
        min_local[i] = ( Pcoordinates[i] )*( params.n_space[i]*params.cell_length[i] );
        max_local[i] = ( ( Pcoordinates[i]+1 ) )*( params.n_space[i]*params.cell_length[i] );
        cell_starting_global_index[i] += Pcoordinates[i]*params.n_space[i];
        cell_starting_global_index[i] -= params.oversize[i];
        center[i] = ( min_local[i]+max_local[i] )*0.5;
        radius += pow( max_local[i] - center[i] + params.cell_length[i], 2 );
    }
    radius = sqrt( radius );

    cell_starting_global_index[0] += n_moved;
    min_local[0] += n_moved*params.cell_length[0];
    max_local[0] += n_moved*params.cell_length[0];
    center   [0] += n_moved*params.cell_length[0];

    //Shift point position by dr/2 for the AM spectral geometry
    //if ( (params.is_spectral) && (params.geometry== "AMcylindrical") ) {
    //    min_local[1] += params.cell_length[1]/2.;
    //    max_local[1] += params.cell_length[1]/2.;
    //    center   [1] += params.cell_length[1]/2.;
    //}
    
}

void Patch::finishCreation( Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition )
{
    // initialize vector of Species (virtual)
    vecSpecies = SpeciesFactory::createVector( params, this );

    // initialize the electromagnetic fields (virtual)
    EMfields   = ElectroMagnFactory::create( params, domain_decomposition, vecSpecies, this );

    // Initialize the collisions
    vecCollisions = CollisionsFactory::create( params, this, vecSpecies );

    // Initialize the particle injector
    particle_injector_vector_ = ParticleInjectorFactory::createVector( params, this, vecSpecies );

    // Initialize the particle walls
    partWalls = new PartWalls( params, this );

    // Initialize the probes
    probes = DiagnosticFactory::createProbes();

    probesInterp = InterpolatorFactory::create( params, this, false );

    if( has_an_MPI_neighbor() ) {
        createType( params );
    }

}


void Patch::finishCloning( Patch *patch, Params &params, SmileiMPI *smpi, unsigned int n_moved, bool with_particles = true )
{
    // clone vector of Species (virtual)
    vecSpecies = SpeciesFactory::cloneVector( patch->vecSpecies, params, this, with_particles );

    // clone the electromagnetic fields (virtual)
    EMfields   = ElectroMagnFactory::clone( patch->EMfields, params, vecSpecies, this, n_moved );

    // clone the collisions
    vecCollisions = CollisionsFactory::clone( patch->vecCollisions );

    // Clone the particle injector
    particle_injector_vector_ = ParticleInjectorFactory::cloneVector( patch->particle_injector_vector_, params, patch);

    // clone the particle walls
    partWalls = new PartWalls( patch->partWalls, this );

    // clone the probes
    probes = DiagnosticFactory::cloneProbes( patch->probes );

    probesInterp = InterpolatorFactory::create( params, this, false );

    if( has_an_MPI_neighbor() ) {
        createType( params );
    }

}

void Patch::finalizeMPIenvironment( Params &params )
{
    int nb_comms( 9 ); // E, B, B_m : min number of comms

    if( params.geometry == "AMcylindrical" ) {
        nb_comms += 9*( params.nmodes - 1 );
    }
    // if envelope is present,
    // add to comms A, A0, Phi, Phi_old, GradPhi (x,y,z components), GradPhi_old (x,y,z components)
    if( params.Laser_Envelope_model ) {
        nb_comms += 10;
    }

    // add comms for species
    nb_comms += 2*vecSpecies.size();

    // Adaptive vectorization:
    if( params.has_adaptive_vectorization ) {
        nb_comms ++;
    }

    // Scalars
    nb_comms ++;

    // Just apply on species & fields to start

    for( unsigned int idiag=0; idiag<EMfields->allFields_avg.size(); idiag++ ) {
        nb_comms += EMfields->allFields_avg[idiag].size();
    }
    nb_comms += EMfields->antennas.size();

    for( unsigned int bcId=0 ; bcId<EMfields->emBoundCond.size() ; bcId++ ) {
        if( EMfields->emBoundCond[bcId] ) {
            for( unsigned int laserId=0 ; laserId < EMfields->emBoundCond[bcId]->vecLaser.size() ; laserId++ ) {
                nb_comms += 4;
            }
        }
        if( EMfields->extFields.size()>0 ) {
            if( dynamic_cast<ElectroMagnBC1D_SM *>( EMfields->emBoundCond[bcId] ) ) {
                nb_comms += 4;
            } else if( dynamic_cast<ElectroMagnBC2D_SM *>( EMfields->emBoundCond[bcId] ) ) {
                nb_comms += 12;
            } else if( dynamic_cast<ElectroMagnBC3D_SM *>( EMfields->emBoundCond[bcId] ) ) {
                nb_comms += 18;
            }
        }
    }
    requests_.resize( nb_comms, MPI_REQUEST_NULL );

}


void Patch::setLocationAndAllocateFields( Params &params, DomainDecomposition *domain_decomposition, VectorPatch &vecPatch )
{
    for( int iDim = 0 ; iDim < nDim_fields_; iDim++ ) {
        oversize[iDim] = params.region_oversize[iDim];
    }

    Pcoordinates.resize( params.nDim_field );


    min_local = vecPatch( 0 )->min_local;
    max_local = vecPatch( 0 )->max_local;
    center   .resize( nDim_fields_, 0. );
    cell_starting_global_index = vecPatch( 0 )->cell_starting_global_index;
    radius = 0.;

    // Constraint to enforce 1 neighboor per side
    //double nppp_root = pow( vecPatch.size(), 1./(double)nDim_fields_ );
    //if ( fabs( (double)(int)nppp_root - nppp_root ) > 0. ) 
    //    ERROR( "Bad choice of decomposition" );

// Coupling v0    for (int i = 0 ; i<nDim_fields_ ; i++) {
// Coupling v0
// Coupling v0        for ( unsigned int ipatch = 0 ; ipatch < vecPatch.size() ; ipatch++  ) {
// Coupling v0            if ( vecPatch(ipatch)->min_local[i] <= min_local[i] ) {
// Coupling v0                min_local[i] = vecPatch(ipatch)->min_local[i];
// Coupling v0                if (vecPatch(ipatch)->MPI_neighbor_[i][0]!=MPI_PROC_NULL)
// Coupling v0                    MPI_neighbor_[i][0] = vecPatch(ipatch)->MPI_neighbor_[i][0];
// Coupling v0                if (vecPatch(ipatch)->neighbor_[i][0]!=MPI_PROC_NULL)
// Coupling v0                    neighbor_[i][0] = (vecPatch(ipatch)->neighbor_[i][0] / vecPatch.size() );
// Coupling v0            }
// Coupling v0            if( vecPatch( ipatch )->max_local[i] >= max_local[i] ) {
// Coupling v0                max_local[i] = vecPatch( ipatch )->max_local[i];
// Coupling v0                if( vecPatch( ipatch )->MPI_neighbor_[i][1]!=MPI_PROC_NULL ) {
// Coupling v0                    MPI_neighbor_[i][1] = vecPatch( ipatch )->MPI_neighbor_[i][1];
// Coupling v0                }
// Coupling v0                if( vecPatch( ipatch )->neighbor_[i][1]!=MPI_PROC_NULL ) {
// Coupling v0                    neighbor_[i][1] = ( vecPatch( ipatch )->neighbor_[i][1] / vecPatch.size() );
// Coupling v0                }
// Coupling v0            }
// Coupling v0            if( vecPatch( ipatch )->cell_starting_global_index[i] <= cell_starting_global_index[i] ) {
// Coupling v0                cell_starting_global_index[i] = vecPatch( ipatch )->cell_starting_global_index[i];
// Coupling v0            }
// Coupling v0        }
// Coupling v0        
// Coupling v0        center[i] = (min_local[i]+max_local[i])*0.5;
// Coupling v0        radius += pow(max_local[i] - center[i] + params.cell_length[i], 2);
// Coupling v0    }

    // New_DD 
    int rk(0);
    MPI_Comm_rank( MPI_COMM_WORLD, &rk );
    int sz(1);
    MPI_Comm_size( MPI_COMM_WORLD, &sz );

    // If current patch is a Domain's patch
    if ( dynamic_cast<RegionDomainDecomposition*>( domain_decomposition ) ) {

        if (nDim_fields_==1) {
            for ( unsigned int xDom = 0 ; xDom < params.number_of_region[0] ; xDom++ )
                for ( unsigned int yDom = 0 ; yDom < params.number_of_region[1] ; yDom++ ) {
                    for ( unsigned int zDom = 0 ; zDom < params.number_of_region[2] ; zDom++ ) {

                        if (params.map_rank[xDom][yDom][zDom] == rk ) {
                            Pcoordinates[0] = xDom;
                            Pcoordinates[1] = yDom;
                            min_local[0] =  params.offset_map[0][xDom]                           * params.cell_length[0];
                            max_local[0] = (params.offset_map[0][xDom]+params.n_space_region[0]) * params.cell_length[0];

                            cell_starting_global_index[0] = params.offset_map[0][xDom];

                            if (xDom>0)
                                MPI_neighbor_[0][0] = params.map_rank[xDom-1][yDom][zDom];
                            else if (params.EM_BCs[0][0]=="periodic") {
                                MPI_neighbor_[0][0] = params.map_rank[xDom-1+params.number_of_region[0]][yDom][zDom];
                            }
                            else 
                                MPI_neighbor_[0][0] = MPI_PROC_NULL;

                            if (xDom<params.number_of_region[0]-1)
                                MPI_neighbor_[0][1] = params.map_rank[xDom+1][yDom][zDom];
                            else if (params.EM_BCs[0][1]=="periodic")
                                MPI_neighbor_[0][1] = params.map_rank[(xDom+1)%params.number_of_region[0]][yDom][zDom];
                            else 
                                MPI_neighbor_[0][1] = MPI_PROC_NULL;
                        }                
                    }
                }

            neighbor_[0][0] = MPI_neighbor_[0][0];
            neighbor_[0][1] = MPI_neighbor_[0][1];

            cell_starting_global_index[0] -= oversize[0];
        } // Fin 1D
    
        if (nDim_fields_==2) {

            for ( unsigned int xDom = 0 ; xDom < params.number_of_region[0] ; xDom++ )
                for ( unsigned int yDom = 0 ; yDom < params.number_of_region[1] ; yDom++ ) {
                    for ( unsigned int zDom = 0 ; zDom < params.number_of_region[2] ; zDom++ ) {
                        
                        if (params.map_rank[xDom][yDom][zDom] == rk ) {
                            
           
                            Pcoordinates[0] =xDom ;
                            Pcoordinates[1] =yDom;

                            //cout << "coords = " << Pcoordinates[0] << " " << Pcoordinates[1] <<endl;

                            min_local[0] =  params.offset_map[0][xDom]                           * params.cell_length[0];
                            max_local[0] = (params.offset_map[0][xDom]+params.n_space_region[0]) * params.cell_length[0];
                            min_local[1] =  params.offset_map[1][yDom]                           * params.cell_length[1];
                            max_local[1] = (params.offset_map[1][yDom]+params.n_space_region[1]) * params.cell_length[1];

                            center[0] = ( min_local[0]+max_local[0] )*0.5;
                            radius += pow( max_local[0] - center[0] + params.cell_length[0], 2 );
                            center[1] = ( min_local[1]+max_local[1] )*0.5;
                            radius += pow( max_local[1] - center[1] + params.cell_length[1], 2 );

                            //Shift point position by dr/2 for the AM spectral geometry
                            //if ( (params.is_spectral) && (params.geometry== "AMcylindrical") ) {
                            //    min_local[1] += params.cell_length[1]/2.;
                            //    max_local[1] += params.cell_length[1]/2.;
                            //    center   [1] += params.cell_length[1]/2.;
                            //}
 
                            cell_starting_global_index[0] = params.offset_map[0][xDom];
                            cell_starting_global_index[1] = params.offset_map[1][yDom];

                            // ---------------- X ----------------
                            if (xDom>0)
                                MPI_neighbor_[0][0] = params.map_rank[xDom-1][yDom][zDom];
                            else if (params.EM_BCs[0][0]=="periodic") {
                                MPI_neighbor_[0][0] = params.map_rank[xDom-1+params.number_of_region[0]][yDom][zDom];
                            }
                            else 
                                MPI_neighbor_[0][0] = MPI_PROC_NULL;

                            if (xDom<params.number_of_region[0]-1)
                                MPI_neighbor_[0][1] = params.map_rank[xDom+1][yDom][zDom];
                            else if (params.EM_BCs[0][1]=="periodic")
                                MPI_neighbor_[0][1] = params.map_rank[(xDom+1)%params.number_of_region[0]][yDom][zDom];
                            else 
                                MPI_neighbor_[0][1] = MPI_PROC_NULL;

                            // ---------------- Y ----------------
                            if (yDom>0)
                                MPI_neighbor_[1][0] = params.map_rank[xDom][yDom-1][zDom];
                            else if (params.EM_BCs[1][0]=="periodic") {
                                MPI_neighbor_[1][0] = params.map_rank[xDom][yDom-1+params.number_of_region[1]][zDom];
                            }
                            else 
                                MPI_neighbor_[1][0] = MPI_PROC_NULL;

                            if (yDom<params.number_of_region[1]-1)
                                MPI_neighbor_[1][1] = params.map_rank[xDom][yDom+1][zDom];
                            else if (params.EM_BCs[1][1]=="periodic")
                                MPI_neighbor_[1][1] = params.map_rank[xDom][(yDom+1)%params.number_of_region[1]][zDom];
                            else 
                                MPI_neighbor_[1][1] = MPI_PROC_NULL;
                        }                
                    }
                }

            //cout << "HERE - " <<Pcoordinates[0] << " " << Pcoordinates[1]<< endl; 
            std::vector<int> xcall( 2, 0 );
            // 1st direction
            xcall[0] = Pcoordinates[0]-1;
            xcall[1] = Pcoordinates[1];
            if( params.EM_BCs[0][0]=="periodic" && xcall[0] < 0 ) {
                xcall[0] += domain_decomposition->ndomain_[0];
            }
            neighbor_[0][0] = domain_decomposition->getDomainId( xcall );
    
            xcall[0] = Pcoordinates[0]+1;
            if( params.EM_BCs[0][0]=="periodic" && xcall[0] >= (int)params.number_of_region[0] ) {
                xcall[0] -= params.number_of_region[0];
            }
            neighbor_[0][1] = domain_decomposition->getDomainId( xcall );
    
            // 2nd direction
            xcall[0] = Pcoordinates[0];
            xcall[1] = Pcoordinates[1]-1;
            if( params.EM_BCs[1][0]=="periodic" && xcall[1] < 0 ) {
                xcall[1] += domain_decomposition->ndomain_[1];
            }
            neighbor_[1][0] = domain_decomposition->getDomainId( xcall );
    
            xcall[1] = Pcoordinates[1]+1;
            if( params.EM_BCs[1][0]=="periodic" && xcall[1] >= (int)params.number_of_region[1] ) {
                xcall[1] -=  params.number_of_region[1];
            }
            neighbor_[1][1] = domain_decomposition->getDomainId( xcall );

            //cout << "\t"<< neighbor_[1][1] << endl;
            //cout << neighbor_[0][0] << "\t h_me \t" << neighbor_[0][1] << endl;
            //cout << "\t"<< neighbor_[1][0] << endl;
            //
            //cout << "\t"<< MPI_neighbor_[1][1] << endl;
            //cout << MPI_neighbor_[0][0] << "\t mpi_me \t" << MPI_neighbor_[0][1] << endl;
            //cout << "\t"<< MPI_neighbor_[1][0] << endl;

            cell_starting_global_index[0] -= oversize[0];
            cell_starting_global_index[1] -= oversize[1];
        } // Fin 2D


        if (nDim_fields_==3) {

            for ( unsigned int xDom = 0 ; xDom < params.number_of_region[0] ; xDom++ )
                for ( unsigned int yDom = 0 ; yDom < params.number_of_region[1] ; yDom++ ) {
                    for ( unsigned int zDom = 0 ; zDom < params.number_of_region[2] ; zDom++ ) {

                        if (params.map_rank[xDom][yDom][zDom] == rk ) {
                            Pcoordinates[0] = xDom;
                            Pcoordinates[1] = yDom;
                            Pcoordinates[2] = zDom;
                            //cout << hindex << " - coords = " << xDom << " " << yDom << " " << zDom << endl;

                            min_local[0] =  params.offset_map[0][xDom]                           * params.cell_length[0];
                            max_local[0] = (params.offset_map[0][xDom]+params.n_space_region[0]) * params.cell_length[0];
                            min_local[1] =  params.offset_map[1][yDom]                           * params.cell_length[1];
                            max_local[1] = (params.offset_map[1][yDom]+params.n_space_region[1]) * params.cell_length[1];
                            min_local[2] =  params.offset_map[2][zDom]                           * params.cell_length[2];
                            max_local[2] = (params.offset_map[2][zDom]+params.n_space_region[2]) * params.cell_length[2];

                            center[0] = ( min_local[0]+max_local[0] )*0.5;
                            radius += pow( max_local[0] - center[0] + params.cell_length[0], 2 );
                            center[1] = ( min_local[1]+max_local[1] )*0.5;
                            radius += pow( max_local[1] - center[1] + params.cell_length[1], 2 );
                            center[2] = ( min_local[2]+max_local[2] )*0.5;
                            radius += pow( max_local[2] - center[2] + params.cell_length[2], 2 );


                            cell_starting_global_index[0] = params.offset_map[0][xDom];
                            cell_starting_global_index[1] = params.offset_map[1][yDom];
                            cell_starting_global_index[2] = params.offset_map[2][zDom];

                            // ---------------- X ----------------
                            if (xDom>0)
                                MPI_neighbor_[0][0] = params.map_rank[xDom-1][yDom][zDom];
                            else if (params.EM_BCs[0][0]=="periodic") {
                                MPI_neighbor_[0][0] = params.map_rank[xDom-1+params.number_of_region[0]][yDom][zDom];
                            }
                            else 
                                MPI_neighbor_[0][0] = MPI_PROC_NULL;

                            if (xDom<params.number_of_region[0]-1)
                                MPI_neighbor_[0][1] = params.map_rank[xDom+1][yDom][zDom];
                            else if (params.EM_BCs[0][1]=="periodic")
                                MPI_neighbor_[0][1] = params.map_rank[(xDom+1)%params.number_of_region[0]][yDom][zDom];
                            else 
                                MPI_neighbor_[0][1] = MPI_PROC_NULL;


                            //cout << MPI_neighbor_[0][0] << " " << " me "  << " "  << MPI_neighbor_[0][1] << endl;
                               
                            // ---------------- Y ----------------
                            if (yDom>0)
                                MPI_neighbor_[1][0] = params.map_rank[xDom][yDom-1][zDom];
                            else if (params.EM_BCs[1][0]=="periodic") {
                                MPI_neighbor_[1][0] = params.map_rank[xDom][yDom-1+params.number_of_region[1]][zDom];
                            }
                            else 
                                MPI_neighbor_[1][0] = MPI_PROC_NULL;

                            if (yDom<params.number_of_region[1]-1)
                                MPI_neighbor_[1][1] = params.map_rank[xDom][yDom+1][zDom];
                            else if (params.EM_BCs[1][1]=="periodic")
                                MPI_neighbor_[1][1] = params.map_rank[xDom][(yDom+1)%params.number_of_region[1]][zDom];
                            else 
                                MPI_neighbor_[1][1] = MPI_PROC_NULL;

                            // ---------------- Z ----------------
                            if (zDom>0)
                                MPI_neighbor_[2][0] = params.map_rank[xDom][yDom][zDom-1];
                            else if (params.EM_BCs[2][0]=="periodic") {
                                MPI_neighbor_[2][0] = params.map_rank[xDom][yDom][zDom-1+params.number_of_region[2]];
                            }
                            else 
                                MPI_neighbor_[2][0] = MPI_PROC_NULL;

                            if (zDom<params.number_of_region[2]-1)
                                MPI_neighbor_[2][1] = params.map_rank[xDom][yDom][zDom+1];
                            else if (params.EM_BCs[2][1]=="periodic")
                                MPI_neighbor_[2][1] = params.map_rank[xDom][yDom][(zDom+1)%params.number_of_region[2]];
                            else 
                                MPI_neighbor_[2][1] = MPI_PROC_NULL;

                        }                
                    }
                }

            //neighbor_[0][0] = MPI_neighbor_[0][0];
            //neighbor_[0][1] = MPI_neighbor_[0][1];
            //neighbor_[1][0] = MPI_neighbor_[1][0];
            //neighbor_[1][1] = MPI_neighbor_[1][1];
            //neighbor_[2][0] = MPI_neighbor_[2][0];
            //neighbor_[2][1] = MPI_neighbor_[2][1];
            std::vector<int> xcall( 3, 0 );

            // 1st direction
            xcall[0] = Pcoordinates[0]-1;
            xcall[1] = Pcoordinates[1];
            xcall[2] = Pcoordinates[2];
            if( params.EM_BCs[0][0]=="periodic" && xcall[0] < 0 ) {
                xcall[0] += domain_decomposition->ndomain_[0];
            }
            neighbor_[0][0] = domain_decomposition->getDomainId( xcall );
            xcall[0] = Pcoordinates[0]+1;
            if( params.EM_BCs[0][0]=="periodic" && xcall[0] >= ( int )domain_decomposition->ndomain_[0] ) {
                xcall[0] -= domain_decomposition->ndomain_[0];
            }
            neighbor_[0][1] = domain_decomposition->getDomainId( xcall );
    
            // 2st direction
            xcall[0] = Pcoordinates[0];
            xcall[1] = Pcoordinates[1]-1;
            xcall[2] = Pcoordinates[2];
            if( params.EM_BCs[1][0]=="periodic" && xcall[1] < 0 ) {
                xcall[1] += domain_decomposition->ndomain_[1];
            }
            neighbor_[1][0] =  domain_decomposition->getDomainId( xcall );
            xcall[1] = Pcoordinates[1]+1;
            if( params.EM_BCs[1][0]=="periodic" && xcall[1] >= ( int )domain_decomposition->ndomain_[1] ) {
                xcall[1] -= domain_decomposition->ndomain_[1];
            }
            neighbor_[1][1] =  domain_decomposition->getDomainId( xcall );
    
            // 3st direction
            xcall[0] = Pcoordinates[0];
            xcall[1] = Pcoordinates[1];
            xcall[2] = Pcoordinates[2]-1;
            if( params.EM_BCs[2][0]=="periodic" && xcall[2] < 0 ) {
                xcall[2] += domain_decomposition->ndomain_[2];
            }
            neighbor_[2][0] =  domain_decomposition->getDomainId( xcall );
            xcall[2] = Pcoordinates[2]+1;
            if( params.EM_BCs[2][0]=="periodic" && xcall[2] >= ( int )domain_decomposition->ndomain_[2] ) {
                xcall[2] -= domain_decomposition->ndomain_[2];
            }
            neighbor_[2][1] =  domain_decomposition->getDomainId( xcall );




            cell_starting_global_index[0] -= oversize[0];
            cell_starting_global_index[1] -= oversize[1];
            cell_starting_global_index[2] -= oversize[2];
        } // Fin 3D

    } // ENd if 
    else if ( domain_decomposition == NULL ) { // If current patch is the global Domain's patch
        if ( params.geometry != "AMcylindrical" ) {
            WARNING ("Global gathering not tested on non AM configuration" ) ;
        }
        for ( int iDim=0 ; iDim<nDim_fields_ ; iDim++ ) {
            Pcoordinates[iDim] = 0;
            min_local[iDim] = 0.;
            max_local[iDim] = params.n_space_global[iDim]*params.cell_length[iDim];

            center[iDim] = ( min_local[iDim]+max_local[iDim] )*0.5;
            radius += pow( max_local[iDim] - center[iDim] + params.cell_length[iDim], 2 );

            cell_starting_global_index[iDim] = -oversize[iDim];

            if (params.EM_BCs[iDim][0]=="periodic") {
                MPI_neighbor_[iDim][0] = rk;
            }
            else
                MPI_neighbor_[iDim][0] = MPI_PROC_NULL;
            if (params.EM_BCs[iDim][1]=="periodic") {
                MPI_neighbor_[iDim][1] = rk;
            }
            else
                MPI_neighbor_[iDim][1] = MPI_PROC_NULL;

        }
        //if ( (params.is_spectral) && (params.geometry== "AMcylindrical") ) {
        //    min_local[1] += params.cell_length[1]/2.;
        //    max_local[1] += params.cell_length[1]/2.;
        //    center   [1] += params.cell_length[1]/2.;
        //}

    }
    else {
        ERROR( "Should not pass here" );
    }

    radius = sqrt(radius);
    

    //cart_updateMPIenv(smpi);

    MPI_me_ = vecPatch( 0 )->MPI_me_;
    /*for( int i = 0 ; i<nDim_fields_ ; i++ ) {
        if( ( MPI_neighbor_[i][0]==MPI_me_ ) && ( params.EM_BCs[i][0]!="periodic" ) ) {
            MPI_neighbor_[i][0] = MPI_PROC_NULL;
        }
        if( ( MPI_neighbor_[i][1]==MPI_me_ ) && ( params.EM_BCs[i][0]!="periodic" ) ) {
            MPI_neighbor_[i][1] = MPI_PROC_NULL;
        }
        if( ( neighbor_[i][0]==MPI_me_ ) && ( params.EM_BCs[i][0]!="periodic" ) ) {
            neighbor_[i][0] = MPI_PROC_NULL;
        }
        if( ( neighbor_[i][1]==MPI_me_ ) && ( params.EM_BCs[i][0]!="periodic" ) ) {
            neighbor_[i][1] = MPI_PROC_NULL;
        }
    }*/
    
    
    //cout << "MPI Nei\t"  << "\t" << MPI_neighbor_[1][1] << endl;
    //cout << "MPI Nei\t"  << MPI_neighbor_[0][0] << "\t" << MPI_me_ << "\t" << MPI_neighbor_[0][1] << endl;
    //cout << "MPI Nei\t"  << "\t" << MPI_neighbor_[1][0] << endl;
    
    EMfields   = ElectroMagnFactory::create( params, domain_decomposition, vecPatch( 0 )->vecSpecies, this );

    vecSpecies.resize( 0 );
    vecCollisions.resize( 0 );
    partWalls = NULL;
    probes.resize( 0 );
    probesInterp = NULL;
    
    if( has_an_MPI_neighbor() ) {
        createType2( params );
    }

    is_small= false;

}


// ---------------------------------------------------------------------------------------------------------------------
// Delete Patch members
// ---------------------------------------------------------------------------------------------------------------------
Patch::~Patch()
{

    if (probesInterp) delete probesInterp;
    
    for( unsigned int i=0; i<probes.size(); i++ ) {
        delete probes[i];
    }

    for( unsigned int i=0; i<vecCollisions.size(); i++ ) {
        delete vecCollisions[i];
    }
    vecCollisions.clear();

    if( partWalls!=NULL ) {
        delete partWalls;
    }

    if( EMfields !=NULL ) {
        delete EMfields;
    }

    for( unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++ ) {
        delete vecSpecies[ispec];
    }
    vecSpecies.clear();
    
    delete rand_;
    
} // END Patch::~Patch


// ---------------------------------------------------------------------------------------------------------------------
// Compute MPI rank of patch neigbors and current patch
// ---------------------------------------------------------------------------------------------------------------------
void Patch::updateTagenv( SmileiMPI *smpi )
{
}
void Patch::updateMPIenv( SmileiMPI *smpi )
{
    // HARDCODED VALUE - to test cartesian decomposition
//    vector<int> npattchpp(2,1); // NDOMAIN PER PROC - WORKS WITH 1 domain per process, ATTENTION PERIDODICITY
//    vector<int> npatchs(2,4);   // NDOMAIN TOTAL
//    vector<int> gCoord(2,0);
    MPI_me_ = smpi->smilei_rk;

    for( int iDim = 0 ; iDim < nDim_fields_ ; iDim++ )
        for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
            MPI_neighbor_[iDim][iNeighbor] = smpi->hrank( neighbor_[iDim][iNeighbor] );

//            gCoord[0] = ( (int)Pcoordinates[0] + (1-iDim) * ((1-iNeighbor)*(-1) + (iNeighbor)) );
//            gCoord[1] = ( (int)Pcoordinates[1] + (iDim  ) * ((1-iNeighbor)*(-1) + (iNeighbor)) );
//            if ( gCoord[0] < 0 )
//                MPI_neighbor_[iDim][iNeighbor] = MPI_PROC_NULL;
//            else if ( gCoord[0] >= npatchs[0] )
//                MPI_neighbor_[iDim][iNeighbor] = MPI_PROC_NULL;
//            else if ( gCoord[1] < 0 )
//                MPI_neighbor_[iDim][iNeighbor] = MPI_PROC_NULL;
//            else if ( gCoord[1] >= npatchs[1] )
//                MPI_neighbor_[iDim][iNeighbor] = MPI_PROC_NULL;
//            else {
//                gCoord[0] = (double)( (int)Pcoordinates[0] + (1-iDim) * ((1-iNeighbor)*(-1) + (iNeighbor)) ) /(double)npattchpp[0];
//                gCoord[1] = (double)( (int)Pcoordinates[1] + (iDim  ) * ((1-iNeighbor)*(-1) + (iNeighbor)) ) /(double)npattchpp[1];
//                MPI_neighbor_[iDim][iNeighbor] = gCoord[0] * npatchs[1] + gCoord[1];
//            }
        }

} // END updateMPIenv

// ---------------------------------------------------------------------------------------------------------------------
// Clean the MPI buffers for communications
// ---------------------------------------------------------------------------------------------------------------------
void Patch::cleanMPIBuffers( int ispec, Params &params )
{
    int ndim = params.nDim_field;

    for( int iDim=0 ; iDim < ndim ; iDim++ ) {
        for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
            vecSpecies[ispec]->MPI_buffer_.partRecv[iDim][iNeighbor].clear();//resize(0,ndim);
            vecSpecies[ispec]->MPI_buffer_.partSend[iDim][iNeighbor].clear();//resize(0,ndim);
            vecSpecies[ispec]->MPI_buffer_.part_index_send[iDim][iNeighbor].clear();
            //vecSpecies[ispec]->MPI_buffer_.part_index_send[iDim][iNeighbor].resize(0);
            vecSpecies[ispec]->MPI_buffer_.part_index_recv_sz[iDim][iNeighbor] = 0;
        }
    }
    vecSpecies[ispec]->indexes_of_particles_to_exchange.clear();
} // cleanMPIBuffers


// ---------------------------------------------------------------------------------------------------------------------
// Split particles Id to send in per direction and per patch neighbor dedicated buffers
// Apply periodicity if necessary
// ---------------------------------------------------------------------------------------------------------------------
void Patch::initExchParticles( SmileiMPI *smpi, int ispec, Params &params )
{
    Particles &cuParticles = ( *vecSpecies[ispec]->particles );
    int ndim = params.nDim_field;
    int idim, check;
    std::vector<int> *indexes_of_particles_to_exchange = &vecSpecies[ispec]->indexes_of_particles_to_exchange;
//    double xmax[3];

    for( int iDim=0 ; iDim < ndim ; iDim++ ) {
        for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
            vecSpecies[ispec]->MPI_buffer_.partRecv[iDim][iNeighbor].clear();//resize(0,ndim);
            vecSpecies[ispec]->MPI_buffer_.partSend[iDim][iNeighbor].clear();//resize(0,ndim);
            vecSpecies[ispec]->MPI_buffer_.part_index_send[iDim][iNeighbor].resize( 0 );
            vecSpecies[ispec]->MPI_buffer_.part_index_recv_sz[iDim][iNeighbor] = 0;
        }
    }

    int n_part_send = indexes_of_particles_to_exchange->size();

    int iPart;

    // Define where particles are going
    //Put particles in the send buffer it belongs to. Priority to lower dimensions.
    if( params.geometry != "AMcylindrical" ) {
        for( int i=0 ; i<n_part_send ; i++ ) {
            iPart = ( *indexes_of_particles_to_exchange )[i];
            check = 0;
            idim = 0;
            //Put indexes of particles in the first direction they will be exchanged and correct their position according to periodicity for the first exchange only.
            while( check == 0 && idim<ndim ) {
                if( cuParticles.position( idim, iPart ) < min_local[idim] ) {
                    if( neighbor_[idim][0]!=MPI_PROC_NULL ) {
                        vecSpecies[ispec]->MPI_buffer_.part_index_send[idim][0].push_back( iPart );
                    }
                    //If particle is outside of the global domain (has no neighbor), it will not be put in a send buffer and will simply be deleted.
                    check = 1;
                } else if( cuParticles.position( idim, iPart ) >= max_local[idim] ) {
                    if( neighbor_[idim][1]!=MPI_PROC_NULL ) {
                        vecSpecies[ispec]->MPI_buffer_.part_index_send[idim][1].push_back( iPart );
                    }
                    check = 1;
                }
                idim++;
            }
        }
    } else { //if (geometry == "AMcylindrical")
        double r_min2, r_max2;
        r_max2 = max_local[1] * max_local[1] ;
        r_min2 = min_local[1] * min_local[1] ;
        for( int i=0 ; i<n_part_send ; i++ ) {
            iPart = ( *indexes_of_particles_to_exchange )[i];
            //Put indexes of particles in the first direction they will be exchanged and correct their position according to periodicity for the first exchange only.
            if( cuParticles.position( 0, iPart ) < min_local[0] ) {
                if( neighbor_[0][0]!=MPI_PROC_NULL ) {
                    if ( (Pcoordinates[0]==0) && ( vecSpecies[ispec]->boundary_conditions[0][0]!="periodic" ) ) {
                        continue;
                    }
                    vecSpecies[ispec]->MPI_buffer_.part_index_send[0][0].push_back( iPart );
                    //MESSAGE("Sending particle to the left x= " << cuParticles.position(0,iPart) <<  " xmin = " <<  min_local[0] );
                }
                //If particle is outside of the global domain (has no neighbor), it will not be put in a send buffer and will simply be deleted.
            } else if( cuParticles.position( 0, iPart ) >= max_local[0] ) {
                if ( (Pcoordinates[0]==params.number_of_patches[0]-1) && ( vecSpecies[ispec]->boundary_conditions[0][1]!="periodic" ) ) {
                    continue;
                }
                if( neighbor_[0][1]!=MPI_PROC_NULL ) {
                    vecSpecies[ispec]->MPI_buffer_.part_index_send[0][1].push_back( iPart );
                    // MESSAGE("Sending particle to the right x= " << cuParticles.position(0,iPart) <<  " xmax = " <<  max_local[0] );
                }
            } else if( cuParticles.distance2ToAxis( iPart ) < r_min2 ) {
                if( neighbor_[1][0]!=MPI_PROC_NULL ) {
                    vecSpecies[ispec]->MPI_buffer_.part_index_send[1][0].push_back( iPart );
                    //MESSAGE("Sending particle to the south r= " << cuParticles.distance2ToAxis(iPart) <<  " rmin2 = " <<  r_min2 );
                }
            } else if( cuParticles.distance2ToAxis( iPart ) >= r_max2 ) {
                if( neighbor_[1][1]!=MPI_PROC_NULL ) {
                    vecSpecies[ispec]->MPI_buffer_.part_index_send[1][1].push_back( iPart );
                    //MESSAGE("Sending particle to the north r= " << cuParticles.distance2ToAxis(iPart) <<  " rmax2 = " <<  r_max2 << " rmin2= " << r_min2 );
                }
            }

        }
    }

} // initExchParticles(... iDim)


// ---------------------------------------------------------------------------------------------------------------------
// For direction iDim, start exchange of number of particles
//   - vecPatch : used for intra-MPI process comm (direct copy using Particels::copyParticles)
//   - smpi     : inhereted from previous SmileiMPI::exchangeParticles()
// ---------------------------------------------------------------------------------------------------------------------
void Patch::exchNbrOfParticles( SmileiMPI *smpi, int ispec, Params &params, int iDim, VectorPatch *vecPatch )
{
    int h0 = ( *vecPatch )( 0 )->hindex;
    /********************************************************************************/
    // Exchange number of particles to exchange to establish or not a communication
    /********************************************************************************/
    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
        if( neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL ) {
            vecSpecies[ispec]->MPI_buffer_.part_index_send_sz[iDim][iNeighbor] = ( vecSpecies[ispec]->MPI_buffer_.part_index_send[iDim][iNeighbor] ).size();

            if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
                //If neighbour is MPI ==> I send him the number of particles I'll send later.
                int local_hindex = hindex - vecPatch->refHindex_;
                int tag = buildtag( local_hindex, iDim+1, iNeighbor+3 );
                MPI_Isend( &( vecSpecies[ispec]->MPI_buffer_.part_index_send_sz[iDim][iNeighbor] ), 1, MPI_INT, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &( vecSpecies[ispec]->MPI_buffer_.srequest[iDim][iNeighbor] ) );
            } else {
                //Else, I directly set the receive size to the correct value.
                ( *vecPatch )( neighbor_[iDim][iNeighbor]- h0 )->vecSpecies[ispec]->MPI_buffer_.part_index_recv_sz[iDim][( iNeighbor+1 )%2] = vecSpecies[ispec]->MPI_buffer_.part_index_send_sz[iDim][iNeighbor];
            }
        } // END of Send

        if( neighbor_[iDim][( iNeighbor+1 )%2]!=MPI_PROC_NULL ) {
            if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
                //If other neighbour is MPI ==> I receive the number of particles I'll receive later.
                int local_hindex = neighbor_[iDim][( iNeighbor+1 )%2] - smpi->patch_refHindexes[ MPI_neighbor_[iDim][( iNeighbor+1 )%2] ];
                int tag = buildtag( local_hindex, iDim+1, iNeighbor+3 );
                MPI_Irecv( &( vecSpecies[ispec]->MPI_buffer_.part_index_recv_sz[iDim][( iNeighbor+1 )%2] ), 1, MPI_INT, MPI_neighbor_[iDim][( iNeighbor+1 )%2], tag, MPI_COMM_WORLD, &( vecSpecies[ispec]->MPI_buffer_.rrequest[iDim][( iNeighbor+1 )%2] ) );
            }
        }
    }//end loop on nb_neighbors.

} // exchNbrOfParticles(... iDim)


void Patch::endNbrOfParticles( SmileiMPI *smpi, int ispec, Params &params, int iDim, VectorPatch *vecPatch )
{
    Particles &cuParticles = ( *vecSpecies[ispec]->particles );

    /********************************************************************************/
    // Wait for end of communications over number of particles
    /********************************************************************************/
    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
        MPI_Status sstat    [2];
        MPI_Status rstat    [2];
        if( neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL ) {
            if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
                MPI_Wait( &( vecSpecies[ispec]->MPI_buffer_.srequest[iDim][iNeighbor] ), &( sstat[iNeighbor] ) );
            }
        }
        if( neighbor_[iDim][( iNeighbor+1 )%2]!=MPI_PROC_NULL ) {
            if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) )  {
                MPI_Wait( &( vecSpecies[ispec]->MPI_buffer_.rrequest[iDim][( iNeighbor+1 )%2] ), &( rstat[( iNeighbor+1 )%2] ) );
                if( vecSpecies[ispec]->MPI_buffer_.part_index_recv_sz[iDim][( iNeighbor+1 )%2]!=0 ) {
                    //If I receive particles over MPI, I initialize my receive buffer with the appropriate size.
                    vecSpecies[ispec]->MPI_buffer_.partRecv[iDim][( iNeighbor+1 )%2].initialize( vecSpecies[ispec]->MPI_buffer_.part_index_recv_sz[iDim][( iNeighbor+1 )%2], cuParticles );
                }
            }
        }
    }

} // END endNbrOfParticles(... iDim)


// ---------------------------------------------------------------------------------------------------------------------
// For direction iDim, finalize receive of number of particles and really send particles
//   - vecPatch : used for intra-MPI process comm (direct copy using Particels::copyParticles)
//   - smpi     : used smpi->periods_
// ---------------------------------------------------------------------------------------------------------------------
void Patch::prepareParticles( SmileiMPI *smpi, int ispec, Params &params, int iDim, VectorPatch *vecPatch )
{
    Particles &cuParticles = ( *vecSpecies[ispec]->particles );

    int n_part_send;
    int h0 = ( *vecPatch )( 0 )->hindex;
    double x_max = params.cell_length[iDim]*( params.n_space_global[iDim] );

    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {

        // n_part_send : number of particles to send to current neighbor
        n_part_send = ( vecSpecies[ispec]->MPI_buffer_.part_index_send[iDim][iNeighbor] ).size();
        if( ( neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL ) && ( n_part_send!=0 ) ) {
            // Enabled periodicity
            if( smpi->periods_[iDim]==1 ) {
                for( int iPart=0 ; iPart<n_part_send ; iPart++ ) {
                    if( ( iNeighbor==0 ) && ( Pcoordinates[iDim] == 0 ) &&( cuParticles.position( iDim, vecSpecies[ispec]->MPI_buffer_.part_index_send[iDim][iNeighbor][iPart] ) < 0. ) ) {
                        cuParticles.position( iDim, vecSpecies[ispec]->MPI_buffer_.part_index_send[iDim][iNeighbor][iPart] )     += x_max;
                    } else if( ( iNeighbor==1 ) && ( Pcoordinates[iDim] == params.number_of_patches[iDim]-1 ) && ( cuParticles.position( iDim, vecSpecies[ispec]->MPI_buffer_.part_index_send[iDim][iNeighbor][iPart] ) >= x_max ) ) {
                        cuParticles.position( iDim, vecSpecies[ispec]->MPI_buffer_.part_index_send[iDim][iNeighbor][iPart] )     -= x_max;
                    }
                }
            }
            // Send particles
            if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
                // If MPI comm, first copy particles in the sendbuffer
                for( int iPart=0 ; iPart<n_part_send ; iPart++ ) {
                    cuParticles.copyParticle( vecSpecies[ispec]->MPI_buffer_.part_index_send[iDim][iNeighbor][iPart], vecSpecies[ispec]->MPI_buffer_.partSend[iDim][iNeighbor] );
                }
            } else {
                //If not MPI comm, copy particles directly in the receive buffer
                for( int iPart=0 ; iPart<n_part_send ; iPart++ ) {
                    cuParticles.copyParticle( vecSpecies[ispec]->MPI_buffer_.part_index_send[iDim][iNeighbor][iPart], ( ( *vecPatch )( neighbor_[iDim][iNeighbor]- h0 )->vecSpecies[ispec]->MPI_buffer_.partRecv[iDim][( iNeighbor+1 )%2] ) );
                }
            }
        } // END of Send

    } // END for iNeighbor

} // END prepareParticles(... iDim)


void Patch::exchParticles( SmileiMPI *smpi, int ispec, Params &params, int iDim, VectorPatch *vecPatch )
{
    int n_part_send, n_part_recv;

    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {

        // n_part_send : number of particles to send to current neighbor
        n_part_send = ( vecSpecies[ispec]->MPI_buffer_.part_index_send[iDim][iNeighbor] ).size();
        if( ( neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL ) && ( n_part_send!=0 ) ) {
            // Send particles
            if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
                // Then send particles
                int local_hindex = hindex - vecPatch->refHindex_;
                int tag = buildtag( local_hindex, iDim+1, iNeighbor+3 );
                vecSpecies[ispec]->typePartSend[( iDim*2 )+iNeighbor] = smpi->createMPIparticles( &( vecSpecies[ispec]->MPI_buffer_.partSend[iDim][iNeighbor] ) );
                MPI_Isend( &( ( vecSpecies[ispec]->MPI_buffer_.partSend[iDim][iNeighbor] ).position( 0, 0 ) ), 1, vecSpecies[ispec]->typePartSend[( iDim*2 )+iNeighbor], MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &( vecSpecies[ispec]->MPI_buffer_.srequest[iDim][iNeighbor] ) );
            }
        } // END of Send

        n_part_recv = vecSpecies[ispec]->MPI_buffer_.part_index_recv_sz[iDim][( iNeighbor+1 )%2];
        if( ( neighbor_[iDim][( iNeighbor+1 )%2]!=MPI_PROC_NULL ) && ( n_part_recv!=0 ) ) {
            if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
                // If MPI comm, receive particles in the recv buffer previously initialized.
                vecSpecies[ispec]->typePartRecv[( iDim*2 )+iNeighbor] = smpi->createMPIparticles( &( vecSpecies[ispec]->MPI_buffer_.partRecv[iDim][( iNeighbor+1 )%2] ) );
                int local_hindex = neighbor_[iDim][( iNeighbor+1 )%2] - smpi->patch_refHindexes[ MPI_neighbor_[iDim][( iNeighbor+1 )%2] ];
                int tag = buildtag( local_hindex, iDim+1, iNeighbor+3 );
                MPI_Irecv( &( ( vecSpecies[ispec]->MPI_buffer_.partRecv[iDim][( iNeighbor+1 )%2] ).position( 0, 0 ) ), 1, vecSpecies[ispec]->typePartRecv[( iDim*2 )+iNeighbor], MPI_neighbor_[iDim][( iNeighbor+1 )%2], tag, MPI_COMM_WORLD, &( vecSpecies[ispec]->MPI_buffer_.rrequest[iDim][( iNeighbor+1 )%2] ) );
            }

        } // END of Recv

    } // END for iNeighbor

} // END exchParticles(... iDim)


// ---------------------------------------------------------------------------------------------------------------------
// For direction iDim, finalize receive of particles, temporary store particles if diagonalParticles
// And store recv particles at their definitive place.
// Call Patch::cleanupSentParticles
//   - vecPatch : used for intra-MPI process comm (direct copy using Particels::copyParticles)
//   - smpi     : used smpi->periods_
// ---------------------------------------------------------------------------------------------------------------------
void Patch::finalizeExchParticles( SmileiMPI *smpi, int ispec, Params &params, int iDim, VectorPatch *vecPatch )
{

    int n_part_send, n_part_recv;

    /********************************************************************************/
    // Wait for end of communications over Particles
    /********************************************************************************/
    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
        MPI_Status sstat    [2];
        MPI_Status rstat    [2];

        n_part_send = vecSpecies[ispec]->MPI_buffer_.part_index_send[iDim][iNeighbor].size();
        n_part_recv = vecSpecies[ispec]->MPI_buffer_.part_index_recv_sz[iDim][( iNeighbor+1 )%2];

        if( ( neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL ) && ( n_part_send!=0 ) ) {
            if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
                MPI_Wait( &( vecSpecies[ispec]->MPI_buffer_.srequest[iDim][iNeighbor] ), &( sstat[iNeighbor] ) );
                MPI_Type_free( &( vecSpecies[ispec]->typePartSend[( iDim*2 )+iNeighbor] ) );
            }
        }
        if( ( neighbor_[iDim][( iNeighbor+1 )%2]!=MPI_PROC_NULL ) && ( n_part_recv!=0 ) ) {
            if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
                MPI_Wait( &( vecSpecies[ispec]->MPI_buffer_.rrequest[iDim][( iNeighbor+1 )%2] ), &( rstat[( iNeighbor+1 )%2] ) );
                MPI_Type_free( &( vecSpecies[ispec]->typePartRecv[( iDim*2 )+iNeighbor] ) );
            }
        }
    }
}

void Patch::cornersParticles( SmileiMPI *smpi, int ispec, Params &params, int iDim, VectorPatch *vecPatch )
{

    int ndim = params.nDim_field;
    int idim, check;

    Particles &cuParticles = ( *vecSpecies[ispec]->particles );

    int n_part_recv;

    /********************************************************************************/
    // Wait for end of communications over Particles
    /********************************************************************************/
    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {

        n_part_recv = vecSpecies[ispec]->MPI_buffer_.part_index_recv_sz[iDim][( iNeighbor+1 )%2];

        if( ( neighbor_[iDim][( iNeighbor+1 )%2]!=MPI_PROC_NULL ) && ( n_part_recv!=0 ) ) {

            // Treat diagonalParticles
            if( iDim < ndim-1 ) { // No need to treat diag particles at last dimension.
                if( params.geometry != "AMcylindrical" ) {
                    for( int iPart=n_part_recv-1 ; iPart>=0; iPart-- ) {
                        check = 0;
                        idim = iDim+1;//We check next dimension
                        while( check == 0 && idim<ndim ) {
                            //If particle not in the domain...
                            if( ( vecSpecies[ispec]->MPI_buffer_.partRecv[iDim][( iNeighbor+1 )%2] ).position( idim, iPart ) < min_local[idim] ) {
                                if( neighbor_[idim][0]!=MPI_PROC_NULL ) { //if neighbour exists
                                    //... copy it at the back of the local particle vector ...
                                    ( vecSpecies[ispec]->MPI_buffer_.partRecv[iDim][( iNeighbor+1 )%2] ).copyParticle( iPart, cuParticles );
                                    //...adjust particles->last_index or cell_keys ...
                                    vecSpecies[ispec]->addSpaceForOneParticle();
                                    //... and add its index to the particles to be sent later...
                                    vecSpecies[ispec]->MPI_buffer_.part_index_send[idim][0].push_back( cuParticles.size()-1 );
                                    //..without forgeting to add it to the list of particles to clean.
                                    vecSpecies[ispec]->addPartInExchList( cuParticles.size()-1 );
                                }
                                //Remove it from receive buffer.
                                ( vecSpecies[ispec]->MPI_buffer_.partRecv[iDim][( iNeighbor+1 )%2] ).eraseParticle( iPart );
                                vecSpecies[ispec]->MPI_buffer_.part_index_recv_sz[iDim][( iNeighbor+1 )%2]--;
                                check = 1;
                            }
                            //Other side of idim
                            else if( ( vecSpecies[ispec]->MPI_buffer_.partRecv[iDim][( iNeighbor+1 )%2] ).position( idim, iPart ) >= max_local[idim] ) {
                                if( neighbor_[idim][1]!=MPI_PROC_NULL ) { //if neighbour exists
                                    ( vecSpecies[ispec]->MPI_buffer_.partRecv[iDim][( iNeighbor+1 )%2] ).copyParticle( iPart, cuParticles );
                                    //...adjust particles->last_index or cell_keys ...
                                    vecSpecies[ispec]->addSpaceForOneParticle();
                                    vecSpecies[ispec]->MPI_buffer_.part_index_send[idim][1].push_back( cuParticles.size()-1 );
                                    vecSpecies[ispec]->addPartInExchList( cuParticles.size()-1 );
                                }
                                ( vecSpecies[ispec]->MPI_buffer_.partRecv[iDim][( iNeighbor+1 )%2] ).eraseParticle( iPart );
                                vecSpecies[ispec]->MPI_buffer_.part_index_recv_sz[iDim][( iNeighbor+1 )%2]--;
                                check = 1;
                            }
                            idim++;
                        }
                    }
                } else { //In AM geometry
                    //In this case, iDim = 0 and idim = iDim + 1 = 1. We only have to check potential comms along R.
                    double r_min2, r_max2;
                    r_min2 = min_local[1]*min_local[1];
                    r_max2 = max_local[1]*max_local[1];
                    for( int iPart=n_part_recv-1 ; iPart>=0; iPart-- ) {
                        //MESSAGE("test particle diag r2 = " << (vecSpecies[ispec]->MPI_buffer_.partRecv[0][(iNeighbor+1)%2]).distance2ToAxis(iPart) << "rmin2 = " << r_min2 << " rmax2 = " << r_max2 );
                        if( ( vecSpecies[ispec]->MPI_buffer_.partRecv[0][( iNeighbor+1 )%2] ).distance2ToAxis( iPart ) < r_min2 ) {
                            if( neighbor_[1][0]!=MPI_PROC_NULL ) { //if neighbour exists
                                //... copy it at the back of the local particle vector ...
                                ( vecSpecies[ispec]->MPI_buffer_.partRecv[0][( iNeighbor+1 )%2] ).copyParticle( iPart, cuParticles );
                                //...adjust particles->last_index or cell_keys ...
                                vecSpecies[ispec]->addSpaceForOneParticle();
                                //... and add its index to the particles to be sent later...
                                vecSpecies[ispec]->MPI_buffer_.part_index_send[1][0].push_back( cuParticles.size()-1 );
                                //..without forgeting to add it to the list of particles to clean.
                                vecSpecies[ispec]->addPartInExchList( cuParticles.size()-1 );
                            }
                            //Remove it from receive buffer.
                            ( vecSpecies[ispec]->MPI_buffer_.partRecv[0][( iNeighbor+1 )%2] ).eraseParticle( iPart );
                            vecSpecies[ispec]->MPI_buffer_.part_index_recv_sz[0][( iNeighbor+1 )%2]--;
                        }
                        //Other side of idim
                        else if( ( vecSpecies[ispec]->MPI_buffer_.partRecv[0][( iNeighbor+1 )%2] ).distance2ToAxis( iPart ) >= r_max2 ) {
                            if( neighbor_[1][1]!=MPI_PROC_NULL ) { //if neighbour exists
                                //MESSAGE("particle diag +R");
                                ( vecSpecies[ispec]->MPI_buffer_.partRecv[0][( iNeighbor+1 )%2] ).copyParticle( iPart, cuParticles );
                                //...adjust particles->last_index or cell_keys ...
                                vecSpecies[ispec]->addSpaceForOneParticle();
                                vecSpecies[ispec]->MPI_buffer_.part_index_send[1][1].push_back( cuParticles.size()-1 );
                                vecSpecies[ispec]->addPartInExchList( cuParticles.size()-1 );
                            }
                            ( vecSpecies[ispec]->MPI_buffer_.partRecv[0][( iNeighbor+1 )%2] ).eraseParticle( iPart );
                            vecSpecies[ispec]->MPI_buffer_.part_index_recv_sz[0][( iNeighbor+1 )%2]--;
                        }
                    }
                }
            }//If not last dim for diagonal particles.
        } //If received something
    } //loop i Neighbor
}

//! Import particles exchanged with surrounding patches/mpi and sort at the same time
void Patch::importAndSortParticles( SmileiMPI *smpi, int ispec, Params &params, VectorPatch *vecPatch )
{

#ifdef  __DETAILED_TIMERS
    double timer;
    timer = MPI_Wtime();
#endif

    vecSpecies[ispec]->sortParticles( params , this);

#ifdef  __DETAILED_TIMERS
    this->patch_timers[13] += MPI_Wtime() - timer;
#endif

} // sortParticles(...)


void Patch::cleanParticlesOverhead( Params &params )
{
    int ndim = params.nDim_field;
    for( unsigned int ispec=0 ; ispec<vecSpecies.size() ; ispec++ ) {
        Particles &cuParticles = ( *vecSpecies[ispec]->particles );

        for( int idim = 0; idim < ndim; idim++ ) {
            for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
                vecSpecies[ispec]->MPI_buffer_.partRecv[idim][iNeighbor].clear();
                vecSpecies[ispec]->MPI_buffer_.partRecv[idim][iNeighbor].shrinkToFit( );
                vecSpecies[ispec]->MPI_buffer_.partSend[idim][iNeighbor].clear();
                vecSpecies[ispec]->MPI_buffer_.partSend[idim][iNeighbor].shrinkToFit( );
                vecSpecies[ispec]->MPI_buffer_.part_index_send[idim][iNeighbor].clear();
                vector<int>( vecSpecies[ispec]->MPI_buffer_.part_index_send[idim][iNeighbor] ).swap( vecSpecies[ispec]->MPI_buffer_.part_index_send[idim][iNeighbor] );
            }
        }

        cuParticles.shrinkToFit(  );
    }

}

// ---------------------------------------------------------------------------------------------------------------------
// Clear vecSpecies[]->indexes_of_particles_to_exchange, suppress particles send and manage memory
// ---------------------------------------------------------------------------------------------------------------------
void Patch::cleanupSentParticles( int ispec, std::vector<int> *indexes_of_particles_to_exchange )
{
    /********************************************************************************/
    // Delete Particles included in the index of particles to exchange. Assumes indexes are sorted.
    /********************************************************************************/
    int ii, iPart;
    std::vector<int> *cufirst_index = &vecSpecies[ispec]->particles->first_index;
    std::vector<int> *culast_index = &vecSpecies[ispec]->particles->last_index;
    Particles &cuParticles = ( *vecSpecies[ispec]->particles );


    // Push lost particles at the end of bins
    for( unsigned int ibin = 0 ; ibin < culast_index->size() ; ibin++ ) {
        ii = indexes_of_particles_to_exchange->size()-1;
        if( ii >= 0 ) { // Push lost particles to the end of the bin
            iPart = ( *indexes_of_particles_to_exchange )[ii];
            while( iPart >= ( *culast_index )[ibin] && ii > 0 ) {
                ii--;
                iPart = ( *indexes_of_particles_to_exchange )[ii];
            }
            while( iPart == ( *culast_index )[ibin]-1 && iPart >= ( *cufirst_index )[ibin] && ii > 0 ) {
                ( *culast_index )[ibin]--;
                ii--;
                iPart = ( *indexes_of_particles_to_exchange )[ii];
            }
            while( iPart >= ( *cufirst_index )[ibin] && ii > 0 ) {
                cuParticles.overwriteParticle( ( *culast_index )[ibin]-1, iPart );
                ( *culast_index )[ibin]--;
                ii--;
                iPart = ( *indexes_of_particles_to_exchange )[ii];
            }
            if( iPart >= ( *cufirst_index )[ibin] && iPart < ( *culast_index )[ibin] ) { //On traite la dernière particule (qui peut aussi etre la premiere)
                cuParticles.overwriteParticle( ( *culast_index )[ibin]-1, iPart );
                ( *culast_index )[ibin]--;
            }
        }
    }


    //Shift the bins in memory
    //Warning: this loop must be executed sequentially. Do not use openMP here.
    for( int unsigned ibin = 1 ; ibin < culast_index->size() ; ibin++ ) { //First bin don't need to be shifted
        ii = ( *cufirst_index )[ibin]-( *culast_index )[ibin-1]; // Shift the bin in memory by ii slots.
        iPart = min( ii, ( *culast_index )[ibin]-( *cufirst_index )[ibin] ); // Number of particles we have to shift = min (Nshift, Nparticle in the bin)
        if( iPart > 0 ) {
            cuParticles.overwriteParticle( ( *culast_index )[ibin]-iPart, ( *culast_index )[ibin-1], iPart );
        }
        ( *culast_index )[ibin] -= ii;
        ( *cufirst_index )[ibin] = ( *culast_index )[ibin-1];
    }

} // END cleanupSentParticles

//Copy positions of all particles of the target species to the positions of the species to update.
//Used for particle initialization on top of another species
void Patch::copyPositions( std::vector<Species *> vecSpecies_to_update )
{
    for( unsigned int i=0; i<vecSpecies_to_update.size(); i++ ) {
        if( vecSpecies_to_update[i]->position_initialization_on_species_==false )
            continue;
        unsigned int target_species = vecSpecies_to_update[i]->position_initialization_on_species_index;
        if( vecSpecies_to_update[i]->getNbrOfParticles() != vecSpecies_to_update[target_species]->getNbrOfParticles() ) {
            ERROR( "Number of particles in species '"<<vecSpecies_to_update[i]->name_<<"' is not equal to the number of particles in species '"<<vecSpecies_to_update[target_species]->name_<<"'." );
        }
        // We copy target_species which is the index of the species, already created, from which species_to_update particles positions are copied.
        vecSpecies_to_update[i]->particles->Position = vecSpecies_to_update[target_species]->particles->Position;
    }
    return;
}
