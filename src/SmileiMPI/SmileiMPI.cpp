#include "SmileiMPI.h"

#include <cmath>
#include <cstring>
#include <algorithm>

#include <iostream>
#include <sstream>
#include <fstream>

#include "Params.h"
#include "Tools.h"

#include "ElectroMagn.h"
#include "ElectroMagnBC1D_SM.h"
#include "ElectroMagnBC2D_SM.h"
#include "ElectroMagnBC3D_SM.h"
#include "ElectroMagnBC2D_PML.h"
#include "ElectroMagnBC3D_PML.h"
#include "ElectroMagnBCAM_PML.h"
#include "EnvelopeBCAM_PML.h"
#include "EnvelopeBC2D_PML.h"
#include "EnvelopeBC3D_PML.h"
#include "Field.h"

#include "Species.h"
#include "PeekAtSpecies.h"
#include "Hilbert_functions.h"
#include "VectorPatch.h"
#include "DomainDecomposition.h"

#include "Diagnostic.h"
#include "DiagnosticScalar.h"
#include "DiagnosticParticleBinning.h"
#include "DiagnosticScreen.h"
#include "DiagnosticRadiationSpectrum.h"
#include "DiagnosticProbes.h"

#include "Laser.h"
#include "LaserEnvelope.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// SmileiMPI constructor :
//     - Call MPI_Init_thread, MPI_THREAD_MULTIPLE required
//     - Set MPI env
// ---------------------------------------------------------------------------------------------------------------------
SmileiMPI::SmileiMPI( int *argc, char ***argv )
{
    test_mode = false;

    // Send information on current simulation
    int mpi_provided;

#ifdef _OPENMP
#ifdef _NO_MPI_TM
    MPI_Init_thread( argc, argv, MPI_THREAD_SERIALIZED, &mpi_provided );
#else
    MPI_Init_thread( argc, argv, MPI_THREAD_MULTIPLE, &mpi_provided );
    if( mpi_provided != MPI_THREAD_MULTIPLE ) {
        ERROR( "MPI_THREAD_MULTIPLE not supported. Compile your MPI library with THREAD_MULTIPLE support." );
    }
#endif
    smilei_omp_max_threads = omp_get_max_threads();
    number_of_cores = min( omp_get_num_procs(), smilei_omp_max_threads );
#else
    MPI_Init( argc, argv );
    smilei_omp_max_threads = 1;
    number_of_cores = 1;
#endif

    world_ = MPI_COMM_WORLD;
    MPI_Comm_size( world_, &smilei_sz );
    MPI_Comm_rank( world_, &smilei_rk );

    MPI_Allreduce( &number_of_cores, &global_number_of_cores, 1, MPI_INT, MPI_SUM, world_ );
} // END SmileiMPI::SmileiMPI


// ---------------------------------------------------------------------------------------------------------------------
// SmileiMPI destructor :
//     - Call MPI_Finalize
// ---------------------------------------------------------------------------------------------------------------------
SmileiMPI::~SmileiMPI()
{
    delete[]periods_;

    MPI_Finalize();

} // END SmileiMPI::~SmileiMPI


// ---------------------------------------------------------------------------------------------------------------------
// Broadcast namelist in world_
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI::bcast( string &val )
{
    int charSize=0;
    if( isMaster() ) {
        charSize = val.size()+1;
    }
    MPI_Bcast( &charSize, 1, MPI_INT, 0, world_ );

    char tmp[charSize];
    if( isMaster() ) {
        strcpy( tmp, val.c_str() );
    }
    MPI_Bcast( tmp, charSize, MPI_CHAR, 0, world_ );

    if( !isMaster() ) {
        val=tmp;
    }

} // END bcast( string )


// ---------------------------------------------------------------------------------------------------------------------
// Broadcast namelist
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI::bcast( int &val )
{
    MPI_Bcast( &val, 1, MPI_INT, 0, world_ );

} // END bcast( int ) in world_


// ---------------------------------------------------------------------------------------------------------------------
// Initialize MPI (per process) environment
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI::init( Params &params, DomainDecomposition *domain_decomposition )
{
    // Initialize patch environment
    patch_count.resize( smilei_sz, 0 );
    capabilities.resize( smilei_sz, 1 );
    Tcapabilities = smilei_sz;

    if( smilei_rk == 0 ) {
        remove( "patch_load.txt" ) ;
    }
    // Initialize patch distribution
    if( !params.restart ) {
        init_patch_count( params, domain_decomposition );
    }

    // Initialize buffers for particles push vectorization
    //     - 1 thread push particles for a unique patch at a given time
    //     - so 1 buffer per thread

    int n_envlaser = PyTools::nComponents( "LaserEnvelope" );

#ifdef _OPENMP
    dynamics_Epart.resize( omp_get_max_threads() );
    dynamics_Bpart.resize( omp_get_max_threads() );
    dynamics_invgf.resize( omp_get_max_threads() );
    dynamics_iold.resize( omp_get_max_threads() );
    dynamics_deltaold.resize( omp_get_max_threads() );
    if( params.geometry == "AMcylindrical" ) {
        dynamics_thetaold.resize( omp_get_max_threads() );
    }

    if( n_envlaser > 0 ) {
        dynamics_GradPHIpart.resize( omp_get_max_threads() );
        dynamics_GradPHI_mpart.resize( omp_get_max_threads() );
        dynamics_PHIpart.resize( omp_get_max_threads() );
        dynamics_PHI_mpart.resize( omp_get_max_threads() );
        dynamics_inv_gamma_ponderomotive.resize( omp_get_max_threads() );
        if (params.envelope_ionization_is_active){
            dynamics_EnvEabs_part.resize( omp_get_max_threads() );
            dynamics_EnvExabs_part.resize( omp_get_max_threads() );
        } else {
            dynamics_EnvEabs_part.clear();
            dynamics_EnvExabs_part.clear();
        }
    }
#else
    dynamics_Epart.resize( 1 );
    dynamics_Bpart.resize( 1 );
    dynamics_invgf.resize( 1 );
    dynamics_iold.resize( 1 );
    dynamics_deltaold.resize( 1 );
    if( params.geometry == "AMcylindrical" ) {
        dynamics_thetaold.resize( 1 );
    }

    if( n_envlaser > 0 ) {
        dynamics_GradPHIpart.resize( 1 );
        dynamics_GradPHI_mpart.resize( 1 );
        dynamics_PHIpart.resize( 1 );
        dynamics_PHI_mpart.resize( 1 );
        dynamics_inv_gamma_ponderomotive.resize( 1 );
        if (params.envelope_ionization_is_active){
            dynamics_EnvEabs_part.resize( 1 );
            dynamics_EnvExabs_part.resize( 1 );
        } else {
            dynamics_EnvEabs_part.clear();
            dynamics_EnvExabs_part.clear();
        }
    }
#endif

    // Set periodicity of the simulated problem
    periods_  = new int[params.nDim_field];
    for( unsigned int i=0 ; i<params.nDim_field ; i++ ) {
        periods_[i] = 0;
        if( params.EM_BCs[i][0]=="periodic" ) {
            periods_[i] = 1;
            MESSAGE( 1, "applied topology for periodic BCs in "<<"xyz"[i]<<"-direction" );
        }
    }

    // Estimate the maximum tag requires by Smilei regarding the current patch distribution
    auto it = max_element(std::begin(patch_count), std::end(patch_count));
    // the maximum tag use the maximum local patch id, iDim=1, iNeghibor=1, 8 for Jx
    int tagmax = buildtag( (*it)-1, 1, 1, 8 );
    // Compare to MPI tag upper bound
    int tagUB = getTagUB();
    if ( tagmax > tagUB ) {
        int ratio = ceil( (double)tagmax/(double)tagUB );
        if (!smilei_rk) {
            ERROR( "The MPI library you are using authorizes as upper bound for a tag : " << tagUB << endl <<
                   "Regarding the number of patches you are using and the number of MPI process, Smilei will at least generate as larger tag : " << tagmax << endl <<
                   "You should use " << ratio << " more MPI process, or " << ratio << " less patches or change of MPI library." << endl <<
                   "This is a trough estimation, which depends of the plasma load imbalance." );
        }
    }
} // END init


// ---------------------------------------------------------------------------------------------------------------------
//  Initialize patch distribution
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI::init_patch_count( Params &params, DomainDecomposition *domain_decomposition )
{

//#ifndef _NOTBALANCED
//    bool use_load_balancing(true);
//    if (!use_load_balancing) {
//        int Npatches = params.number_of_patches[0];
//        for (unsigned int i = 1; i < params.nDim_field; i++)
//            Npatches *=  params.number_of_patches[i]; // Total number of patches.
//        if (Npatches!=smilei_sz) ERROR("number of patches abd MPI processes");
//        for (unsigned int i=0; i<smilei_sz; i++) patch_count[i] = 1;
//        return;
//    }
//#endif

    std::vector<unsigned int> Pcoordinates( 3, 0 );
    unsigned int Npatches, r, Ncur, tot_ncells_perpatch;
    double Tload, Tcur, Lcur, total_load=0, local_load, above_target, below_target;

    unsigned int tot_species_number = PyTools::nComponents( "Species" );

    // Define capabilities here if not default.
    //Capabilities of devices hosting the different mpi processes. All capabilities are assumed to be equal for the moment.
    //Compute total capability: Tcapabilities. Uncomment if cpability != 1 per MPI rank
    //Tcapabilities = 0;
    //for (unsigned int i = 0; i < smilei_sz; i++)
    //    Tcapabilities += capabilities[i];

    //Compute target load: Tload = Total load * local capability / Total capability.

    // Some initialization of the box parameters
    Npatches = params.tot_number_of_patches;
    tot_ncells_perpatch = 1;
    vector<double> x_cell( 3, 0. );
    for( unsigned int i = 0; i < params.nDim_field; i++ ) {
        tot_ncells_perpatch *= params.n_space[i]+2*params.oversize[i];
    }

    // First, distribute all patches evenly
    unsigned int Npatches_local = Npatches / smilei_sz, FirstPatch_local;
    int remainder = Npatches % smilei_sz;
    if( smilei_rk < remainder ) {
        Npatches_local++;
        FirstPatch_local = Npatches_local * smilei_rk;
    } else {
        FirstPatch_local = Npatches_local * smilei_rk + remainder;
    }
//    // Test
//    int tot, loc=Npatches_local;
//    MPI_Allreduce( &loc, &tot, 1, MPI_INT, MPI_SUM, world_ );
//    if( tot != Npatches ) ERROR("Npatches should be "<<Npatches<<" but it is "<<tot);

    // Second, prepare the profiles for each species
    vector<PeekAtSpecies *> peek;
    peek.reserve( tot_species_number );
    for( unsigned int ispecies = 0; ispecies < tot_species_number; ispecies++ ) {
        peek.push_back( new PeekAtSpecies( params, ispecies ) );
    }

    // Third, loop over local patches to obtain their approximate load
    vector<double> PatchLoad( Npatches_local, 1. );
    if( !( params.has_load_balancing && params.initial_balance ) ) {
        total_load = Npatches_local; //We don't balance the simulation, all patches have a load of 1.
    } else {
        for( unsigned int ipatch=0; ipatch<Npatches_local; ipatch++ ) {
            // Get patch coordinates
            unsigned int hindex = FirstPatch_local + ipatch;
            Pcoordinates = domain_decomposition->getDomainCoordinates( hindex );
            for( unsigned int i=0 ; i<params.nDim_field ; i++ ) {
                x_cell[i] = ( Pcoordinates[i]+0.5 )*params.patch_dimensions[i];
            }
            //Accumulate particles load of the current patch
            for( unsigned int ispecies = 0; ispecies < tot_species_number; ispecies++ ) {
                local_load = peek[ispecies]->numberOfParticlesInPatch( x_cell );
                // Consider whether this species is frozen
                double time_frozen( 0. );
                PyTools::extract( "time_frozen", time_frozen, "Species", ispecies );
                if( time_frozen > 0. ) {
                    local_load *= params.frozen_particle_load;
                }
                // Add the load of the species to the current patch load
                PatchLoad[ipatch] += local_load;
            }
            //Add grid contribution to the load.
            PatchLoad[ipatch] += tot_ncells_perpatch*params.cell_load-1; //-1 to compensate the initialization to 1.
            total_load += PatchLoad[ipatch];
        }
    }
    // Clean the temporary species profiles
    for( unsigned int i=0 ; i<tot_species_number ; i++ ) {
        delete peek[i];
    }
    peek.clear();

    // Fourth, the arrangement of patches is balanced

    // Initialize loads
    MPI_Reduce( &total_load, &Tload, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    Tload /= Tcapabilities; //Target load for each mpi process.
    Tcur = Tload * capabilities[0];  //Init.
    r = 0;  //Start by finding work for rank 0.
    Ncur = 0; // Number of patches assigned to current rank r.
    Lcur = 0.; //Load assigned to current rank r.

    int res_distributed( 0 );
    // MPI master loops patches and figures the best arrangement
    if( smilei_rk==0 ) {
        int rk = 0;
        MPI_Status status;
        while( true ) { // loop cpu ranks
            unsigned int hindex = 0;
            for( unsigned int ipatch=0; ipatch < Npatches_local; ipatch++ ) {
                local_load = PatchLoad[ipatch];
                Lcur += local_load; //Add grid contribution to the load.
                Ncur++; // Try to assign current patch to rank r.

                if( r < ( unsigned int )smilei_sz-1 ) {

                    if( Lcur > Tcur || smilei_sz-r >= Npatches-hindex ) { //Load target is exceeded or we have as many patches as procs left.
                        above_target = Lcur - Tcur;  //Including current patch, we exceed target by that much.
                        below_target = Tcur - ( Lcur-local_load ); // Excluding current patch, we mis the target by that much.
                        if( ( above_target > below_target ) && ( Ncur!=1 ) ) { // If we're closer to target without the current patch...
                            patch_count[r] = Ncur-1;      // ... include patches up to current one.
                            Ncur = 1;
                            //Lcur = local_load;
                        } else {                          //Else ...
                            patch_count[r] = Ncur;        //...assign patches including the current one.
                            Ncur = 0;
                            //Lcur = 0.;
                        }
                        res_distributed += patch_count[r];
                        if( Npatches - res_distributed <= smilei_sz-1-( r+1 ) ) {
                            // look for the last rank with more than one patch
                            int lrwmtop = r;
                            while( patch_count[lrwmtop]<=1 ) {
                                lrwmtop--;
                            }
                            patch_count[lrwmtop]--;
                            res_distributed--;
                            Ncur++;
                        }

                        r++; //Move on to the next rank.
                        //Tcur = Tload * capabilities[r];  //Target load for current rank r.
                        Tcur += Tload * capabilities[r];  //Target load for current rank r.
                    }
                }// End if on r.
                hindex++;
            }// End loop on patches for rank rk
            patch_count[smilei_sz-1] = Ncur; // the last MPI process takes what's left.

            // Go to next rank
            rk++;
            if( rk >= smilei_sz ) {
                break;
            }

            // Get the load of patches pre-calculated by the next rank
            if( rk == remainder ) {
                Npatches_local--;
                PatchLoad.resize( Npatches_local );
            }
            MPI_Recv( &PatchLoad[0], Npatches_local, MPI_DOUBLE, rk, rk, world_, &status );
        }

        // The master cpu also writes the patch count to the file
        ofstream fout;
        fout.open( "patch_load.txt" );
        fout << "Target load = " << Tload << endl;
        for( rk=0; rk<smilei_sz; rk++ ) {
            fout << "patch count = " << patch_count[rk]<<endl;
        }
        fout.close();

        // The other MPIs send their pre-calculated information
    } else {
        MPI_Send( &PatchLoad[0], Npatches_local, MPI_DOUBLE, 0, smilei_rk, world_ );
    }

    // Lastly, the patch count is broadcast to all ranks
    MPI_Bcast( &patch_count[0], smilei_sz, MPI_INT, 0, world_ );

    patch_refHindexes.resize( patch_count.size(), 0 );
    patch_refHindexes[0] = 0;
    for( int rk=1 ; rk<smilei_sz ; rk++ ) {
        patch_refHindexes[rk] = patch_refHindexes[rk-1] + patch_count[rk-1];
    }

} // END init_patch_count


// ---------------------------------------------------------------------------------------------------------------------
//  Recompute patch distribution
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI::recompute_patch_count( Params &params, VectorPatch &vecpatches, double time_dual )
{

    unsigned int ncells_perpatch, j;
    int Ncur;
    double Tload, Tload_loc, Tcur, cells_load, target, Tscan, largest_patch_loc, largest_patch;
    bool recompute_tload = true;
    //Load of a cell = cell_load*load of a particle.
    //Load of a frozen particle = frozen_particle_load*load of a particle.
    std::vector<double> Lp, Lp_left, Lp_right;
    ofstream fout;

    if( isMaster() ) {
        fout.open( "patch_load.txt", std::ofstream::out | std::ofstream::app );
    }

    MPI_Status status, status0, status1;
    MPI_Request request0, request1;

    ncells_perpatch = params.n_space[0]+2*params.oversize[0]; //Initialization
    for( unsigned int idim = 1; idim < params.nDim_field; idim++ ) {
        ncells_perpatch *= params.n_space[idim]+2*params.oversize[idim];
    }

    unsigned int tot_species_number = vecpatches( 0 )->vecSpecies.size();
    cells_load = ncells_perpatch*params.cell_load ;

    Lp.resize( patch_count[smilei_rk] );
    if( smilei_rk > 0 ) {
        Lp_left.resize( patch_count[smilei_rk-1] );
    }
    if( smilei_rk < smilei_sz-1 ) {
        Lp_right.resize( patch_count[smilei_rk+1] );
    }



    while( recompute_tload ) {

        Tload_loc = 0.;
        Ncur = 0; // Variation of the number of patches assigned to current rank r.
        for( unsigned int ipatch=0; ipatch < ( unsigned int )patch_count[smilei_rk]; ipatch++ ) {
            Lp[ipatch] =  cells_load ;
        }

        //Compute particle contribution to Local Loads of each Patch (Lp)
        for( unsigned int ipatch=0; ipatch < ( unsigned int )patch_count[smilei_rk]; ipatch++ ) {
            for( unsigned int ispecies = 0; ispecies < tot_species_number; ispecies++ ) {
                Lp[ipatch] += vecpatches( ipatch )->vecSpecies[ispecies]->getNbrOfParticles()*( 1+( params.frozen_particle_load-1 )*( time_dual < vecpatches( ipatch )->vecSpecies[ispecies]->time_frozen_ ) ) ;
            }
            Tload_loc += Lp[ipatch];
        }

        largest_patch_loc = *max_element( Lp.begin(), Lp.end() );

        //Tscan = total load carried by previous ranks and me
        MPI_Scan( &Tload_loc, &Tscan, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        //Tload = total load carried by all ranks
        MPI_Allreduce( &Tload_loc, &Tload, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        //Evaluate largest patch of the simulation
        MPI_Allreduce( &largest_patch_loc, &largest_patch, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );

        Tload /= Tcapabilities; //Target load for each mpi process.

        //This algorithm does not support single patches having a load larger than the target load per MPI rank.
        //If this happens, the code multiplies the cell load coefficient in order to be able to continue.
        if( largest_patch >= Tload ) {
            params.cell_load *= 2.;
            cells_load = ncells_perpatch*params.cell_load ;
            WARNING( "Dynamic Load balancing had to increase cell load coefficient because of an overloaded patch with respect to the target load per MPI rank. Try using smaller patches or less MPI ranks." );
        } else {
            recompute_tload = false;
        }
    }

    //Communicate the detail of the load of each patch to neighbouring MPI ranks
    if( smilei_rk < smilei_sz-1 ) {
        MPI_Isend( &( Lp[0] ), patch_count[smilei_rk], MPI_DOUBLE, smilei_rk+1, 0, MPI_COMM_WORLD, &request0 );
    }
    if( smilei_rk > 0 ) {
        MPI_Isend( &( Lp[0] ), patch_count[smilei_rk], MPI_DOUBLE, smilei_rk-1, 1, MPI_COMM_WORLD, &request1 );
        MPI_Recv( &( Lp_left[0] ), patch_count[smilei_rk-1], MPI_DOUBLE, smilei_rk-1, 0, MPI_COMM_WORLD, &status0 );
    }
    if( smilei_rk < smilei_sz-1 ) {
        MPI_Recv( &( Lp_right[0] ), patch_count[smilei_rk+1], MPI_DOUBLE, smilei_rk+1, 1, MPI_COMM_WORLD, &status1 );
    }


    if( smilei_rk > 0 ) {
        MPI_Wait( &request1, &status );
    }
    if( smilei_rk < smilei_sz-1 ) {
        MPI_Wait( &request0, &status );
    }

    if( smilei_rk > 0 ) {
        //Tcur is now initialized as the total load currently carried by previous ranks.
        Tcur = Tscan - Tload_loc;
        //Check if my rank should start with additional patches from left neighbour.
        target = smilei_rk*Tload; //target here points at the optimal begining for current rank
        if( Tcur > target ) {
            j = Lp_left.size()-1;
            while( abs( Tcur-target ) > abs( Tcur-Lp_left[j] - target ) && j>0 ) { //Leave at least 1 patch to my neighbour.
                Tcur -= Lp_left[j];
                j--;
                Ncur++;
            }
        } else {
            //  Check if some of my patches should be given to my left neighbour.
            j = 0;
            while( ( abs( Tcur-target ) > abs( Tcur+Lp[j]-target ) ) && ( j < ( unsigned int )patch_count[smilei_rk]-1 ) ) { //Keep at least 1 patch from my original set of patches
                Tcur += Lp[j];
                j++;
                Ncur --;
            }
        }
    }

    if( smilei_rk < smilei_sz-1 ) {
        //Tcur is now initialized as the total load carried by previous ranks + my load.
        Tcur = Tscan;
        target = ( smilei_rk+1 )*Tload;

        //Check if my rank should start with additional patches from right neighbour ...
        if( Tcur < target ) {
            j = 0;
            while( ( abs( Tcur-target ) > abs( Tcur+Lp_right[j] - target ) ) && ( j<( unsigned int )patch_count[smilei_rk+1] - 1 ) ) { //Leave at least 1 patch to my neighbour
                Tcur += Lp_right[j];
                j++;
                Ncur++;
            }

        } else {
            //  Check if some of my patches should be given to my right neighbour.
            j = patch_count[smilei_rk]-1;
            while( abs( Tcur-target ) > abs( Tcur-Lp[j]-target ) && j > 0 ) { //Keep at least 1 patch from my original set of patches
                Tcur -= Lp[j];
                j--;
                Ncur --;
            }
        }
    }

    //Ncur is the variation of number of patches owned by current rank.
    //Stores in Ncur the final patch count of this rank
    Ncur += patch_count[smilei_rk] ;

    //Ncur now has to be gathered to all as target_patch_count[smilei_rk]
    MPI_Allgather( &Ncur, 1, MPI_INT, &patch_count[0], 1, MPI_INT, MPI_COMM_WORLD );

    patch_refHindexes[0] = 0;
    for( int rk=1 ; rk<smilei_sz ; rk++ ) {
        patch_refHindexes[rk] = patch_refHindexes[rk-1] + patch_count[rk-1];
    }

    //Write patch_load.txt
    if( smilei_rk==0 ) {
        fout << "\tt = " << time_dual << endl;
        for( int irk=0; irk<smilei_sz; irk++ ) {
            fout << " patch_count[" << irk << "] = " << patch_count[irk] << endl;
        }
        fout.close();
    }

    return;

} // END recompute_patch_count


// ----------------------------------------------------------------------
// Returns the rank of the MPI process currently owning patch h.
// ----------------------------------------------------------------------
int SmileiMPI::hrank( int h )
{
    if( h == MPI_PROC_NULL ) {
        return MPI_PROC_NULL;
    }

    int patch_counter, rank;
    rank=0;
    patch_counter = patch_count[0];
    while( h >= patch_counter ) {
        rank++;
        patch_counter += patch_count[rank];
    }
    return rank;
} // END hrank


// ----------------------------------------------------------------------
// Create MPI type to exchange all particles properties of particles
// ----------------------------------------------------------------------
MPI_Datatype SmileiMPI::createMPIparticles( Particles *particles )
{
    int nbrOfProp = particles->double_prop_.size() + particles->short_prop_.size() + particles->uint64_prop_.size();

    MPI_Aint address[nbrOfProp];
    for( unsigned int iprop=0 ; iprop<particles->double_prop_.size() ; iprop++ ) {
        MPI_Get_address( &( ( *( particles->double_prop_[iprop] ) )[0] ), &( address[iprop] ) );
    }
    for( unsigned int iprop=0 ; iprop<particles->short_prop_.size() ; iprop++ ) {
        MPI_Get_address( &( ( *( particles->short_prop_[iprop] ) )[0] ), &( address[particles->double_prop_.size()+iprop] ) );
    }
    for( unsigned int iprop=0 ; iprop<particles->uint64_prop_.size() ; iprop++ ) {
        MPI_Get_address( &( ( *( particles->uint64_prop_[iprop] ) )[0] ), &( address[particles->double_prop_.size()+particles->short_prop_.size()+iprop] ) );
    }

    int nbr_parts[nbrOfProp];
    // number of elements per property
    for( int i=0 ; i<nbrOfProp ; i++ ) {
        nbr_parts[i] = particles->size();
    }

    MPI_Aint disp[nbrOfProp];
    // displacement between 2 properties
    disp[0] = 0;
    for( int i=1 ; i<nbrOfProp ; i++ ) {
        disp[i] = address[i] - address[0];
    }

    MPI_Datatype partDataType[nbrOfProp];
    // define MPI type of each property, default is DOUBLE
    for( unsigned int i=0 ; i<particles->double_prop_.size() ; i++ ) {
        partDataType[i] = MPI_DOUBLE;
    }
    for( unsigned int iprop=0 ; iprop<particles->short_prop_.size() ; iprop++ ) {
        partDataType[ particles->double_prop_.size()+iprop] = MPI_SHORT;
    }
    for( unsigned int iprop=0 ; iprop<particles->uint64_prop_.size() ; iprop++ ) {
        partDataType[ particles->double_prop_.size()+particles->short_prop_.size()+iprop] = MPI_UNSIGNED_LONG_LONG;
    }

    MPI_Datatype typeParticlesMPI;
    MPI_Type_create_struct( nbrOfProp, &( nbr_parts[0] ), &( disp[0] ), &( partDataType[0] ), &typeParticlesMPI );
    MPI_Type_commit( &typeParticlesMPI );

    return typeParticlesMPI;

} // END createMPIparticles


// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// -----------------------------------------       PATCH SEND / RECV METHODS        ------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI::isend( Patch *patch, int to, int tag, Params &params, bool send_xmax_bc )
{
    //MPI_Request request;

    // Count number max of comms :
    int irequest = 0;

    isend_species( patch, to, irequest, tag, params );
    isend_fields ( patch, to, irequest, tag, params, send_xmax_bc );

} // END isend( Patch )


void SmileiMPI::isend_species( Patch *patch, int to, int &irequest, int tag, Params &params )
{

    // number of species
    unsigned int nspec = patch->vecSpecies.size();

    // Adaptive vectorization:
    // In the case of the adaptive mixed sort Vectorization,
    // we communicate the operator state (vectorized_operators variable)
    // to deduce the bin number (particles->last_index.size())
    // In both adaptive cases :
    //   - a reconfiguration of operators is done after patch exchange (DLB and MW)
    //   - default values of the bin number is defined by the vectorized conf
    if( params.vectorization_mode == "adaptive_mixed_sort" ) {
        // Parameter vectorized_operators
        patch->buffer_vecto.resize( nspec );
        for( unsigned int ispec=0; ispec<nspec; ispec++ ) {
            patch->buffer_vecto[ispec] = patch->vecSpecies[ispec]->vectorized_operators;
        }
        MPI_Isend( &patch->buffer_vecto[0], nspec, MPI_INT, to, tag+irequest, MPI_COMM_WORLD, &patch->requests_[irequest] );
        irequest ++;
    }

    // For the particles
    for( unsigned int ispec=0; ispec<nspec; ispec++ ) {
        isend( &( patch->vecSpecies[ispec]->particles->last_index ), to, tag+irequest+2*ispec+1, patch->requests_[irequest+2*ispec] );
        if( patch->vecSpecies[ispec]->getNbrOfParticles() > 0 ) {
            patch->vecSpecies[ispec]->exchangePatch = createMPIparticles( patch->vecSpecies[ispec]->particles );
            isend( patch->vecSpecies[ispec]->particles, to, tag+irequest+2*ispec, patch->vecSpecies[ispec]->exchangePatch, patch->requests_[irequest+2*ispec+1] );
        }
    }
    irequest += 2*nspec;

    // Send some scalars
    unsigned int nscalars = 4 + ( params.hasMCRadiation || params.hasLLRadiation || params.hasNielRadiation );
    patch->buffer_scalars_particles.resize( nscalars*nspec );
    for( unsigned int ispec=0; ispec<nspec; ispec++ ) {
        unsigned int i = ispec*nscalars;
        patch->buffer_scalars_particles[i+0] = patch->vecSpecies[ispec]->nrj_bc_lost; // lost at boundaries
        patch->buffer_scalars_particles[i+1] = patch->vecSpecies[ispec]->nrj_new_part_; // injected
        patch->buffer_scalars_particles[i+2] = patch->vecSpecies[ispec]->nrj_mw_out; // lost by moving window
        patch->buffer_scalars_particles[i+3] = patch->vecSpecies[ispec]->nrj_mw_inj; // gained by moving window
        if( params.hasMCRadiation || params.hasLLRadiation || params.hasNielRadiation ) {
            patch->buffer_scalars_particles[i+4] = patch->vecSpecies[ispec]->nrj_radiated_; // radiated energy
        }
    }
    MPI_Isend( &patch->buffer_scalars_particles[0], patch->buffer_scalars_particles.size(), MPI_DOUBLE, to, tag + irequest, world_, &patch->requests_[irequest] );
    irequest ++;
}

void SmileiMPI::isend_fields( Patch *patch, int to, int &irequest, int tag, Params &params, bool send_xmax_bc )
{
    // Send fields
    if( params.geometry != "AMcylindrical" ) {
        isend( patch->EMfields, to, irequest, patch->requests_, tag, send_xmax_bc );
    } else {
        isend( patch->EMfields, to, irequest, patch->requests_, tag, static_cast<ElectroMagnAM *>( patch->EMfields )->El_.size(), send_xmax_bc );
    }

    // Send some scalars
    unsigned int nscalars = 2 + 2*params.nDim_field;
    patch->buffer_scalars_fields.resize( nscalars );
    patch->buffer_scalars_fields[0] = patch->EMfields->nrj_mw_out; // lost by moving window
    patch->buffer_scalars_fields[1] = patch->EMfields->nrj_mw_inj; // lost by moving window
    for( unsigned int jp=0; jp<2; jp++ ) { //directions (xmin/xmax, ymin/ymax, zmin/zmax)
        for( unsigned int i=0 ; i<params.nDim_field ; i++ ) { //axis 0=x, 1=y, 2=z
            patch->buffer_scalars_fields[2+i*2+jp] = patch->EMfields->poynting[jp][i];
        }
    }
    MPI_Isend( &patch->buffer_scalars_fields[0], patch->buffer_scalars_fields.size(), MPI_DOUBLE, to, tag + irequest, world_, &patch->requests_[irequest] );
    irequest ++;
} // END isend( Patch )


void SmileiMPI::waitall( Patch *patch )
{

    for( unsigned int ireq=0; ireq<patch->requests_.size() ; ireq++ ) {
        MPI_Status status;
        if( patch->requests_[ireq] != MPI_REQUEST_NULL ) {
            MPI_Wait( &( patch->requests_[ireq] ), &status );
        }
        //AB: This operation is done in MPI_Wait already.
        //patch->requests_[ireq] = MPI_REQUEST_NULL;
    }

    for( int ispec=0 ; ispec<( int )patch->vecSpecies.size() ; ispec++ ) {
        if( patch->vecSpecies[ispec]->getNbrOfParticles() > 0 ) {
            if( patch->vecSpecies[ispec]->exchangePatch != MPI_DATATYPE_NULL ) {
                MPI_Type_free( &( patch->vecSpecies[ispec]->exchangePatch ) );
                patch->vecSpecies[ispec]->exchangePatch = MPI_DATATYPE_NULL;
            }
        }
    }

}

void SmileiMPI::recv( Patch *patch, int from, int tag, Params &params, bool recv_xmin_bc )
{
    // Receive species
    recv_species( patch, from, tag, params );

    // Receive EM fields
    recv_fields( patch, from, tag, params, recv_xmin_bc );

} // END recv ( Patch )


void SmileiMPI::recv_species( Patch *patch, int from, int &tag, Params &params )
{
    MPI_Datatype recvParts;
    int nbrOfPartsRecv;

    // number of species
    unsigned int nspec = patch->vecSpecies.size();

    // Adaptive vectorization:
    // In the case of the adaptive mixed sort Vectorization,
    // we communicate the operator state (vectorized_operators variable)
    // to deduce the bin number (particles->last_index.size())
    // In both adaptive cases :
    //   - a reconfiguration of operators is done after patch exchange (DLB and MW)
    //   - default values of the bin number is defined by the vectorized conf
    if( params.vectorization_mode == "adaptive_mixed_sort" ) {
        // Parameter vectorized_operators
        MPI_Status status;
        patch->buffer_vecto.resize( nspec );
        MPI_Recv( &patch->buffer_vecto[0], nspec, MPI_INT, from, tag, MPI_COMM_WORLD, &status );
        tag ++;
        for( unsigned int ispec=0; ispec<nspec; ispec++ ) {
            patch->vecSpecies[ispec]->vectorized_operators = patch->buffer_vecto[ispec];
            if( ! patch->buffer_vecto[ispec] ) {
                patch->vecSpecies[ispec]->particles->last_index.resize( 1 );
                patch->vecSpecies[ispec]->particles->first_index.resize( 1 );
            }
        }
    }

    for( unsigned int ispec=0; ispec<nspec; ispec++ ) {
        //Receive last_index
        recv( &patch->vecSpecies[ispec]->particles->last_index, from, tag+2*ispec+1 );
        //Reconstruct first_index from last_index
        memcpy( &( patch->vecSpecies[ispec]->particles->first_index[1] ), &( patch->vecSpecies[ispec]->particles->last_index[0] ), ( patch->vecSpecies[ispec]->particles->last_index.size()-1 )*sizeof( int ) );
        patch->vecSpecies[ispec]->particles->first_index[0]=0;
        //Prepare patch for receiving particles
        nbrOfPartsRecv = patch->vecSpecies[ispec]->particles->last_index.back();
        patch->vecSpecies[ispec]->particles->initialize( nbrOfPartsRecv, params.nDim_particle, params.keep_position_old );
        //Receive particles
        if( nbrOfPartsRecv > 0 ) {
            recvParts = createMPIparticles( patch->vecSpecies[ispec]->particles );
            recv( patch->vecSpecies[ispec]->particles, from, tag+2*ispec, recvParts );
            MPI_Type_free( &( recvParts ) );
        }
        /*std::cerr << "Species: " << ispec
                  << " particles->last_index: " <<  patch->vecSpecies[ispec]->particles->last_index[0]
                  << " Number of particles: " << patch->vecSpecies[ispec]->particles->size() <<'\n';*/
    }
    tag += 2*nspec;

    // Receive some scalars
    unsigned int nscalars = 4 + ( params.hasMCRadiation || params.hasLLRadiation || params.hasNielRadiation );
    patch->buffer_scalars_particles.resize( nscalars*nspec );
    MPI_Status status;
    MPI_Recv( &patch->buffer_scalars_particles[0], patch->buffer_scalars_particles.size(), MPI_DOUBLE, from, tag, world_, &status );
    tag++;
    for( unsigned int ispec=0; ispec<nspec; ispec++ ) {
        unsigned int i = ispec*nscalars;
        patch->vecSpecies[ispec]->nrj_bc_lost   = patch->buffer_scalars_particles[i+0];
        patch->vecSpecies[ispec]->nrj_new_part_ = patch->buffer_scalars_particles[i+1];
        patch->vecSpecies[ispec]->nrj_mw_out    = patch->buffer_scalars_particles[i+2];
        patch->vecSpecies[ispec]->nrj_mw_inj    = patch->buffer_scalars_particles[i+3];
        if( params.hasMCRadiation || params.hasLLRadiation || params.hasNielRadiation ) {
            patch->vecSpecies[ispec]->nrj_radiated_ = patch->buffer_scalars_particles[i+4];
        }
    }
}

void SmileiMPI::recv_fields( Patch *patch, int from, int &tag, Params &params, bool recv_xmin_bc )
{
    // Receive EM fields
    patch->EMfields->initAntennas( patch, params );
    if( params.geometry != "AMcylindrical" ) {
        recv( patch->EMfields, from, tag, recv_xmin_bc );
    } else {
        recv( patch->EMfields, from, tag, static_cast<ElectroMagnAM *>( patch->EMfields )->El_.size(), recv_xmin_bc );
    }

    // Receive some scalars
    unsigned int nscalars = 2 + 2*params.nDim_field;
    patch->buffer_scalars_fields.resize( nscalars );
    MPI_Status status;
    MPI_Recv( &patch->buffer_scalars_fields[0], patch->buffer_scalars_fields.size(), MPI_DOUBLE, from, tag, world_, &status );
    tag++;
    patch->EMfields->nrj_mw_out =  patch->buffer_scalars_fields[0] ;
    patch->EMfields->nrj_mw_inj =  patch->buffer_scalars_fields[1] ;
    for( unsigned int jp=0; jp<2; jp++ ) { //directions (xmin/xmax, ymin/ymax, zmin/zmax)
        for( unsigned int i=0 ; i<params.nDim_field ; i++ ) { //axis 0=x, 1=y, 2=z
            patch->EMfields->poynting[jp][i] = patch->buffer_scalars_fields[2+i*2+jp];
        }
    }
} // END recv ( Patch )


void SmileiMPI::isend( Particles *particles, int to, int tag, MPI_Datatype typePartSend, MPI_Request &request )
{
    MPI_Isend( &( particles->position( 0, 0 ) ), 1, typePartSend, to, tag, MPI_COMM_WORLD, &request );

} // END isend( Particles )


void SmileiMPI::recv( Particles *particles, int to, int tag, MPI_Datatype typePartRecv )
{
    MPI_Status status;
    MPI_Recv( &( particles->position( 0, 0 ) ), 1, typePartRecv, to, tag, MPI_COMM_WORLD, &status );

} // END recv( Particles )


// Assuming vec.size() is known (number of species). Asynchronous.
void SmileiMPI::isend( std::vector<int> *vec, int to, int tag, MPI_Request &request )
{
    MPI_Isend( &( ( *vec )[0] ), vec->size(), MPI_INT, to, tag, MPI_COMM_WORLD, &request );

} // End isend

void SmileiMPI::recv( std::vector<int> *vec, int from, int tag )
{
    MPI_Status status;
    MPI_Recv( &( ( *vec )[0] ), vec->size(), MPI_INT, from, tag, MPI_COMM_WORLD, &status );

} // End recv

// Assuming vec.size() is known (number of species). Asynchronous.
void SmileiMPI::isend( std::vector<double> *vec, int to, int tag, MPI_Request &request )
{
    MPI_Isend( &( ( *vec )[0] ), vec->size(), MPI_DOUBLE, to, tag, MPI_COMM_WORLD, &request );

} // End isend

void SmileiMPI::recv( std::vector<double> *vec, int from, int tag )
{
    MPI_Status status;
    MPI_Recv( &( ( *vec )[0] ), vec->size(), MPI_DOUBLE, from, tag, MPI_COMM_WORLD, &status );

} // End recv

template <typename Tpml> //ElectroMagnBC2D_PML or ElectroMagnBC3D_PML
void  SmileiMPI::send_PML(ElectroMagn *EM, Tpml embc, int bcId, int to, int &irequest, vector<MPI_Request> &requests, int tag, bool send_xmax_bc){
    if (embc->Hx_) {
        if(!send_xmax_bc && bcId>1 && EM->isXmax){ //When not sending xmax_bc (MovingWindow), the size of the non longitudinal pml sent must be tailored on corner cases.

            isend( embc->Hx_, to, tag+irequest, requests[irequest], (EM->dimPrim[0]+0));// * (ncells_pml_domain_yminmax+1) );
            irequest++;
            isend( embc->Hy_, to, tag+irequest, requests[irequest], (EM->dimPrim[0]+1));// * (ncells_pml_domain_yminmax+0)  );
            irequest++;
            isend( embc->Hz_, to, tag+irequest, requests[irequest], (EM->dimPrim[0]+1));// * (ncells_pml_domain_yminmax+1)  );
            irequest++;
            isend( embc->Bx_, to, tag+irequest, requests[irequest], (EM->dimPrim[0]+0));// * (ncells_pml_domain_yminmax+1)  );
            irequest++;
            isend( embc->By_, to, tag+irequest, requests[irequest], (EM->dimPrim[0]+1));// * (ncells_pml_domain_yminmax+0)  );
            irequest++;
            isend( embc->Bz_, to, tag+irequest, requests[irequest], (EM->dimPrim[0]+1));// * (ncells_pml_domain_yminmax+1)  );
            irequest++;
            isend( embc->Ex_, to, tag+irequest, requests[irequest], (EM->dimPrim[0]+1));// * (ncells_pml_domain_yminmax+0)  );
            irequest++;
            isend( embc->Ey_, to, tag+irequest, requests[irequest], (EM->dimPrim[0]+0));// * (ncells_pml_domain_yminmax+1)  );
            irequest++;
            isend( embc->Ez_, to, tag+irequest, requests[irequest], (EM->dimPrim[0]+0));// * (ncells_pml_domain_yminmax+0)  );
            irequest++;
            isend( embc->Dx_, to, tag+irequest, requests[irequest], (EM->dimPrim[0]+1));// * (ncells_pml_domain_yminmax+0)  );
            irequest++;
            isend( embc->Dy_, to, tag+irequest, requests[irequest], (EM->dimPrim[0]+0));// * (ncells_pml_domain_yminmax+1)  );
            irequest++;
            isend( embc->Dz_, to, tag+irequest, requests[irequest], (EM->dimPrim[0]+0));// * (ncells_pml_domain_yminmax+0)  );
            irequest++;
        } else {
            isend( embc->Hx_, to, tag + irequest, requests[irequest] );
            irequest++;
            isend( embc->Hy_, to, tag + irequest, requests[irequest] );
            irequest++;
            isend( embc->Hz_, to, tag + irequest, requests[irequest] );
            irequest++;
            isend( embc->Bx_, to, tag + irequest, requests[irequest] );
            irequest++;
            isend( embc->By_, to, tag + irequest, requests[irequest] );
            irequest++;
            isend( embc->Bz_, to, tag + irequest, requests[irequest] );
            irequest++;
            isend( embc->Ex_, to, tag + irequest, requests[irequest] );
            irequest++;
            isend( embc->Ey_, to, tag + irequest, requests[irequest] );
            irequest++;
            isend( embc->Ez_, to, tag + irequest, requests[irequest] );
            irequest++;
            isend( embc->Dx_, to, tag + irequest, requests[irequest] );
            irequest++;
            isend( embc->Dy_, to, tag + irequest, requests[irequest] );
            irequest++;
            isend( embc->Dz_, to, tag + irequest, requests[irequest] );
            irequest++;
        }
    }
    // Envelope communication
    if (EM->envelope!=NULL) {
        if (static_cast<EnvelopeBC2D_PML *>( EM->envelope->EnvBoundCond[bcId] )){
            EnvelopeBC2D_PML *embcenv = static_cast<EnvelopeBC2D_PML *>( EM->envelope->EnvBoundCond[bcId] );

            if (embcenv->A_n_ ) {
                if(!send_xmax_bc && bcId>1 && EM->isXmax){ //When not sending xmax_bc (MovingWindow), the size of the non longitudinal pml sent must be tailored on corner cases.
                    isendComplex( embcenv->A_n_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                    irequest++;
                    isendComplex( embcenv->A_nm1_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                    irequest++;
                    isendComplex( embcenv->u1_nm1_x_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                    irequest++;
                    isendComplex( embcenv->u2_nm1_x_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                    irequest++;
                    isendComplex( embcenv->u3_nm1_x_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                    irequest++;
                    isendComplex( embcenv->u1_nm1_y_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                    irequest++;
                    isendComplex( embcenv->u2_nm1_y_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                    irequest++;
                    isendComplex( embcenv->u3_nm1_y_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                    irequest++;
                } else {
                    isendComplex( embcenv->A_n_, to, tag+irequest, requests[irequest] );
                    irequest++;
                    isendComplex( embcenv->A_nm1_, to, tag+irequest, requests[irequest] );
                    irequest++;
                    isendComplex( embcenv->u1_nm1_x_, to, tag+irequest, requests[irequest] );
                    irequest++;
                    isendComplex( embcenv->u2_nm1_x_, to, tag+irequest, requests[irequest] );
                    irequest++;
                    isendComplex( embcenv->u3_nm1_x_, to, tag+irequest, requests[irequest] );
                    irequest++;
                    if (bcId > 1) { //Sending Transverse PML
                        isendComplex( embcenv->u1_nm1_y_, to, tag+irequest, requests[irequest] );
                        irequest++;
                        isendComplex( embcenv->u2_nm1_y_, to, tag+irequest, requests[irequest] );
                        irequest++;
                        isendComplex( embcenv->u3_nm1_y_, to, tag+irequest, requests[irequest] );
                        irequest++;
                    }
                }
            }
        } else if (static_cast<EnvelopeBC3D_PML *>( EM->envelope->EnvBoundCond[bcId] )){ 
            EnvelopeBC3D_PML *embcenv = static_cast<EnvelopeBC3D_PML *>( EM->envelope->EnvBoundCond[bcId] );
            if (embcenv->A_n_ ) {
                if(!send_xmax_bc && bcId>1 && EM->isXmax){ //When not sending xmax_bc (MovingWindow), the size of the non longitudinal pml sent must be tailored on corner cases.
                    isendComplex( embcenv->A_n_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                    irequest++;
                    isendComplex( embcenv->A_nm1_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                    irequest++;
                    isendComplex( embcenv->u1_nm1_x_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                    irequest++;
                    isendComplex( embcenv->u2_nm1_x_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                    irequest++;
                    isendComplex( embcenv->u3_nm1_x_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                    irequest++;
                    isendComplex( embcenv->u1_nm1_y_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                    irequest++;
                    isendComplex( embcenv->u2_nm1_y_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                    irequest++;
                    isendComplex( embcenv->u3_nm1_y_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                    irequest++;
                    if (bcId > 3){
                        isendComplex( embcenv->u1_nm1_z_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                        irequest++;
                        isendComplex( embcenv->u2_nm1_z_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                        irequest++;
                        isendComplex( embcenv->u3_nm1_z_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                        irequest++;
                    }
                } else {
                    isendComplex( embcenv->A_n_, to, tag+irequest, requests[irequest] );
                    irequest++;
                    isendComplex( embcenv->A_nm1_, to, tag+irequest, requests[irequest] );
                    irequest++;
                    isendComplex( embcenv->u1_nm1_x_, to, tag+irequest, requests[irequest] );
                    irequest++;
                    isendComplex( embcenv->u2_nm1_x_, to, tag+irequest, requests[irequest] );
                    irequest++;
                    isendComplex( embcenv->u3_nm1_x_, to, tag+irequest, requests[irequest] );
                    irequest++;
                    if (bcId > 1 ) { //Sending Transverse PML
                        isendComplex( embcenv->u1_nm1_y_, to, tag+irequest, requests[irequest] );
                        irequest++;
                        isendComplex( embcenv->u2_nm1_y_, to, tag+irequest, requests[irequest] );
                        irequest++;
                        isendComplex( embcenv->u3_nm1_y_, to, tag+irequest, requests[irequest] );
                        irequest++;
                        if (bcId > 3) { //Sending Transverse PML
                            isendComplex( embcenv->u1_nm1_z_, to, tag+irequest, requests[irequest] );
                            irequest++;
                            isendComplex( embcenv->u2_nm1_z_, to, tag+irequest, requests[irequest] );
                            irequest++;
                            isendComplex( embcenv->u3_nm1_z_, to, tag+irequest, requests[irequest] );
                            irequest++;
                        }
                    }
                }
            }
        }
    }
}

void SmileiMPI::isend( ElectroMagn *EM, int to, int &irequest, vector<MPI_Request> &requests, int tag, bool send_xmax_bc )
{

    isend( EM->Ex_, to, tag+irequest, requests[irequest] );
    irequest++;
    isend( EM->Ey_, to, tag+irequest, requests[irequest] );
    irequest++;
    isend( EM->Ez_, to, tag+irequest, requests[irequest] );
    irequest++;
    isend( EM->Bx_, to, tag+irequest, requests[irequest] );
    irequest++;
    isend( EM->By_, to, tag+irequest, requests[irequest] );
    irequest++;
    isend( EM->Bz_, to, tag+irequest, requests[irequest] );
    irequest++;
    isend( EM->Bx_m, to, tag+irequest, requests[irequest] );
    irequest++;
    isend( EM->By_m, to, tag+irequest, requests[irequest] );
    irequest++;
    isend( EM->Bz_m, to, tag+irequest, requests[irequest] );
    irequest++;

    // if laser envelope is present, send it
    // send also Phi, Phi_m, GradPhi, GradPhi_m
    if( EM->envelope!=NULL ) {
        isendComplex( EM->envelope->A_, to, tag+irequest, requests[irequest] );
        irequest++;
        isendComplex( EM->envelope->A0_, to, tag+irequest, requests[irequest] );
        irequest++;
        isend( EM->envelope->Phi_, to, tag+irequest, requests[irequest] );
        irequest++;
        isend( EM->envelope->Phi_m, to, tag+irequest, requests[irequest] );
        irequest++;
        isend( EM->envelope->GradPhix_, to, tag+irequest, requests[irequest] );
        irequest++;
        isend( EM->envelope->GradPhix_m, to, tag+irequest, requests[irequest] );
        irequest++;
        isend( EM->envelope->GradPhiy_, to, tag+irequest, requests[irequest] );
        irequest++;
        isend( EM->envelope->GradPhiy_m, to, tag+irequest, requests[irequest] );
        irequest++;
        isend( EM->envelope->GradPhiz_, to, tag+irequest, requests[irequest] );
        irequest++;
        isend( EM->envelope->GradPhiz_m, to, tag+irequest, requests[irequest] );
        irequest++;

    }

    for( unsigned int idiag=0; idiag<EM->allFields_avg.size(); idiag++ ) {
        for( unsigned int ifield=0; ifield<EM->allFields_avg[idiag].size(); ifield++ ) {
            isend( EM->allFields_avg[idiag][ifield], to, tag+irequest, requests[irequest] );
            irequest++;
        }
    }

    for( unsigned int antennaId=0 ; antennaId<EM->antennas.size() ; antennaId++ ) {
        isend( EM->antennas[antennaId].field, to, tag+irequest, requests[irequest] );
        irequest++;
    }

    for( unsigned int bcId=0 ; bcId<EM->emBoundCond.size() ; bcId++ ) {
        if( ! EM->emBoundCond[bcId] ) {
            continue;
        }

        for( unsigned int laserId=0 ; laserId < EM->emBoundCond[bcId]->vecLaser.size() ; laserId++ ) {

            Laser *laser = EM->emBoundCond[bcId]->vecLaser[laserId];
            if( !( laser->spacetime[0] ) && !( laser->spacetime[1] ) ) {
                LaserProfileSeparable *profile;
                profile = static_cast<LaserProfileSeparable *>( laser->profiles[0] );
                if( ! profile->space_envelope ) {
                    continue;
                }
                isend( profile->space_envelope, to, tag+irequest, requests[irequest] );
                irequest++;
                isend( profile->phase, to, tag+irequest, requests[irequest] );
                irequest++;
                profile = static_cast<LaserProfileSeparable *>( laser->profiles[1] );
                isend( profile->space_envelope, to, tag+irequest, requests[irequest] );
                irequest++;
                isend( profile->phase, to, tag+irequest, requests[irequest] );
                irequest++;
            }
        }

        if( EM->extFields.size()>0 ) {

            if( dynamic_cast<ElectroMagnBC1D_SM *>( EM->emBoundCond[bcId] ) ) {
                ElectroMagnBC1D_SM *embc = static_cast<ElectroMagnBC1D_SM *>( EM->emBoundCond[bcId] );
                MPI_Isend( &( embc->By_val ), 1, MPI_DOUBLE, to, tag+irequest, MPI_COMM_WORLD, &requests[irequest] );
                irequest++;
                MPI_Isend( &( embc->Bz_val ), 1, MPI_DOUBLE, to, tag+irequest, MPI_COMM_WORLD, &requests[irequest] );
                irequest++;
            } else if( dynamic_cast<ElectroMagnBC2D_SM *>( EM->emBoundCond[bcId] ) ) {
                // BCs at the x-border
                ElectroMagnBC2D_SM *embc = static_cast<ElectroMagnBC2D_SM *>( EM->emBoundCond[bcId] );

                if( embc->B_val[0].size() ) {
                    isend( &embc->B_val[0], to, tag+irequest, requests[irequest] );
                    irequest++;
                }
                if( embc->B_val[1].size() ) {
                    isend( &embc->B_val[1], to, tag+irequest, requests[irequest] );
                    irequest++;
                }
                if( embc->B_val[2].size() ) {
                    isend( &embc->B_val[2], to, tag+irequest, requests[irequest] );
                    irequest++;
                }

            } else if( dynamic_cast<ElectroMagnBC3D_SM *>( EM->emBoundCond[bcId] ) ) {
                ElectroMagnBC3D_SM *embc = static_cast<ElectroMagnBC3D_SM *>( EM->emBoundCond[bcId] );

                // BCs at the border
                if( embc->B_val[0] ) {
                    isend( embc->B_val[0], to, tag+irequest, requests[irequest] );
                    irequest++;
                }
                if( embc->B_val[1] ) {
                    isend( embc->B_val[1], to, tag+irequest, requests[irequest] );
                    irequest++;
                }
                if( embc->B_val[2] ) {
                    isend( embc->B_val[2], to, tag+irequest, requests[irequest] );
                    irequest++;
                }

            }
        }
        if( (dynamic_cast<ElectroMagnBC2D_PML *>( EM->emBoundCond[bcId] ) || dynamic_cast<ElectroMagnBC3D_PML *>( EM->emBoundCond[bcId] )) && (bcId != 1 || send_xmax_bc) ){
            if( dynamic_cast<ElectroMagnBC2D_PML *>( EM->emBoundCond[bcId] )){
                ElectroMagnBC2D_PML *embc = static_cast<ElectroMagnBC2D_PML *>( EM->emBoundCond[bcId] );
                send_PML(EM, embc, bcId, to, irequest,requests, tag, send_xmax_bc);
            } else {
                ElectroMagnBC3D_PML *embc = static_cast<ElectroMagnBC3D_PML *>( EM->emBoundCond[bcId] );
                send_PML(EM, embc, bcId, to, irequest,requests, tag, send_xmax_bc);
            }
        }
    }
} // End isend ( ElectroMagn )

void SmileiMPI::isend( ElectroMagn *EM, int to, int &irequest, vector<MPI_Request> &requests, int tag, unsigned int nmodes, bool send_xmax_bc )
{

    ElectroMagnAM *EMAM = static_cast<ElectroMagnAM *>( EM );
    for( unsigned int imode =0; imode < nmodes; imode++ ) {
        isendComplex( EMAM->El_[imode], to, tag+irequest, requests[irequest] );
        irequest++;
        isendComplex( EMAM->Er_[imode], to, tag+irequest, requests[irequest] );
        irequest++;
        isendComplex( EMAM->Et_[imode], to, tag+irequest, requests[irequest] );
        irequest++;
        isendComplex( EMAM->Bl_[imode], to, tag+irequest, requests[irequest] );
        irequest++;
        isendComplex( EMAM->Br_[imode], to, tag+irequest, requests[irequest] );
        irequest++;
        isendComplex( EMAM->Bt_[imode], to, tag+irequest, requests[irequest] );
        irequest++;
        isendComplex( EMAM->Bl_m[imode], to, tag+irequest, requests[irequest] );
        irequest++;
        isendComplex( EMAM->Br_m[imode], to, tag+irequest, requests[irequest] );
        irequest++;
        isendComplex( EMAM->Bt_m[imode], to, tag+irequest, requests[irequest] );
        irequest++;
    }

    // if laser envelope is present, send it
    // send also Phi, Phi_m, GradPhi, GradPhi_m
    if( EM->envelope!=NULL ) {
        isendComplex( EM->envelope->A_, to, tag+irequest, requests[irequest] );
        irequest++;
        isendComplex( EM->envelope->A0_, to, tag+irequest, requests[irequest] );
        irequest++;
        isend( EM->envelope->Phi_, to, tag+irequest, requests[irequest] );
        irequest++;
        isend( EM->envelope->Phi_m, to, tag+irequest, requests[irequest] );
        irequest++;
        isend( EM->envelope->GradPhil_, to, tag+irequest, requests[irequest] );
        irequest++;
        isend( EM->envelope->GradPhil_m, to, tag+irequest, requests[irequest] );
        irequest++;
        isend( EM->envelope->GradPhir_, to, tag+irequest, requests[irequest] );
        irequest++;
        isend( EM->envelope->GradPhir_m, to, tag+irequest, requests[irequest] );
        irequest++;

    }

    for( unsigned int idiag=0; idiag<EM->allFields_avg.size(); idiag++ ) {
        for( unsigned int ifield=0; ifield<EM->allFields_avg[idiag].size(); ifield++ ) {
            isend( EM->allFields_avg[idiag][ifield], to, tag+irequest, requests[irequest] );
            irequest++;
        }
    }

    for( unsigned int antennaId=0 ; antennaId<EM->antennas.size() ; antennaId++ ) {
        isend( EM->antennas[antennaId].field, to, tag+irequest, requests[irequest] );
        irequest++;
    }

    for( unsigned int bcId=0 ; bcId<EM->emBoundCond.size() ; bcId++ ) {
        if( ! EM->emBoundCond[bcId] ) {
            continue;
        }

        for( unsigned int laserId=0 ; laserId < EM->emBoundCond[bcId]->vecLaser.size() ; laserId++ ) {

            Laser *laser = EM->emBoundCond[bcId]->vecLaser[laserId];
            if( !( laser->spacetime[0] ) && !( laser->spacetime[1] ) ) {
                LaserProfileSeparable *profile;
                profile = static_cast<LaserProfileSeparable *>( laser->profiles[0] );
                if( ! profile->space_envelope ) {
                    continue;
                }
                isend( profile->space_envelope, to, tag+irequest, requests[irequest] );
                irequest++;
                isend( profile->phase, to, tag+irequest, requests[irequest] );
                irequest++;
                profile = static_cast<LaserProfileSeparable *>( laser->profiles[1] );
                isend( profile->space_envelope, to, tag+irequest, requests[irequest] );
                irequest++;
                isend( profile->phase, to, tag+irequest, requests[irequest] );
                irequest++;
            }
        }


        if( dynamic_cast<ElectroMagnBCAM_PML *>( EM->emBoundCond[bcId] ) && (bcId != 1 || send_xmax_bc)  ){
            ElectroMagnBCAM_PML *embc = static_cast<ElectroMagnBCAM_PML *>( EM->emBoundCond[bcId] );
            // if I have PML && if the receiver also have PMLs <=> the hindex I send to touches the same boundary I am dealing with now
            if (embc->Hl_[0] ) {
                if(!send_xmax_bc && bcId>1 && EMAM->isXmax){ //When not sending xmax_bc (MovingWindow), the size of the non longitudinal pml sent must be tailored on corner cases.
                    for( unsigned int imode =0; imode < nmodes; imode++ ) {
                        isendComplex( embc->Hl_[imode], to, tag+irequest, requests[irequest], (EM->dimPrim[0]+0) );
                        irequest++;
                        isendComplex( embc->Hr_[imode], to, tag+irequest, requests[irequest], (EM->dimPrim[0]+1) );
                        irequest++;
                        isendComplex( embc->Ht_[imode], to, tag+irequest, requests[irequest], (EM->dimPrim[0]+1) );
                        irequest++;
                        isendComplex( embc->Bl_[imode], to, tag+irequest, requests[irequest], (EM->dimPrim[0]+0) );
                        irequest++;
                        isendComplex( embc->Br_[imode], to, tag+irequest, requests[irequest], (EM->dimPrim[0]+1) );
                        irequest++;
                        isendComplex( embc->Bt_[imode], to, tag+irequest, requests[irequest], (EM->dimPrim[0]+1) );
                        irequest++;
                        isendComplex( embc->El_[imode], to, tag+irequest, requests[irequest], (EM->dimPrim[0]+1) );
                        irequest++;
                        isendComplex( embc->Er_[imode], to, tag+irequest, requests[irequest], (EM->dimPrim[0]+0) );
                        irequest++;
                        isendComplex( embc->Et_[imode], to, tag+irequest, requests[irequest], (EM->dimPrim[0]+0) );
                        irequest++;
                        isendComplex( embc->Dl_[imode], to, tag+irequest, requests[irequest], (EM->dimPrim[0]+1) );
                        irequest++;
                        isendComplex( embc->Dr_[imode], to, tag+irequest, requests[irequest], (EM->dimPrim[0]+0) );
                        irequest++;
                        isendComplex( embc->Dt_[imode], to, tag+irequest, requests[irequest], (EM->dimPrim[0]+0) );
                        irequest++;
                    }
                } else {
                    for( unsigned int imode =0; imode < nmodes; imode++ ) {
                        isendComplex( embc->Hl_[imode], to, tag+irequest, requests[irequest] );
                        irequest++;
                        isendComplex( embc->Hr_[imode], to, tag+irequest, requests[irequest] );
                        irequest++;
                        isendComplex( embc->Ht_[imode], to, tag+irequest, requests[irequest] );
                        irequest++;
                        isendComplex( embc->Bl_[imode], to, tag+irequest, requests[irequest] );
                        irequest++;
                        isendComplex( embc->Br_[imode], to, tag+irequest, requests[irequest] );
                        irequest++;
                        isendComplex( embc->Bt_[imode], to, tag+irequest, requests[irequest] );
                        irequest++;
                        isendComplex( embc->El_[imode], to, tag+irequest, requests[irequest] );
                        irequest++;
                        isendComplex( embc->Er_[imode], to, tag+irequest, requests[irequest] );
                        irequest++;
                        isendComplex( embc->Et_[imode], to, tag+irequest, requests[irequest] );
                        irequest++;
                        isendComplex( embc->Dl_[imode], to, tag+irequest, requests[irequest] );
                        irequest++;
                        isendComplex( embc->Dr_[imode], to, tag+irequest, requests[irequest] );
                        irequest++;
                        isendComplex( embc->Dt_[imode], to, tag+irequest, requests[irequest] );
                        irequest++;
                    }
                }
                if (EM->envelope!=NULL){
                    EnvelopeBCAM_PML *embcenv = static_cast<EnvelopeBCAM_PML *>( EM->envelope->EnvBoundCond[bcId] );
                    if (embcenv->A_n_ ) {
                        if(!send_xmax_bc && bcId==3 && EMAM->isXmax){ //When not sending xmax_bc (MovingWindow), the size of the non longitudinal pml sent must be tailored on corner cases.
                            isendComplex( embcenv->A_n_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                            irequest++;
                            isendComplex( embcenv->A_nm1_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                            irequest++;
                            isendComplex( embcenv->G_n_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                            irequest++;
                            isendComplex( embcenv->G_nm1_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                            irequest++;
                            isendComplex( embcenv->u1_nm1_l_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                            irequest++;
                            isendComplex( embcenv->u2_nm1_l_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                            irequest++;
                            isendComplex( embcenv->u3_nm1_l_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                            irequest++;
                            isendComplex( embcenv->u1_nm1_r_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                            irequest++;
                            isendComplex( embcenv->u2_nm1_r_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                            irequest++;
                            isendComplex( embcenv->u3_nm1_r_, to, tag+irequest, requests[irequest], EM->dimPrim[0] );
                            irequest++;
                        } else {
                            isendComplex( embcenv->A_n_, to, tag+irequest, requests[irequest] );
                            irequest++;
                            isendComplex( embcenv->A_nm1_, to, tag+irequest, requests[irequest] );
                            irequest++;
                            isendComplex( embcenv->G_n_, to, tag+irequest, requests[irequest] );
                            irequest++;
                            isendComplex( embcenv->G_nm1_, to, tag+irequest, requests[irequest] );
                            irequest++;
                            isendComplex( embcenv->u1_nm1_l_, to, tag+irequest, requests[irequest] );
                            irequest++;
                            isendComplex( embcenv->u2_nm1_l_, to, tag+irequest, requests[irequest] );
                            irequest++;
                            isendComplex( embcenv->u3_nm1_l_, to, tag+irequest, requests[irequest] );
                            irequest++;
                            if (bcId == 3) { //Sending Radial PML
                                isendComplex( embcenv->u1_nm1_r_, to, tag+irequest, requests[irequest] );
                                irequest++;
                                isendComplex( embcenv->u2_nm1_r_, to, tag+irequest, requests[irequest] );
                                irequest++;
                                isendComplex( embcenv->u3_nm1_r_, to, tag+irequest, requests[irequest] );
                                irequest++;
                            }
                        }
                    }
                }
            }
        }
    }
} // End isend ( ElectroMagn LRT )

template <typename Tpml>
int  SmileiMPI::recv_PML(ElectroMagn *EM, Tpml embc, int bcId, int from, int tag, bool recv_xmin_bc){

    if (embc->Hx_) {
            if(!recv_xmin_bc && bcId>1 && EM->isXmin){ //When not receiving xmin_bc (MovingWindow), the position of the pml received must be tailored on xmin for non longitudinal PML corner array.
                recvShifted( embc->Hx_, from, tag, embc->ncells_pml_xmin );
                tag++;
                recvShifted( embc->Hy_, from, tag, embc->ncells_pml_xmin );
                tag++;
                recvShifted( embc->Hz_, from, tag, embc->ncells_pml_xmin );
                tag++;
                recvShifted( embc->Bx_, from, tag, embc->ncells_pml_xmin );
                tag++;
                recvShifted( embc->By_, from, tag, embc->ncells_pml_xmin );
                tag++;
                recvShifted( embc->Bz_, from, tag, embc->ncells_pml_xmin );
                tag++;
                recvShifted( embc->Ex_, from, tag, embc->ncells_pml_xmin );
                tag++;
                recvShifted( embc->Ey_, from, tag, embc->ncells_pml_xmin );
                tag++;
                recvShifted( embc->Ez_, from, tag, embc->ncells_pml_xmin );
                tag++;
                recvShifted( embc->Dx_, from, tag, embc->ncells_pml_xmin );
                tag++;
                recvShifted( embc->Dy_, from, tag, embc->ncells_pml_xmin );
                tag++;
                recvShifted( embc->Dz_, from, tag, embc->ncells_pml_xmin );
                tag++;
            } else {
                recv( embc->Hx_, from, tag );
                tag++;
                recv( embc->Hy_, from, tag );
                tag++;
                recv( embc->Hz_, from, tag );
                tag++;
                recv( embc->Bx_, from, tag );
                tag++;
                recv( embc->By_, from, tag );
                tag++;
                recv( embc->Bz_, from, tag );
                tag++;
                recv( embc->Ex_, from, tag );
                tag++;
                recv( embc->Ey_, from, tag );
                tag++;
                recv( embc->Ez_, from, tag );
                tag++;
                recv( embc->Dx_, from, tag );
                tag++;
                recv( embc->Dy_, from, tag );
                tag++;
                recv( embc->Dz_, from, tag );
                tag++;
            }
    }
    // Receive envelope communication
    if (EM->envelope!=NULL){
        if (static_cast<EnvelopeBC2D_PML *>( EM->envelope->EnvBoundCond[bcId] )){
            EnvelopeBC2D_PML *embcenv = static_cast<EnvelopeBC2D_PML *>( EM->envelope->EnvBoundCond[bcId] );
            if (embcenv->A_n_) {
                if(!recv_xmin_bc && bcId>1 && EM->isXmin){ //When not receiving xmin_bc (MovingWindow), the position of the pml received must be tailored on xmin for non longitudinal PML corner array.
                    recvComplexShifted( embcenv->A_n_, from, tag, embc->ncells_pml_xmin );
                    tag++;
                    recvComplexShifted( embcenv->A_nm1_, from, tag, embc->ncells_pml_xmin );
                    tag++;
                    recvComplexShifted( embcenv->u1_nm1_x_, from, tag, embc->ncells_pml_xmin );
                    tag++;
                    recvComplexShifted( embcenv->u2_nm1_x_, from, tag, embc->ncells_pml_xmin );
                    tag++;
                    recvComplexShifted( embcenv->u3_nm1_x_, from, tag, embc->ncells_pml_xmin );
                    tag++;
                    recvComplexShifted( embcenv->u1_nm1_y_, from, tag, embc->ncells_pml_xmin );
                    tag++;
                    recvComplexShifted( embcenv->u2_nm1_y_, from, tag, embc->ncells_pml_xmin );
                    tag++;
                    recvComplexShifted( embcenv->u3_nm1_y_, from, tag, embc->ncells_pml_xmin );
                    tag++;

                } else {
                    recvComplex( embcenv->A_n_, from, tag );
                    tag++;
                    recvComplex( embcenv->A_nm1_, from, tag );
                    tag++;
                    recvComplex( embcenv->u1_nm1_x_, from, tag );
                    tag++;
                    recvComplex( embcenv->u2_nm1_x_, from, tag );
                    tag++;
                    recvComplex( embcenv->u3_nm1_x_, from, tag );
                    tag++;
                    if (bcId > 1) { //Receiving Transverse PML
                        recvComplex( embcenv->u1_nm1_y_, from, tag );
                        tag++;
                        recvComplex( embcenv->u2_nm1_y_, from, tag );
                        tag++;
                        recvComplex( embcenv->u3_nm1_y_, from, tag );
                        tag++;
                    }
                }
            }
        } else if (static_cast<EnvelopeBC3D_PML *>( EM->envelope->EnvBoundCond[bcId] )){
            EnvelopeBC3D_PML *embcenv = static_cast<EnvelopeBC3D_PML *>( EM->envelope->EnvBoundCond[bcId] );
            if (embcenv->A_n_) {
                if(!recv_xmin_bc && bcId>1 && EM->isXmin){ //When not receiving xmin_bc (MovingWindow), the position of the pml received must be tailored on xmin for non longitudinal PML corner array.
                    recvComplexShifted( embcenv->A_n_, from, tag, embc->ncells_pml_xmin );
                    tag++;
                    recvComplexShifted( embcenv->A_nm1_, from, tag, embc->ncells_pml_xmin );
                    tag++;
                    recvComplexShifted( embcenv->u1_nm1_x_, from, tag, embc->ncells_pml_xmin );
                    tag++;
                    recvComplexShifted( embcenv->u2_nm1_x_, from, tag, embc->ncells_pml_xmin );
                    tag++;
                    recvComplexShifted( embcenv->u3_nm1_x_, from, tag, embc->ncells_pml_xmin );
                    tag++;
                    recvComplexShifted( embcenv->u1_nm1_y_, from, tag, embc->ncells_pml_xmin );
                    tag++;
                    recvComplexShifted( embcenv->u2_nm1_y_, from, tag, embc->ncells_pml_xmin );
                    tag++;
                    recvComplexShifted( embcenv->u3_nm1_y_, from, tag, embc->ncells_pml_xmin );
                    tag++;
                    if (bcId >3){
                        recvComplexShifted( embcenv->u1_nm1_z_, from, tag, embc->ncells_pml_xmin );
                        tag++;
                        recvComplexShifted( embcenv->u2_nm1_z_, from, tag, embc->ncells_pml_xmin );
                        tag++;
                        recvComplexShifted( embcenv->u3_nm1_z_, from, tag, embc->ncells_pml_xmin );
                        tag++;
                    }
                } else {
                    recvComplex( embcenv->A_n_, from, tag );
                    tag++;
                    recvComplex( embcenv->A_nm1_, from, tag );
                    tag++;
                    recvComplex( embcenv->u1_nm1_x_, from, tag );
                    tag++;
                    recvComplex( embcenv->u2_nm1_x_, from, tag );
                    tag++;
                    recvComplex( embcenv->u3_nm1_x_, from, tag );
                    tag++;
                    if (bcId > 1) { //Receiving Transverse PML
                        recvComplex( embcenv->u1_nm1_y_, from, tag );
                        tag++;
                        recvComplex( embcenv->u2_nm1_y_, from, tag );
                        tag++;
                        recvComplex( embcenv->u3_nm1_y_, from, tag );
                        tag++;
                        if (bcId > 3) { //Receiving Transverse PML
                            recvComplex( embcenv->u1_nm1_z_, from, tag );
                            tag++;
                            recvComplex( embcenv->u2_nm1_z_, from, tag );
                            tag++;
                            recvComplex( embcenv->u3_nm1_z_, from, tag );
                            tag++;
                        }
                    }
                }
            }

        }
    }
    return tag;
}

void SmileiMPI::recv( ElectroMagn *EM, int from, int &tag, bool recv_xmin_bc )
{
    recv( EM->Ex_, from, tag );
    tag++;
    recv( EM->Ey_, from, tag );
    tag++;
    recv( EM->Ez_, from, tag );
    tag++;
    recv( EM->Bx_, from, tag );
    tag++;
    recv( EM->By_, from, tag );
    tag++;
    recv( EM->Bz_, from, tag );
    tag++;
    recv( EM->Bx_m, from, tag );
    tag++;
    recv( EM->By_m, from, tag );
    tag++;
    recv( EM->Bz_m, from, tag );
    tag++;

    if( EM->envelope!=NULL ) {
        recvComplex( EM->envelope->A_, from, tag );
        tag++;
        recvComplex( EM->envelope->A0_, from, tag );
        tag++;
        recv( EM->envelope->Phi_, from, tag );
        tag++;
        recv( EM->envelope->Phi_m, from, tag );
        tag++;
        recv( EM->envelope->GradPhix_, from, tag );
        tag++;
        recv( EM->envelope->GradPhix_m, from, tag );
        tag++;
        recv( EM->envelope->GradPhiy_, from, tag );
        tag++;
        recv( EM->envelope->GradPhiy_m, from, tag );
        tag++;
        recv( EM->envelope->GradPhiz_, from, tag );
        tag++;
        recv( EM->envelope->GradPhiz_m, from, tag );
        tag++;
    }

    for( unsigned int idiag=0; idiag<EM->allFields_avg.size(); idiag++ ) {
        for( unsigned int ifield=0; ifield<EM->allFields_avg[idiag].size(); ifield++ ) {
            recv( EM->allFields_avg[idiag][ifield], from, tag );
            tag++;
        }
    }

    for( int antennaId=0 ; antennaId<( int )EM->antennas.size() ; antennaId++ ) {
        recv( EM->antennas[antennaId].field, from, tag );
        tag++;
    }

    for( unsigned int bcId=0 ; bcId<EM->emBoundCond.size() ; bcId++ ) {
        if( ! EM->emBoundCond[bcId] ) {
            continue;
        }

        for( unsigned int laserId=0 ; laserId<EM->emBoundCond[bcId]->vecLaser.size() ; laserId++ ) {
            Laser *laser = EM->emBoundCond[bcId]->vecLaser[laserId];
            if( !( laser->spacetime[0] ) && !( laser->spacetime[1] ) ) {
                LaserProfileSeparable *profile;
                profile = static_cast<LaserProfileSeparable *>( laser->profiles[0] );
                if( ! profile->space_envelope ) {
                    continue;
                }
                recv( profile->space_envelope, from, tag );
                tag++;
                recv( profile->phase, from, tag );
                tag++;
                profile = static_cast<LaserProfileSeparable *>( laser->profiles[1] );
                recv( profile->space_envelope, from, tag );
                tag++;
                recv( profile->phase, from, tag );
                tag++;
            }
        }

        if( EM->extFields.size()>0 ) {

            if( dynamic_cast<ElectroMagnBC1D_SM *>( EM->emBoundCond[bcId] ) ) {
                ElectroMagnBC1D_SM *embc = static_cast<ElectroMagnBC1D_SM *>( EM->emBoundCond[bcId] );
                MPI_Status status;
                MPI_Recv( &( embc->By_val ), 1, MPI_DOUBLE, from, tag, MPI_COMM_WORLD, &status );
                tag++;
                MPI_Recv( &( embc->Bz_val ), 1, MPI_DOUBLE, from, tag, MPI_COMM_WORLD, &status );
                tag++;
            } else if( dynamic_cast<ElectroMagnBC2D_SM *>( EM->emBoundCond[bcId] ) ) {
                // BCs at the x-border
                ElectroMagnBC2D_SM *embc = static_cast<ElectroMagnBC2D_SM *>( EM->emBoundCond[bcId] );

                if( embc->B_val[0].size() ) {
                    recv( &embc->B_val[0], from, tag );
                    tag++;
                }
                if( embc->B_val[1].size() ) {
                    recv( &embc->B_val[1], from, tag );
                    tag++;
                }
                if( embc->B_val[2].size() ) {
                    recv( &embc->B_val[2], from, tag );
                    tag++;
                }

            } else if( dynamic_cast<ElectroMagnBC3D_SM *>( EM->emBoundCond[bcId] ) ) {
                ElectroMagnBC3D_SM *embc = static_cast<ElectroMagnBC3D_SM *>( EM->emBoundCond[bcId] );

                // BCs at the border
                if( embc->B_val[0] ) {
                    recv( embc->B_val[0], from, tag );
                    tag++;
                }
                if( embc->B_val[1] ) {
                    recv( embc->B_val[1], from, tag );
                    tag++;
                }
                if( embc->B_val[2] ) {
                    recv( embc->B_val[2], from, tag );
                    tag++;
                }

            }
        }
        if( (dynamic_cast<ElectroMagnBC2D_PML *>( EM->emBoundCond[bcId] ) || dynamic_cast<ElectroMagnBC3D_PML *>( EM->emBoundCond[bcId] )) && (bcId != 0 || recv_xmin_bc) ){
            if( dynamic_cast<ElectroMagnBC2D_PML *>( EM->emBoundCond[bcId] )){
                ElectroMagnBC2D_PML *embc = static_cast<ElectroMagnBC2D_PML *>( EM->emBoundCond[bcId] );
                tag = recv_PML(EM, embc, bcId, from, tag, recv_xmin_bc);
            } else {
                ElectroMagnBC3D_PML *embc = static_cast<ElectroMagnBC3D_PML *>( EM->emBoundCond[bcId] );
                tag = recv_PML(EM, embc, bcId, from, tag, recv_xmin_bc);
            }
        }
    }

} // End recv ( ElectroMagn )



void SmileiMPI::recv( ElectroMagn *EM, int from, int &tag, unsigned int nmodes, bool recv_xmin_bc )
{
    ElectroMagnAM *EMAM = static_cast<ElectroMagnAM *>( EM );
    for( unsigned int imode =0; imode < nmodes; imode++ ) {
        recvComplex( EMAM->El_[imode], from, tag );
        tag++;
        recvComplex( EMAM->Er_[imode], from, tag );
        tag++;
        recvComplex( EMAM->Et_[imode], from, tag );
        tag++;
        recvComplex( EMAM->Bl_[imode], from, tag );
        tag++;
        recvComplex( EMAM->Br_[imode], from, tag );
        tag++;
        recvComplex( EMAM->Bt_[imode], from, tag );
        tag++;
        recvComplex( EMAM->Bl_m[imode], from, tag );
        tag++;
        recvComplex( EMAM->Br_m[imode], from, tag );
        tag++;
        recvComplex( EMAM->Bt_m[imode], from, tag );
        tag++;
    }

    if( EM->envelope!=NULL ) {
        recvComplex( EM->envelope->A_, from, tag );
        tag++;
        recvComplex( EM->envelope->A0_, from, tag );
        tag++;
        recv( EM->envelope->Phi_, from, tag );
        tag++;
        recv( EM->envelope->Phi_m, from, tag );
        tag++;
        recv( EM->envelope->GradPhil_, from, tag );
        tag++;
        recv( EM->envelope->GradPhil_m, from, tag );
        tag++;
        recv( EM->envelope->GradPhir_, from, tag );
        tag++;
        recv( EM->envelope->GradPhir_m, from, tag );
        tag++;
    }

    for( unsigned int idiag=0; idiag<EM->allFields_avg.size(); idiag++ ) {
        for( unsigned int ifield=0; ifield<EM->allFields_avg[idiag].size(); ifield++ ) {
            recv( EM->allFields_avg[idiag][ifield], from, tag );
            tag++;
        }
    }

    for( int antennaId=0 ; antennaId<( int )EM->antennas.size() ; antennaId++ ) {
        recv( EM->antennas[antennaId].field, from, tag );
        tag++;
    }

    for( unsigned int bcId=0 ; bcId<EM->emBoundCond.size() ; bcId++ ) {
        if( ! EM->emBoundCond[bcId] ) {
            continue;
        }

        for( unsigned int laserId=0 ; laserId<EM->emBoundCond[bcId]->vecLaser.size() ; laserId++ ) {
            Laser *laser = EM->emBoundCond[bcId]->vecLaser[laserId];
            if( !( laser->spacetime[0] ) && !( laser->spacetime[1] ) ) {
                LaserProfileSeparable *profile;
                profile = static_cast<LaserProfileSeparable *>( laser->profiles[0] );
                if( ! profile->space_envelope ) {
                    continue;
                }
                recv( profile->space_envelope, from, tag );
                tag++;
                recv( profile->phase, from, tag );
                tag++;
                profile = static_cast<LaserProfileSeparable *>( laser->profiles[1] );
                recv( profile->space_envelope, from, tag );
                tag++;
                recv( profile->phase, from, tag );
                tag++;
            }
        }

        if( dynamic_cast<ElectroMagnBCAM_PML *>( EM->emBoundCond[bcId] ) && (bcId != 0 || recv_xmin_bc)){
            ElectroMagnBCAM_PML *embc = static_cast<ElectroMagnBCAM_PML *>( EM->emBoundCond[bcId] );
            if (embc->Hl_[0]) {
                if(!recv_xmin_bc && bcId>1 && EMAM->isXmin){ //When not receiving xmin_bc (MovingWindow), the position of the pml received must be tailored on xmin for non longitudinal PML corner array.
                    for( unsigned int imode =0; imode < nmodes; imode++ ) {
                        recvComplexShifted( embc->Hl_[imode], from, tag, embc->ncells_pml_lmin);
                        tag++;
                        recvComplexShifted( embc->Hr_[imode], from, tag, embc->ncells_pml_lmin);
                        tag++;
                        recvComplexShifted( embc->Ht_[imode], from, tag, embc->ncells_pml_lmin);
                        tag++;
                        recvComplexShifted( embc->Bl_[imode], from, tag, embc->ncells_pml_lmin);
                        tag++;
                        recvComplexShifted( embc->Br_[imode], from, tag, embc->ncells_pml_lmin);
                        tag++;
                        recvComplexShifted( embc->Bt_[imode], from, tag, embc->ncells_pml_lmin);
                        tag++;
                        recvComplexShifted( embc->El_[imode], from, tag, embc->ncells_pml_lmin);
                        tag++;
                        recvComplexShifted( embc->Er_[imode], from, tag, embc->ncells_pml_lmin);
                        tag++;
                        recvComplexShifted( embc->Et_[imode], from, tag, embc->ncells_pml_lmin);
                        tag++;
                        recvComplexShifted( embc->Dl_[imode], from, tag, embc->ncells_pml_lmin);
                        tag++;
                        recvComplexShifted( embc->Dr_[imode], from, tag, embc->ncells_pml_lmin);
                        tag++;
                        recvComplexShifted( embc->Dt_[imode], from, tag, embc->ncells_pml_lmin);
                        tag++;
                    }
                } else {
                    for( unsigned int imode =0; imode < nmodes; imode++ ) {
                        recvComplex( embc->Hl_[imode], from, tag );
                        tag++;
                        recvComplex( embc->Hr_[imode], from, tag );
                        tag++;
                        recvComplex( embc->Ht_[imode], from, tag );
                        tag++;
                        recvComplex( embc->Bl_[imode], from, tag );
                        tag++;
                        recvComplex( embc->Br_[imode], from, tag );
                        tag++;
                        recvComplex( embc->Bt_[imode], from, tag );
                        tag++;
                        recvComplex( embc->El_[imode], from, tag );
                        tag++;
                        recvComplex( embc->Er_[imode], from, tag );
                        tag++;
                        recvComplex( embc->Et_[imode], from, tag );
                        tag++;
                        recvComplex( embc->Dl_[imode], from, tag );
                        tag++;
                        recvComplex( embc->Dr_[imode], from, tag );
                        tag++;
                        recvComplex( embc->Dt_[imode], from, tag );
                        tag++;
                    }
                }
                if (EM->envelope!=NULL) {
                    EnvelopeBCAM_PML *embcenv = static_cast<EnvelopeBCAM_PML *>( EM->envelope->EnvBoundCond[bcId] );
                    if (embcenv->A_n_) {
                        if(!recv_xmin_bc && bcId==3 && EMAM->isXmin){ //When not receiving xmin_bc (MovingWindow), the position of the pml received must be tailored on xmin for non longitudinal PML corner array.
                            recvComplexShifted( embcenv->A_n_, from, tag, embc->ncells_pml_lmin );
                            tag++;
                            recvComplexShifted( embcenv->A_nm1_, from, tag, embc->ncells_pml_lmin );
                            tag++;
                            recvComplexShifted( embcenv->G_n_, from, tag, embc->ncells_pml_lmin );
                            tag++;
                            recvComplexShifted( embcenv->G_nm1_, from, tag, embc->ncells_pml_lmin );
                            tag++;
                            recvComplexShifted( embcenv->u1_nm1_l_, from, tag, embc->ncells_pml_lmin );
                            tag++;
                            recvComplexShifted( embcenv->u2_nm1_l_, from, tag, embc->ncells_pml_lmin );
                            tag++;
                            recvComplexShifted( embcenv->u3_nm1_l_, from, tag, embc->ncells_pml_lmin );
                            tag++;
                            recvComplexShifted( embcenv->u1_nm1_r_, from, tag, embc->ncells_pml_lmin );
                            tag++;
                            recvComplexShifted( embcenv->u2_nm1_r_, from, tag, embc->ncells_pml_lmin );
                            tag++;
                            recvComplexShifted( embcenv->u3_nm1_r_, from, tag, embc->ncells_pml_lmin );
                            tag++;

                        } else {
                            recvComplex( embcenv->A_n_, from, tag );
                            tag++;
                            recvComplex( embcenv->A_nm1_, from, tag );
                            tag++;
                            recvComplex( embcenv->G_n_, from, tag );
                            tag++;
                            recvComplex( embcenv->G_nm1_, from, tag );
                            tag++;
                            recvComplex( embcenv->u1_nm1_l_, from, tag );
                            tag++;
                            recvComplex( embcenv->u2_nm1_l_, from, tag );
                            tag++;
                            recvComplex( embcenv->u3_nm1_l_, from, tag );
                            tag++;
                            if (bcId == 3) { //Receiving Radial PML
                                recvComplex( embcenv->u1_nm1_r_, from, tag );
                                tag++;
                                recvComplex( embcenv->u2_nm1_r_, from, tag );
                                tag++;
                                recvComplex( embcenv->u3_nm1_r_, from, tag );
                                tag++;
                            }
                        }
                    }
                }

            }
        }

    }

} // End recv ( ElectroMagn LRT )

void SmileiMPI::isend( Field *field, int to, int tag, MPI_Request &request )
{
    //This version of isend(Field) sends the whole array
    MPI_Isend( &( ( *field )( 0 ) ), field->globalDims_, MPI_DOUBLE, to, tag, MPI_COMM_WORLD, &request );

} // End isend ( Field )

void SmileiMPI::isend( Field *field, int to, int tag, MPI_Request &request, int x_first )
{
    //This version of isendComplex(Field) sends only the first "x_first" columns of the array
    int data_size = x_first;
    for (unsigned int idim=1; idim < field->dims_.size(); idim++){
        data_size *= field->dims_[idim];
    }
    MPI_Isend( &( ( *field )( 0 ) ), data_size, MPI_DOUBLE, to, tag, MPI_COMM_WORLD, &request );

} // End isend ( Field )

void SmileiMPI::isendComplex( Field *field, int to, int tag, MPI_Request &request )
{
    cField *cf = static_cast<cField *>( field );
    //This version of isendComplex(Field) sends the whole array
    MPI_Isend( &( ( *cf )( 0 ) ), 2*field->globalDims_, MPI_DOUBLE, to, tag, MPI_COMM_WORLD, &request );

}
void SmileiMPI::isendComplex( Field *field, int to, int tag, MPI_Request &request, int x_first )
{
    cField *cf = static_cast<cField *>( field );
    int data_size = x_first;
    for (unsigned int idim=1; idim < field->dims_.size(); idim++){
        data_size *= field->dims_[idim];
    }
    //This version of isendComplex(Field) sends only the first "x_first" columns of the array
    MPI_Isend( &( ( *cf )( 0 ) ), 2*data_size, MPI_DOUBLE, to, tag, MPI_COMM_WORLD, &request );

} // End isendComplex ( Field )

void SmileiMPI::sendComplex( Field *field, int to, int tag )
{
    cField *cf = static_cast<cField *>( field );
    MPI_Send( &( ( *cf )( 0 ) ), 2*field->globalDims_, MPI_DOUBLE, to, tag, MPI_COMM_WORLD );

} // End isendComplex ( Field )


void SmileiMPI::send(Field* field, int to, int tag)
{
    MPI_Send( &((*field)(0)),field->globalDims_, MPI_DOUBLE, to, tag, MPI_COMM_WORLD );

} // End isend ( Field )


void SmileiMPI::recv( Field *field, int from, int tag)
{
    MPI_Status status;
    //origin shifts the reception position in the array and reduces the received buffer size.
    MPI_Recv( &( ( *field )( 0 ) ), field->globalDims_, MPI_DOUBLE, from, tag, MPI_COMM_WORLD, &status );

} // End recv ( Field )

void SmileiMPI::recvShifted( Field *field, int from, int tag, int xshift )
{
    MPI_Status status;
    int data_shift = xshift;
    for (unsigned int idim=1; idim < field->dims_.size(); idim++){
        data_shift *= field->dims_[idim];
    }
    //Shifts the reception position in the array along the x dimension and reduces the received buffer size.
    MPI_Recv( &( ( *field )( data_shift ) ), field->globalDims_ - data_shift, MPI_DOUBLE, from, tag, MPI_COMM_WORLD, &status );

} // End recv ( Field )

void SmileiMPI::recvComplex( Field *field, int from, int tag )
{
    MPI_Status status;
    cField *cf = static_cast<cField *>( field );
    MPI_Recv( &( ( *cf )( 0 ) ), 2*(field->globalDims_ ), MPI_DOUBLE, from, tag, MPI_COMM_WORLD, &status );

} // End recv ( Field )

void SmileiMPI::recvComplexShifted( Field *field, int from, int tag, int xshift )
{
    MPI_Status status;
    cField *cf = static_cast<cField *>( field );
    int data_shift = xshift;
    for (unsigned int idim=1; idim < field->dims_.size(); idim++){
        data_shift *= field->dims_[idim];
    }
    //Shifts the reception position in the array along the x dimension and reduces the received buffer size.
    MPI_Recv( &( ( *cf )( data_shift ) ), 2*(field->globalDims_ - data_shift), MPI_DOUBLE, from, tag, MPI_COMM_WORLD, &status );

} // End recv ( Field )
void SmileiMPI::irecvComplex( Field *field, int from, int tag, MPI_Request &request )
{
    cField *cf = static_cast<cField *>( field );
    MPI_Irecv( &( ( *cf )( 0 ) ), 2*field->globalDims_, MPI_DOUBLE, from, tag, MPI_COMM_WORLD, &request );

} // End recv ( Field )

void SmileiMPI::irecv(Field* field, int from, int tag, MPI_Request& request)
{
    MPI_Irecv( &((*field)(0)),2*field->globalDims_, MPI_DOUBLE, from, tag, MPI_COMM_WORLD, &request );

} // End recv ( Field )


void SmileiMPI::isend( ProbeParticles *probe, int to, int tag, unsigned int nDim_particles )
{
    MPI_Request request;
    // send offset
    MPI_Isend( &( probe->offset_in_file ), 1, MPI_INT, to, tag, MPI_COMM_WORLD, &request );
    // send number of particles
    int nPart = probe->particles.size();
    MPI_Isend( &nPart, 1, MPI_INT, to, tag+1, MPI_COMM_WORLD, &request );
    // send particles
    if( nPart>0 )
        for( unsigned int i=0; i<nDim_particles; i++ ) {
            MPI_Isend( &( probe->particles.Position[i][0] ), nPart, MPI_DOUBLE, to, tag+1+i, MPI_COMM_WORLD, &request );
        }

} // End isend ( probes )


void SmileiMPI::recv( ProbeParticles *probe, int from, int tag, unsigned int nDim_particles )
{
    MPI_Status status;
    // receive offset
    MPI_Recv( &( probe->offset_in_file ), 1, MPI_INT, from, tag, MPI_COMM_WORLD, &status );
    // receive number of particles
    int nPart;
    MPI_Recv( &nPart, 1, MPI_INT, from, tag+1, MPI_COMM_WORLD, &status );
    // Resize particles
    probe->particles.initialize( nPart, nDim_particles, false );
    // receive particles
    if( nPart>0 )
        for( unsigned int i=0; i<nDim_particles; i++ ) {
            MPI_Recv( &( probe->particles.Position[i][0] ), nPart, MPI_DOUBLE, from, tag+1+i, MPI_COMM_WORLD, &status );
        }

} // End recv ( probes )

//! Wrapper for integer MPI communication
void SmileiMPI::isend( int *integer, int to, int tag, unsigned int nDim_particles, MPI_Request &request )
{
    MPI_Isend( &integer, 1, MPI_INT, to, tag, MPI_COMM_WORLD, &request );
} // End isend ( integer )

//! Wrapper for integer MPI communication
void SmileiMPI::recv( int *integer, int from, int tag, unsigned int nDim_particles )
{
    MPI_Status status;
    MPI_Recv( &integer, 1, MPI_INT, from, tag, MPI_COMM_WORLD, &status );
} // End recv ( integer )

// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ------------------------------------------      DIAGS MPI SYNC     --------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------


// ---------------------------------------------------------------------------------------------------------------------
// Wrapper of MPI synchronization of all computing diags
//   - concerns    : scalars, phasespace, particles
//   - not concern : probes, fields, track particles (each patch write its own data)
//   - called in VectorPatch::runAllDiags(...)
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI::computeGlobalDiags( Diagnostic *diag, int timestep )
{
    if ( DiagnosticScalar* scalar = dynamic_cast<DiagnosticScalar*>( diag ) ) {
        computeGlobalDiags(scalar, timestep);
    } else if (DiagnosticScreen* screen = dynamic_cast<DiagnosticScreen*>( diag )) {
        computeGlobalDiags(screen, timestep);
    } else if (DiagnosticRadiationSpectrum* rad = dynamic_cast<DiagnosticRadiationSpectrum*>( diag )) {
        computeGlobalDiags(rad, timestep);
    } else if (DiagnosticParticleBinning* particles = dynamic_cast<DiagnosticParticleBinning*>( diag )) {
        computeGlobalDiags(particles, timestep);
    }
}


// ---------------------------------------------------------------------------------------------------------------------
// MPI synchronization of scalars diags
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI::computeGlobalDiags( DiagnosticScalar *scalars, int itime )
{

    if( !scalars->timeSelection->theTimeIsNow( itime ) ) {
        return;
    }

    // Reduce all scalars that should be summed
    int n_sum = scalars->values_SUM.size();
    double *d_sum = &scalars->values_SUM[0];
    MPI_Reduce( isMaster()?MPI_IN_PLACE:d_sum, d_sum, n_sum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

    if( scalars->necessary_fieldMinMax_any ) {
        // Reduce all scalars that are a "min" and its location
        int n_min = scalars->values_MINLOC.size();
        val_index *d_min = &scalars->values_MINLOC[0];
        MPI_Reduce( isMaster()?MPI_IN_PLACE:d_min, d_min, n_min, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD );

        // Reduce all scalars that are a "max" and its location
        int n_max = scalars->values_MAXLOC.size();
        val_index *d_max = &scalars->values_MAXLOC[0];
        MPI_Reduce( isMaster()?MPI_IN_PLACE:d_max, d_max, n_max, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD );
    }

    // Complete the computation of the scalars after all reductions
    if( isMaster() ) {

        // Calculate average Z
        for( unsigned int ispec=0; ispec<scalars->sDens.size(); ispec++ )
            if( scalars->sDens[ispec] && scalars->necessary_species[ispec] ) {
                *scalars->sZavg[ispec] = ( double )*scalars->sZavg[ispec] / ( double )*scalars->sDens[ispec];
            }

        // total energy in the simulation
        if( scalars->necessary_Utot ) {
            double Ukin = *scalars->Ukin;
            double Uelm = *scalars->Uelm;
            double Urad = *scalars->Urad;
            *scalars->Utot = Ukin + Uelm + Urad;
        }

        // expected total energy
        if( scalars->necessary_Uexp ) {
            // total energy at time 0
            if( itime==0 ) {
                scalars->Energy_time_zero = *scalars->Utot;
            }
            // Global kinetic energy, and BC losses/gains
            double Ukin_bnd     = *scalars->Ukin_bnd    ;
            double Ukin_out_mvw = *scalars->Ukin_out_mvw;
            double Ukin_inj_mvw = *scalars->Ukin_inj_mvw;
            double Ukin_new     = *scalars->Ukin_new    ;
            // Global elm energy, and BC losses/gains
            double Uelm_bnd     = *scalars->Uelm_bnd    ;
            double Uelm_out_mvw = *scalars->Uelm_out_mvw;
            double Uelm_inj_mvw = *scalars->Uelm_inj_mvw;
            // Global radiated energy
            //double Urad = *scalars->Urad;
            // expected total energy
            double Uexp = scalars->Energy_time_zero
                + Uelm_bnd + Ukin_inj_mvw + Uelm_inj_mvw + Ukin_new
                - ( Ukin_bnd + Ukin_out_mvw + Uelm_out_mvw );
            *scalars->Uexp = Uexp;
        }

        if( scalars->necessary_Ubal ) {
            // energy balance
            double Ubal = ( double )*scalars->Utot - ( double )*scalars->Uexp;
            *scalars->Ubal = Ubal;

            if( scalars->necessary_Ubal_norm ) {
                // the normalized energy balanced is normalized with respect to the current energy
                scalars->EnergyUsedForNorm = *scalars->Utot;
                // normalized energy balance
                double Ubal_norm( 0. );
                if( scalars->EnergyUsedForNorm>0. ) {
                    Ubal_norm = Ubal / scalars->EnergyUsedForNorm;
                }

                *scalars->Ubal_norm = Ubal_norm;
            }
        }

    }
} // END computeGlobalDiags(DiagnosticScalar& scalars ...)


// ---------------------------------------------------------------------------------------------------------------------
// MPI synchronization of diags particle binning
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI::computeGlobalDiags( DiagnosticParticleBinning *diagParticles, int itime )
{
    if( itime - diagParticles->timeSelection->previousTime() == diagParticles->time_average-1 ) {
        MPI_Reduce( diagParticles->filename.size()?MPI_IN_PLACE:&diagParticles->data_sum[0], &diagParticles->data_sum[0], diagParticles->output_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

        if( !isMaster() ) {
            diagParticles->clear();
        }
    }
} // END computeGlobalDiags(DiagnosticParticleBinning* diagParticles ...)

// ---------------------------------------------------------------------------------------------------------------------
// MPI synchronization of diags screen
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI::computeGlobalDiags( DiagnosticScreen *diagScreen, int itime )
{
    if( diagScreen->timeSelection->theTimeIsNow( itime ) ) {
        MPI_Reduce( diagScreen->filename.size()?MPI_IN_PLACE:&diagScreen->data_sum[0], &diagScreen->data_sum[0], diagScreen->output_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

        if( !isMaster() ) {
            diagScreen->clear();
        }
    }
} // END computeGlobalDiags(DiagnosticScreen* diagScreen ...)

// ---------------------------------------------------------------------------------------------------------------------
// MPI synchronization of diags radiation
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI::computeGlobalDiags(DiagnosticRadiationSpectrum* diagRad, int itime)
{
    if (itime - diagRad->timeSelection->previousTime() == diagRad->time_average-1) {
        MPI_Reduce( diagRad->filename.size()?MPI_IN_PLACE:&diagRad->data_sum[0], &diagRad->data_sum[0], diagRad->output_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

        if( !isMaster() ) {
            diagRad->clear();
        }
    }
} // END computeGlobalDiags(DiagnosticRadiationSpectrum*  ...)
