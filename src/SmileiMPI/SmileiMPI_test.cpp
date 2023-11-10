#include "SmileiMPI_test.h"

#include <cmath>
#include <cstring>

#include <iostream>
#include <sstream>
#include <fstream>

#include "Params.h"

using namespace std;


// Obtain an integer argument (next in the argument list)
int get_integer_argument( int *argc, char ***argv )
{
    if( *argc < 2 ) {
        return -1;
    }
    std::string s = ( *argv )[1];
    for( unsigned int i=0; i<s.size(); i++ )
        if( !std::isdigit( s[i] ) ) {
            return -1;
        }
    ( *argv )++;
    ( *argc )--;
    return std::stoi( s.c_str() );
}

// ---------------------------------------------------------------------------------------------------------------------
// SmileiMPI_test constructor
// ---------------------------------------------------------------------------------------------------------------------
SmileiMPI_test::SmileiMPI_test( int *argc, char ***argv )
{
    test_mode = true;

    // If first argument is a number, interpret as the number of MPIs
    int nMPI = get_integer_argument( argc, argv );
    if( nMPI < 0 ) {
        nMPI = 1;
    }
    if( nMPI == 0 ) {
        ERROR( "The first argument in test mode (number of MPIs) cannot be zero" );
    }
    // If first argument is a number, interpret as the number of OMPs
    int nOMP = get_integer_argument( argc, argv );
    if( nOMP < 0 ) {
        nOMP = 1;
    }
    if( nOMP == 0 ) {
        ERROR( "The second argument in test mode (number of OMPs) cannot be zero" );
    }
    // Return a test-mode SmileiMPI

    int mpi_provided;
#ifdef _OPENMP
    MPI_Init_thread( argc, argv, MPI_THREAD_MULTIPLE, &mpi_provided );
#ifndef _NO_MPI_TM
    if( mpi_provided != MPI_THREAD_MULTIPLE ) {
        ERROR( "MPI_THREAD_MULTIPLE not supported. Compile your MPI library with THREAD_MULTIPLE support." );
    }
#endif
    omp_set_num_threads( 1 );
#else
    MPI_Init( argc, argv );
#endif

    world_ = MPI_COMM_WORLD;
    MPI_Comm_size( world_, &smilei_sz );
    MPI_Comm_rank( world_, &smilei_rk );

    if( smilei_sz > 1 ) {
        ERROR( "Test mode cannot be run with several MPI processes. Instead, indicate the MPIxOMP intended partition after the -T argument." );
    }

    smilei_sz = nMPI;
    smilei_rk = 0;
    smilei_omp_max_threads = nOMP;
    number_of_cores = nOMP;
    global_number_of_cores = nMPI * nOMP;

    MESSAGE( "    ----- TEST MODE WITH INTENDED PARTITION "<<nMPI<<"x"<<nOMP<<"-----" );

} // END SmileiMPI_test::SmileiMPI_test


// ---------------------------------------------------------------------------------------------------------------------
// SmileiMPI_test destructor :
// ---------------------------------------------------------------------------------------------------------------------
SmileiMPI_test::~SmileiMPI_test()
{
    remove( "smilei.py" );
} // END SmileiMPI_test::~SmileiMPI_test


// ---------------------------------------------------------------------------------------------------------------------
// Initialize MPI (per process) environment
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI_test::init( Params &params, DomainDecomposition *domain_decomposition )
{
    // Initialize patch environment
    patch_count.resize( smilei_sz, 0 );
    Tcapabilities = smilei_sz;
    use_BTIS3 = params.use_BTIS3;

    remove( "patch_load.txt" );

    // Initialize patch distribution
    if( !params.restart ) {
        init_patch_count( params, domain_decomposition );
    }

    // Initialize buffers for particles push vectorization
    dynamics_Epart.resize( 1 );
    dynamics_Bpart.resize( 1 );
    dynamics_invgf.resize( 1 );
    dynamics_iold.resize( 1 );
    dynamics_deltaold.resize( 1 );

    // Set periodicity of the simulated problem
    periods_  = new int[params.nDim_field];
    for( unsigned int i=0 ; i<params.nDim_field ; i++ ) {
        periods_[i] = 0;
        if( params.EM_BCs[i][0]=="periodic" ) {
            periods_[i] = 1;
            MESSAGE( 1, "applied topology for periodic BCs in "<<"xyz"[i]<<"-direction" );
        }
    }
} // END init


// ---------------------------------------------------------------------------------------------------------------------
//  Initialize patch distribution
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI_test::init_patch_count( Params &params, DomainDecomposition */*domain_decomposition*/ )
{

    unsigned int Npatches;
    //unsigned int ncells_perpatch;

    //Compute target load: Tload = Total load * local capability / Total capability.

    // Some initialization of the box parameters
    Npatches = params.tot_number_of_patches;
    //ncells_perpatch = 1;
    // vector<double> cell_xmin( 3, 0. ), cell_xmax( 3, 1. ), cell_dx( 3, 2. ), x_cell( 3, 0 );
    // for( unsigned int i = 0; i < params.nDim_field; i++ ) {
    //     ncells_perpatch *= params.patch_size_[i]+2*params.oversize[i];
    //     if( params.cell_length[i]!=0. ) {
    //         cell_dx[i] = params.cell_length[i];
    //     }
    // }

    // First, distribute all patches evenly
    unsigned int Npatches_local = Npatches / smilei_sz;
    int remainder = Npatches % smilei_sz;
    for( int rk=0; rk<smilei_sz; rk++ ) {
        if( rk < remainder ) {
            patch_count[rk] = Npatches_local + 1;
        } else {
            patch_count[rk] = Npatches_local;
        }
    }

} // END init_patch_count
