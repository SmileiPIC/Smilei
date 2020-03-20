#include <algorithm>
#include <cmath>
#include <ctime>
#include <iomanip>

#define SMILEI_IMPORT_ARRAY

#include "PyTools.h"
#include "Params.h"
#include "Species.h"
#include "Tools.h"
#include "SmileiMPI.h"
#include "H5.h"
#include "LaserPropagator.h"

#include "pyinit.pyh"
#include "pyprofiles.pyh"
#include "pycontrol.pyh"

using namespace std;

namespace Rand
{
std::random_device device;
std::mt19937 gen( device() );

std::uniform_real_distribution<double> uniform_distribution( 0., 1. );
double uniform()
{
    return uniform_distribution( gen );
}

std::uniform_real_distribution<double> uniform_distribution1( 0., 1.-1e-11 );
double uniform1()
{
    return uniform_distribution1( gen );
}

std::uniform_real_distribution<double> uniform_distribution2( -1., 1. );
double uniform2()
{
    return uniform_distribution2( gen );
}
double normal( double stddev )
{
    std::normal_distribution<double> normal_distribution( 0., stddev );
    return normal_distribution( gen );
}
}

#define DO_EXPAND(VAL)  VAL ## 1
#define EXPAND(VAL)     DO_EXPAND(VAL)
#ifdef SMILEI_USE_NUMPY
#if !defined(NUMPY_IMPORT_ARRAY_RETVAL) || (EXPAND(NUMPY_IMPORT_ARRAY_RETVAL) == 1)
void smilei_import_array()   // python 2
{
    import_array();
}
#else
void *smilei_import_array()   // hack for python3
{
    import_array();
    return NULL;
}
#endif
#endif


// ---------------------------------------------------------------------------------------------------------------------
// Params : open & parse the input data file, test that parameters are coherent
// ---------------------------------------------------------------------------------------------------------------------
Params::Params( SmileiMPI *smpi, std::vector<std::string> namelistsFiles ) :
    namelist( "" )
{

    MESSAGE( "HDF5 version "<<H5_VERS_MAJOR << "." << H5_VERS_MINOR << "." << H5_VERS_RELEASE );

    if( ( H5_VERS_MAJOR< 1 ) ||
            ( ( H5_VERS_MAJOR==1 ) && ( H5_VERS_MINOR< 8 ) ) ||
            ( ( H5_VERS_MAJOR==1 ) && ( H5_VERS_MINOR==8 ) && ( H5_VERS_RELEASE<16 ) ) ) {
        WARNING( "Smilei suggests using HDF5 version 1.8.16 or newer" );
    }

    if( namelistsFiles.size()==0 ) {
        ERROR( "No namelists given!" );
    }

    //string commandLineStr("");
    //for (unsigned int i=0;i<namelistsFiles.size();i++) commandLineStr+="\""+namelistsFiles[i]+"\" ";
    //MESSAGE(1,commandLineStr);

    //init Python
    PyTools::openPython();
#ifdef SMILEI_USE_NUMPY
    smilei_import_array();
    // Workaround some numpy multithreading bug
    // https://github.com/numpy/numpy/issues/5856
    // We basically call the command numpy.seterr(all="ignore")
    PyObject *numpy = PyImport_ImportModule( "numpy" );
    string seterr( "seterr" );
    string sChar( "s" );
    Py_DECREF( PyObject_CallMethod( numpy, &seterr[0], &sChar[0], "ignore" ) );
    Py_DECREF( numpy );
#endif

    // Print python version
    MESSAGE( 1, "Python version "<<PyTools::python_version() );

    // First, we tell python to filter the ctrl-C kill command (or it would prevent to kill the code execution).
    // This is done separately from other scripts because we don't want it in the concatenated python namelist.
    PyTools::checkPyError();
    string command = "import signal\nsignal.signal(signal.SIGINT, signal.SIG_DFL)";
    if( !PyRun_SimpleString( command.c_str() ) ) {
        PyTools::checkPyError();
    }

    PyObject *Py_main = PyImport_AddModule( "__main__" );
    PyObject *globals = PyModule_GetDict( Py_main );

    // Running pyinit.py
    runScript( string( reinterpret_cast<const char *>( pyinit_py ), pyinit_py_len ), "pyinit.py", globals );

    runScript( Tools::merge( "smilei_version='", string( __VERSION ), "'\n" ), string( __VERSION ), globals );

    // Set the _test_mode to False
    PyObject_SetAttrString( Py_main, "_test_mode", Py_False );
    PyTools::checkPyError();

    // here we add the rank, in case some script need it
    PyModule_AddIntConstant( Py_main, "smilei_mpi_rank", smpi->getRank() );

    // here we add the MPI size, in case some script need it
    PyModule_AddIntConstant( Py_main, "smilei_mpi_size", smpi->getSize() );
    namelist += string( "smilei_mpi_size = " ) + to_string( smpi->getSize() ) + "\n";

    // here we add the larget int, important to get a valid seed for randomization
    PyModule_AddIntConstant( Py_main, "smilei_rand_max", RAND_MAX );
    namelist += string( "smilei_rand_max = " ) + to_string( RAND_MAX ) + "\n\n";

    // Running pyprofiles.py
    runScript( string( reinterpret_cast<const char *>( pyprofiles_py ), pyprofiles_py_len ), "pyprofiles.py", globals );

    namelist += "\n\n";
    namelist += "\"\"\"\n";
    namelist += "-----------------------------------------------------------------------\n";
    namelist += "BEGINNING OF THE USER NAMELIST\n";
    namelist += "\"\"\"\n\n";

    // Running the namelists
    for( vector<string>::iterator it=namelistsFiles.begin(); it!=namelistsFiles.end(); it++ ) {
        string strNamelist="";
        if( smpi->isMaster() ) {
            ifstream istr( it->c_str() );
            // If file
            if( istr.is_open() ) {
                std::stringstream buffer;
                buffer << istr.rdbuf();
                strNamelist+=buffer.str();
                // If command
            } else {
                command = *it;
                // Remove quotes
                unsigned int s = command.size();
                if( s>1 && command.substr( 0, 1 )=="\"" && command.substr( s-1, 1 )=="\"" ) {
                    command = command.substr( 1, s - 2 );
                }
                // Add to namelist
                strNamelist = Tools::merge( "# Smilei:) From command line:\n", command );
            }
            strNamelist +="\n";
        }
        smpi->bcast( strNamelist );
        runScript( strNamelist, ( *it ), globals );
    }

    // Running pycontrol.py
    runScript( string( reinterpret_cast<const char *>( pycontrol_py ), pycontrol_py_len ), "pycontrol.py", globals );

    // Run custom pre-processing function
    MESSAGE( 1, "Check for function preprocess()" );
    PyTools::runPyFunction( "preprocess" );
    PyErr_Clear();

    smpi->barrier();

    // Error if no block Main() exists
    if( PyTools::nComponents( "Main" ) == 0 ) {
        ERROR( "Block Main() not defined" );
    }

    // CHECK namelist on python side
    PyTools::runPyFunction( "_smilei_check" );
    smpi->barrier();

    // Python makes the checkpoint dir tree
    if( ! smpi->test_mode ) {
        PyTools::runPyFunction( "_prepare_checkpoint_dir" );
        smpi->barrier();
    }

    // Now the string "namelist" contains all the python files concatenated
    // It is written as a file: smilei.py
    if( smpi->isMaster() ) {
        ofstream out_namelist( "smilei.py" );
        if( out_namelist.is_open() ) {
            out_namelist << "# coding: utf-8" << endl << endl ;
            out_namelist << namelist;
            out_namelist.close();
        }
    }

    // random seed
    if( PyTools::extract( "random_seed", random_seed, "Main" ) ) {
        // Init of the seed for the vectorized C++ random generator recommended by Intel
        // See https://software.intel.com/en-us/articles/random-number-function-vectorization
        srand48( random_seed );
        // Init of the seed for the C++ random generator
        Rand::gen = std::mt19937( random_seed );
    }

    // communication pattern initialized as partial B exchange
    full_B_exchange = false;
    // communication pattern initialized as partial A, Phi exchange for envelope simulations
    full_Envelope_exchange = false;

    // --------------
    // Stop & Restart
    // --------------

    restart = false;
    std::vector<std::string> _unused_restart_files;
    if( PyTools::nComponents( "Checkpoints" )>0 && PyTools::extract( "restart_files", _unused_restart_files, "Checkpoints" ) ) {
        MESSAGE( 1, "Code will restart" );
        restart=true;
    }

    // ---------------------
    // Normalisation & units
    // ---------------------

    reference_angular_frequency_SI = 0.;
    PyTools::extract( "reference_angular_frequency_SI", reference_angular_frequency_SI, "Main" );


    // -------------------
    // Simulation box info
    // -------------------

    // geometry of the simulation
    geometry = "";
    if( !PyTools::extract( "geometry", geometry, "Main" ) ) {
        ERROR( "Parameter Main.geometry is required" );
    }
    if( geometry!="1Dcartesian" && geometry!="2Dcartesian" && geometry!="3Dcartesian" && geometry!="AMcylindrical" ) {
        ERROR( "Main.geometry `" << geometry << "` invalid" );
    }
    setDimensions();

    // interpolation order
    PyTools::extract( "interpolation_order", interpolation_order, "Main" );
    if( interpolation_order!=2 && interpolation_order!=4 ) {
        ERROR( "Main.interpolation_order " << interpolation_order << " not defined" );
    }

    //!\todo (MG to JD) Please check if this parameter should still appear here
    // Disabled, not compatible for now with particles sort
    // if ( !PyTools::extract("exchange_particles_each", exchange_particles_each) )
    exchange_particles_each = 1;

    PyTools::extract( "every_clean_particles_overhead", every_clean_particles_overhead, "Main" );

    // TIME & SPACE RESOLUTION/TIME-STEPS

    // reads timestep & cell_length
    PyTools::extract( "timestep", timestep, "Main" );
    res_time = 1.0/timestep;

    PyTools::extract( "cell_length", cell_length, "Main" );
    if( cell_length.size()!=nDim_field ) {
        ERROR( "Dimension of cell_length ("<< cell_length.size() << ") != " << nDim_field << " for geometry " << geometry );
    }
    res_space.resize( nDim_field );
    for( unsigned int i=0; i<nDim_field; i++ ) {
        res_space[i] = 1.0/cell_length[i];
    }
    // Number of modes in AMcylindrical geometry
    PyTools::extract( "number_of_AM", nmodes, "Main" );

    nmodes_rel_field_init = 1;

    // Number of modes in AMcylindrical geometry for relativistic field initialization
    // if not specified, it will be equal to the number of modes of the simulation
    PyTools::extract( "number_of_AM_relativistic_field_initialization", nmodes_rel_field_init, "Main" );
    if (nmodes_rel_field_init>nmodes){
        ERROR( "The number of AM modes computed in relativistic field initialization must be lower or equal than the number of modes of the simulation" );
    }

    
    // simulation duration & length
    PyTools::extract( "simulation_time", simulation_time, "Main" );

    PyTools::extract( "grid_length", grid_length, "Main" );
    if( grid_length.size()!=nDim_field ) {
        ERROR( "Dimension of grid_length ("<< grid_length.size() << ") != " << nDim_field << " for geometry " << geometry );
    }


    //! Boundary conditions for ElectroMagnetic Fields
    if( !PyTools::extract( "EM_boundary_conditions", EM_BCs, "Main" ) ) {
        ERROR( "Electromagnetic boundary conditions (EM_boundary_conditions) not defined" );
    }

    if( EM_BCs.size() == 0 ) {
        ERROR( "EM_boundary_conditions cannot be empty" );
    } else if( EM_BCs.size() == 1 ) {
        while( EM_BCs.size() < nDim_field ) {
            EM_BCs.push_back( EM_BCs[0] );
        }
    } else if( EM_BCs.size() != nDim_field ) {
        ERROR( "EM_boundary_conditions must be the same size as the number of dimensions" );
    }

    for( unsigned int iDim=0; iDim<nDim_field; iDim++ ) {
        if( EM_BCs[iDim].size() == 1 ) { // if just one type is specified, then take the same bc type in a given dimension
            EM_BCs[iDim].push_back( EM_BCs[iDim][0] );
        } else if( ( EM_BCs[iDim][0] != EM_BCs[iDim][1] ) && ( EM_BCs[iDim][0] == "periodic" || EM_BCs[iDim][1] == "periodic" ) ) {
            ERROR( "EM_boundary_conditions along dimension "<<"012"[iDim]<<" cannot be periodic only on one side" );
        }
    }

    int n_envlaser = PyTools::nComponents( "LaserEnvelope" );
    if( n_envlaser >=1 ) {
        Laser_Envelope_model = true;
        //! Boundary conditions for Envelope Field
        if( !PyTools::extract( "Envelope_boundary_conditions", Env_BCs, "LaserEnvelope" ) ) {
            ERROR( "Envelope_boundary_conditions not defined" );
        }

        if( Env_BCs.size() == 0 ) {
            ERROR( "Envelope_boundary_conditions cannot be empty" );
        } else if( Env_BCs.size() == 1 ) {
            while( Env_BCs.size() < nDim_field ) {
                Env_BCs.push_back( Env_BCs[0] );
            }
        } else if( Env_BCs.size() != nDim_field ) {
            ERROR( "Envelope_boundary_conditions must be the same size as the number of dimensions" );
        }

        for( unsigned int iDim=0; iDim<nDim_field; iDim++ ) {
            if( Env_BCs[iDim].size() == 1 ) { // if just one type is specified, then take the same bc type in a given dimension
                Env_BCs[iDim].push_back( Env_BCs[iDim][0] );
            }
            //    else if ( (Env_BCs[iDim][0] != Env_BCs[iDim][1]) &&  (Env_BCs[iDim][0] == "periodic" || Env_BCs[iDim][1] == "periodic") )
            //        ERROR("Envelope_boundary_conditions along "<<"xyz"[iDim]<<" cannot be periodic only on one side");
        }
        
        if (!PyTools::extract( "envelope_solver", envelope_solver, "LaserEnvelope" )){
            ERROR("envelope_solver not defined in the LaserEnvelope block");
        }
        if ( (envelope_solver != "explicit") && (envelope_solver != "explicit_reduced_dispersion") ){
            ERROR("Unknown envelope_solver - only 'explicit' and 'explicit_reduced_dispersion' are available. ");
        }
        if (envelope_solver == "explicit_reduced_dispersion"){
            full_Envelope_exchange = true;
        }

    }

    for( unsigned int iDim = 0 ; iDim < nDim_field; iDim++ ) {
        if( EM_BCs[iDim][0] == "buneman" || EM_BCs[iDim][1] == "buneman" ) {
            full_B_exchange = true;
            open_boundaries = true;
        }
        if( EM_BCs[iDim][0] == "silver-muller" || EM_BCs[iDim][1] == "silver-muller" ) {
            open_boundaries = true;
        }
    }

    PyTools::extract( "EM_boundary_conditions_k", EM_BCs_k, "Main" );
    if( EM_BCs_k.size() == 0 ) {
        //Gives default value
        for( unsigned int iDim=0; iDim<nDim_field; iDim++ ) {
            std::vector<double> temp_k;

            for( unsigned int iiDim=0; iiDim<iDim; iiDim++ ) {
                temp_k.push_back( 0. );
            }
            temp_k.push_back( 1. );
            for( unsigned int iiDim=iDim+1; iiDim<nDim_field; iiDim++ ) {
                temp_k.push_back( 0. );
            }
            EM_BCs_k.push_back( temp_k );
            for( unsigned int iiDim=0; iiDim<nDim_field; iiDim++ ) {
                temp_k[iiDim] *= -1. ;
            }
            EM_BCs_k.push_back( temp_k );
        }
    }

    //Complete with zeros if not defined
    if( EM_BCs_k.size() == 1 ) {
        while( EM_BCs_k.size() < nDim_field*2 ) {
            EM_BCs_k.push_back( EM_BCs_k[0] );
        }
    } else if( EM_BCs_k.size() != nDim_field*2 ) {
        ERROR( "EM_boundary_conditions_k must be the same size as the number of faces." );
    }
    for( unsigned int iDim=0; iDim<nDim_field*2; iDim++ ) {
        if( EM_BCs_k[iDim].size() != nDim_field ) {
            ERROR( "EM_boundary_conditions_k must have exactly " << nDim_field << " elements along dimension "<<"-+"[iDim%2]<<"012"[iDim/2] );
        }
        if( EM_BCs_k[iDim][iDim/2] == 0. ) {
            ERROR( "EM_boundary_conditions_k must have a non zero normal component along dimension "<<"-+"[iDim%2]<<"012"[iDim/2] );
        }

    }
    save_magnectic_fields_for_SM = true;
    PyTools::extract( "save_magnectic_fields_for_SM", save_magnectic_fields_for_SM, "Main" );

    // -----------------------------------
    // MAXWELL SOLVERS & FILTERING OPTIONS
    // -----------------------------------

    time_fields_frozen=0.0;
    PyTools::extract( "time_fields_frozen", time_fields_frozen, "Main" );

    // Poisson Solver
    PyTools::extract( "solve_poisson", solve_poisson, "Main" );
    PyTools::extract( "poisson_max_iteration", poisson_max_iteration, "Main" );
    PyTools::extract( "poisson_max_error", poisson_max_error, "Main" );
    // Relativistic Poisson Solver
    PyTools::extract( "solve_relativistic_poisson", solve_relativistic_poisson, "Main" );
    PyTools::extract( "relativistic_poisson_max_iteration", relativistic_poisson_max_iteration, "Main" );
    PyTools::extract( "relativistic_poisson_max_error", relativistic_poisson_max_error, "Main" );

    // PXR parameters
    PyTools::extract( "is_spectral", is_spectral, "Main" );
    if( is_spectral ) {
        full_B_exchange=true;
    }
    PyTools::extract( "is_pxr", is_pxr, "Main" );

    // Maxwell Solver
    PyTools::extract( "maxwell_solver", maxwell_sol, "Main" );
    if( maxwell_sol == "Lehe" ) {
        full_B_exchange=true;
    }

    // Current filter properties
    currentFilter_passes = 0;
    int nCurrentFilter = PyTools::nComponents( "CurrentFilter" );
    for( int ifilt = 0; ifilt < nCurrentFilter; ifilt++ ) {
        string model;
        PyTools::extract( "model", model, "CurrentFilter", ifilt );
        if( model != "binomial" ) {
            ERROR( "Currently, only the `binomial` model is available in CurrentFilter()" );
        }
        PyTools::extract( "passes", currentFilter_passes, "CurrentFilter", ifilt );
    }

    // Field filter properties
    Friedman_filter = false;
    Friedman_theta = 0;
    int nFieldFilter = PyTools::nComponents( "FieldFilter" );
    for( int ifilt = 0; ifilt < nFieldFilter; ifilt++ ) {
        string model;
        PyTools::extract( "model", model, "FieldFilter", ifilt );
        if( model != "Friedman" ) {
            ERROR( "Currently, only the `Friedman` model is available in FieldFilter()" );
        }
        if( geometry != "2Dcartesian" ) {
            ERROR( "Currently, the `Friedman` field filter is only availble in `2Dcartesian` geometry" );
        }
        Friedman_filter = true;
        PyTools::extract( "theta", Friedman_theta, "FieldFilter", ifilt );
        if( Friedman_filter && ( Friedman_theta==0. ) ) {
            WARNING( "Friedman filter is applied but parameter theta is set to zero" );
        }
        if( ( Friedman_theta<0. ) || ( Friedman_theta>1. ) ) {
            ERROR( "Friedman filter theta = " << Friedman_theta << " must be between 0 and 1" );
        }
    }


    // testing the CFL condition
    //!\todo (MG) CFL cond. depends on the Maxwell solv. ==> HERE JUST DONE FOR YEE!!!
    double res_space2=0;
    for( unsigned int i=0; i<nDim_field; i++ ) {
        res_space2 += res_space[i]*res_space[i];
    }
    if( geometry == "AMcylindrical" ) {
        res_space2 += ( ( nmodes-1 )*( nmodes-1 )-1 )*res_space[1]*res_space[1];
    }
    dtCFL=1.0/sqrt( res_space2 );
    if( timestep>dtCFL ) {
        WARNING( "CFL problem: timestep=" << timestep << " should be smaller than " << dtCFL );
    }



    // clrw
    PyTools::extract( "clrw", clrw, "Main" );



    // --------------------
    // Number of patches
    // --------------------
    if( !PyTools::extract( "number_of_patches", number_of_patches, "Main" ) ) {
        ERROR( "The parameter `number_of_patches` must be defined as a list of integers" );
    }

    tot_number_of_patches = 1;
    for( unsigned int iDim=0 ; iDim<nDim_field ; iDim++ ) {
        tot_number_of_patches *= number_of_patches[iDim];
    }

    if( tot_number_of_patches == ( unsigned int )( smpi->getSize() ) ) {
        one_patch_per_MPI = true;
    } else {
        one_patch_per_MPI = false;
        if( tot_number_of_patches < ( unsigned int )( smpi->getSize() ) ) {
            ERROR( "The total number of patches "<<tot_number_of_patches<<" must be greater or equal to the number of MPI processes "<<smpi->getSize() );
        }
    }
#ifdef _OPENMP
    if( tot_number_of_patches < ( unsigned int )( smpi->getSize()*smpi->getOMPMaxThreads() ) ) {
        WARNING( "Resources allocated "<<( smpi->getSize()*smpi->getOMPMaxThreads() )<<" underloaded regarding the total number of patches "<<tot_number_of_patches );
    }
#endif

    global_factor.resize( nDim_field, 1 );
    PyTools::extract( "global_factor", global_factor, "Main" );
    norder.resize( nDim_field, 1 );
    norder.resize( nDim_field, 1 );
    PyTools::extract( "norder", norder, "Main" );
    //norderx=norder[0];
    //nordery=norder[1];
    //norderz=norder[2];


    if( PyTools::extract( "patch_arrangement", patch_arrangement, "Main" ) ) {
        WARNING( "Change patches distribution to " << patch_arrangement );
    } else {
        patch_arrangement = "hilbertian";
        WARNING( "Use default distribution : " << patch_arrangement );
    }


    int total_number_of_hilbert_patches = 1;
    if( patch_arrangement == "hilbertian" ) {
        for( unsigned int iDim=0 ; iDim<nDim_field ; iDim++ ) {
            total_number_of_hilbert_patches *= number_of_patches[iDim];
            if( ( number_of_patches[iDim] & ( number_of_patches[iDim]-1 ) ) != 0 ) {
                ERROR( "Number of patches in each direction must be a power of 2" );
            }
        }
    }


    if( PyTools::nComponents( "LoadBalancing" )>0 ) {
        // get parameter "every" which describes a timestep selection
        load_balancing_time_selection = new TimeSelection(
            PyTools::extract_py( "every", "LoadBalancing" ), "Load balancing"
        );
        PyTools::extract( "cell_load", cell_load, "LoadBalancing" );
        PyTools::extract( "frozen_particle_load", frozen_particle_load, "LoadBalancing" );
        PyTools::extract( "initial_balance", initial_balance, "LoadBalancing" );
    } else {
        load_balancing_time_selection = new TimeSelection();
    }

    has_load_balancing = ( smpi->getSize()>1 )  && ( ! load_balancing_time_selection->isEmpty() );

    if( has_load_balancing && patch_arrangement != "hilbertian" ) {
        ERROR( "Dynamic load balancing is only available for Hilbert decomposition" );
    }
    if( has_load_balancing && total_number_of_hilbert_patches < 2*smpi->getSize() ) {
        ERROR( "Dynamic load balancing requires to use at least 2 patches per MPI process." );
    }

    mi.resize( 3, 0 );
    while( ( number_of_patches[0] >> mi[0] ) >1 ) {
        mi[0]++ ;
    }
    if( number_of_patches.size()>1 ) {
        while( ( number_of_patches[1] >> mi[1] ) >1 ) {
            mi[1]++ ;
        }
        if( number_of_patches.size()>2 )
            while( ( number_of_patches[2] >> mi[2] ) >1 ) {
                mi[2]++ ;
            }
    }

    // Activation of the vectorized subroutines
    vectorization_mode = "off";
    has_adaptive_vectorization = false;
    adaptive_vecto_time_selection = nullptr;

    if( PyTools::nComponents( "Vectorization" )>0 ) {
        // Extraction of the vectorization mode
        PyTools::extract( "mode", vectorization_mode, "Vectorization" );
        if( !( vectorization_mode == "off" ||
                vectorization_mode == "on" ||
                vectorization_mode == "adaptive" ||
                vectorization_mode == "adaptive_mixed_sort" ) ) {
            ERROR( "In block `Vectorization`, parameter `mode` must be `off`, `on`, `adaptive`" );
        } else if( vectorization_mode == "adaptive_mixed_sort" || vectorization_mode == "adaptive" ) {
            has_adaptive_vectorization = true;
        }

        // Check that we are in 3D, adaptive mode not possible in 2d
        if( vectorization_mode == "adaptive_mixed_sort" || vectorization_mode == "adaptive" ) {
            if (nDim_particle != 3) {
                ERROR("In block `Vectorization`, `adaptive` mode only available in 3D")
            }
        }

        // Default mode for the adaptive mode
        PyTools::extract( "initial_mode", adaptive_default_mode, "Vectorization" );
        if( !( adaptive_default_mode == "off" ||
                adaptive_default_mode == "on" ) ) {
            ERROR( "In block `Vectorization`, parameter `default` must be `off` or `on`" );
        }

        // get parameter "every" which describes a timestep selection
        if( ! adaptive_vecto_time_selection )
            adaptive_vecto_time_selection = new TimeSelection(
                PyTools::extract_py( "reconfigure_every", "Vectorization" ), "Adaptive vectorization"
            );
    }

    ;
    PyTools::extract( "cell_sorting", cell_sorting, "Main" );
    //MESSAGE("Sorting per cell : " << cell_sorting );
    //if (cell_sorting)
    //    vectorization_mode = "on";
    
    // In case of collisions, ensure particle sort per cell
    if( PyTools::nComponents( "Collisions" ) > 0 ) {

        if( geometry!="1Dcartesian"
                && geometry!="2Dcartesian"
                && geometry!="3Dcartesian" ) {
            ERROR( "Collisions only valid for cartesian geometries for the moment" )
        }
        
        // collisions need sorting per cell
        if( vectorization_mode == "adaptive_mixed_sort" ) {
            ERROR( "Collisions are incompatible with the vectorization mode 'adaptive_mixed_sort'." )
        }
        
        if( vectorization_mode == "off" ) {
            cell_sorting = true;
        }
            
    }

    // Read the "print_every" parameter
    print_every = ( int )( simulation_time/timestep )/10;
    PyTools::extract( "print_every", print_every, "Main" );
    if( !print_every ) {
        print_every = 1;
    }

    // Read the "print_expected_disk_usage" parameter
    if( ! PyTools::extract( "print_expected_disk_usage", print_expected_disk_usage, "Main" ) ) {
        ERROR( "The parameter `Main.print_expected_disk_usage` must be True or False" );
    }

    // -------------------------------------------------------
    // Checking species order
    // -------------------------------------------------------
    // read from python namelist the number of species
    unsigned int tot_species_number = PyTools::nComponents( "Species" );

    double mass, mass2=0;

    for( unsigned int ispec = 0; ispec < tot_species_number; ispec++ ) {
        PyTools::extract( "mass", mass, "Species", ispec );
        if( mass == 0 ) {
            for( unsigned int ispec2 = ispec+1; ispec2 < tot_species_number; ispec2++ ) {
                PyTools::extract( "mass", mass2, "Species", ispec2 );
                if( mass2 > 0 ) {
                    ERROR( "the photon species (mass==0) should be defined after the particle species (mass>0)" );
                }
            }
        }
    }

    // -------------------------------------------------------
    // Parameters for the synchrotron-like radiation losses
    // -------------------------------------------------------
    hasMCRadiation = false ;// Default value
    hasLLRadiation = false ;// Default value
    hasNielRadiation = false ;// Default value


    // Loop over all species to check if the radiation losses are activated
    std::string radiation_model = "none";
    for( unsigned int ispec = 0; ispec < tot_species_number; ispec++ ) {

        PyTools::extract( "radiation_model", radiation_model, "Species", ispec );

        // Cancelation of the letter case for `radiation_model`
        std::transform( radiation_model.begin(), radiation_model.end(), radiation_model.begin(), ::tolower );

        if( radiation_model=="monte-carlo" || radiation_model=="mc" ) {
            this->hasMCRadiation = true;
        } else if( radiation_model=="landau-lifshitz"
                   || radiation_model=="ll"
                   || radiation_model=="corrected-landau-lifshitz"
                   || radiation_model=="cll" ) {
            this->hasLLRadiation = true;
        } else if( radiation_model=="niel" ) {
            this->hasNielRadiation = true;
        }
    }

    // -------------------------------------------------------
    // Parameters for the mutliphoton Breit-Wheeler pair decay
    // -------------------------------------------------------
    this->hasMultiphotonBreitWheeler = false ;// Default value

    std::vector<std::string> multiphoton_Breit_Wheeler( 2 );
    for( unsigned int ispec = 0; ispec < tot_species_number; ispec++ ) {

        if( PyTools::extract( "multiphoton_Breit_Wheeler", multiphoton_Breit_Wheeler, "Species", ispec ) ) {
            this->hasMultiphotonBreitWheeler = true;
        }
    }

    // -------------------------------------------------------
    // Compute useful quantities and introduce normalizations
    // also defines defaults values for the species lengths
    // -------------------------------------------------------
    compute();

    // -------------------------------------------------------
    // Print
    // -------------------------------------------------------
    smpi->barrier();
    if( smpi->isMaster() ) {
        print_init();
    }
    smpi->barrier();

    // -------------------------------------------------------
    // Handle the pre-processing of LaserOffset
    // -------------------------------------------------------
    unsigned int n_laser = PyTools::nComponents( "Laser" );
    unsigned int n_laser_offset = 0;
    LaserPropagator propagateX;
    
    for( unsigned int i_laser=0; i_laser<n_laser; i_laser++ ) {
        double offset = 0.;
        
        // If this laser has the hidden _offset attribute
        if( PyTools::extract( "_offset", offset, "Laser", i_laser ) ) {
            
            bool propagate = false;
            PyTools::extract( "_propagate", propagate, "Laser", i_laser );
            if( ! propagate ) continue;
            
            // Extract the angle
            double angle_z = 0.;
            PyTools::extract( "_angle", angle_z, "Laser", i_laser );
            
            // Prepare propagator
            if( n_laser_offset == 0 ) {
                TITLE( "Pre-processing LaserOffset" );
                propagateX.init( this, smpi, 0 );
            }
            
            MESSAGE( 1, "LaserOffset #"<< n_laser_offset );
            
            // Extract the file name
            string file( "" );
            PyTools::extract( "file", file, "Laser", i_laser );
            
            // Extract the list of profiles and verify their content
            PyObject *p = PyTools::extract_py( "_profiles", "Laser", i_laser );
            vector<PyObject *> profiles;
            vector<int> profiles_n = {1, 2};
            if( ! PyTools::convert( p, profiles ) ) {
                ERROR( "For LaserOffset #" << n_laser_offset << ": space_time_profile must be a list of 2 profiles" );
            }
            Py_DECREF( p );
            if( profiles.size()!=2 ) {
                ERROR( "For LaserOffset #" << n_laser_offset << ": space_time_profile needs 2 profiles." );
            }
            if( profiles[1] == Py_None ) {
                profiles  .pop_back();
                profiles_n.pop_back();
            }
            if( profiles[0] == Py_None ) {
                profiles  .erase( profiles  .begin() );
                profiles_n.erase( profiles_n.begin() );
            }
            if( profiles.size() == 0 ) {
                ERROR( "For LaserOffset #" << n_laser_offset << ": space_time_profile cannot be [None, None]" );
            }
            for( unsigned int i=0; i<profiles.size(); i++ ) {
                int nargs = PyTools::function_nargs( profiles[i] );
                if( nargs<0 ) {
                    ERROR( "For LaserOffset #" << n_laser_offset << ": space_time_profile["<<i<<"] not callable" );
                }
                if( nargs != ( int ) nDim_field ) {
                    ERROR( "For LaserOffset #" << n_laser_offset << ": space_time_profile["<<i<<"] requires " << nDim_field << " arguments but has " << nargs );
                }
            }
            
            // Extract _keep_n_strongest_modes
            int keep_n_strongest_modes=0;
            if( !PyTools::extract( "_keep_n_strongest_modes", keep_n_strongest_modes, "Laser", i_laser ) || keep_n_strongest_modes<1 ) {
                ERROR( "For LaserOffset #" << n_laser_offset << ": keep_n_strongest_modes must be a positive integer" );
            }
            
            // Make the propagation happen and write out the file
            if( ! smpi->test_mode && ! restart ) {
                propagateX( profiles, profiles_n, offset, file, keep_n_strongest_modes, angle_z );
            }
            
            n_laser_offset ++;
        }
    }
    
    check_consistency();
}

Params::~Params()
{
    if( load_balancing_time_selection ) {
        delete load_balancing_time_selection;
    }
    if( adaptive_vecto_time_selection ) {
        delete adaptive_vecto_time_selection;
    }
    PyTools::closePython();
}

// ---------------------------------------------------------------------------------------------------------------------
// Compute useful values (normalisation, time/space step, etc...)
// ---------------------------------------------------------------------------------------------------------------------
void Params::compute()
{
    // time-related parameters
    // -----------------------

    // number of time-steps
    n_time   = ( int )( simulation_time/timestep );

    // simulation time & time-step value
    double entered_simulation_time = simulation_time;
    simulation_time = ( double )( n_time ) * timestep;
    if( simulation_time!=entered_simulation_time )
        WARNING( "simulation_time has been redefined from " << entered_simulation_time
                 << " to " << simulation_time << " to match timestep." );


    // grid/cell-related parameters
    // ----------------------------
    n_space.resize( 3, 1 );
    cell_length.resize( 3 );
    n_space_global.resize( 3, 1 ); //! \todo{3 but not real size !!! Pbs in Species::Species}
    oversize.resize( 3, 0 );
    patch_dimensions.resize( 3, 0. );
    cell_volume=1.0;
    n_cell_per_patch = 1;

    // compute number of cells & normalized lengths
    for( unsigned int i=0; i<nDim_field; i++ ) {
        n_space[i] = round( grid_length[i]/cell_length[i] );
        double entered_grid_length = grid_length[i];
        grid_length[i] = ( double )( n_space[i] )*cell_length[i]; // ensure that nspace = grid_length/cell_length
        if( grid_length[i]!=entered_grid_length ) {
            WARNING( "grid_length[" << i << "] has been redefined from " << entered_grid_length << " to " << grid_length[i] << " to match n x cell_length (" << scientific << setprecision( 4 ) << grid_length[i]-entered_grid_length <<")" );
        }
        cell_volume *= cell_length[i];
    }
    if( geometry == "AMcylindrical" ) {
        cell_volume *= 2 * M_PI;
    }
    // create a 3d equivalent of n_space & cell_length
    for( unsigned int i=nDim_field; i<3; i++ ) {
        cell_length[i]=0.0;
    }

    for( unsigned int i=0; i<nDim_field; i++ ) {
        oversize[i]  = max( interpolation_order, ( unsigned int )( norder[i]/2+1 ) ) + ( exchange_particles_each-1 );;
        n_space_global[i] = n_space[i];
        n_space[i] /= number_of_patches[i];
        if( n_space_global[i]%number_of_patches[i] !=0 ) {
            ERROR( "ERROR in dimension " << i <<". Number of patches = " << number_of_patches[i] << " must divide n_space_global = " << n_space_global[i] );
        }
        if( n_space[i] <= 2*oversize[i]+1 ) {
            ERROR( "ERROR in dimension " << i <<". Patches length = "<<n_space[i] << " cells must be at least " << 2*oversize[i] +2 << " cells long. Increase number of cells or reduce number of patches in this direction. " );
        }
        patch_dimensions[i] = n_space[i] * cell_length[i];
        n_cell_per_patch *= n_space[i];
    }

    // Set clrw if not set by the user
    if( clrw == -1 ) {

        // default value
        clrw = n_space[0];

        // check cache issue for interpolation/projection
        int cache_threshold( 3200 ); // sizeof( L2, Sandy Bridge-HASWELL ) / ( 10 * sizeof(double) )
        // Compute the "transversal bin size"
        int bin_size( 1 );
        for( unsigned int idim = 1 ; idim < nDim_field ; idim++ ) {
            bin_size *= ( n_space[idim]+1+2*oversize[idim] );
        }

        // IF Ionize or pair generation : clrw = n_space_x_pp ?
        if( ( clrw+1+2*oversize[0] ) * bin_size > ( unsigned int ) cache_threshold ) {
            int clrw_max = cache_threshold / bin_size - 1 - 2*oversize[0];
            if( clrw_max > 0 ) {
                for( clrw=clrw_max ; clrw > 0 ; clrw-- )
                    if( ( ( clrw+1+2*oversize[0] ) * bin_size <= ( unsigned int ) cache_threshold ) && ( n_space[0]%clrw==0 ) ) {
                        break;
                    }
            } else {
                clrw = 1;
            }
            WARNING( "Particles cluster width set to : " << clrw );
        }

    }

    // clrw != n_space[0] is not compatible
    // with the adaptive vectorization for the moment
    if( vectorization_mode == "adaptive_mixed_sort" || vectorization_mode == "adaptive" ) {
        if( clrw != ( int )( n_space[0] ) ) {
            clrw = ( int )( n_space[0] );
            WARNING( "Particles cluster width set to: " << clrw << " for the adaptive vectorization mode" );
        }
    }

    // Verify that clrw divides n_space[0]
    if( n_space[0]%clrw != 0 ) {
        ERROR( "The parameter clrw must divide the number of cells in one patch (in dimension x)" );
    }

}


void Params::check_consistency()
{
    if( vectorization_mode != "off" ) {

        if( ( geometry=="1Dcartesian" ) ) {
            ERROR( "Vectorized algorithms not implemented for this geometry" );
        }

        if( ( geometry=="2Dcartesian" ) && ( interpolation_order==4 ) ) {
            ERROR( "4th order vectorized algorithms not implemented in 2D" );
        }


        if( hasMultiphotonBreitWheeler ) {
            WARNING( "Performances of advanced physical processes which generates new particles could be degraded for the moment !" );
            WARNING( "\t The improvment of their integration in vectorized algorithm is in progress." );
        }

    }

}


// ---------------------------------------------------------------------------------------------------------------------
// Set dimensions according to geometry
// ---------------------------------------------------------------------------------------------------------------------
void Params::setDimensions()
{
    if( geometry=="1Dcartesian" ) {
        nDim_particle=1;
        nDim_field=1;
    } else if( geometry=="2Dcartesian" ) {
        nDim_particle=2;
        nDim_field=2;
    } else if( geometry=="3Dcartesian" ) {
        nDim_particle=3;
        nDim_field=3;
    } else if( geometry=="AMcylindrical" ) {
        nDim_particle=3;
        nDim_field=2;
    } else {
        ERROR( "Geometry: " << geometry << " not defined" );
    }
}



// ---------------------------------------------------------------------------------------------------------------------
// Printing out the data at initialisation
// ---------------------------------------------------------------------------------------------------------------------
void Params::print_init()
{
    TITLE( "Geometry: " << geometry );
    MESSAGE( 1, "Interpolation order : " <<  interpolation_order );
    MESSAGE( 1, "Maxwell solver : " <<  maxwell_sol );
    MESSAGE( 1, "(Time resolution, Total simulation time) : (" << res_time << ", " << simulation_time << ")" );
    MESSAGE( 1, "(Total number of iterations,   timestep) : (" << n_time << ", " << timestep << ")" );
    MESSAGE( 1, "           timestep  = " << timestep/dtCFL << " * CFL" );

    for( unsigned int i=0 ; i<grid_length.size() ; i++ ) {
        MESSAGE( 1, "dimension " << i << " - (Spatial resolution, Grid length) : (" << res_space[i] << ", " << grid_length[i] << ")" );
        MESSAGE( 1, "            - (Number of cells,    Cell length)  : " << "(" << n_space_global[i] << ", " << cell_length[i] << ")" );
        MESSAGE( 1, "            - Electromagnetic boundary conditions: " << "(" << EM_BCs[i][0] << ", " << EM_BCs[i][1] << ")" );
        if( open_boundaries ) {
            cout << setprecision( 2 );
            cout << "                     - Electromagnetic boundary conditions k    : " << "( [" << EM_BCs_k[2*i][0] ;
            for( unsigned int ii=1 ; ii<grid_length.size() ; ii++ ) {
                cout << ", " << EM_BCs_k[2*i][ii] ;
            }
            cout << "] , [" << EM_BCs_k[2*i+1][0] ;
            for( unsigned int ii=1 ; ii<grid_length.size() ; ii++ ) {
                cout << ", " << EM_BCs_k[2*i+1][ii] ;
            }
            cout << "] )" << endl;
        }
    }

    if( currentFilter_passes > 0 ) {
        MESSAGE( 1, "Binomial current filtering : "<< currentFilter_passes << " passes" );
    }
    if( Friedman_filter ) {
        MESSAGE( 1, "Friedman field filtering : theta = " << Friedman_theta );
    }
    if( full_B_exchange ) {
        MESSAGE( 1, "All components of B are exchanged at synchronization" );
    }

    if( has_load_balancing ) {
        TITLE( "Load Balancing: " );
        if( initial_balance ) {
            MESSAGE( 1, "Computational load is initially balanced between MPI ranks. (initial_balance = true) " );
        } else {
            MESSAGE( 1, "Patches are initially homogeneously distributed between MPI ranks. (initial_balance = false) " );
        }
        MESSAGE( 1, "Happens: " << load_balancing_time_selection->info() );
        MESSAGE( 1, "Cell load coefficient = " << cell_load );
        MESSAGE( 1, "Frozen particle load coefficient = " << frozen_particle_load );
    }

    TITLE( "Vectorization: " );
    MESSAGE( 1, "Mode: " << vectorization_mode );
    if( vectorization_mode == "adaptive_mixed_sort" || vectorization_mode == "adaptive" ) {
        MESSAGE( 1, "Default mode: " << adaptive_default_mode );
        MESSAGE( 1, "Time selection: " << adaptive_vecto_time_selection->info() );
    }

}

// ---------------------------------------------------------------------------------------------------------------------
// Printing out some data at a given timestep
// ---------------------------------------------------------------------------------------------------------------------
void Params::print_timestep( unsigned int itime, double time_dual, Timer &timer )
{
    double before = timer.getTime();
    timer.update();
    double now = timer.getTime();
    ostringstream my_msg;
    my_msg << "  " << setw( timestep_width ) << itime << "/" << n_time << " "
           << "  " << scientific << setprecision( 4 ) << setw( 12 ) << time_dual << " "
           << "  " << scientific << setprecision( 4 ) << setw( 12 ) << now << " "
           << "  " << "(" << scientific << setprecision( 4 ) << setw( 12 ) << now - before << " )"
           ;
    #pragma omp master
    MESSAGE( my_msg.str() );
    #pragma omp barrier
}

void Params::print_timestep_headers()
{
    timestep_width = log10( n_time ) + 1;
    if( timestep_width<3 ) {
        timestep_width = 3;
    }
    ostringstream my_msg;
    my_msg << setw( timestep_width*2+4 ) << " timestep "
           << setw( 15 ) << "sim time "
           << setw( 15 ) << "cpu time [s] "
           << "  (" << setw( 12 ) << "diff [s]" << " )"
           ;
    MESSAGE( my_msg.str() );
}


// Print information about the MPI aspects
void Params::print_parallelism_params( SmileiMPI *smpi )
{
    if( smpi->isMaster() ) {
#ifndef _NO_MPI_TM
        MESSAGE( 1, "MPI_THREAD_MULTIPLE enabled" );
#else
        MESSAGE( 1, "MPI_THREAD_MULTIPLE not enabled" );
#endif
        MESSAGE( 1, "Number of MPI process : " << smpi->getSize() );
        MESSAGE( 1, "Number of patches : " );
        for( unsigned int iDim=0 ; iDim<nDim_field ; iDim++ ) {
            MESSAGE( 2, "dimension " << iDim << " - number_of_patches : " << number_of_patches[iDim] );
        }

        MESSAGE( 1, "Patch size :" );
        for( unsigned int iDim=0 ; iDim<nDim_field ; iDim++ ) {
            MESSAGE( 2, "dimension " << iDim << " - n_space : " << n_space[iDim] << " cells." );
        }

        MESSAGE( 1, "Dynamic load balancing: " << load_balancing_time_selection->info() );
    }

    if( smpi->isMaster() ) {
        TITLE( "OpenMP" );
#ifdef _OPENMP
//    int nthds(0);
//#pragma omp parallel shared(nthds)
//    {
//        nthds = omp_get_num_threads();
//    }
        MESSAGE( 1, "Number of thread per MPI process : " << smpi->getOMPMaxThreads() );
#else
        MESSAGE( "Disabled" );
#endif
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// Finds requested species in the list of existing species.
// Returns an array of the numbers of the requested species.
// Note that there might be several species that have the same "name" or "type"
//  so that we have to search for all possibilities.
vector<unsigned int> Params::FindSpecies( vector<Species *> &vecSpecies, vector<string> requested_species )
{
    bool species_found;
    vector<unsigned int> result;
    unsigned int i;
    vector<string> existing_species;

    // Make an array of the existing species names
    existing_species.resize( 0 );
    for( unsigned int ispec=0 ; ispec<vecSpecies.size() ; ispec++ ) {
        existing_species.push_back( vecSpecies[ispec]->name_ );
    }

    // Loop over group of requested species
    for( unsigned int rs=0 ; rs<requested_species.size() ; rs++ ) {
        species_found = false;
        // Loop over existing species
        for( unsigned int es=0 ; es<existing_species.size() ; es++ ) {
            if( requested_species[rs] == existing_species[es] ) { // if found
                species_found = true;
                // Add to the list and sort
                for( i=0 ; i<result.size() ; i++ ) {
                    if( es == result[i] ) {
                        break;    // skip if duplicate
                    }
                    if( es <  result[i] ) {
                        result.insert( result.begin()+i, es ); // insert at the right place
                        break;
                    }
                }
                // Put at the end if not put earlier
                if( i == result.size() ) {
                    result.push_back( es );
                }
            }
        }
        if( !species_found ) {
            ERROR( "Species `" << requested_species[rs] << "` was not found." );
        }
    }

    return result;
}


//! Run string as python script and add to namelist
void Params::runScript( string command, string name, PyObject *scope )
{
    PyTools::checkPyError();
    namelist+=command;
    if( name.size()>0 ) {
        MESSAGE( 1, "Parsing " << name );
    }
    PyObject *result = PyRun_String( command.c_str(), Py_file_input, scope, scope );
    PyTools::checkPyError();
    if( !result ) {
        ERROR( "error parsing "<< name );
    }
    Py_DECREF( result );
}

//! run the python functions cleanup (user defined) and _keep_python_running (in pycontrol.py)
void Params::cleanup( SmileiMPI *smpi )
{
    // call cleanup function from the user namelist (it can be used to free some memory
    // from the python side) while keeping the interpreter running
    MESSAGE( 1, "Checking for cleanup() function:" );
    PyTools::runPyFunction( "cleanup" );
    // this will reset error in python in case cleanup doesn't exists
    PyErr_Clear();

    smpi->barrier();

    // this function is defined in the Python/pyontrol.py file and should return false if we can close
    // the python interpreter
    MESSAGE( 1, "Calling python _keep_python_running() :" );
    if( PyTools::runPyFunction<bool>( "_keep_python_running" ) ) {
        MESSAGE( 2, "Keeping Python interpreter alive" );
    } else {
        MESSAGE( 2, "Closing Python" );
        PyErr_Print();
        Py_Finalize();
    }
    smpi->barrier();
}

bool Params::isSpeciesField( string field_name )
{
    if( geometry!="AMcylindrical" ) {
        if( ( field_name.at( 0 )=='J' && field_name.length()>2 )
                || ( field_name.at( 0 )=='R' && field_name.length()>3 ) ) {
            return true;
        }
    } else {
        if( ( ( field_name.at( 0 )=='J' && field_name.length()>8 )
                && ( field_name.substr( 2, 6 )!="_mode_" || field_name.find( "mode_" ) != field_name.rfind( "mode_" ) ) )
                || ( ( field_name.at( 0 )=='R' && field_name.length()>9 )
                     && ( field_name.substr( 3, 6 )!="_mode_" || field_name.find( "mode_" ) != field_name.rfind( "mode_" ) ) ) ) {
            return true;
        }
    }
    return false;
}
