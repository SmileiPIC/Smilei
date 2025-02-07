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
#if PY_MAJOR_VERSION < 3
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
        ERROR_NAMELIST( "No namelist (input file) given!","" );
    }

    //string commandLineStr("");
    //for (unsigned int i=0;i<namelistsFiles.size();i++) commandLineStr+="\""+namelistsFiles[i]+"\" ";
    //MESSAGE(1,commandLineStr);

    //init Python
    PyTools::openPython();
    // Print python version
    MESSAGE( "Python version "<<PyTools::python_version() );

#ifdef SMILEI_USE_NUMPY
    smilei_import_array();
    // Workaround some numpy multithreading bug
    // https://github.com/numpy/numpy/issues/5856
    // We basically call the command numpy.seterr(all="ignore")
    PyObject *numpy = PyImport_ImportModule( "numpy" );
    string seterr( "seterr" );
    string sChar( "s" );
    Py_DECREF( PyObject_CallMethod( numpy, &seterr[0], &sChar[0], "ignore" ) );
    string numpy_version = "";
    PyTools::getAttr( numpy, "__version__", numpy_version );
    MESSAGE( "Numpy version " << numpy_version );
    Py_DECREF( numpy );
#else
    WARNING("Numpy not found. Some options will not be available");
#endif


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

    // we add the rank, in case some script needs it
    PyModule_AddIntConstant( Py_main, "smilei_mpi_rank", smpi->getRank() );

    // we add the MPI size, in case some script needs it
    PyModule_AddIntConstant( Py_main, "smilei_mpi_size", smpi->getSize() );
    namelist += string( "smilei_mpi_size = " ) + to_string( smpi->getSize() ) + "\n";

    // we add the openMP size, in case some script needs it
    PyModule_AddIntConstant( Py_main, "smilei_omp_threads", smpi->getOMPMaxThreads() );
    namelist += string( "smilei_omp_threads = " ) + to_string( smpi->getOMPMaxThreads() ) + "\n";

    // we add the total number of cores, in case some script needs it
    PyModule_AddIntConstant( Py_main, "smilei_total_cores", smpi->getGlobalNumCores() );
    namelist += string( "smilei_total_cores = " ) + to_string( smpi->getGlobalNumCores() ) + "\n";

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
        ERROR_NAMELIST( "Block Main() not defined",LINK_NAMELIST + std::string("#main-variables") );
    }

    // CHECK namelist on python side
    PyTools::runPyFunction( "_smilei_check" );
    smpi->barrier();

    // Python makes the checkpoint dir tree
    if( ! smpi->test_mode ) {
        PyTools::runPyFunction( "_prepare_checkpoint_dir" );
        smpi->barrier();
    }

    // Call python function _keep_python_running (see pyontrol.py)
    // Return false if we can close the python interpreter
    MESSAGE( 1, "Calling python _keep_python_running() :" );
    keep_python_running_ = PyTools::runPyFunction<bool>( "_keep_python_running" );

    // random seed
    random_seed = 0;
    if( PyTools::extractOrNone( "random_seed", random_seed, "Main" ) ) {
        // Init of the seed for the vectorized C++ random generator recommended by Intel
        // See https://software.intel.com/en-us/articles/random-number-function-vectorization
        srand48( random_seed );
        // Init of the seed for the C++ random generator
        Rand::gen = std::mt19937( random_seed );
    }

    // communication pattern initialized as partial B exchange
    full_B_exchange = false;
    // communication pattern initialized as partial A, Phi exchange for envelope simulations
    full_Envelope_exchange = true;

    // --------------
    // Stop & Restart
    // --------------

    restart = false;
    std::vector<std::string> _unused_restart_files;
    if( PyTools::nComponents( "Checkpoints" )>0 && PyTools::extractV( "restart_files", _unused_restart_files, "Checkpoints" ) ) {
        MESSAGE( 1, "Code will restart" );
        restart=true;
    }

    // ---------------------
    // Normalisation & units
    // ---------------------

    reference_angular_frequency_SI = 0.;
    PyTools::extract( "reference_angular_frequency_SI", reference_angular_frequency_SI, "Main"   );

    // -------------------
    // Simulation box info
    // -------------------

    // geometry of the simulation
    geometry = "";
    PyTools::extract( "geometry", geometry, "Main"  );
    if( geometry!="1Dcartesian" && geometry!="2Dcartesian" && geometry!="3Dcartesian" && geometry!="AMcylindrical" ) {
        ERROR_NAMELIST( "Main.geometry `" << geometry << "` invalid", LINK_NAMELIST + std::string("#main-variables") );
    }
    setDimensions();

    // Maxwell Solver
    PyTools::extract( "maxwell_solver", maxwell_sol, "Main"   );
    is_spectral = false;
    is_pxr = false;
    if( maxwell_sol == "Lehe" || maxwell_sol == "Bouchard" || maxwell_sol == "M4" ) {
        full_B_exchange=true;
    } else if( maxwell_sol == "spectral" ) {
        is_spectral = true;
        is_pxr = true;
        full_B_exchange = true;
    } else if( maxwell_sol == "picsar" ) {
        is_pxr = true;
    }

#ifndef _PICSAR
    if (is_pxr) {
        ERROR_NAMELIST( "Smilei not linked with picsar, use make config=picsar", "https://smileipic.github.io/Smilei/install_PICSAR.html" );
    }
#endif

    // interpolation order
    PyTools::extract( "interpolation_order", interpolation_order, "Main"  );
    if( geometry=="AMcylindrical") {
        if( interpolation_order != 1 && is_spectral ) {
            ERROR_NAMELIST( "Main.interpolation_order " << interpolation_order << " should be 1 for PSATD solver",
            LINK_NAMELIST + std::string("#main-variables") );
        }
    } else if( interpolation_order!=2 && interpolation_order!=4 && !is_spectral ) {
        ERROR_NAMELIST( "Main.interpolation_order " << interpolation_order << " should be 2 or 4",
        LINK_NAMELIST + std::string("#main-variables"));
    }

    // Interpolation scheme
    PyTools::extract( "interpolator", interpolator_, "Main"  );

    // Cancelation of the letter case
    std::transform( interpolator_.begin(), interpolator_.end(), interpolator_.begin(), ::tolower );

    if (interpolator_ != "wt" && interpolator_ != "momentum-conserving") {
        ERROR_NAMELIST( "Parameter `Main.interpolator` should be `momentum-conserving` or `wt`.",
        LINK_NAMELIST + std::string("#main-variables"));
    }

    if( ( interpolator_  == "wt") &&
        (geometry != "1Dcartesian")                &&
        (geometry != "2Dcartesian")                &&
        (geometry != "3Dcartesian")               ) {
        ERROR_NAMELIST( "Interpolator `wt` not implemented for geometry: " << geometry << ".",
        LINK_NAMELIST + std::string("#main-variables") );
    }



    //!\todo (MG to JD) Please check if this parameter should still appear here
    // Disabled, not compatible for now with particles sort
    // if ( !PyTools::extract("exchange_particles_each", exchange_particles_each) )
    exchange_particles_each = 1;

    PyTools::extract( "every_clean_particles_overhead", every_clean_particles_overhead, "Main"   );

    // TIME & SPACE RESOLUTION/TIME-STEPS

    // reads timestep & cell_length
    PyTools::extract( "timestep", timestep, "Main"   );
    res_time = 1.0/timestep;

    PyTools::extractV( "cell_length", cell_length, "Main" );
    if( cell_length.size()!=nDim_field ) {
        ERROR_NAMELIST( "Dimension of cell_length ("<< cell_length.size() << ") != " << nDim_field << " for geometry " << geometry,
        LINK_NAMELIST + std::string("#main-variables"));
    }
    res_space.resize( nDim_field );
    for( unsigned int i=0; i<nDim_field; i++ ) {
        res_space[i] = 1.0/cell_length[i];
    }
    // Number of modes in AMcylindrical geometry
    PyTools::extract( "number_of_AM", nmodes, "Main"   );

    nmodes_rel_field_init = 1; // default value

    // Number of modes in AMcylindrical geometry for relativistic field initialization
    // if not specified, it will be equal to 1
    PyTools::extract( "number_of_AM_relativistic_field_initialization", nmodes_rel_field_init, "Main"   );
    if (nmodes_rel_field_init>nmodes){
        ERROR_NAMELIST( "The number of AM modes computed in relativistic field initialization must be lower or equal than the number of modes of the simulation",
                    LINK_NAMELIST + std::string("#main-variables") );
    }

    nmodes_classical_Poisson_field_init = 1; // default value

    // Number of modes in AMcylindrical geometry for non relativistic field initialization with Poisson solver
    // if not specified, it will be equal to 1
    PyTools::extract( "number_of_AM_classical_Poisson_solver", nmodes_classical_Poisson_field_init, "Main"   );
    if (nmodes_classical_Poisson_field_init>nmodes){
        ERROR_NAMELIST( "The number of AM modes computed in classical Poisson solver must be lower or equal than the number of modes of the simulation",
                        LINK_NAMELIST + std::string("#main-variables") );
    }


    // simulation duration & length
    PyTools::extract( "simulation_time", simulation_time, "Main"   );

    PyTools::extractV( "grid_length", grid_length, "Main" );
    if( grid_length.size()!=nDim_field ) {
        ERROR_NAMELIST( "Dimension of grid_length ("<< grid_length.size() << ") != " << nDim_field << " for geometry " << geometry,
        LINK_NAMELIST + std::string("#main-variables"));
    }


    //! Boundary conditions for ElectroMagnetic Fields
    if( !PyTools::extractVV( "EM_boundary_conditions", EM_BCs, "Main" ) ) {
        ERROR_NAMELIST( "Electromagnetic boundary conditions (EM_boundary_conditions) not defined",
                        LINK_NAMELIST + std::string("#main-variables"));
    }

    if( EM_BCs.size() == 0 ) {
        ERROR_NAMELIST( "EM_boundary_conditions cannot be empty",
               LINK_NAMELIST + std::string("#main-variables") );
    } else if( EM_BCs.size() == 1 ) {
        while( EM_BCs.size() < nDim_field ) {
            EM_BCs.push_back( EM_BCs[0] );
        }
    } else if( EM_BCs.size() != nDim_field ) {
        ERROR_NAMELIST( "EM_boundary_conditions must be the same size as the number of dimensions",
                        LINK_NAMELIST + std::string("#main-variables") );
    }

    bool use_pml = false;
    for( unsigned int iDim=0; iDim<nDim_field; iDim++ ) {
        if( EM_BCs[iDim].size() == 1 ) { // if just one type is specified, then take the same bc type in a given dimension
            EM_BCs[iDim].push_back( EM_BCs[iDim][0] );
        } else if( EM_BCs[iDim][0] != EM_BCs[iDim][1] && ( EM_BCs[iDim][0] == "periodic" || EM_BCs[iDim][1] == "periodic" ) ) {
            ERROR_NAMELIST( "EM_boundary_conditions along "<<"xyz"[iDim]<<" cannot be periodic only on one side.",
                            LINK_NAMELIST + std::string("#main-variables") );
        } else if( EM_BCs[iDim][0] != EM_BCs[iDim][1] && ( EM_BCs[iDim][0] == "PML" || EM_BCs[iDim][1] == "PML" ) ) {
            ERROR_NAMELIST( "EM_boundary_conditions along "<<"xyz"[iDim]<<" cannot be PML only on one side for the moment.",
                            LINK_NAMELIST + std::string("#main-variables") );
        }
        if( is_spectral && geometry != "AMcylindrical" && ( EM_BCs[iDim][0] != "periodic" || EM_BCs[iDim][1] != "periodic" ) ) {
            ERROR_NAMELIST( "EM_boundary_conditions along "<<"xyz"[iDim]<<" must be periodic for spectral solver in cartesian geometry.",
                            LINK_NAMELIST + std::string("#main-variables") );
        }
        if (EM_BCs[iDim][0] == "PML" || EM_BCs[iDim][1] == "PML"){ use_pml = true; }
    }

    int n_envlaser = PyTools::nComponents( "LaserEnvelope" );
    if( n_envlaser >=1 ) {
        Laser_Envelope_model = true;
        //! Boundary conditions for Envelope Field
        if( !PyTools::extractVV( "Envelope_boundary_conditions", Env_BCs, "LaserEnvelope" ) ) {
            ERROR_NAMELIST( "Envelope_boundary_conditions not defined",
                            LINK_NAMELIST + std::string("#laser-envelope-model") );
        }

        if( Env_BCs.size() == 0 ) {
            ERROR_NAMELIST( "Envelope_boundary_conditions cannot be empty",
                            LINK_NAMELIST + std::string("#laser-envelope-model"));
        } else if( Env_BCs.size() == 1 ) {
            while( Env_BCs.size() < nDim_field ) {
                Env_BCs.push_back( Env_BCs[0] );
            }
        } else if( Env_BCs.size() != nDim_field ) {
            ERROR_NAMELIST( "Envelope_boundary_conditions must be the same size as the number of dimensions",
                           LINK_NAMELIST + std::string("#laser-envelope-model") );
        }

        for( unsigned int iDim=0; iDim<nDim_field; iDim++ ) {
            if( Env_BCs[iDim].size() == 1 ) { // if just one type is specified, then take the same bc type in a given dimension
                Env_BCs[iDim].push_back( Env_BCs[iDim][0] );
            }
            //    else if ( (Env_BCs[iDim][0] != Env_BCs[iDim][1]) &&  (Env_BCs[iDim][0] == "periodic" || Env_BCs[iDim][1] == "periodic") )
            //        ERROR("Envelope_boundary_conditions along "<<"xyz"[iDim]<<" cannot be periodic only on one side");
        }


        // Find if at least one species is ionized by envelope
        int n_species = PyTools::nComponents( "Species" );
        std:: string ionization_model;
        for (int i_species=0;i_species<n_species;i_species++){
            PyTools::extract( "ionization_model", ionization_model, "Species", i_species );
            if (ionization_model=="tunnel_envelope" || ionization_model=="tunnel_envelope_averaged" ){
                envelope_ionization_is_active = true;
                break;
            }
        }

        // Read envelope solver for the envelope equation
        PyTools::extract( "envelope_solver", envelope_solver, "LaserEnvelope" );
        if ( (envelope_solver != "explicit") && (envelope_solver != "explicit_reduced_dispersion") ){
            ERROR_NAMELIST("Unknown envelope_solver - only 'explicit' and 'explicit_reduced_dispersion' are available. ",
                           LINK_NAMELIST + std::string("#laser-envelope-model"));
        }
        if (geometry=="1Dcartesian"){
            full_Envelope_exchange = false;
        }

        PyTools::extractVV( "Env_pml_sigma_parameters", envelope_pml_sigma_parameters, "LaserEnvelope" );
        PyTools::extractVV( "Env_pml_kappa_parameters", envelope_pml_kappa_parameters, "LaserEnvelope" );
        PyTools::extractVV( "Env_pml_alpha_parameters", envelope_pml_alpha_parameters, "LaserEnvelope" );
    }

    open_boundaries.resize( nDim_field );
    for( unsigned int iDim = 0 ; iDim < nDim_field; iDim++ ) {
        open_boundaries[iDim].resize( 2, false );
        for( unsigned int j = 0; j < 2; j++ ) {
            if( EM_BCs[iDim][j] == "buneman" ) {
                full_B_exchange = true;
                open_boundaries[iDim][j] = true;
            } else if( EM_BCs[iDim][j] == "silver-muller" ) {
                open_boundaries[iDim][j] = true;
            }
        }
    }

    PyTools::extractVV( "EM_boundary_conditions_k", EM_BCs_k, "Main" );
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
        ERROR_NAMELIST( "EM_boundary_conditions_k must be the same size as the number of faces.",
                        LINK_NAMELIST + std::string("#main-variables"));
    }
    for( unsigned int iDim=0; iDim<nDim_field*2; iDim++ ) {
        if( EM_BCs_k[iDim].size() != nDim_field ) {
            ERROR_NAMELIST( "EM_boundary_conditions_k must have exactly " << nDim_field << " elements along dimension "<<"-+"[iDim%2]<<"012"[iDim/2],
                   LINK_NAMELIST + std::string("#main-variables") );
        }
        if( EM_BCs_k[iDim][iDim/2] == 0. ) {
            ERROR_NAMELIST( "EM_boundary_conditions_k must have a non zero normal component along dimension "<<"-+"[iDim%2]<<"012"[iDim/2],
                           LINK_NAMELIST + std::string("#main-variables") );
        }

    }
    save_magnectic_fields_for_SM = true;
    PyTools::extract( "save_magnectic_fields_for_SM", save_magnectic_fields_for_SM, "Main"   );

    if (use_pml){
        PyTools::extractVV( "number_of_pml_cells", number_of_pml_cells, "Main" );
        if( number_of_pml_cells.size() == 1 ) {
            while( number_of_pml_cells.size() < nDim_field ) {
                number_of_pml_cells.push_back( number_of_pml_cells[0] );
            }
        } else if( number_of_pml_cells.size() != nDim_field ) {
            ERROR_NAMELIST( "number_of_pml_cells must be the same size as the number of dimensions",
                            LINK_NAMELIST + std::string("#main-variables") );
        }
        for( unsigned int iDim=0; iDim<nDim_field; iDim++ ) {
            if( number_of_pml_cells[iDim].size() == 1 ) { // if just one type is specified, then take the same bc type in a given dimension
                number_of_pml_cells[iDim].push_back( number_of_pml_cells[iDim][0] );
            }
        }
    }

    // -----------------------------------
    // POISSON & FILTERING OPTIONS
    // -----------------------------------

    time_fields_frozen=0.0;
    PyTools::extract( "time_fields_frozen", time_fields_frozen, "Main"   );

    // Poisson Solver
    PyTools::extract( "solve_poisson", solve_poisson, "Main"   );
    PyTools::extract( "poisson_max_iteration", poisson_max_iteration, "Main"   );
    PyTools::extract( "poisson_max_error", poisson_max_error, "Main"   );
    // Relativistic Poisson Solver
    PyTools::extract( "solve_relativistic_poisson", solve_relativistic_poisson, "Main"   );
    PyTools::extract( "relativistic_poisson_max_iteration", relativistic_poisson_max_iteration, "Main"   );
    PyTools::extract( "relativistic_poisson_max_error", relativistic_poisson_max_error, "Main"   );

    // Use BTIS3 interpolation method to reduce the effects of numerical Cherenkov radiation
    // This method is detailed in P.-L. Bourgeois and X. Davoine (2023) https://doi.org/10.1017/S0022377823000223
    use_BTIS3 = false;
    PyTools::extract( "use_BTIS3_interpolation", use_BTIS3, "Main"   );
    if (use_BTIS3 && interpolation_order != 2 && (interpolation_order != 1 || geometry != "AMcylindrical" || is_spectral==true )){
        if (geometry=="AMcylindrical"){
            ERROR("B-TIS3 interpolation is not implemented for PSATD solver.");
        } else {
            ERROR("B-TIS3 interpolation is implemented only at order 2 for Cartesian geometries.");
        }
    }
    
    // Current filter properties
    int nCurrentFilter = PyTools::nComponents( "CurrentFilter" );
    for( int ifilt = 0; ifilt < nCurrentFilter; ifilt++ ) {
        PyTools::extract( "model", currentFilter_model, "CurrentFilter", ifilt );
        if( (currentFilter_model != "binomial")&&(currentFilter_model != "customFIR") ) {
            ERROR_NAMELIST( "Currently, only the `binomial` and `customFIR` model is available in CurrentFilter()",
            LINK_NAMELIST + std::string("#current-filtering") );
        }

        if(currentFilter_model == "customFIR") {
            PyTools::extractV( "kernelFIR", currentFilter_kernelFIR, "CurrentFilter", ifilt );
            if( currentFilter_kernelFIR.size() < 3 ) {
                ERROR_NAMELIST( "Kernel have to measure 3 taps at least. For example the binomial FIR kernel on three tapis [0.25,0.50,0.25]",
                                LINK_NAMELIST + std::string("#current-filtering"));
            }
        }

        PyTools::extractV( "passes", currentFilter_passes, "CurrentFilter", ifilt );  //test list

        if( currentFilter_passes.size() == 0 ) {
            ERROR_NAMELIST( "passes in block 'CurrentFilter' cannot be empty",  LINK_NAMELIST + std::string("#current-filtering"));
        } else if( currentFilter_passes.size() == 1 ) {
            while( currentFilter_passes.size() < nDim_field ) {
                currentFilter_passes.push_back( currentFilter_passes[0] );
            }
        } else if( currentFilter_passes.size() != nDim_field ) {
            ERROR_NAMELIST( "passes in block 'CurrentFilter' must be the same size as the number of field dimensions",  LINK_NAMELIST + std::string("#current-filtering") );
        }
    }

    // Field filter properties
    Friedman_filter = false;
    Friedman_theta = 0;
    int nFieldFilter = PyTools::nComponents( "FieldFilter" );
    for( int ifilt = 0; ifilt < nFieldFilter; ifilt++ ) {
        string model;
        PyTools::extract( "model", model, "FieldFilter", ifilt );
        if( model != "Friedman" ) {
            ERROR_NAMELIST( "Currently, only the `Friedman` model is available in FieldFilter()",  LINK_NAMELIST + std::string("#field-filtering"));
        }
        Friedman_filter = true;
        PyTools::extract( "theta", Friedman_theta, "FieldFilter", ifilt );
        if( Friedman_filter && ( Friedman_theta==0. ) ) {
            CAREFUL(0, "Friedman filter is applied but parameter theta is set to zero" );
        }
        if( ( Friedman_theta<0. ) || ( Friedman_theta>1. ) ) {
            ERROR_NAMELIST( "Friedman filter theta = " << Friedman_theta << " must be between 0 and 1",  LINK_NAMELIST + std::string("#field-filtering") );
        }
        if (geometry=="3Dcartesian"){
            ERROR("Friedman filter is not yet supported for `3Dcartesian geometry`");
        }
    }


    // testing the CFL condition
    //!\todo (MG) CFL cond. depends on the Maxwell solv. ==> HERE JUST DONE FOR YEE!!!
    double res_space2=0;
    for( unsigned int i=0; i<nDim_field; i++ ) {
        res_space2 += res_space[i]*res_space[i];
    }
    if( geometry == "AMcylindrical" ) {
        if( !is_spectral ){
            double alpha;
            switch (nmodes){
                case 1: {
                    alpha = 0.210486;
                    break;
                }
                case 2: {
                    alpha = 0.591305;
                    break;
                }
                case 3: {
                    alpha = 3.5234;
                    break;
                }
                case 4: {
                    alpha = 8.51041;
                    break;
                }
                case 5: {
                    alpha = 15.5059;
                    break;
                }
                default:
                    alpha = (nmodes-1)*(nmodes-1)-1.;
            }
            res_space2 += alpha * res_space[1]*res_space[1];
        } else {
            res_space2 = max(res_space[0], res_space[1]) * max(res_space[0], res_space[1]);
            if( timestep != min(cell_length[0], cell_length[1]) ) {
                CAREFUL( 0," timestep=" << timestep << " is not equal to optimal timestep for this solver = " << min(cell_length[0], cell_length[1])  );
            }
        }
    }
    dtCFL=1.0/sqrt( res_space2 );
    if( timestep>dtCFL && !is_spectral ) {
        WARNING( "CFL problem: timestep=" << timestep << " should be smaller than " << dtCFL );
    }

    // cluster_width_
    PyTools::extract( "cluster_width", cluster_width_, "Main"   );


    // --------------------
    // Number of patches
    // --------------------
    if( !PyTools::extractV( "number_of_patches", number_of_patches, "Main" ) ) {
        ERROR_NAMELIST( "The parameter `number_of_patches` must be defined as a list of integers",
        LINK_NAMELIST + std::string("#main-variables")  );
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
            ERROR_NAMELIST( "The total number of patches "<<tot_number_of_patches
                          <<" must be greater or equal to the number of MPI processes "
                          <<smpi->getSize(), LINK_NAMELIST + std::string("#main-variables"));
        }
    }
#ifdef _OPENMP
    if( tot_number_of_patches < ( unsigned int )( smpi->getSize()*smpi->getOMPMaxThreads() ) ) {
        WARNING( "Resources allocated "<<( smpi->getSize()*smpi->getOMPMaxThreads() )
               <<" underloaded regarding the total number of patches "<<tot_number_of_patches );
    }
#endif


    spectral_solver_order.resize( nDim_field, 1 );
    PyTools::extractV( "spectral_solver_order", spectral_solver_order, "Main" );

    initial_rotational_cleaning = false;
    if( is_spectral && geometry == "AMcylindrical" ) {
        PyTools::extract( "initial_rotational_cleaning", initial_rotational_cleaning, "Main" );
        if( initial_rotational_cleaning && smpi->getSize() > 1 ) {
            CAREFUL(0,"Rotational cleaning (laser initialization) is not parallelized for now and may use a large amount of memory.");
        }
    }

    PyTools::extract( "patch_arrangement", patch_arrangement, "Main"  );
    CAREFUL( 0,"Patches distribution: " << patch_arrangement );

    int total_number_of_hilbert_patches = 1;
    if( patch_arrangement == "hilbertian" ) {
        for( unsigned int iDim=0 ; iDim<nDim_field ; iDim++ ) {
            total_number_of_hilbert_patches *= number_of_patches[iDim];
            if( ( number_of_patches[iDim] & ( number_of_patches[iDim]-1 ) ) != 0 ) {
                ERROR_NAMELIST( "Number of patches in each direction must be a power of 2",  LINK_NAMELIST + std::string("#main-variables")  );
            }
        }
    }




    if( PyTools::nComponents( "LoadBalancing" )>0 ) {
        // get parameter "every" which describes a timestep selection
        load_balancing_time_selection = new TimeSelection(
            PyTools::extract_py( "every", "LoadBalancing" ), "Load balancing"
        );
        PyTools::extract( "cell_load", cell_load, "LoadBalancing"   );
        PyTools::extract( "frozen_particle_load", frozen_particle_load, "LoadBalancing"   );
        PyTools::extract( "initial_balance", initial_balance, "LoadBalancing"   );
    } else {
        load_balancing_time_selection = new TimeSelection();
    }

    has_load_balancing = ( smpi->getSize()>1 )  && ( ! load_balancing_time_selection->isEmpty() );

    if( has_load_balancing && patch_arrangement != "hilbertian" ) {
        ERROR_NAMELIST( "Dynamic load balancing is only available for Hilbert decomposition",  LINK_NAMELIST + std::string("#main-variables") );
    }
    if( has_load_balancing && total_number_of_hilbert_patches < 2*smpi->getSize() ) {
        ERROR_NAMELIST( "Dynamic load balancing requires to use at least 2 patches per MPI process.",  LINK_NAMELIST + std::string("#main-variables") );
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

    bool defined_cell_sort = true;
    if (!PyTools::extractOrNone( "cell_sorting", cell_sorting_, "Main"  )){
    //cell_sorting is undefined by the user
        defined_cell_sort = false;
        cell_sorting_ = false;
    }

    // Activation of the vectorized subroutines
    vectorization_mode = "off";
    has_adaptive_vectorization = false;
    adaptive_vecto_time_selection = nullptr;

    if( PyTools::nComponents( "Vectorization" )>0 ) {
        // Extraction of the vectorization mode
        PyTools::extract( "mode", vectorization_mode, "Vectorization"   );
        if( !( vectorization_mode == "off" ||
                vectorization_mode == "on" ||
                vectorization_mode == "adaptive" ||
                vectorization_mode == "adaptive_mixed_sort" ) ) {
            ERROR_NAMELIST( "In block `Vectorization`, parameter `mode` must be `off`, `on`, `adaptive`",  LINK_NAMELIST + std::string("#vectorization") );
        } else if( vectorization_mode == "adaptive_mixed_sort" || vectorization_mode == "adaptive" ) {
            has_adaptive_vectorization = true;
        }

        if (geometry=="1Dcartesian" && vectorization_mode != "off") {
            vectorization_mode = "off";
            has_adaptive_vectorization = false;
            WARNING("In 1D, the vectorization block does not apply. `vectorization back to `off`.")
        }

        if (use_BTIS3 && vectorization_mode != "off") {
            ERROR("B-TIS3 interpolator not yet implemented in vectorized mode.")
        }
        
        // Cell sorting not defined by the user
        if (!defined_cell_sort) {
            if (vectorization_mode == "off") {
                cell_sorting_ = false;
            } else {
                cell_sorting_ = true;
            }
        }

        // Cell sorting explicitely defined by the user
	    if (defined_cell_sort){
            // cell sorting explicitely set on
            if (cell_sorting_) {
                if (vectorization_mode == "off") {
                    WARNING(" Cell sorting `cell_sorting` cannot be used when vectorization is off for the moment. Vectorization is automatically activated.")
                }
            // cell sorting explicitely set off
            } else {
                if (!( vectorization_mode == "off")) {
                    ERROR_NAMELIST(" Cell sorting `cell_sorting` must be allowed in order to use vectorization.",
                        LINK_NAMELIST + std::string("#vectorization"))
                }
            }
        }


        // Default mode for the adaptive mode
        PyTools::extract( "initial_mode", adaptive_default_mode, "Vectorization"   );
        if( !( adaptive_default_mode == "off" ||
                adaptive_default_mode == "on" ) ) {
            ERROR_NAMELIST( "In block `Vectorization`, parameter `initial_mode` must be `off` or `on`",  LINK_NAMELIST + std::string("#vectorization") );
        }

        // get parameter "every" which describes a timestep selection
        if( ! adaptive_vecto_time_selection )
            adaptive_vecto_time_selection = new TimeSelection(
                PyTools::extract_py( "reconfigure_every", "Vectorization" ), "Adaptive vectorization"
            );
    }

    PyTools::extract( "gpu_computing", gpu_computing, "Main" );
    if( gpu_computing ) {
#if( defined( SMILEI_ACCELERATOR_GPU_OACC ) && defined( _OPENACC ) ) || defined( SMILEI_ACCELERATOR_GPU_OMP )
        // If compiled for GPU and asking for GPU
        MESSAGE( 1, "Smilei will run on GPU devices" );
#else
        // If compiled for CPU and asking for GPU
        ERROR( "Smilei is not compiled for GPU" );
#endif
    } else {
#if defined( _OPENACC ) || defined( SMILEI_ACCELERATOR_GPU_OMP )
        // If compiled for GPU and asking for CPU
        ERROR( "Smilei needs to be executed on GPU, set Main.gpu_computing = True" );
#else
        // If compiled for CPU and asking for CPU
        MESSAGE( 1, "Smilei will run on CPU devices" );
#endif
    }

    // In case of collisions, ensure particle sort per cell
    if( PyTools::nComponents( "Collisions" ) > 0 ) {

        // collisions need sorting per cell
        if (defined_cell_sort && cell_sorting_ == false){
            ERROR_NAMELIST(" Cell sorting or vectorization must be allowed in order to use collisions.",  LINK_NAMELIST + std::string("#collisions-reactions"));
        }

        // Force cell sorting and later the adaptive vectorization mode in scalar mode
        cell_sorting_ = true;

        if( geometry!="1Dcartesian"
                && geometry!="2Dcartesian"
                && geometry!="3Dcartesian" ) {
            //ERROR_NAMELIST( "Collisions only valid for cartesian geometries for the moment",  LINK_NAMELIST + std::string("#collisions-reactions") );
            WARNING( "Collisions in AM geometry is experimental and valid only with a single mode" );
        }

    }

    // Read the "print_every" parameter
    print_every = ( int )( simulation_time/timestep )/10;
    PyTools::extractOrNone( "print_every", print_every, "Main" );
    if( !print_every ) {
        print_every = 1;
    }

    // Read the "print_expected_disk_usage" parameter
    PyTools::extract( "print_expected_disk_usage", print_expected_disk_usage, "Main"   );

    // Decide when necessary to keep position_old
    keep_position_old = false;
    DEBUGEXEC( keep_position_old = true );

    // -------------------------------------------------------
    // Checking species order
    // -------------------------------------------------------
    // read from python namelist the number of species
    unsigned int tot_species_number = PyTools::nComponents( "Species" );

    double mass, mass2=0;
    std::string merging_method;

    for( unsigned int ispec = 0; ispec < tot_species_number; ispec++ ) {
        PyTools::extract( "mass", mass, "Species", ispec );
        if( mass == 0 ) {
            for( unsigned int ispec2 = ispec+1; ispec2 < tot_species_number; ispec2++ ) {
                PyTools::extract( "mass", mass2, "Species", ispec2 );
                if( mass2 > 0 ) {
                    ERROR_NAMELIST( "the photon species (mass==0) should be defined after the particle species (mass>0)",  LINK_NAMELIST + std::string("#species") );
                }
            }
        }
        //Use cell sorting if merge is used.
        PyTools::extract( "merging_method", merging_method, "Species", ispec );
        if (merging_method != "none"){

            if (defined_cell_sort && !cell_sorting_){
                ERROR_NAMELIST(" Cell sorting or vectorization must be allowed in order to use particle merging.",  LINK_NAMELIST + std::string("#collisions-reactions"));
            }

            // Force cell sorting and later the adaptive vectorization mode in scalar mode
            cell_sorting_ = true;

        }
    }

    // Force adaptive vectorization in scalar mode if cell_sorting requested
    if ( cell_sorting_ ) {

        if( vectorization_mode == "adaptive_mixed_sort" ) {
            ERROR_NAMELIST( "Cell sorting (required by Collision or Merging) is incompatible with the vectorization mode 'adaptive_mixed_sort'.",  LINK_NAMELIST + std::string("#vectorization") );
        } else if ( vectorization_mode == "off" ) {
            vectorization_mode            = "adaptive";
            has_adaptive_vectorization    = true;
            adaptive_default_mode         = "off";
            adaptive_vecto_time_selection = new TimeSelection();
        }
    }

    // -------------------------------------------------------
    // Parameters for the synchrotron-like radiation losses
    // -------------------------------------------------------
    has_MC_radiation_ = false ;// Default value
    has_LL_radiation_ = false ;// Default value
    has_Niel_radiation_ = false ;// Default value
    has_diag_radiation_spectrum_ = false; // Default value

    // Loop over all species to check if the radiation losses are activated
    std::string radiation_model = "none";
    for( unsigned int ispec = 0; ispec < tot_species_number; ispec++ ) {

        PyTools::extract( "radiation_model", radiation_model, "Species", ispec );

        // Cancelation of the letter case for `radiation_model`
        std::transform( radiation_model.begin(), radiation_model.end(), radiation_model.begin(), ::tolower );

        if( radiation_model=="monte-carlo" || radiation_model=="mc" ) {
            has_MC_radiation_ = true;
        } else if( radiation_model=="landau-lifshitz"
                   || radiation_model=="ll"
                   || radiation_model=="corrected-landau-lifshitz"
                   || radiation_model=="cll" ) {
            has_LL_radiation_ = true;
        } else if( radiation_model=="niel" ) {
            has_Niel_radiation_ = true;
        }
        else if (radiation_model=="diagradiationspectrum")
        {
            has_diag_radiation_spectrum_ = true;
        }
    }

    // -------------------------------------------------------
    // Parameters for the multiphoton Breit-Wheeler pair decay
    // -------------------------------------------------------
    has_multiphoton_Breit_Wheeler_ = false ;// Default value
    std::vector<std::string> multiphoton_Breit_Wheeler( 2 );
    for( unsigned int ispec = 0; ispec < tot_species_number; ispec++ ) {
        if( PyTools::extractV( "multiphoton_Breit_Wheeler", multiphoton_Breit_Wheeler, "Species", ispec ) ) {
            has_multiphoton_Breit_Wheeler_ = true;
        }
    }

    // -------------------------------------------------------
    // Compute useful quantities and introduce normalizations
    // also defines defaults values for the species lengths
    // -------------------------------------------------------
    compute();

    // add the read or computed value of cluster_width_ to the content of smilei.py
    namelist += string( "Main.cluster_width= " ) + to_string( cluster_width_ ) + "\n";

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

    for( unsigned int i_laser=0; i_laser<n_laser; i_laser++ ) {
        double offset = 0.;

        // If this laser has the hidden _offset attribute
        if( PyTools::extractOrNone( "_offset", offset, "Laser", i_laser ) ) {

            bool propagate = false;
            PyTools::extract( "_propagate", propagate, "Laser", i_laser );
            if( ! propagate ) continue;

            // Extract the angle
            double angle_z = 0.;
            PyTools::extract( "_angle", angle_z, "Laser", i_laser );

            // Extract _fft_time_window
            double fft_time_window = 0.;
            PyTools::extract( "_fft_time_window", fft_time_window, "Laser", i_laser );

            // Extract _fft_time_step
            double fft_time_step = 0.;
            PyTools::extract( "_fft_time_step", fft_time_step, "Laser", i_laser );

            // Extract _number_of_processes
            int number_of_processes = 0;
            MPI_Comm comm;
            if( PyTools::extractOrNone( "_number_of_processes", number_of_processes, "Laser", i_laser ) ) {
                int color = smpi->getRank() < number_of_processes ? 1 : MPI_UNDEFINED;
                MPI_Comm_split( smpi->world(), color, smpi->getRank(), &comm );
            } else {
                number_of_processes = smpi->getSize();
                MPI_Comm_dup( smpi->world(), &comm );
            }

            // Extract the file name
            string file( "" );
            PyTools::extract( "file", file, "Laser", i_laser );

            // Extract the list of profiles and verify their content
            PyObject *p = PyTools::extract_py( "_profiles", "Laser", i_laser );
            vector<PyObject *> profiles;
            if( ! PyTools::py2pyvector( p, profiles ) ) {
                ERROR_NAMELIST( "For LaserOffset #" << n_laser_offset << ": space_time_profile must be a list of 2 profiles",  LINK_NAMELIST + std::string("#lasers") );
            }
            Py_DECREF( p );
            if( profiles.size() != 2 ) {
                ERROR_NAMELIST( "For LaserOffset #" << n_laser_offset << ": space_time_profile needs 2 profiles.",  LINK_NAMELIST + std::string("#lasers") );
            }
            vector<int> profiles_n;
            vector<PyObject *> profiles_kept;
            for( unsigned int i = 0; i < 2; i++ ) {
                if( profiles[i] != Py_None ) {
                    profiles_kept.push_back( profiles[i] );
                    profiles_n.push_back( i + 1 );
                }
            }
            if( profiles_kept.size() == 0 ) {
                ERROR_NAMELIST( "For LaserOffset #" << n_laser_offset << ": space_time_profile cannot be [None, None]", LINK_NAMELIST + std::string("#lasers") );
            }
            for( unsigned int i=0; i<profiles_kept.size(); i++ ) {
                int nargs = PyTools::function_nargs( profiles_kept[i] );
                if( nargs == -2 ) {
                    ERROR_NAMELIST( "For LaserOffset #" << n_laser_offset << ": space_time_profile["<<i<<"] not callable", LINK_NAMELIST + std::string("#lasers") );
                }
                if( nargs >= 0 && nargs != ( int ) nDim_field ) {
                    ERROR_NAMELIST( "For LaserOffset #" << n_laser_offset << ": space_time_profile["<<i<<"] requires " << nDim_field << " arguments but has " << nargs, LINK_NAMELIST + std::string("#lasers") );
                }
            }

            // Extract _keep_n_strongest_modes
            int keep_n_strongest_modes=0;
            PyTools::extract( "_keep_n_strongest_modes", keep_n_strongest_modes, "Laser", i_laser );
            if( keep_n_strongest_modes<1 ) {
                ERROR_NAMELIST( "For LaserOffset #" << n_laser_offset << ": keep_n_strongest_modes must be a positive integer", LINK_NAMELIST + std::string("#lasers") );
            }

            // Extract box_side
            string box_side = "";
            PyTools::extract( "box_side", box_side, "Laser", i_laser );
            unsigned int normal_axis = 0;
            if( box_side == "xmin" || box_side == "xmax" ) {
                normal_axis = 0;
            } else if( box_side == "ymin" || box_side == "ymax" ) {
                if( geometry != "2Dcartesian" && geometry != "3Dcartesian" ) {
                    ERROR_NAMELIST( "For LaserOffset #" << n_laser_offset << ": box_side `ymin` or `ymax` requires 2D or 3D geometry", LINK_NAMELIST + std::string("#lasers") );
                }
                normal_axis = 1;
            } else if( box_side == "zmin" || box_side == "zmax" ) {
                if( geometry != "3Dcartesian" ) {
                    ERROR_NAMELIST( "For LaserOffset #" << n_laser_offset << ": box_side `zmin` or `zmax` requires 3D geometry", LINK_NAMELIST + std::string("#lasers") );
                }
                normal_axis = 2;
            } else {
                ERROR_NAMELIST( "For LaserOffset #" << n_laser_offset << ": box_side must be `xmin`, `xmax`, `ymin`, `ymax`, `zmin` or `zmax`", LINK_NAMELIST + std::string("#lasers") );
            }

            if( smpi->getRank() < number_of_processes && ! restart ) {
                if( n_laser_offset == 0 ) {
                    TITLE( "Calculate LaserOffset" );
                }
                // Prepare propagator
                MESSAGE( 1, "LaserOffset #"<< n_laser_offset );
                LaserPropagator propagateX( this, normal_axis, fft_time_window, fft_time_step, comm );

                // Make the propagation happen and write out the file
                if( ! smpi->test_mode ) {
                    propagateX( profiles_kept, profiles_n, offset, file, keep_n_strongest_modes, angle_z );
                }
            }
            
            for( auto p: profiles ) {
                Py_DECREF( p );
            }
            
            n_laser_offset ++;
        }
    }

    check_consistency();
    
    
    // Run the _writeInfo function that creates a small pickle file with basic info
    if( ! smpi->test_mode && smpi->isMaster() ) {
        PyTools::runPyFunction( "writeInfo" );
        PyTools::checkPyError();
    }
}

Params::~Params()
{
    if( load_balancing_time_selection ) {
        delete load_balancing_time_selection;
    }
    if( adaptive_vecto_time_selection ) {
        delete adaptive_vecto_time_selection;
    }
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
    patch_size_.resize( 3, 1 );
    cell_length.resize( 3 );
    global_size_.resize( 3, 1 );
    oversize.resize( 3, 0 );
    region_oversize.resize( 3, 0 );
    patch_dimensions.resize( 3, 0. );
    cell_volume=1.0;
    n_cell_per_patch = 1;

    multiple_decomposition = PyTools::nComponents( "MultipleDecomposition" )>0;

    // compute number of cells & normalized lengths
    for( unsigned int i=0; i<nDim_field; i++ ) {
        patch_size_[i] = round( grid_length[i]/cell_length[i] );
        double entered_grid_length = grid_length[i];
        grid_length[i] = ( double )( patch_size_[i] )*cell_length[i]; // ensure that nspace = grid_length/cell_length
        if( grid_length[i]!=entered_grid_length ) {
            WARNING( "grid_length[" << i << "] has been redefined from " << entered_grid_length << " to " << grid_length[i] << " to match n x cell_length (" << scientific << setprecision( 4 ) << grid_length[i]-entered_grid_length <<")" );
        }
        cell_volume *= cell_length[i];
    }
    if( geometry == "AMcylindrical" ) {
        cell_volume *= 2 * M_PI;
    }
    // create a 3d equivalent of patch_size_ & cell_length
    for( unsigned int i=nDim_field; i<3; i++ ) {
        cell_length[i]=0.0;
    }

    //Define number of cells per patch and number of ghost cells
    for( unsigned int i=0; i<nDim_field; i++ ) {
        PyTools::extract( "custom_oversize", custom_oversize, "Main"  );
        if (maxwell_sol == "Bouchard" && custom_oversize < 4 ) {
             ERROR_NAMELIST( "With `Bouchard` solver the oversize have to be greater than 4", LINK_NAMELIST + std::string("#main-variables") );
        }
        if( ! multiple_decomposition ) {
            oversize[i]  = std::max( interpolation_order, std::max( ( unsigned int )( spectral_solver_order[i]/2+1 ),custom_oversize ) ) + ( exchange_particles_each-1 );
            if( currentFilter_model == "customFIR" && oversize[i] < (currentFilter_kernelFIR.size()-1)/2 ) {
                ERROR_NAMELIST( "With the `customFIR` current filter model, the ghost cell number (oversize) = " << oversize[i] << " have to be >= " << (currentFilter_kernelFIR.size()-1)/2 << ", the (kernelFIR size - 1)/2", LINK_NAMELIST + std::string("#current-filtering")  );
            }
        } else {
            oversize[i] = interpolation_order + ( exchange_particles_each-1 );
        }
        global_size_[i] = patch_size_[i];
        patch_size_[i] /= number_of_patches[i];
        if( global_size_[i]%number_of_patches[i] !=0 ) {
            ERROR_NAMELIST( "ERROR in dimension " << i <<". Number of patches = " << number_of_patches[i] << " must divide global_size_ = " << global_size_[i], LINK_NAMELIST + std::string("#main-variables") );
        }
        if( patch_size_[i] <= 2*oversize[i]+1 ) {
            ERROR_NAMELIST( "ERROR in dimension " << i <<". Patches length = "<<patch_size_[i] << " cells must be at least " << 2*oversize[i] +2 << " cells long. Increase number of cells or reduce number of patches in this direction. ",LINK_NAMELIST + std::string("#main-variables") );
        }
        patch_dimensions[i] = patch_size_[i] * cell_length[i];
        n_cell_per_patch *= patch_size_[i];
    }

    // Set cluster_width_ if not set by the user
    if( cluster_width_ == -1 ) {
#if defined( SMILEI_ACCELERATOR_GPU )
        cluster_width_ = patch_size_[0];
        // On GPU, dont do the CPU automatic cluster_width computation, only one
        // bin is expected.
        // NOTE: In OMP GPU offloading and 2D, the true number of cluster is
        // redefined in nvidiaParticles::prepareBinIndex.
#else
        // default value
        cluster_width_ = patch_size_[0];

        // check cache issue for interpolation/projection
        int cache_threshold( 3200 ); // sizeof( L2, Sandy Bridge-HASWELL ) / ( 10 * sizeof(double) )
        // Compute the "transversal bin size"
        int bin_size( 1 );
        for( unsigned int idim = 1 ; idim < nDim_field ; idim++ ) {
            bin_size *= ( patch_size_[idim]+1+2*oversize[idim] );
        }

        // IF Ionize or pair generation : cluster_width_ = n_space_x_pp ?
        if( ( cluster_width_+1+2*oversize[0] ) * bin_size > ( unsigned int ) cache_threshold ) {
            int clrw_max = cache_threshold / bin_size - 1 - 2*oversize[0];
            if( clrw_max > 0 ) {
                for( cluster_width_=clrw_max ; cluster_width_ > 0 ; cluster_width_-- )
                    if( ( ( cluster_width_+1+2*oversize[0] ) * bin_size <= ( unsigned int ) cache_threshold ) && ( patch_size_[0]%cluster_width_==0 ) ) {
                        break;
                    }
            } else {
                cluster_width_ = 1;
            }
            WARNING( "Particles cluster width `cluster_width` set to : " << cluster_width_ );
        }
#endif
    } else {
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
        ERROR( "Cluster width size is not specifiable in OpenMP GPU mode! " );
#endif
    }

    // cluster_width_ != patch_size_[0] is not compatible
    // with the adaptive vectorization for the moment
    if( vectorization_mode == "adaptive_mixed_sort" || vectorization_mode == "adaptive" ) {
        if( cluster_width_ != ( int )( patch_size_[0] ) ) {
            cluster_width_ = ( int )( patch_size_[0] );
            WARNING( "Particles cluster width set to: " << cluster_width_ << " for the adaptive vectorization mode" );
        }
    }


    // Verify that cluster_width_ divides patch_size_[0] or patch_size_[n] in GPU mode
#if defined( SMILEI_ACCELERATOR_GPU_OMP ) || defined( SMILEI_ACCELERATOR_GPU_OACC )
    const int kClusterWidth = getGPUClusterWidth();

    if( kClusterWidth < 0 ) {
        // getGPUClusterWidth failing means that binning is not supported for
        // nDim_particle space dimensions. cluster_width_ shall default to
        // patch_size_[0] which is guarenteed to divide itself, no need to check.
    } else {
        for( std::size_t dimension_id = 0; dimension_id < nDim_particle; ++dimension_id ) {
            if( ( patch_size_[dimension_id] % kClusterWidth ) != 0 ) {
                ERROR_NAMELIST( "The parameter `cluster_width`==" << kClusterWidth << " must divide the number of cells in a patch, in all dimensions.",
                                LINK_NAMELIST + std::string( "#main-variables" ) );
            }
        }
    }
#else
    if( ( patch_size_[0] % cluster_width_ ) != 0 ) {
        ERROR_NAMELIST( "The parameter `cluster_width` must divide the number of cells in one patch (in dimension x)",
                        LINK_NAMELIST + std::string( "#main-variables" ) );
    }
#endif

    // Define domain decomposition if double grids are used for particles and fields
    region_size_ = patch_size_; // by default, region sizes are set to those of the patch
    region_oversize = oversize;
    if( multiple_decomposition ) {
        multiple_decompose();
        full_B_exchange = true;
    }
}


void Params::check_consistency()
{
    if( vectorization_mode != "off" ) {

        if( ( geometry=="1Dcartesian" ) ) {
            WARNING( "Vectorized and scalar algorithms are the same in 1D Cartesian geometry." );
        }

        // if( ( geometry=="2Dcartesian" ) && ( interpolation_order==4 ) ) {
        //     ERROR( "4th order vectorized algorithms not implemented in 2D" );
        // }
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
        ERROR_NAMELIST( "Geometry: " << geometry << " not defined", LINK_NAMELIST + std::string("#main-variables") );
    }
}



// ---------------------------------------------------------------------------------------------------------------------
// Printing out the data at initialisation
// ---------------------------------------------------------------------------------------------------------------------
void Params::print_init()
{
    if( full_B_exchange ) {
        MESSAGE( 1, "All components of B are exchanged at synchronization" );
    }

    TITLE( "Geometry: " << geometry );
    MESSAGE( 1, "Interpolation order : " <<  interpolation_order );
    if (use_BTIS3){
        MESSAGE(1, "B-TIS3 interpolation scheme activated")
    }
    MESSAGE( 1, "Maxwell solver : " <<  maxwell_sol );
    MESSAGE( 1, "simulation duration = " << simulation_time <<",   total number of iterations = " << n_time);
    MESSAGE( 1, "timestep = " << timestep << " = " << timestep/dtCFL << " x CFL,   time resolution = " << res_time);

    ostringstream gl;
    gl << "Grid length: ";
    for( unsigned int i=0 ; i<grid_length.size() ; i++ ) {
        gl << grid_length[i] << ( i<grid_length.size()-1 ? ", " : "" );
    }
    MESSAGE( 1, gl.str() );

    ostringstream cl;
    cl << "Cell length: ";
    for( unsigned int i=0 ; i<cell_length.size() ; i++ ) {
        cl << cell_length[i] << ( i<cell_length.size()-1 ? ", " : "" );
    }
    MESSAGE( 1, cl.str() );

    ostringstream nc;
    nc << "Number of cells: " ;
    for( unsigned int i=0 ; i<nDim_field ; i++ ) {
        nc << global_size_[i] << ( i<nDim_field-1 ? ", " : "" );
    }
    MESSAGE( 1, nc.str() );

    ostringstream sr;
    sr << "Spatial resolution: ";
    for( unsigned int i=0 ; i<nDim_field ; i++ ) {
        sr << res_space[i] << ( i<nDim_field-1 ? ", " : "" );
    }
    MESSAGE( 1, sr.str() );

    if (cell_sorting_) {
        MESSAGE( 1, "Cell sorting: activated" );
    }

    TITLE( "Electromagnetic boundary conditions" );
    string xyz = geometry=="AMcylindrical" ? "xr" : "xyz";
    for( unsigned int i=0 ; i<grid_length.size() ; i++ ) {
        for( unsigned int j=0 ; j<2 ; j++ ) {
            ostringstream bc( "" );
            bc << xyz[i] << (  j==0 ? "min " : "max ") << EM_BCs[i][j];
            if( open_boundaries[i][j] ) {
                bc << setprecision( 2 ) << ", absorbing vector " << "[" << EM_BCs_k[2*i+j][0];
                for( unsigned int ii=1 ; ii<grid_length.size() ; ii++ ) {
                    bc << ", " << EM_BCs_k[2*i+j][ii] ;
                }
                bc << "]";
            }
            if (EM_BCs[i][j] == "PML"){
               bc << "  " << number_of_pml_cells[i][j] << " cells.";
            }
            MESSAGE( 1, bc.str() );
        }
    }

    if ( Laser_Envelope_model ) {
        TITLE( "Laser Envelope parameters" );
        ostringstream info( "" );
        info << "\tpolarization angle : " << envelope_polarization_phi << endl;
        info << "\t\tellipticity        : " << envelope_ellipticity << endl;
        info << "\t\tEnvelope solver    : " << envelope_solver << endl;
        MESSAGE( 1, info.str() );
        for( unsigned int i=0 ; i<grid_length.size() ; i++ ) {
            MESSAGE( 1, "\tdimension " << i );
            MESSAGE( 1, "\t- Envelope boundary conditions: " << "(" << Env_BCs[i][0] << ", " << Env_BCs[i][1] << ")" );
        }
    }


    if( currentFilter_passes.size() > 0 ){
        TITLE( "Current filtering" );
        if( *std::max_element(std::begin(currentFilter_passes), std::end(currentFilter_passes)) > 0 ) {
            for( unsigned int idim=0 ; idim < nDim_field ; idim++ ){
                std::string strpass = (currentFilter_passes[idim] > 1 ? "passes" : "pass");
                MESSAGE( 1, currentFilter_model << " current filtering: " << currentFilter_passes[idim] << " " << strpass << " along dimension " << idim );
            }
        }
    }
    if( Friedman_filter ) {
        TITLE( "Field filtering" );
        MESSAGE( 1, "Friedman field filtering : theta = " << Friedman_theta );
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
void Params::print_timestep( SmileiMPI *smpi, unsigned int itime, double time_dual, Timer &timer, double npart )
{
    if( smpi->isMaster() ) {
        double before = timer.getTime();
        timer.update();
        double now = timer.getTime();

        ostringstream push_time;
        if( npart * (double) print_every > 0 ) {
            push_time << setw( 14 ) << (int)( 1e9 * (now - before) * (double) smpi->getGlobalNumCores() / ( npart * (double) print_every ) );
        } else {
            push_time << "  ??";
        }

        #pragma omp master
        MESSAGE(
            "  " << setw( timestep_width ) << itime << "/" << n_time << " "
            << "  " << scientific << setprecision( 4 ) << setw( 12 ) << time_dual << " "
            << "  " << scientific << setprecision( 4 ) << setw( 12 ) << now << " "
            << "  " << "(" << scientific << setprecision( 4 ) << setw( 12 ) << now - before << " )"
            << "  " << push_time.str() << " "
        );
        #pragma omp barrier
    }
}

void Params::print_timestep_headers( SmileiMPI *smpi )
{
    timestep_width = log10( n_time ) + 1;
    if( timestep_width<3 ) {
        timestep_width = 3;
    }
    CAREFUL(0,"The following `push time` assumes a global number of "<<smpi->getGlobalNumCores()<<" cores (hyperthreading is unknown)" );
    MESSAGE(
        setw( timestep_width*2+4 ) << " timestep "
        << setw( 15 ) << "sim time "
        << setw( 15 ) << "cpu time [s] "
        << "  (" << setw( 12 ) << "diff [s]" << " )"
        << setw( 17 ) << "   push time [ns]"
    );
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

        MESSAGE( 1, "Number of MPI processes: " << smpi->getSize() );

#ifdef _OPENMP
        MESSAGE( 1, "Number of threads per MPI process : " << smpi->getOMPMaxThreads() );
#else
        MESSAGE( 1, "OpenMP disabled" );
#endif

        ostringstream np;
        np << "Number of patches: " << number_of_patches[0];
        for( unsigned int iDim=1 ; iDim<nDim_field ; iDim++ ) {
            np << " x " << number_of_patches[iDim];
        }
        MESSAGE( 1, np.str() );

        ostringstream ps;
        ps << "Number of cells in one patch: " << patch_size_[0];
        for( unsigned int iDim=1 ; iDim<nDim_field ; iDim++ ) {
            ps << " x " << patch_size_[iDim];
        }
        MESSAGE( 1, ps.str() );

        MESSAGE( 1, "Dynamic load balancing: " << load_balancing_time_selection->info() );
    }

}

// ---------------------------------------------------------------------------------------------------------------------
// Finds requested species in the list of existing species.
// Returns an array of the numbers of the requested species.
// Note that there might be several species that have the same "name"
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
        ERROR( "error parsing "<< name << "\n Check out the namelist options: " << LINK_NAMELIST );
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

    if( keep_python_running_ ) {
        MESSAGE( 1, "Keeping Python interpreter alive" );
    } else {
        MESSAGE( 1, "Closing Python" );
        PyErr_Print();
        Py_Finalize();
    }
    smpi->barrier();
}


void Params::multiple_decompose()
{
    // Compute the oversize of the region
    if( is_spectral ) {
        for( unsigned int i=0; i<nDim_field; i++ ){
            region_oversize[i]  = std::max( interpolation_order, ( unsigned int )( spectral_solver_order[i]/2+1 ) ) + ( exchange_particles_each-1 );
        }
    } else {
        for( unsigned int i=0; i<nDim_field; i++ ){
            region_oversize[i]  = interpolation_order + ( exchange_particles_each-1 );
        }
    }
    PyTools::extract( "region_ghost_cells", region_ghost_cells, "MultipleDecomposition" );
    for( unsigned int i=0; i<nDim_field; i++ ) {
        region_oversize[i] = std::max( region_oversize[i], region_ghost_cells );
    }
    if( is_spectral && geometry == "AMcylindrical" )  {
        //Force ghost cells number in L when spectral
        region_oversize[0] = region_ghost_cells;
        //Force zero ghost cells in R when spectral
        WARNING("Forcing region ghost-cell size along r from " << region_oversize[1] << " to " <<  oversize[1])
        region_oversize[1] = oversize[1];
    }

    number_of_region.resize( 3, 1 );

    int rk(0);
    MPI_Comm_rank( MPI_COMM_WORLD, &rk );

    if( nDim_field==1 ) {
        multiple_decompose_1D();
    } else if( nDim_field==2 ) {
        multiple_decompose_2D();
    } else if( nDim_field==3 ) {
        multiple_decompose_3D();
    }

    map_rank.resize( number_of_region[0] );
    for( unsigned int ix = 0 ; ix < number_of_region[0] ; ix++ ) {
        map_rank[ix].resize( number_of_region[1] );
        for( unsigned int iy = 0 ; iy < number_of_region[1] ; iy++ ) {
            map_rank[ix][iy].resize( number_of_region[2] );
        }
    }

    int new_rk(0);
    unsigned int ijk[3];
    region_coordinates.resize( nDim_field );
    // - Build the map of MPI ranks in 3D
    // - Set the coordinates of the current region

    for(ijk[0] = 0 ; ijk[0] < number_of_region[0] ; ijk[0]++ ) {
        for( ijk[1] = 0 ; ijk[1] < number_of_region[1] ; ijk[1]++ ) {
            for( ijk[2] = 0 ; ijk[2] < number_of_region[2] ; ijk[2]++ ) {
                map_rank[ijk[0]][ijk[1]][ijk[2]] = new_rk;
                if( new_rk == rk ) {
                    for( unsigned int iDim = 0; iDim < nDim_field; iDim++ ) {
                        region_coordinates[iDim] = ijk[iDim];
                    }
                }
                new_rk++;
            }
        }
    }

    // Build the map of offset, contains offset for each domain, expressed in number of cells
    offset_map.resize( nDim_field );
    for( unsigned int iDim = 0 ; iDim < nDim_field ; iDim++ ) {
        offset_map[iDim].resize( number_of_region[iDim] );
        int nlocal_i = number_of_patches[iDim] / number_of_region[iDim];
        for( unsigned int iDom = 0 ; iDom < number_of_region[iDim] ; iDom++ ) {
            offset_map[iDim][iDom] = iDom * nlocal_i * patch_size_[iDim];
        }
    }

    // Compute size of local domain
    for( unsigned int iDim = 0 ; iDim < nDim_field ; iDim++ ) {
        if( region_coordinates[iDim] != (int)number_of_region[iDim]-1 ) {
            region_size_[iDim] = offset_map[iDim][region_coordinates[iDim]+1] - offset_map[iDim][region_coordinates[iDim]];
        }
        else {
            region_size_[iDim] = global_size_[iDim] - offset_map[iDim][region_coordinates[iDim]];
        }
    }

    print_multiple_decomposition_params();
}


void Params::print_multiple_decomposition_params()
{
    int rk(0);
    int sz(1);

    MPI_Comm_rank( MPI_COMM_WORLD, &rk );
    MPI_Comm_size( MPI_COMM_WORLD, &sz );

    if (rk==0) {
        cout << "Number of regions : ";
        for ( unsigned int iDim  = 0 ; iDim < nDim_field ; iDim++ )
            cout << number_of_region[iDim] << " ";
        cout << endl;
    }
    MPI_Barrier( MPI_COMM_WORLD );
    std::cout << std::flush;

    for ( int irk = 0 ; irk < sz ; irk++ )  {
        if ( irk == rk) {
            cout << " MPI_rank = " << rk << endl;
            cout << "\tcoords = ";
            for ( unsigned int iDim  = 0 ; iDim < nDim_field ; iDim++ ) cout << region_coordinates[iDim] << " ";
            cout << endl;
            cout << "\tsize :  ";
            for ( unsigned int iDim  = 0 ; iDim < nDim_field ; iDim++ ) cout << region_size_[iDim] << " ";
            cout << endl;
        }
        MPI_Barrier( MPI_COMM_WORLD );
        std::cout << std::flush;
    }
}

void Params::multiple_decompose_1D()
{
    int rk(0);
    int sz(1);

    MPI_Comm_rank( MPI_COMM_WORLD, &rk );
    MPI_Comm_size( MPI_COMM_WORLD, &sz );

    // Number of domain in 1D
    number_of_region[0] = sz;

}


void Params::multiple_decompose_2D()
{
    int rk(0);
    int sz(1);

    MPI_Comm_rank( MPI_COMM_WORLD, &rk );
    MPI_Comm_size( MPI_COMM_WORLD, &sz );


    if ( geometry != "AMcylindrical" || !is_spectral ) {
        // Number of domain in 2D
        double tmp(0.);
        tmp  = (double) number_of_patches[0] / (double) number_of_patches[1];

        number_of_region[0] = min( sz, max(1, (int)sqrt ( (double)sz*tmp) ) );
        number_of_region[1] = (int)(sz / number_of_region[0]);

        while ( number_of_region[0]*number_of_region[1] != (unsigned int) sz ) {
            if (number_of_region[0]>=number_of_region[1] ) {
                number_of_region[0]++;
                number_of_region[1] = (int)(sz / number_of_region[0]);
            }
            else {
                number_of_region[1]++;
                number_of_region[0] = (int)(sz / number_of_region[1]);
            }
        }
    }
    else { // AM and spectral
        number_of_region[0] = sz;
        number_of_region[1] = 1;
        if (number_of_patches[0]<(unsigned int)sz) {
            ERROR( "In AM, the number of patches in dimension 0, here " << number_of_patches[0]
                   << ",  must be at least equal to the number of MPI process which is here " << sz );
        }
    }
    //cout << "ndomain : " << number_of_region[0] << " " << number_of_region[1] << " " << number_of_region[2] << endl;

}


void Params::multiple_decompose_3D()
{
    int rk(0);
    int sz(1);

    MPI_Comm_rank( MPI_COMM_WORLD, &rk );
    MPI_Comm_size( MPI_COMM_WORLD, &sz );

    // Number of domain in 3D
    // Decomposition in 2 times, X and larger side
    double tmp = (double)(number_of_patches[0]*number_of_patches[0]) / (double)(number_of_patches[1]*number_of_patches[2]);
    number_of_region[0] = min( sz, max(1, (int) (cbrt (sz*tmp)) ) );

    int rest = (int)(sz / number_of_region[0]);
    while ( (int)number_of_region[0]*rest != sz ) {
        if ((int)number_of_region[0]>=rest ) {
            number_of_region[0]++;
            rest = (int)(sz / number_of_region[0]);
        }
        else {
            rest++;
            number_of_region[0] = (int)(sz / rest);
        }
    }
    // then the 2 last sides
    double tmp2 = number_of_patches[1] / number_of_patches[2];
    number_of_region[1] = min( rest, max(1, (int)sqrt ( (double)rest*tmp2 ) ) );
    number_of_region[2] = (int)( (double)rest / (double)number_of_region[1] );
    while ( number_of_region[1]*number_of_region[2] != (unsigned int)rest ) {
        if (number_of_region[1]>=number_of_region[2] ) {
            number_of_region[1]++;
            number_of_region[2] = (int)(rest / number_of_region[1]);
        }
        else {
            number_of_region[2]++;
            number_of_region[1] = (int)(rest / number_of_region[2]);
        }
    }
    if ( (number_of_region[0]*number_of_region[1]*number_of_region[2] != (unsigned int)sz ) && (!rk) ) {
        ERROR( "The total number of regions ("<< number_of_region[0]*number_of_region[1]*number_of_region[2]
                << ") is not equal to the number of MPI processes ("
                << number_of_region[0] << "*" << number_of_region[1] << "*" << number_of_region[2]
                << " != " << sz << ")" );
    }

}

string Params::speciesField( string field_name )
{
    if( geometry != "AMcylindrical" ) {
        size_t i1 = field_name.find( "_" );
        size_t l = field_name.length();
        if( i1 != string::npos && l-i1 > 2 ) {
            string field = field_name.substr( 0, i1 );
            if( field == "Jx" || field == "Jy" || field == "Jz" || field == "Rho" ) {
                return field_name.substr( i1+1, l-i1-1 );
            }
        }
    } else {
        size_t i1 = field_name.find( "_" );
        size_t i2 = field_name.rfind( "_mode_" );
        if( i1 != string::npos && i2 != string::npos && i2-i1 > 2 ) {
            return field_name.substr( i1+1, i2-i1-1 );
        }
    }
    return "";
}

#if defined( SMILEI_ACCELERATOR_GPU_OMP ) || defined( SMILEI_ACCELERATOR_GPU_OACC )

bool Params::isGPUParticleBinningAvailable() const
{
    return getGPUClusterWidth() != -1 &&
           getGPUClusterGhostCellBorderWidth() != -1;
}

#else

bool Params::isGPUParticleBinningAvailable() const
{
    return false;
}

#endif

#if defined( SMILEI_ACCELERATOR_GPU_OMP ) || defined( SMILEI_ACCELERATOR_GPU_OACC )

int Params::getGPUClusterWidth() const
{
    return getGPUClusterWidth( nDim_particle );
}

int Params::getGPUClusterGhostCellBorderWidth() const
{
    return getGPUClusterGhostCellBorderWidth( interpolation_order );
}

int Params::getGPUClusterCellVolume() const
{
    // Compute pow(getGPUClusterWidth(), nDim_particle)
    static const int kClusterCellVolume = getGPUClusterWidth() *
                                          ( nDim_particle >= 2 ? getGPUClusterWidth() : 1 ) *
                                          ( nDim_particle >= 3 ? getGPUClusterWidth() : 1 );
    return isGPUParticleBinningAvailable() ?
               kClusterCellVolume :
               -1; // Propagate the error if the dimension is not supported
}

int Params::getGPUInterpolationClusterCellVolume() const
{
    return isGPUParticleBinningAvailable() ?
               getGPUInterpolationClusterCellVolume( nDim_particle, interpolation_order ) :
               -1; // Propagate the error if the dimension is not supported
}

int Params::getGPUBinCount( int dimension_id ) const
{
    const int cells_in_dimension = patch_size_[dimension_id - 1];

    const int kGPUBinCount = cells_in_dimension / getGPUClusterWidth();

    return kGPUBinCount;
}

int Params::getGPUBinCount() const
{
    const int cells_in_cluster_volume = getGPUClusterCellVolume();

    if( cells_in_cluster_volume < 0 ) {
        // Unsupported dimension
        return -1;
    }

    const int kGPUBinCount = n_cell_per_patch / cells_in_cluster_volume;

    SMILEI_ASSERT( ( kGPUBinCount * cells_in_cluster_volume ) == n_cell_per_patch );

    return kGPUBinCount;
}

#endif
