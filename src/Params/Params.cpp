#include "PyTools.h"
#include "Params.h"
#include "Species.h"
#include <cmath>
#include <iomanip>
#include "Tools.h"
#include "SmileiMPI.h"

#include "pyinit.pyh"
#include "pyprofiles.pyh"
#include "pycontrol.pyh"

#include <algorithm>

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Params : open & parse the input data file, test that parameters are coherent
// ---------------------------------------------------------------------------------------------------------------------
Params::Params(SmileiMPI* smpi, std::vector<std::string> namelistsFiles) :
namelist("")
{
    string commandLineStr("");
    for (int i=0;i<namelistsFiles.size();i++) commandLineStr+=namelistsFiles[i]+" ";
    MESSAGE(1,commandLineStr);
    
    //init Python
    PyTools::openPython();
    
    // First, we tell python to filter the ctrl-C kill command (or it would prevent to kill the code execution).
    // This is done separately from other scripts because we don't want it in the concatenated python namelist.
    PyTools::checkPyError();
    string command = "import signal\nsignal.signal(signal.SIGINT, signal.SIG_DFL)";
    if( !PyRun_SimpleString(command.c_str()) ) PyTools::checkPyError();
    
    // Running pyinit.py
    runScript(string(reinterpret_cast<const char*>(pyinit_py), pyinit_py_len), "pyinit.py");
    // here we add the rank, in case some script need it
    PyModule_AddIntConstant(PyImport_AddModule("__main__"), "smilei_mpi_rank", smpi->getRank());
    
    // Running pyfunctons.py
    runScript(string(reinterpret_cast<const char*>(pyprofiles_py), pyprofiles_py_len), "pyprofiles.py");
    
    // Running the namelists
    runScript("############### BEGIN USER NAMELISTS/COMMANDS ###############\n");
    for (vector<string>::iterator it=namelistsFiles.begin(); it!=namelistsFiles.end(); it++) {
        string strNamelist="";
        if (smpi->isMaster()) {
            ifstream istr(it->c_str());
            if (istr.is_open()) {
                MESSAGE(1,"Reading file " << *it);
                std::stringstream buffer;
                buffer << istr.rdbuf();
                strNamelist+=buffer.str();
            } else {
                strNamelist = "# Smilei:) From command line :\n" + (*it);
            }
            strNamelist +="\n";
        }
        smpi->bcast(strNamelist);
        runScript(strNamelist,(*it));
    }
    runScript("################ END USER NAMELISTS/COMMANDS  ################\n");
    // Running pycontrol.py
    runScript(string(reinterpret_cast<const char*>(pycontrol_py), pycontrol_py_len),"pycontrol.py");
    
    PyTools::runPyFunction("_smilei_check");
    
    
    // Now the string "namelist" contains all the python files concatenated
    // It is written as a file: smilei.py
    if (smpi->isMaster()) {
        ofstream out_namelist("smilei.py");
        if (out_namelist.is_open()) {
            out_namelist << namelist;
            out_namelist.close();
        }
    }
    
    unsigned int random_seed=0;
    if (!PyTools::extract("random_seed", random_seed)) {
        random_seed = time(NULL);
    }
    srand(random_seed);
    
    // --------------
    // Stop & Restart
    // --------------   
    
    restart=false;
    PyTools::extract("restart", restart);
    if (restart) MESSAGE("Code running from restart"); //! \todo Give info on restart properties
        
    
    // ---------------------
    // Normalisation & units
    // ---------------------
    
    wavelength_SI = 0.;
    PyTools::extract("wavelength_SI",wavelength_SI);
    
    
    // -------------------
    // Simulation box info
    // -------------------
    
    // geometry of the simulation
    PyTools::extract("dim", geometry);
    if (geometry!="1d3v" && geometry!="2d3v") {
        ERROR("Geometry " << geometry << " does not exist");
    }
    setDimensions();
    
    // interpolation order
    PyTools::extract("interpolation_order", interpolation_order);
    if (interpolation_order!=2 && interpolation_order!=4) {
        ERROR("Interpolation/projection order " << interpolation_order << " not defined");
    }
    if (geometry=="2d3v" && interpolation_order==4) {
        ERROR("Interpolation/projection order " << interpolation_order << " not yet defined in 2D");
    }
    
    //!\todo (MG to JD) Please check if this parameter should still appear here
    // Disabled, not compatible for now with particles sort
    // if ( !PyTools::extract("exchange_particles_each", exchange_particles_each) )
    exchange_particles_each = 1;
    
    
    // TIME & SPACE RESOLUTION/TIME-STEPS
    
    // reads timestep & cell_length
    PyTools::extract("timestep", timestep);
    res_time = 1.0/timestep;
    PyTools::extract("cell_length",cell_length);
    if (cell_length.size()!=nDim_field) {
        ERROR("Dimension of cell_length ("<< cell_length.size() << ") != " << nDim_field << " for geometry " << geometry);
    }
    res_space.resize(nDim_field);
    for (unsigned int i=0;i<nDim_field;i++){
        res_space[i] = 1.0/cell_length[i];
    }
    
    time_fields_frozen=0.0;
    PyTools::extract("time_fields_frozen", time_fields_frozen);
    
    // testing the CFL condition
    //!\todo (MG) CFL cond. depends on the Maxwell solv. ==> Move this computation to the ElectroMagn Solver
    double res_space2=0;
    for (unsigned int i=0; i<nDim_field; i++) {
        res_space2 += res_space[i]*res_space[i];
    }
    dtCFL=1.0/sqrt(res_space2);
    if ( timestep>dtCFL ) {
        WARNING("CFL problem: timestep=" << timestep << " should be smaller than " << dtCFL);
    }
    
    
    // simulation duration & length
    PyTools::extract("sim_time", sim_time);
    
    PyTools::extract("sim_length",sim_length);
    if (sim_length.size()!=nDim_field) {
        ERROR("Dimension of sim_length ("<< sim_length.size() << ") != " << nDim_field << " for geometry " << geometry);
    }
    
    
    //! Boundary conditions for ElectroMagnetic Fields
    if ( !PyTools::extract("bc_em_type_x", bc_em_type_x)  ) {
        ERROR("Electromagnetic boundary condition type (bc_em_type_x) not defined" );
    }
    if (bc_em_type_x.size()==1) { // if just one type is specified, then take the same bc type in a given dimension
        bc_em_type_x.resize(2); bc_em_type_x[1]=bc_em_type_x[0];
    }
    if ( geometry == "2d3v" || geometry == "3d3v" ) {
        if ( !PyTools::extract("bc_em_type_y", bc_em_type_y) )
            ERROR("Electromagnetic boundary condition type (bc_em_type_y) not defined" );
        if (bc_em_type_y.size()==1) { // if just one type is specified, then take the same bc type in a given dimension
            bc_em_type_y.resize(2); bc_em_type_y[1]=bc_em_type_y[0];
        }
    }
    if ( geometry == "3d3v" ) {
        if ( !PyTools::extract("bc_em_type_z", bc_em_type_z) )
            ERROR("Electromagnetic boundary condition type (bc_em_type_z) not defined" );
        if (bc_em_type_z.size()==1) { // if just one type is specified, then take the same bc type in a given dimension
            bc_em_type_z.resize(2); bc_em_type_z[1]=bc_em_type_z[0];
        }
    }
    
    // Maxwell Solver 
	PyTools::extract("maxwell_sol", maxwell_sol);


    // ------------------------
    // Moving window parameters
    // ------------------------
    if (!PyTools::extract("nspace_win_x",nspace_win_x)) {
        nspace_win_x = 0;
    }
    
    if (!PyTools::extract("t_move_win",t_move_win)) {
        t_move_win = 0.0;
    }
    
    if (!PyTools::extract("vx_win",vx_win)) {
        vx_win = 1.;
    }
    
    if (!PyTools::extract("clrw",clrw)) {
        clrw = 1;
    }
    
    global_every=0;
    
    PyTools::extract("every",global_every);
    
    // --------------------
    // Number of processors
    // --------------------
    if ( !PyTools::extract("number_of_procs", number_of_procs) )
        number_of_procs.resize(nDim_field, 0);
    
    // -------------------------------------------------------
    // Compute usefull quantities and introduce normalizations
    // also defines defaults values for the species lengths
    // -------------------------------------------------------
    compute();
    
}

Params::~Params() {
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
    n_time   = (int)(sim_time/timestep);
    
    // simulation time & time-step value
    double entered_sim_time = sim_time;
    sim_time = (double)(n_time) * timestep;
    if (sim_time!=entered_sim_time)
        WARNING("sim_time has been redefined from " << entered_sim_time << " to " << sim_time << " to match nxtimestep (" << scientific << setprecision(4) << sim_time - entered_sim_time<< ")" );
    
    
    // grid/cell-related parameters
    // ----------------------------
    n_space.resize(3);
    cell_length.resize(3);
    cell_volume=1.0;
    if (nDim_field==res_space.size() && nDim_field==sim_length.size()) {
        
        // compute number of cells & normalized lengths
        for (unsigned int i=0; i<nDim_field; i++) {
            n_space[i]         = round(sim_length[i]/cell_length[i]);
            double entered_sim_length = sim_length[i];
            sim_length[i]      = (double)(n_space[i])*cell_length[i]; // ensure that nspace = sim_length/cell_length
            if (sim_length[i]!=entered_sim_length)
                WARNING("sim_length[" << i << "] has been redefined from " << entered_sim_length << " to " << sim_length[i] << " to match n x cell_length (" << scientific << setprecision(4) << sim_length[i]-entered_sim_length <<")");
            cell_volume   *= cell_length[i];
        }
        // create a 3d equivalent of n_space & cell_length
        for (unsigned int i=nDim_field; i<3; i++) {
            n_space[i]=1;
            cell_length[i]=0.0;
        }
        // compute number of cells per cluster
        n_cell_per_cluster = clrw * n_space[1] * n_space[2];
        
    } else {
        ERROR("Problem with the definition of nDim_field");
    }
    
    //!\todo (MG to JD) Are these 2 lines really necessary ? It seems to me it has just been done before
    n_space.resize(3, 1);
    cell_length.resize(3, 0.);	    //! \todo{3 but not real size !!! Pbs in Species::Species}
    
    n_space_global.resize(3, 1);	//! \todo{3 but not real size !!! Pbs in Species::Species}
    oversize.resize(3, 0);
    
}


// ---------------------------------------------------------------------------------------------------------------------
// Set dimensions according to geometry
// ---------------------------------------------------------------------------------------------------------------------
void Params::setDimensions()
{
    if (geometry=="1d3v") {
        nDim_particle=1;
        nDim_field=1;
    } else if (geometry=="2d3v") {
        nDim_particle=2;
        nDim_field=2;
    } else if (geometry=="3d3v") {
        nDim_particle=3;
        nDim_field=3;
    } else if (geometry=="2drz") {
        nDim_particle=3;
        nDim_field=2;
    } else {
        ERROR("Geometry: " << geometry << " not defined");
    }
}



// ---------------------------------------------------------------------------------------------------------------------
// Printing out the data at initialisation
// ---------------------------------------------------------------------------------------------------------------------
void Params::print()
{
    // Numerical parameters
    // ---------------------
    TITLE("Geometry: " << geometry);
    MESSAGE(1,"(nDim_particle, nDim_field) : (" << nDim_particle << ", "<< nDim_field << ")");
    MESSAGE(1,"Interpolation_order : " <<  interpolation_order);
    MESSAGE(1,"(res_time, sim_time) : (" << res_time << ", " << sim_time << ")");
    MESSAGE(1,"(n_time,   timestep) : (" << n_time << ", " << timestep << ")");
    MESSAGE(1,"           timestep  = " << timestep/dtCFL << " * CFL");
    
    for ( unsigned int i=0 ; i<sim_length.size() ; i++ ){
        MESSAGE(1,"dimension " << i << " - (res_space, sim_length) : (" << res_space[i] << ", " << sim_length[i] << ")");
        MESSAGE(1,"            - (n_space,  cell_length) : " << "(" << n_space[i] << ", " << cell_length[i] << ")");
    }
    
    
}

// ---------------------------------------------------------------------------------------------------------------------
// Finds requested species in the list of existing species.
// Returns an array of the numbers of the requested species.
// Note that there might be several species that have the same "name" or "type"
//  so that we have to search for all possibilities.
vector<unsigned int> Params::FindSpecies(vector<Species*>& vecSpecies, vector<string> requested_species)
{
    bool species_found;
    vector<unsigned int> result;
    unsigned int i;
    vector<string> existing_species;
    
    // Make an array of the existing species names
    existing_species.resize(0);
    for (unsigned int ispec=0 ; ispec<vecSpecies.size() ; ispec++) {
        existing_species.push_back( vecSpecies[ispec]->species_type );
    }
    
    // Loop over group of requested species
    for (unsigned int rs=0 ; rs<requested_species.size() ; rs++) {
        species_found = false;
        // Loop over existing species
        for (unsigned int es=0 ; es<existing_species.size() ; es++) {
            if (requested_species[rs] == existing_species[es]) { // if found
                species_found = true;
                // Add to the list and sort
                for (i=0 ; i<result.size() ; i++) {
                    if (es == result[i]) break; // skip if duplicate
                    if (es <  result[i]) {
                        result.insert(result.begin()+i,es); // insert at the right place
                        break;
                    }
                }
                // Put at the end if not put earlier
                if (i == result.size()) result.push_back(es);
            }
        }
        if (!species_found)
            ERROR("Species `" << requested_species[rs] << "` was not found.");
    }
    
    return result;
}


//! Run string as python script and add to namelist
void Params::runScript(string command, string name) {
    PyTools::checkPyError();
    namelist+=command;
    if (name.size()>0)  MESSAGE(1,"Passing to python " << name);
    int retval=PyRun_SimpleString(command.c_str());
    if (retval==-1) {
        ERROR("error parsing "<< name);
        PyTools::checkPyError();
    }
}

//! run the python functions cleanup (user defined) and _keep_python_running (in pycontrol.py)
void Params::cleanup(SmileiMPI* smpi) {
    // call cleanup function from the user namelist (it can be used to free some memory 
    // from the python side) while keeping the interpreter running
    MESSAGE(1,"Checking for cleanup() function:");
    PyTools::runPyFunction("cleanup");
    // this will reset error in python in case cleanup doesn't exists
    PyErr_Clear();
    
    smpi->barrier();
    
    // this function is defined in the Python/pyontrol.py file and should return false if we can close
    // the python interpreter
    MESSAGE(1,"Calling python _keep_python_running() :");    
    if (PyTools::runPyFunction<bool>("_keep_python_running")) {
        MESSAGE(2,"Keeping Python interpreter alive");
    } else {
        MESSAGE(2,"Closing Python");
        PyErr_Print();
        Py_Finalize();
    }
    smpi->barrier();
}


