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

#include "H5.h"

#include <algorithm>

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Params : open & parse the input data file, test that parameters are coherent
// ---------------------------------------------------------------------------------------------------------------------
Params::Params(SmileiMPI* smpi, std::vector<std::string> namelistsFiles) :
namelist("")
{
    
    if((((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>16)) || \
        ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
        (H5_VERS_MAJOR>1))) {
        WARNING("Smilei suggests using hdf5 version 1.8.16 you're using "<<H5_VERS_MAJOR << "." << H5_VERS_MINOR << "." << H5_VERS_RELEASE);
        WARNING("Newer version are not tested and may cause the code to behave incorrectly");
        WARNING("See http://hdf-forum.184993.n3.nabble.com/Segmentation-fault-using-H5Dset-extent-in-parallel-td4029082.html");
    }
    

    if (namelistsFiles.size()==0) ERROR("No namelists given!");

    //string commandLineStr("");
    //for (unsigned int i=0;i<namelistsFiles.size();i++) commandLineStr+="\""+namelistsFiles[i]+"\" ";
    //MESSAGE(1,commandLineStr);
    
    //init Python
    PyTools::openPython();
    
    // Print python version
    MESSAGE(1, "Python version "<<PyTools::python_version());
    
    // First, we tell python to filter the ctrl-C kill command (or it would prevent to kill the code execution).
    // This is done separately from other scripts because we don't want it in the concatenated python namelist.
    PyTools::checkPyError();
    string command = "import signal\nsignal.signal(signal.SIGINT, signal.SIG_DFL)";
    if( !PyRun_SimpleString(command.c_str()) ) PyTools::checkPyError();
    
    // Running pyinit.py
    runScript(string(reinterpret_cast<const char*>(pyinit_py), pyinit_py_len), "pyinit.py");

    // Running pyprofiles.py
    runScript(string(reinterpret_cast<const char*>(pyprofiles_py), pyprofiles_py_len), "pyprofiles.py");
    
    // here we add the rank, in case some script need it
    PyModule_AddIntConstant(PyImport_AddModule("__main__"), "smilei_mpi_rank", smpi->getRank());
    
    // here we add the MPI size, in case some script need it
    PyModule_AddIntConstant(PyImport_AddModule("__main__"), "smilei_mpi_size", smpi->getSize());
    
    // here we add the larget int, important to get a valid seed for randomization
    PyModule_AddIntConstant(PyImport_AddModule("__main__"), "smilei_rand_max", RAND_MAX);
    
    // Running the namelists
    runScript("############### BEGIN USER NAMELISTS/COMMANDS ###############\n");
    for (vector<string>::iterator it=namelistsFiles.begin(); it!=namelistsFiles.end(); it++) {
        string strNamelist="";
        if (smpi->isMaster()) {
            ifstream istr(it->c_str());
            if (istr.is_open()) {
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
    
    smpi->barrier();
    
    // Error if no block Main() exists
    if( PyTools::nComponents("Main") == 0 )
        ERROR("Block Main() not defined");
    
    // CHECK namelist on python side
    PyTools::runPyFunction("_smilei_check");
    smpi->barrier();
    
    // output dir: we force this to be the same on all mpi nodes
    string output_dir("");
    PyTools::extract("output_dir", output_dir, "Main");
    PyTools::checkPyError();
    if (!output_dir.empty()) {
        if (chdir(output_dir.c_str()) != 0) {
            WARNING("Could not chdir to output_dir = " << output_dir);
        }
    }
    
    // Now the string "namelist" contains all the python files concatenated
    // It is written as a file: smilei.py
    if (smpi->isMaster()) {
        ofstream out_namelist("smilei.py");
        if (out_namelist.is_open()) {
            out_namelist << "# coding: utf-8" << endl << endl ;
            out_namelist << namelist;
            out_namelist.close();
        }
    }
    
    // random seed
    unsigned int random_seed=0;
    if (!PyTools::extract("random_seed", random_seed, "Main")) {
        random_seed = time(NULL);
    }
    srand(random_seed);
    
    // --------------
    // Stop & Restart
    // --------------
    
    restart = false;
    restart_dir = "";
    if( PyTools::nComponents("DumpRestart")>0 && PyTools::extract("restart_dir", restart_dir, "DumpRestart") ) {
        restart = true;
        if( restart_dir.at(restart_dir.length()-1)!='/' ) restart_dir+="/";
        MESSAGE("Code running from restart in directory "<<restart_dir);
    }
    
    // ---------------------
    // Normalisation & units
    // ---------------------
    
    referenceAngularFrequency_SI = 0.;
    PyTools::extract("referenceAngularFrequency_SI",referenceAngularFrequency_SI, "Main");
    
    
    // -------------------
    // Simulation box info
    // -------------------
    
    // geometry of the simulation
    geometry = "";
    if( !PyTools::extract("geometry", geometry, "Main") )
        ERROR("Parameter Main.geometry is required");
    if (geometry!="1d3v" && geometry!="2d3v")
        ERROR("Main.geometry `" << geometry << "` invalid");
    setDimensions();
    
    // interpolation order
    PyTools::extract("interpolation_order", interpolation_order, "Main");
    if (interpolation_order!=2 && interpolation_order!=4) {
        ERROR("Main.interpolation_order " << interpolation_order << " not defined");
    }
    if (geometry=="2d3v" && interpolation_order==4) {
        ERROR("Main.interpolation_order = 4 " << interpolation_order << " not yet available in 2D");
    }
    
    //!\todo (MG to JD) Please check if this parameter should still appear here
    // Disabled, not compatible for now with particles sort
    // if ( !PyTools::extract("exchange_particles_each", exchange_particles_each) )
    exchange_particles_each = 1;
    
    
    // TIME & SPACE RESOLUTION/TIME-STEPS
    
    // reads timestep & cell_length
    PyTools::extract("timestep", timestep, "Main");
    res_time = 1.0/timestep;
    
    PyTools::extract("cell_length",cell_length, "Main");
    if (cell_length.size()!=nDim_field) {
        ERROR("Dimension of cell_length ("<< cell_length.size() << ") != " << nDim_field << " for geometry " << geometry);
    }
    res_space.resize(nDim_field);
    for (unsigned int i=0;i<nDim_field;i++){
        res_space[i] = 1.0/cell_length[i];
    }
    
    time_fields_frozen=0.0;
    PyTools::extract("time_fields_frozen", time_fields_frozen, "Main");
    
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
    PyTools::extract("sim_time", sim_time, "Main");
    
    PyTools::extract("sim_length",sim_length, "Main");
    if (sim_length.size()!=nDim_field) {
        ERROR("Dimension of sim_length ("<< sim_length.size() << ") != " << nDim_field << " for geometry " << geometry);
    }
    
    
    //! Boundary conditions for ElectroMagnetic Fields
    if ( !PyTools::extract("bc_em_type_x", bc_em_type_x, "Main")  ) {
        ERROR("Electromagnetic boundary condition type (bc_em_type_x) not defined" );
    }
    if (bc_em_type_x.size()==1) { // if just one type is specified, then take the same bc type in a given dimension
        bc_em_type_x.resize(2); bc_em_type_x[1]=bc_em_type_x[0];
    }
    if ( geometry == "2d3v" || geometry == "3d3v" ) {
        if ( !PyTools::extract("bc_em_type_y", bc_em_type_y, "Main") )
            ERROR("Electromagnetic boundary condition type (bc_em_type_y) not defined" );
        if (bc_em_type_y.size()==1) { // if just one type is specified, then take the same bc type in a given dimension
            bc_em_type_y.resize(2); bc_em_type_y[1]=bc_em_type_y[0];
        }
    }
    if ( geometry == "3d3v" ) {
        if ( !PyTools::extract("bc_em_type_z", bc_em_type_z, "Main") )
            ERROR("Electromagnetic boundary condition type (bc_em_type_z) not defined" );
        if (bc_em_type_z.size()==1) { // if just one type is specified, then take the same bc type in a given dimension
            bc_em_type_z.resize(2); bc_em_type_z[1]=bc_em_type_z[0];
        }
    }
    
    // Maxwell Solver 
    PyTools::extract("maxwell_sol", maxwell_sol, "Main");
    
    
    if (!PyTools::extract("clrw",clrw, "Main")) {
        clrw = 1;
    }
        
    // --------------------
    // Number of patches
    // --------------------
    if ( !PyTools::extract("number_of_patches", number_of_patches, "Main") ) {
        ERROR("The parameter `number_of_patches` must be defined as a list of integers");
    }
    for ( unsigned int iDim=0 ; iDim<nDim_field ; iDim++ )
        if( (number_of_patches[iDim] & (number_of_patches[iDim]-1)) != 0)
            ERROR("Number of patches in each direction must be a power of 2");
    
    tot_number_of_patches = 1;
    for ( unsigned int iDim=0 ; iDim<nDim_field ; iDim++ )
        tot_number_of_patches *= number_of_patches[iDim];
    
    if ( tot_number_of_patches == smpi->getSize() ){
        one_patch_per_MPI = true;
    } else {
        one_patch_per_MPI = false;
        if (tot_number_of_patches < smpi->getSize())
            ERROR("The total number of patches must be greater or equal to the number of MPI processes"); 
    }
#ifdef _OPENMP
    if ( tot_number_of_patches < smpi->getSize()*omp_get_max_threads() )
        WARNING( "Resources allocated underloaded regarding the total number of patches" );
#endif
    
    
    balancing_every = 150;
    coef_cell = 1.;
    coef_frozen = 0.1;
    if( PyTools::nComponents("LoadBalancing")>0 ) {
        PyTools::extract("every"      , balancing_every, "LoadBalancing");
        PyTools::extract("coef_cell"  , coef_cell      , "LoadBalancing");
        PyTools::extract("coef_frozen", coef_frozen    , "LoadBalancing");
    }
    
    //mi.resize(nDim_field, 0);
    mi.resize(3, 0);
    while ((number_of_patches[0] >> mi[0]) >1) mi[0]++ ;
    if (number_of_patches.size()>1){
        while ((number_of_patches[1] >> mi[1]) >1) mi[1]++ ;
        if (number_of_patches.size()>2)
            while ((number_of_patches[2] >> mi[2]) >1) mi[2]++ ;
    }
    
    // -------------------------------------------------------
    // Compute usefull quantities and introduce normalizations
    // also defines defaults values for the species lengths
    // -------------------------------------------------------
    compute();
    
    // Print 
    smpi->barrier();
    if ( smpi->isMaster() ) print();
    smpi->barrier();
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
        
    } else {
        ERROR("Problem with the definition of nDim_field");
    }
    
    //!\todo (MG to JD) Are these 2 lines really necessary ? It seems to me it has just been done before
    n_space.resize(3, 1);
    cell_length.resize(3, 0.);            //! \todo{3 but not real size !!! Pbs in Species::Species}
    
    n_space_global.resize(3, 1);        //! \todo{3 but not real size !!! Pbs in Species::Species}
    oversize.resize(3, 0);

    //n_space_global.resize(nDim_field, 0);
    for (unsigned int i=0; i<nDim_field; i++){
        oversize[i]  = interpolation_order + (exchange_particles_each-1);;
        n_space_global[i] = n_space[i];
        n_space[i] /= number_of_patches[i];
        if(n_space_global[i]%number_of_patches[i] !=0) ERROR("ERROR in dimension " << i <<". Number of patches = " << number_of_patches[i] << " must divide n_space_global = " << n_space_global[i]);
        if ( n_space[i] <= 2*oversize[i] ) ERROR ( "ERROR in dimension " << i <<". Patches length = "<<n_space[i] << " cells must be at least " << 2*oversize[i] +1 << " cells long. Increase number of cells or reduce number of patches in this direction. " );
    }
    
    // compute number of cells per patch
    n_cell_per_patch = n_space[0] * n_space[1] * n_space[2];
    
    // Verify that clrw divides n_space[0]
    if( n_space[0]%clrw != 0 )
        ERROR("The parameter clrw must divide the number of cells in one patch (in dimension x)");

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
        MESSAGE(1,"            - (n_space_global,  cell_length) : " << "(" << n_space_global[i] << ", " << cell_length[i] << ")");
    }

    TITLE("Load Balancing: ");
    MESSAGE(1,"Load balancing every " << balancing_every << " iterations.");
    MESSAGE(1,"Cell load coefficient = " << coef_cell );
    MESSAGE(1,"Frozen particle load coefficient = " << coef_frozen );
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
    if (name.size()>0)  MESSAGE(1,"Parsing " << name);
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


