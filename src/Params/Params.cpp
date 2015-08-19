#include "Params.h"
#include <cmath>
#include "Tools.h"
#include "PyTools.h"
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
    //init Python
    PyTools::openPython();
    initPython(smpi,namelistsFiles);
    
    unsigned int random_seed=0;
    if (!PyTools::extract("random_seed", random_seed)) {
        random_seed = time(NULL);
    }
    srand(random_seed);
    
    
    dump_step=0;
    PyTools::extract("dump_step", dump_step);
    
    dump_minutes=0.0;
    PyTools::extract("dump_minutes", dump_minutes);
    
    exit_after_dump=true;
    PyTools::extract("exit_after_dump", exit_after_dump);
    
    // --------------
    // Stop & Restart
    // --------------   
    
    restart=false;
    PyTools::extract("restart", restart);
    if (restart) MESSAGE("Code running from restart"); //! \todo Give info on restart properties
    
    //!\todo MG is this still used ?? I cannot find it anywhere
    check_stop_file=false;
    PyTools::extract("check_stop_file", check_stop_file);
    
    dump_file_sequence=2;
    PyTools::extract("dump_file_sequence", dump_file_sequence);
    dump_file_sequence=std::max((unsigned int)1,dump_file_sequence);
    
    
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
        ERROR("CFL problem: timestep=" << timestep << " should be smaller than " << dtCFL);
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
    
    
    // ------------------
    // Species properties
    // ------------------
    readSpecies();
    
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
    computeSpecies();
    
}

Params::~Params() {
    PyTools::closePython();
}
void Params::initPython(SmileiMPI *smpi, std::vector<std::string> namelistsFiles){
    // here we add the rank, in case some script need it
    PyModule_AddIntConstant(PyImport_AddModule("__main__"), "smilei_mpi_rank", smpi->getRank());
    
    // First, we tell python to filter the ctrl-C kill command (or it would prevent to kill the code execution).
    // This is done separately from other scripts because we don't want it in the concatenated python namelist.
    PyTools::checkPyError();
    string command = "import signal\nsignal.signal(signal.SIGINT, signal.SIG_DFL)";
    if( !PyRun_SimpleString(command.c_str()) ) PyTools::checkPyError();
    
    // Running pyinit.py
    pyRunScript(string(reinterpret_cast<const char*>(Python_pyinit_py), Python_pyinit_py_len), "pyinit.py");
    
    // Running pyfunctons.py
    pyRunScript(string(reinterpret_cast<const char*>(Python_pyprofiles_py), Python_pyprofiles_py_len), "pyprofiles.py");
    
    // Running the namelists
    pyRunScript("############### BEGIN USER NAMELISTS ###############\n");
    for (vector<string>::iterator it=namelistsFiles.begin(); it!=namelistsFiles.end(); it++) {
        MESSAGE(1,"Reading file " << *it);
        string strNamelist="";
        if (smpi->isMaster()) {
            ifstream istr(it->c_str());
            if (istr.is_open()) {
                string oneLine;
                while (getline(istr, oneLine)) {
                    strNamelist += oneLine + "\n";
                }
            } else {
                ERROR("File " << (*it) << " does not exists");
            }
            strNamelist +="\n";
        }
        smpi->bcast(strNamelist);
        pyRunScript(strNamelist,(*it));
    }
    pyRunScript("################ END USER NAMELISTS ################\n");    
    // Running pycontrol.py
    pyRunScript(string(reinterpret_cast<const char*>(Python_pycontrol_py), Python_pycontrol_py_len),"pycontrol.py");
    
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
}


void Params::readSpecies() {
    bool ok;
    for (unsigned int ispec = 0; ispec < PyTools::nComponents("Species"); ispec++) {
        SpeciesStructure tmpSpec;
        PyTools::extract("species_type",tmpSpec.species_type,"Species",ispec);
        if(tmpSpec.species_type.empty()) {
            ERROR("For species #" << ispec << " empty species_type");
        }
        PyTools::extract("initPosition_type",tmpSpec.initPosition_type ,"Species",ispec);
        if (tmpSpec.initPosition_type.empty()) {
            ERROR("For species #" << ispec << " empty initPosition_type");
        } else if ( (tmpSpec.initPosition_type!="regular")&&(tmpSpec.initPosition_type!="random") ) {
            ERROR("For species #" << ispec << " bad definition of initPosition_type " << tmpSpec.initPosition_type);
        }
        
        PyTools::extract("initMomentum_type",tmpSpec.initMomentum_type ,"Species",ispec);
        if ( (tmpSpec.initMomentum_type=="mj") || (tmpSpec.initMomentum_type=="maxj") ) {
            tmpSpec.initMomentum_type="maxwell-juettner";
        }
        if (   (tmpSpec.initMomentum_type!="cold")
            && (tmpSpec.initMomentum_type!="maxwell-juettner")
            && (tmpSpec.initMomentum_type!="rectangular") ) {
            ERROR("For species #" << ispec << " bad definition of initMomentum_type");
        }
        
        tmpSpec.c_part_max = 1.0;// default value
        PyTools::extract("c_part_max",tmpSpec.c_part_max,"Species",ispec);
        
        if( !PyTools::extract("mass",tmpSpec.mass ,"Species",ispec) ) {
            ERROR("For species #" << ispec << ", mass not defined.");
        }
        
        tmpSpec.dynamics_type = "norm"; // default value
        if (!PyTools::extract("dynamics_type",tmpSpec.dynamics_type ,"Species",ispec) )
            WARNING("For species #" << ispec << ", dynamics_type not defined: assumed = 'norm'.");
        if (tmpSpec.dynamics_type!="norm"){
            ERROR("dynamics_type different than norm not yet implemented");
        }
        
        tmpSpec.time_frozen = 0.0; // default value
        PyTools::extract("time_frozen",tmpSpec.time_frozen ,"Species",ispec);
        if (tmpSpec.time_frozen > 0 && \
            tmpSpec.initMomentum_type!="cold") {
            WARNING("For species #" << ispec << " possible conflict between time-frozen & not cold initialization");
        }
        
        tmpSpec.radiating = false; // default value
        PyTools::extract("radiating",tmpSpec.radiating ,"Species",ispec);
        if (tmpSpec.dynamics_type=="rrll" && (!tmpSpec.radiating)) {
            WARNING("For species #" << ispec << ", dynamics_type='rrll' forcing radiating=True");
            tmpSpec.radiating=true;
        }
        
        if (!PyTools::extract("bc_part_type_west",tmpSpec.bc_part_type_west,"Species",ispec) )
            ERROR("For species #" << ispec << ", bc_part_type_west not defined");
        if (!PyTools::extract("bc_part_type_east",tmpSpec.bc_part_type_east,"Species",ispec) )
            ERROR("For species #" << ispec << ", bc_part_type_east not defined");
        
        if (nDim_particle>1) {
            if (!PyTools::extract("bc_part_type_south",tmpSpec.bc_part_type_south,"Species",ispec) )
                ERROR("For species #" << ispec << ", bc_part_type_south not defined");
            if (!PyTools::extract("bc_part_type_north",tmpSpec.bc_part_type_north,"Species",ispec) )
                ERROR("For species #" << ispec << ", bc_part_type_north not defined");
        }
        
        // for thermalizing BCs on particles check if thermT is correctly defined
        bool thermTisDefined=false;
        if ( (tmpSpec.bc_part_type_west=="thermalize") || (tmpSpec.bc_part_type_east=="thermalize") ){
            thermTisDefined=PyTools::extract("thermT",tmpSpec.thermT,"Species",ispec);
            if (!thermTisDefined) ERROR("thermT needs to be defined for species " <<ispec<< " due to x-BC thermalize");
        }
        if ( (nDim_particle==2) && (!thermTisDefined) &&
             (tmpSpec.bc_part_type_south=="thermalize" || tmpSpec.bc_part_type_north=="thermalize") ) {
            thermTisDefined=PyTools::extract("thermT",tmpSpec.thermT,"Species",ispec);
            if (!thermTisDefined) ERROR("thermT needs to be defined for species " <<ispec<< " due to y-BC thermalize");
        }
        if (thermTisDefined) {
            if (tmpSpec.thermT.size()==1) {
                tmpSpec.thermT.resize(3);
                for (unsigned int i=1; i<3;i++)
                    tmpSpec.thermT[i]=tmpSpec.thermT[0];
            }
        } else {
            tmpSpec.thermT.resize(3);
            for (unsigned int i=0; i<3;i++)
                tmpSpec.thermT[i]=0.0;
        }
        
        
        tmpSpec.ionization_model = "none"; // default value
        PyTools::extract("ionization_model", tmpSpec.ionization_model, "Species",ispec);
        
        ok = PyTools::extract("atomic_number", tmpSpec.atomic_number, "Species",ispec);
        if( !ok && tmpSpec.ionization_model!="none" ) {
            ERROR("For species #" << ispec << ", `atomic_number` not found => required for the ionization model .");
        }
        
        tmpSpec.isTest = false; // default value
        PyTools::extract("isTest",tmpSpec.isTest ,"Species",ispec);
        if (tmpSpec.ionization_model!="none" && (!tmpSpec.isTest)) {
            ERROR("Disabled for now : test & ionized");
        }
        
        // Species geometry
        // ----------------
        
        // Density
        bool ok1, ok2;
        ok1 = extractOneProfile("nb_density"    , tmpSpec.dens_profile, ispec);
        ok2 = extractOneProfile("charge_density", tmpSpec.dens_profile, ispec);
        if(  ok1 &&  ok2 ) ERROR("For species #" << ispec << ", cannot define both `nb_density` and `charge_density`.");
        if( !ok1 && !ok2 ) ERROR("For species #" << ispec << ", must define `nb_density` or `charge_density`.");
        if( ok1 ) tmpSpec.density_type = "nb";
        if( ok2 ) tmpSpec.density_type = "charge";
        
        // Number of particles per cell
        if( !extractOneProfile("n_part_per_cell", tmpSpec.ppc_profile, ispec) )
            ERROR("For species #" << ispec << ", n_part_per_cell not found or not understood");
        
        // Charge
        if( !extractOneProfile("charge", tmpSpec.charge_profile, ispec) )
            ERROR("For species #" << ispec << ", charge not found or not understood");
        
        // Mean velocity
        vector<ProfileStructure*> vecMvel;
        extractVectorOfProfiles("mean_velocity", vecMvel, ispec);
        tmpSpec.mvel_x_profile = *(vecMvel[0]);
        tmpSpec.mvel_y_profile = *(vecMvel[1]);
        tmpSpec.mvel_z_profile = *(vecMvel[2]);
        
        // Temperature
        vector<ProfileStructure*> vecTemp;
        extractVectorOfProfiles("temperature", vecTemp, ispec);
        tmpSpec.temp_x_profile = *(vecTemp[0]);
        tmpSpec.temp_y_profile = *(vecTemp[1]);
        tmpSpec.temp_z_profile = *(vecTemp[2]);
        
        
        // Save the Species params
        // -----------------------
        species_param.push_back(tmpSpec);
    }
}

bool Params::extractProfile(PyObject *mypy, ProfileStructure &P)
{
    double val;
    // If the profile is only a double, then convert to a constant function
    if( PyTools::convert(mypy, val) ) {
        // Extract the function "constant"
        PyObject* constantFunction = PyTools::extract_py("constant");
        // Create the argument which has the value of the profile
        PyObject* arg = PyTuple_New(1);
        PyTuple_SET_ITEM(arg, 0, PyFloat_FromDouble(val));
        // Create the constant anonymous function
        PyObject * tmp = PyObject_Call(constantFunction, arg, NULL);
        P.py_profile = tmp;
        return true;
    } else if (mypy && PyCallable_Check(mypy)) {
        P.py_profile=mypy;
        return true;
    }
    return false;
}

bool Params::extractOneProfile(string varname, ProfileStructure &P, int ispec) {
    PyObject *mypy = PyTools::extract_py(varname, "Species", ispec);
    if( !extractProfile(mypy, P) ) return false;
    return true;
}

void Params::extractVectorOfProfiles(string varname, vector<ProfileStructure*> &Pvec, int ispec)
{
    Pvec.resize(3);
    vector<PyObject*> pvec = PyTools::extract_pyVec(varname, "Species", ispec);
    int len = pvec.size();
    if( len==3 ) {
        for(int i=0; i<len; i++) {
            Pvec[i] = new ProfileStructure();
            if( !extractProfile(pvec[i], *(Pvec[i])) )
                ERROR("For species #" << ispec << ", "<<varname<<"["<<i<<"] not understood");
        }
    } else if ( len==1 ) {
        Pvec[0] = new ProfileStructure();
        if( !extractProfile(pvec[0], *(Pvec[0])) )
            ERROR("For species #" << ispec << ", "<<varname<<" not understood");
        Pvec[1] = Pvec[0];
        Pvec[2] = Pvec[0];
    } else {
        ERROR("For species #" << ispec << ", "<<varname<<" needs 1 or 3 components.");
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
    n_time   = (int)(res_time*sim_time);
    
    //!\todo MG Is this really necessary now ?
    // simulation time & time-step value
    timestep = 1.0/res_time;
    sim_time = (double)(n_time) * timestep;
    
    
    
    // grid/cell-related parameters
    // ----------------------------
    n_space.resize(3);
    cell_length.resize(3);
    cell_volume=1.0;
    if (nDim_field==res_space.size() && nDim_field==sim_length.size()) {
        
        // compute number of cells & normalized lengths
        for (unsigned int i=0; i<nDim_field; i++) {
            //!\todo MG Is this really necessary now ?
            cell_length[i] = 1.0/res_space[i];
            n_space[i]     = round(sim_length[i]/cell_length[i]);
            sim_length[i]  = (double)(n_space[i])*cell_length[i]; // ensure that nspace = sim_length/cell_length
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
// Compute useful values for Species-related quantities
// ---------------------------------------------------------------------------------------------------------------------
void Params::computeSpecies()
{
    // Loop on all species
    for (unsigned int ispec=0; ispec< species_param.size(); ispec++) {
        
        // here I save the dimension of the pb (to use in BoundaryConditionType.h)
        species_param[ispec].nDim_fields = nDim_field;
        
        // define thermal velocity as \sqrt{T/m}
        species_param[ispec].thermalVelocity.resize(3);
        species_param[ispec].thermalMomentum.resize(3);
        for (unsigned int i=0; i<3; i++) {
            species_param[ispec].thermalVelocity[i]=sqrt(2.0*species_param[ispec].thermT[i]/species_param[ispec].mass);
            species_param[ispec].thermalMomentum[i]=species_param[ispec].mass * species_param[ispec].thermalVelocity[i];
        }
        
    }//end loop on all species (ispec)
    
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
    MESSAGE("Numerical parameters");
    MESSAGE(1,"Geometry : " << geometry)
    MESSAGE(1,"(nDim_particle, nDim_field) : (" << nDim_particle << ", "<< nDim_field << ")");
    MESSAGE(1,"Interpolation_order : " <<  interpolation_order);
    MESSAGE(1,"(res_time, sim_time) : (" << res_time << ", " << sim_time << ")");
    MESSAGE(1,"(n_time,   timestep) : (" << n_time << ", " << timestep << ")");
    MESSAGE(1,"           timestep  = " << timestep/dtCFL << " * CFL");
    
    for ( unsigned int i=0 ; i<sim_length.size() ; i++ ){
        MESSAGE(1,"dimension " << i << " - (res_space, sim_length) : (" << res_space[i] << ", " << sim_length[i] << ")");
        MESSAGE(1,"            - (n_space,  cell_length) : " << "(" << n_space[i] << ", " << cell_length[i] << ")");
    }
    
    // Plasma related parameters
    // -------------------------
    MESSAGE("Plasma related parameters");
    MESSAGE(1,"n_species       : " << species_param.size());
    for ( unsigned int i=0 ; i<species_param.size() ; i++ ) {
        MESSAGE(1,"            species_type : "<< species_param[i].species_type);
    }
    
    
}

// Finds requested species in the list of existing species.
// Returns an array of the numbers of the requested species.
// Note that there might be several species that have the same "name" or "type"
//  so that we have to search for all possibilities.
vector<unsigned int> Params::FindSpecies( vector<string> requested_species)
{
    bool species_found;
    vector<unsigned int> result;
    unsigned int i;
    vector<string> existing_species;
    
    // Make an array of the existing species names
    existing_species.resize(0);
    for (unsigned int ispec=0 ; ispec<species_param.size() ; ispec++) {
        existing_species.push_back( species_param[ispec].species_type );
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
void Params::pyRunScript(string command, string name) {
    PyTools::checkPyError();
    namelist+=command;
    if (name.size()>0)  MESSAGE(1,"Passing to python " << name);
    DEBUG(">>>>>>>>>>>>>>> passing this to python:\n" <<command);
    int retval=PyRun_SimpleString(command.c_str());
    DEBUG("<<<<<<<<<<<<<<< from " << name);
    if (retval==-1) {
        ERROR("error parsing "<< name);
        PyTools::checkPyError();
    }
}

//! run the python functions cleanup (user defined) and _keep_python_running (in pycontrol.py)
void Params::cleanup() {
    
    // call cleanup function from the user namelist (it can be used to free some memory 
    // from the python side) while keeping the interpreter running
    MESSAGE(1,"Checking for cleanup() function:");
    PyTools::runPyFunction("cleanup");
    // this will reset error in python in case cleanup doesn't exists
    PyErr_Clear();
    
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
}


