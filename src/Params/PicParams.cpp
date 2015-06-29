#include "PicParams.h"
#include <cmath>
#include "Tools.h"
#include "PyTools.h"
#include "InputData.h"

#include <algorithm>

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// PicParams : open & parse the input data file, test that parameters are coherent
// ---------------------------------------------------------------------------------------------------------------------
PicParams::PicParams(InputData &ifile) {
    
    
    // --------------
    // Stop & Restart
    // --------------   
    dump_step=0;
    ifile.extract("dump_step", dump_step);
    
    dump_minutes=0.0;
    ifile.extract("dump_minutes", dump_minutes);
    
    exit_after_dump=true;
    ifile.extract("exit_after_dump", exit_after_dump);
    
    restart=false;
    ifile.extract("restart", restart);
    if (restart) MESSAGE("Code running from restart"); //! \todo Give info on restart properties
    
    check_stop_file=false;
    ifile.extract("check_stop_file", check_stop_file);
    
    dump_file_sequence=2;
    ifile.extract("dump_file_sequence", dump_file_sequence);
    dump_file_sequence=std::max((unsigned int)1,dump_file_sequence);
    
    
    // ---------------------
    // Normalisation & units
    // ---------------------
    
/*MG150609    ifile.extract("sim_units",sim_units);
    if (sim_units == "normalized") {
        conv_fac = 1.0;
    }
    else if (sim_units == "wavelength") {
        conv_fac = 2.0*M_PI;
        WARNING("Wavelength-related units are used for code entries but not for outputs (apart from log file)");
    }
    else {
        ERROR("Simulation units sim_units" << sim_units << " not specified or inexisting");
    }
 */
    
    wavelength_SI = 0.;
    ifile.extract("wavelength_SI",wavelength_SI);
    
    
    // -------------------
    // Simulation box info
    // -------------------
    
    // geometry of the simulation
    ifile.extract("dim", geometry);
    if (geometry!="1d3v" && geometry!="2d3v") {
        ERROR("Geometry " << geometry << " does not exist");
    }
    setDimensions();
    
    // interpolation order
    ifile.extract("interpolation_order", interpolation_order);
    if (interpolation_order!=2 && interpolation_order!=4) {
        ERROR("Interpolation/projection order " << interpolation_order << " not defined");
    }
    if (geometry=="2d3v" && interpolation_order==4) {
        ERROR("Interpolation/projection order " << interpolation_order << " not yet defined in 2D");
    }
    
    //!\todo (MG to JD) Please check if this parameter should still appear here
    // Disabled, not compatible for now with particles sort
    // if ( !ifile.extract("exchange_particles_each", exchange_particles_each) )
    exchange_particles_each = 1;
    
    
    // TIME & SPACE RESOLUTION/TIME-STEPS
    
/*MG150609    // definition or res_time & res_space
    bool defbyRes = ifile.extract("res_time", res_time);
    ifile.extract("res_space",res_space);
    if ( (res_space.size()!=0) && (res_space.size()!=nDim_field) ) {
        ERROR("Dimension of res_space ("<< res_space.size() << ") != " << nDim_field << " for geometry " << geometry);
    }
    
    // definition of time_step & cell_length (if res_time & res_space are not defined)
    if (!defbyRes) {
        ifile.extract("timestep", timestep);
        res_time = 1.0/timestep;
        ifile.extract("cell_length",cell_length);
        if (cell_length.size()!=nDim_field) {
            ERROR("Dimension of cell_length ("<< cell_length.size() << ") != " << nDim_field << " for geometry " << geometry);
        }
        res_space.resize(nDim_field);
        for (unsigned int i=0;i<nDim_field;i++){
            res_space[i] = 1.0/cell_length[i];
        }
    }
*/
    // reads timestep & cell_length
    ifile.extract("timestep", timestep);
    res_time = 1.0/timestep;
    ifile.extract("cell_length",cell_length);
    if (cell_length.size()!=nDim_field) {
        ERROR("Dimension of cell_length ("<< cell_length.size() << ") != " << nDim_field << " for geometry " << geometry);
    }
    res_space.resize(nDim_field);
    for (unsigned int i=0;i<nDim_field;i++){
        res_space[i] = 1.0/cell_length[i];
    }

/*MG150609    // check that res_space has the good dimension
    if (res_space.size()!=nDim_field) {
        ERROR("Dimension of res_space: "<< res_space.size() << " != " << nDim_field << " for geometry " << geometry);
    }*/
    
    time_fields_frozen=0.0;
    ifile.extract("time_fields_frozen", time_fields_frozen);
    
    // testing the CFL condition
    //!\todo (MG) CFL cond. depends on the Maxwell solv. ==> Move this computation to the ElectroMagn Solver
    double res_space2=0;
    for (unsigned int i=0; i<nDim_field; i++) {
        res_space2 += res_space[i]*res_space[i];
    }
    dtCFL=1.0/sqrt(res_space2);
    if ( timestep>dtCFL ) {
        ERROR("Possible CFL problem: timestep=" << timestep << " should be smaller than " << dtCFL);
    }
    
    
    // simulation duration & length
    ifile.extract("sim_time", sim_time);
    
    ifile.extract("sim_length",sim_length);
    if (sim_length.size()!=nDim_field) {
        ERROR("Dimension of sim_length ("<< sim_length.size() << ") != " << nDim_field << " for geometry " << geometry);
    }
    
    
    //! Boundary conditions for ElectroMagnetic Fields
    if ( !ifile.extract("bc_em_type_x", bc_em_type_x)  ) {
        ERROR("Electromagnetic boundary condition type (bc_em_type_x) not defined" );
    }
    if (bc_em_type_x.size()==1) { // if just one type is specified, then take the same bc type in a given dimension
        bc_em_type_x.resize(2); bc_em_type_x[1]=bc_em_type_x[0];
    }
    if ( geometry == "2d3v" || geometry == "3d3v" ) {
        if ( !ifile.extract("bc_em_type_y", bc_em_type_y) )
            ERROR("Electromagnetic boundary condition type (bc_em_type_y) not defined" );
        if (bc_em_type_y.size()==1) { // if just one type is specified, then take the same bc type in a given dimension
            bc_em_type_y.resize(2); bc_em_type_y[1]=bc_em_type_y[0];
        }
    }
    if ( geometry == "3d3v" ) {
        if ( !ifile.extract("bc_em_type_z", bc_em_type_z) )
            ERROR("Electromagnetic boundary condition type (bc_em_type_z) not defined" );
        if (bc_em_type_z.size()==1) { // if just one type is specified, then take the same bc type in a given dimension
            bc_em_type_z.resize(2); bc_em_type_z[1]=bc_em_type_z[0];
        }
    }

    
    // ------------------------
    // Moving window parameters
    // ------------------------
    if (!ifile.extract("nspace_win_x",nspace_win_x)) {
        nspace_win_x = 0;
    }
    
    if (!ifile.extract("t_move_win",t_move_win)) {
        t_move_win = 0.0;
    }
    
    if (!ifile.extract("vx_win",vx_win)) {
        vx_win = 1.;
    }
    
    if (!ifile.extract("clrw",clrw)) {
        clrw = 1;
    }
    
    
    // ------------------
    // Species properties
    // ------------------
    readSpecies(ifile);
    
    global_every=0;
    
    ifile.extract("every",global_every);
    
    // --------------------
    // Number of processors
    // --------------------
    if ( !ifile.extract("number_of_procs", number_of_procs) )
        number_of_procs.resize(nDim_field, 0);
    
    // -------------------------------------------------------
    // Compute usefull quantities and introduce normalizations
    // also defines defaults values for the species lengths
    // -------------------------------------------------------
    compute();
    computeSpecies();
    
}

void PicParams::readSpecies(InputData &ifile) {
    bool ok;
    n_species=ifile.nComponents("Species");
    for (unsigned int ispec = 0; ispec < n_species; ispec++) {
        SpeciesStructure tmpSpec;

        ifile.extract("species_type",tmpSpec.species_type,"Species",ispec);
        if(tmpSpec.species_type.empty()) {
            ERROR("For species #" << ispec << " empty species_type");
        }
        ifile.extract("initPosition_type",tmpSpec.initPosition_type ,"Species",ispec);
        if (tmpSpec.initPosition_type.empty()) {
            ERROR("For species #" << ispec << " empty initPosition_type");
        } else if ( (tmpSpec.initPosition_type!="regular")&&(tmpSpec.initPosition_type!="random") ) {
            ERROR("For species #" << ispec << " bad definition of initPosition_type " << tmpSpec.initPosition_type);
        }
        
        ifile.extract("initMomentum_type",tmpSpec.initMomentum_type ,"Species",ispec);
        if ( (tmpSpec.initMomentum_type=="mj") || (tmpSpec.initMomentum_type=="maxj") ) {
            tmpSpec.initMomentum_type="maxwell-juettner";
        }
        if (   (tmpSpec.initMomentum_type!="cold")
            && (tmpSpec.initMomentum_type!="maxwell-juettner")
            && (tmpSpec.initMomentum_type!="rectangular") ) {
            ERROR("For species #" << ispec << " bad definition of initMomentum_type");
        }
        
        tmpSpec.c_part_max = 1.0;// default value
        ifile.extract("c_part_max",tmpSpec.c_part_max,"Species",ispec);
        
        if( !ifile.extract("mass",tmpSpec.mass ,"Species",ispec) ) {
            ERROR("For species #" << ispec << ", mass not defined.");
        }
        
        tmpSpec.dynamics_type = "norm"; // default value
        if (!ifile.extract("dynamics_type",tmpSpec.dynamics_type ,"Species",ispec) )
            WARNING("For species #" << ispec << ", dynamics_type not defined: assumed = 'norm'.");
        if (tmpSpec.dynamics_type!="norm"){
            ERROR("dynamics_type different than norm not yet implemented");
        }
        
        tmpSpec.time_frozen = 0.0; // default value
        ifile.extract("time_frozen",tmpSpec.time_frozen ,"Species",ispec);
        if (tmpSpec.time_frozen > 0 && \
            tmpSpec.initMomentum_type!="cold") {
            WARNING("For species #" << ispec << " possible conflict between time-frozen & not cold initialization");
        }
        
        tmpSpec.radiating = false; // default value
        ifile.extract("radiating",tmpSpec.radiating ,"Species",ispec);
        if (tmpSpec.dynamics_type=="rrll" && (!tmpSpec.radiating)) {
            WARNING("For species #" << ispec << ", dynamics_type='rrll' forcing radiating=True");
            tmpSpec.radiating=true;
        }
        
        if (!ifile.extract("bc_part_type_west",tmpSpec.bc_part_type_west,"Species",ispec) )
            ERROR("For species #" << ispec << ", bc_part_type_west not defined");
        if (!ifile.extract("bc_part_type_east",tmpSpec.bc_part_type_east,"Species",ispec) )
            ERROR("For species #" << ispec << ", bc_part_type_east not defined");
        
        if (nDim_particle>1) {
            if (!ifile.extract("bc_part_type_south",tmpSpec.bc_part_type_south,"Species",ispec) )
                ERROR("For species #" << ispec << ", bc_part_type_south not defined");
            if (!ifile.extract("bc_part_type_north",tmpSpec.bc_part_type_north,"Species",ispec) )
                ERROR("For species #" << ispec << ", bc_part_type_north not defined");
        }
        
        tmpSpec.ionization_model = "none"; // default value
        ifile.extract("ionization_model", tmpSpec.ionization_model, "Species",ispec);
        
        ok = ifile.extract("atomic_number", tmpSpec.atomic_number, "Species",ispec);
        if( !ok && tmpSpec.ionization_model!="none" ) {
            ERROR("For species #" << ispec << ", `atomic_number` not found => required for the ionization model .");
        }
        
        
        // Species geometry
        // ----------------
        
        // Density
        bool ok1, ok2;
        ok1 = extractOneProfile(ifile, "nb_density"    , tmpSpec.dens_profile, ispec);
        ok2 = extractOneProfile(ifile, "charge_density", tmpSpec.dens_profile, ispec);
        if(  ok1 &&  ok2 ) ERROR("For species #" << ispec << ", cannot define both `nb_density` and `charge_density`.");
        if( !ok1 && !ok2 ) ERROR("For species #" << ispec << ", must define `nb_density` or `charge_density`.");
        if( ok1 ) tmpSpec.density_type = "nb";
        if( ok2 ) tmpSpec.density_type = "charge";
        // Number of particles per cell
        if( !extractOneProfile(ifile, "n_part_per_cell", tmpSpec.ppc_profile, ispec) )
            ERROR("For species #" << ispec << ", n_part_per_cell not found or not understood");
        // Charge
        if( !extractOneProfile(ifile, "charge", tmpSpec.charge_profile, ispec) )
            ERROR("For species #" << ispec << ", charge not found or not understood");
        // Mean velocity
        vector<ProfileStructure*> vecMvel;
        extractVectorOfProfiles(ifile, "mean_velocity", vecMvel, ispec);
        tmpSpec.mvel_x_profile = *(vecMvel[0]);
        tmpSpec.mvel_y_profile = *(vecMvel[1]);
        tmpSpec.mvel_z_profile = *(vecMvel[2]);
        // Temperature
        vector<ProfileStructure*> vecTemp;
        extractVectorOfProfiles(ifile, "temperature", vecTemp, ispec);
        tmpSpec.temp_x_profile = *(vecTemp[0]);
        tmpSpec.temp_y_profile = *(vecTemp[1]);
        tmpSpec.temp_z_profile = *(vecTemp[2]);
        
        species_param.push_back(tmpSpec);
    }
}

bool PicParams::extractProfile(InputData &ifile, PyObject *mypy, ProfileStructure &P)
{
    double val;
    // If the profile is only a double, then convert to a constant function
    if( PyTools::convert(mypy, val) ) {
        // Extract the function "constant"
        PyObject* constantFunction = ifile.extract_py("constant");
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

bool PicParams::extractOneProfile(InputData &ifile, string varname, ProfileStructure &P, int ispec) {
    PyObject *mypy = ifile.extract_py(varname, "Species", ispec);
    if( !extractProfile(ifile, mypy, P) ) return false;
    return true;
}

void PicParams::extractVectorOfProfiles(InputData &ifile, string varname, vector<ProfileStructure*> &Pvec, int ispec)
{
    Pvec.resize(3);
    vector<PyObject*> pvec = ifile.extract_pyVec(varname, "Species", ispec);
    int len = pvec.size();
    if( len==3 ) {
        for(int i=0; i<len; i++) {
            Pvec[i] = new ProfileStructure();
            if( !extractProfile(ifile, pvec[i], *(Pvec[i])) )
                ERROR("For species #" << ispec << ", "<<varname<<"["<<i<<"] not understood");
        }
    } else if ( len==1 ) {
        Pvec[0] = new ProfileStructure();
        if( !extractProfile(ifile, pvec[0], *(Pvec[0])) )
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
void PicParams::compute()
{
    // time-related parameters
    // -----------------------
    
    // number of time-steps
    n_time   = (int)(res_time*sim_time);
    
    // simulation time & time-step value
    timestep = 1.0/res_time;//*MG150609 conv_fac/res_time;
    sim_time = (double)(n_time) * timestep;
    
/*MG150609   // time during which Maxwell's eqs. are not solved (cst fields)
    time_fields_frozen *= conv_fac;
    
    // time after which the moving-window is turned on
    t_move_win *= conv_fac;
 */
    
    
    // grid/cell-related parameters
    // ----------------------------
    n_space.resize(3);
    cell_length.resize(3);
    cell_volume=1.0;
    if (nDim_field==res_space.size() && nDim_field==sim_length.size()) {
        
        // compute number of cells & normalized lengths
        for (unsigned int i=0; i<nDim_field; i++) {
            /*MG150609cell_length[i] = conv_fac/res_space[i];
            sim_length[i] *= conv_fac;
            n_space[i]     = round(sim_length[i]/cell_length[i]);
            sim_length[i]  = (double)(n_space[i])*cell_length[i]; // ensure that nspace = sim_length/cell_length
            cell_volume   *= cell_length[i];*/
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
void PicParams::computeSpecies()
{
    
    // Loop on all species
    for (unsigned int ispec=0; ispec< species_param.size(); ispec++) {
        
        // here I save the dimension of the pb (to use in BoundaryConditionType.h)
        if (geometry=="1d3v")
            species_param[ispec].nDim_fields = 1;
        else if (geometry=="2d3v")
            species_param[ispec].nDim_fields = 2;
        else if (geometry=="3d3v")
            species_param[ispec].nDim_fields = 3;
        
        // define thermal velocity as \sqrt{T/m}
        species_param[ispec].thermalVelocity.resize(3);
        species_param[ispec].thermalMomentum.resize(3);
        for (unsigned int i=0; i<3; i++) {
            // \fixme: This line is problematic because it uses temperature, which is a python function
            // species_param[ispec].thermalVelocity[i] = sqrt( 2.0 *species_param[ispec].temperature[i]/species_param[ispec].mass );
            species_param[ispec].thermalVelocity[i] = 0.;
            species_param[ispec].thermalMomentum[i] = species_param[ispec].mass * species_param[ispec].thermalVelocity[i];
        }
        
    }//end loop on all species (ispec)
    
}


// ---------------------------------------------------------------------------------------------------------------------
// Set dimensions according to geometry
// ---------------------------------------------------------------------------------------------------------------------
void PicParams::setDimensions()
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
void PicParams::print()
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
    MESSAGE(1,"n_species       : " << n_species);
    for ( unsigned int i=0 ; i<n_species ; i++ ) {
        MESSAGE(1,"            species_type : "<< species_param[i].species_type);
    }
    
    
}

