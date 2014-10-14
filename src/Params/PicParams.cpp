#include "PicParams.h"
#include <cmath>
#include "Tools.h"

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
    
    ifile.extract("sim_units",sim_units);
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
    
    
    ifile.extract("wavelength_SI",wavelength_SI);
    
    
    // -------------------
    // Simulation box info
    // -------------------
    
    ifile.extract("dim", geometry);
    if (geometry!="1d3v" && geometry!="2d3v") {
        ERROR("Geometry " << geometry << " does not exist");
    }
    setDimensions();
    
    
    ifile.extract("interpolation_order", interpolation_order);
    if (interpolation_order!=2 && interpolation_order!=4) {
        ERROR("Interpolation/projection order " << interpolation_order << " not defined");
    }
    if (geometry=="2d3v" && interpolation_order==4) {
        ERROR("Interpolation/projection order " << interpolation_order << " not yet defined in 2D");
    }
    
    
    // Disabled, not compatible for now with particles sort
    // if ( !ifile.extract("exchange_particles_each", exchange_particles_each) )
    //!\todo (MG to JD) Please check if this parameter should still appear here
    exchange_particles_each = 1;
    
    
    // definition or res_time & res_space
    bool defbyRes = ifile.extract("res_time", res_time);
    ifile.extract("res_space",res_space);
    if ((res_space.size()!=0)&&(res_space.size()!=nDim_field)) {
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
    
    // check that res_space has the good dimension
    if (res_space.size()!=nDim_field) {
        ERROR("Dimension of res_space: "<< res_space.size() << " != " << nDim_field << " for geometry " << geometry);
    }
    
    
    // testing the CFL condition
    double res_space2 = 0.0;
    for (unsigned int i=0; i<res_space.size(); i++) {
        res_space2 += (res_space[i]*res_space[i]);
    }
    if ( (sqrt(res_space2) > res_time) || (res_time < *min_element(res_space.begin(),res_space.end())) ) {
        WARNING("Possible CFL problem: res_time = "<<res_time<<" < "<<*min_element(res_space.begin(),res_space.end()));
    }
    
    
    // simulation duration & length
    ifile.extract("sim_time", sim_time);
    
    ifile.extract("sim_length",sim_length);
    if (sim_length.size()!=nDim_field) {
        ERROR("Dimension of sim_length ("<< sim_length.size() << ") != " << nDim_field << " for geometry " << geometry);
    }
    
    
    //! Boundary conditions for ElectroMagnetic Fields
    if ( !ifile.extract("bc_em_type_long", bc_em_type_long)  ) {
        ERROR("bc_em_type_long not defined" );
    }
    if ( geometry == "2d3v" ) {
        if ( !ifile.extract("bc_em_type_trans", bc_em_type_trans) )
            ERROR("bc_em_type_trans not defined" );
    }

    
    // ------------------------
    // Moving window parameters
    // ------------------------
    if (!ifile.extract("res_space_win_x",res_space_win_x)) {
        res_space_win_x = 0;
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
    n_species=0;
    
    while (ifile.existGroup("species",n_species)) {
        SpeciesStructure tmpSpec;
        
        ifile.extract("species_type",tmpSpec.species_type,"species",0,n_species);
        
        ifile.extract("initialization_type",tmpSpec.initialization_type ,"species",0,n_species);
        if (tmpSpec.initialization_type=="mj" || tmpSpec.initialization_type=="m-j") {
            tmpSpec.initialization_type="maxwell-juettner";
        }
        
        ifile.extract("n_part_per_cell",tmpSpec.n_part_per_cell,"species",0,n_species);
        ifile.extract("c_part_max",tmpSpec.c_part_max,"species",0,n_species);
        ifile.extract("mass",tmpSpec.mass ,"species",0,n_species);
        ifile.extract("charge",tmpSpec.charge ,"species",0,n_species);
        ifile.extract("density",tmpSpec.density ,"species",0,n_species);
        
        ifile.extract("mean_velocity",tmpSpec.mean_velocity ,"species",0,n_species);
        if (tmpSpec.mean_velocity.size()!=3) {
            WARNING("mean_velocity should be of size 3 : it is put to zero");
            tmpSpec.mean_velocity.resize(3);
            tmpSpec.mean_velocity[0]=tmpSpec.mean_velocity[1]=tmpSpec.mean_velocity[2]=0.0;
        }
        ifile.extract("temperature",tmpSpec.temperature ,"species",0,n_species);
        if (tmpSpec.temperature.size()==1) {
            tmpSpec.temperature.resize(3);
            tmpSpec.temperature[1]=tmpSpec.temperature[2]=tmpSpec.temperature[0];
            WARNING("Isotropic temperature T="<< tmpSpec.temperature[0] << " for species " << n_species);
        }
        
        ifile.extract("dynamics_type",tmpSpec.dynamics_type ,"species",0,n_species);
        if (tmpSpec.dynamics_type!="norm"){
            ERROR("dynamics_type different than norm not yet implemented");
        }
        
        ifile.extract("time_frozen",tmpSpec.time_frozen ,"species",0,n_species);
        if (tmpSpec.time_frozen > 0 && \
            tmpSpec.initialization_type=="maxwell-juettner") {
            WARNING("For species " << n_species << " possible conflict in Maxwell-Juettner initialization");
        }
        ifile.extract("radiating",tmpSpec.radiating ,"species",0,n_species);
        if (tmpSpec.dynamics_type=="rrll" && (!tmpSpec.radiating)) {
            WARNING("dynamics_type rrll forcing radiating true");
            tmpSpec.radiating=true;
        }
        
        if (!ifile.extract("bc_part_type_long",tmpSpec.bc_part_type_long,"species",0,n_species) )
            ERROR("bc_part_type_long not defined for species " << n_species );
        
        if (nDim_particle>1)
            if (!ifile.extract("bc_part_type_trans ",tmpSpec.bc_part_type_trans,"species",0,n_species) )
                ERROR("bc_part_type_trans not defined for species " << n_species );
        
        if ( !ifile.extract("ionization_model", tmpSpec.ionization_model, "species",0,n_species) )
            tmpSpec.ionization_model = "none";
        
        ifile.extract("atomic_number", tmpSpec.atomic_number, "species",0,n_species);
        
        
        // Species geometry
        // ----------------
        ifile.extract("species_geometry", tmpSpec.species_geometry,"species",0,n_species);
        if( (tmpSpec.species_geometry!="constant") && (tmpSpec.species_geometry!="trapezoidal")
           && (tmpSpec.species_geometry!="gaussian") && (tmpSpec.species_geometry!="sine") ) {
            ERROR("Species_geometry: " << tmpSpec.species_geometry << " not defined for species " << n_species);
        }
        
        // getting vacuum_length & defining default values
        bool vacuum_length_isDefined = ifile.extract("vacuum_length", tmpSpec.vacuum_length,"species",0,n_species);
        if (!vacuum_length_isDefined) {
            tmpSpec.vacuum_length.resize(1);
            tmpSpec.vacuum_length[0] = 0.0;
            WARNING("No vacuum length defined in x-direction, automatically put to 0 for species " << n_species);
        }
        if ( (geometry=="2d3v") || (geometry=="3d3v") ) {
            if (tmpSpec.vacuum_length.size()<2) {
                tmpSpec.vacuum_length.resize(2);
                tmpSpec.vacuum_length[1] = 0.0;
                WARNING("No vacuum length defined in y-direction, automatically put to 0 for species " << n_species);
            }
        }
        if (geometry=="3d3v") {
            if (tmpSpec.vacuum_length.size()<3) {
                tmpSpec.vacuum_length.resize(3);
                tmpSpec.vacuum_length[2] = 0.0;
                WARNING("No vacuum length defined in z-direction, automatically put to 0 for species " << n_species);
            }
        }

        // getting dens_length_{x,y,z} & defining default values
        bool dens_length_x_isDefined = ifile.extract("dens_length_x", tmpSpec.dens_length_x,"species",0,n_species);
        if (!dens_length_x_isDefined) {
            tmpSpec.dens_length_x.resize(1);
            tmpSpec.dens_length_x[0] = sim_length[0] - tmpSpec.vacuum_length[0];
            WARNING("No dens_length_x defined, automatically put to " << tmpSpec.dens_length_x[0] << " for species " << n_species);
        }
        
        if ( (geometry=="2d3v") || (geometry=="3d3v") ) {
            bool dens_length_y_isDefined = ifile.extract("dens_length_y", tmpSpec.dens_length_y,"species",0,n_species);
            if (!dens_length_y_isDefined) {
                tmpSpec.dens_length_y.resize(1);
                tmpSpec.dens_length_y[0] = sim_length[1] - tmpSpec.vacuum_length[1];
                WARNING("No dens_length_y defined, automatically put to " << tmpSpec.dens_length_y[0] << " for species " << n_species);
            }
        }
        
        if ( geometry=="3d3v" ) {
            bool dens_length_z_isDefined = ifile.extract("dens_length_z", tmpSpec.dens_length_z,"species",0,n_species);
            if (!dens_length_z_isDefined) {
                tmpSpec.dens_length_z.resize(1);
                tmpSpec.dens_length_z[0] = sim_length[2] - tmpSpec.vacuum_length[2];
                WARNING("No dens_length_z defined, automatically put to " << tmpSpec.dens_length_z[0] << " for species " << n_species);
            }
        }
        
        // getting additional parameters for the density profile (check DensityProfile for definitions)
        ifile.extract("dens_dbl_params", tmpSpec.dens_dbl_params,"species",0,n_species);
        ifile.extract("dens_int_params", tmpSpec.dens_int_params,"species",0,n_species);
        
        
        species_param.push_back(tmpSpec);
 
        n_species++;
    }
    
    global_every=0;
    
    ifile.extract("every",global_every);
        
    // --------------------
    // Number of processors
    // --------------------
    if ( !ifile.extract("number_of_procs", number_of_procs) )
        number_of_procs.resize(nDim_field, 0);
    
    compute();
    
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
    timestep = conv_fac/res_time;
    sim_time = (double)(n_time) * timestep;
    
    // time after which the moving-window is turned on
    t_move_win *= conv_fac;

    // time during which particles are frozen
    for (unsigned int i=0; i<n_species; i++) {
        species_param[i].time_frozen *= conv_fac;
    }
    
    
    // grid/cell-related parameters
    // ----------------------------
    n_space.resize(3);
    cell_length.resize(3);
    cell_volume=1.0;
    if (nDim_field==res_space.size() && nDim_field==sim_length.size()) {
        
        // compute number of cells & normalized lengths
        for (unsigned int i=0; i<nDim_field; i++) {
            cell_length[i] = conv_fac/res_space[i];
            sim_length[i] *= conv_fac;//(double)(n_space[i]) * cell_length[i];
            n_space[i]     = (int)(sim_length[i]/cell_length[i]);//(int)(res_space[i]*sim_length[i]);
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
    cell_length.resize(3, 0.);	    //! \todo{3 but not real size !!! Pbs in Species::Species}
    
    n_space_global.resize(3, 1);	//! \todo{3 but not real size !!! Pbs in Species::Species}
    oversize.resize(3, 0);

    
    // species-related length normalization
    // ------------------------------------
    for (unsigned int ispec=0; ispec< species_param.size(); ispec++) {
        
        // normalizing the vacuum lengths
        for (unsigned int i=0; i<species_param[ispec].vacuum_length.size(); i++)
            species_param[ispec].vacuum_length[i] *= conv_fac;
        
        // normalizing the density-related lengths
        for (unsigned int i=0; i<species_param[ispec].dens_length_x.size(); i++)
            species_param[ispec].dens_length_x[i]   *= conv_fac;
        
        if ( (geometry=="2d3v") || (geometry=="3d3v") ) {
            for (unsigned int i=0; i<species_param[ispec].dens_length_y.size(); i++)
                species_param[ispec].dens_length_y[i]   *= conv_fac;
        }
        
        if ( geometry=="3d3v" ) {
            for (unsigned int i=0; i<species_param[ispec].dens_length_z.size(); i++)
                species_param[ispec].dens_length_z[i]   *= conv_fac;
        }
        
    }
    
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
    
    for ( unsigned int i=0 ; i<sim_length.size() ; i++ ){
        MESSAGE(1,"dimension " << i << " - (res_space, sim_length) : (" << res_space[i] << ", " << sim_length[i] << ")");
        MESSAGE(1,"            - (n_space,  cell_length) : " << "(" << n_space[i] << ", " << cell_length[i] << ")");
    }

    // Plasma related parameters
    // -------------------------
    MESSAGE("Plasma related parameters");
    MESSAGE(1,"n_species       : " << n_species);
    for ( unsigned int i=0 ; i<n_species ; i++ ) {
        MESSAGE(1,"species_geometry : " << species_param[i].species_geometry);
        MESSAGE(1,"            (species_type, number of particles/cell) : ("<< species_param[i].species_type
                << ", " << species_param[i].n_part_per_cell << ")");
    }
    

}

