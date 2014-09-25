#include "PicParams.h"
#include <cmath>
#include "Tools.h"

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
    if (restart) MESSAGE("Code run from restart"); //! \todo Give info on restart properties
	
    check_stop_file=false;
    ifile.extract("check_stop_file", check_stop_file);
	
    dump_file_sequence=2;
    ifile.extract("dump_file_sequence", dump_file_sequence);
    dump_file_sequence=std::max((unsigned int)1,dump_file_sequence);
	
    
    // -------------------
    // Simulation box info
    // -------------------
    
    ifile.extract("sim_units",sim_units);
    if (sim_units == "physicist") {
        //! \todo Change units to code units
    }
    
    ifile.extract("wavelength_SI",wavelength_SI);
    
    ifile.extract("dim", geometry);
    if (geometry=="1d3v" && geometry=="2d3v") {
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
    
    //if ( !ifile.extract("exchange_particles_each", exchange_particles_each) )
    exchange_particles_each = 1;
    
    if ( !ifile.extract("use_transverse_periodic", use_transverse_periodic) ) {
        use_transverse_periodic = true;
    }
    /*else if (!use_transverse_periodic) {
     for (unsigned int i=0; i<n_species; i++)
     species_param[i].bc_part_type = "stop";
     }*/
    
    ifile.extract("res_time", res_time);
    
    ifile.extract("sim_time", sim_time);
    
    ifile.extract("res_space",res_space);
    if (res_space.size()!=nDim_field) {
        ERROR("Dimension of res_space ("<< res_space.size() << ") != " << nDim_field << " for geometry " << geometry);
    }
    double Dx2 = 0.0;
    for (short int i=0; i<res_space.size(); i++) {
        Dx2 += 1.0/(res_space[i]*res_space[i]);
    }
    if (sqrt(Dx2) < 1.0/res_time) {
        WARNING("Possible CFL problem: time step = " << 1.0/res_time << " > Dx = " << sqrt(Dx2) );
    }
    
    ifile.extract("sim_length",sim_length);
    if (sim_length.size()!=nDim_field) {
        ERROR("Dimension of sim_length ("<< sim_length.size() << ") != " << nDim_field << " for geometry " << geometry);
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
    
    
    // ---------------
    // Plasma geometry
    // ---------------
    //!  \todo Move plasma geometry properties to species-related properties & clean all this
    
    ifile.extract("plasma_geometry", plasma_geometry);
    if ( (plasma_geometry=="constant") || (plasma_geometry=="crossx") || (plasma_geometry=="crossy")) {
        ifile.extract("plasma_length", plasma_length);
        ifile.extract("vacuum_length", vacuum_length);
        if (plasma_length.size()!=nDim_field || vacuum_length.size()!=nDim_field) {
            ERROR("plasma_length and vacuum_length dimension should be " << nDim_field);
        }
        for (unsigned int i=0; i<nDim_field; i++) {
            if (vacuum_length[i]+plasma_length[i] > sim_length[i])
                WARNING("plasma + vacuum  dimension " << i << " > " << sim_length[i]);
        }       
    } else if (plasma_geometry=="trap") {
        ifile.extract("plasma_length", plasma_length);
        ifile.extract("vacuum_length", vacuum_length);
        ifile.extract("slope_length",  slope_length);
        
        for (unsigned int i=0; i<nDim_field; i++) {
            if (vacuum_length[i]+plasma_length[i] > sim_length[i])
                WARNING("plasma + vacuum " << i << " > " << sim_length[i]);
        }
        
        
        //symmetric density profile
        if (slope_length.size()!=0){
            if (plasma_length.size()!=nDim_field || vacuum_length.size()!=nDim_field
                || slope_length.size()!=nDim_field) {
                ERROR("plasma_length, vacuum_length and slope_length dimension should be " << nDim_field);
            }
            
        }
        //not symmetric density profile
        else{
            ifile.extract("left_slope_length",left_slope_length);
            ifile.extract("right_slope_length",right_slope_length);
            if (plasma_length.size()!=nDim_field || vacuum_length.size()!=nDim_field
                || left_slope_length.size()!=nDim_field|| right_slope_length.size()!=nDim_field) {
                ERROR("plasma_length, vacuum_length and slope_length dimension should be " << nDim_field);
            }
        }
        
    }else if(plasma_geometry=="triangular"){
        ifile.extract("plasma_length", plasma_length);
        ifile.extract("vacuum_length", vacuum_length);
        ifile.extract("slope_length", left_slope_length);
        right_slope_length.resize(left_slope_length.size());
        for(unsigned int i=0;i<nDim_field;i++) {
            right_slope_length[i]=plasma_length[i]-left_slope_length[i];
        }
        
    }
    else if(plasma_geometry=="gaussian"){
        ifile.extract("plasma_length", plasma_length);
        ifile.extract("vacuum_length", vacuum_length);
        ifile.extract("cut",cut);
        ifile.extract("plateau",plateau);
        sigma.resize(nDim_field);
        if(cut.size()==0){
            cut.resize(nDim_field);
            for(unsigned int i=0;i<nDim_field;i++) cut[i]=3.0;
        }
        if(plateau.size()==0){
            plateau.resize(nDim_field);
            for(unsigned int i=0;i<nDim_field;i++) plateau[i]=0.0;
        }
        if(plasma_length.size()!=0){
            for(unsigned int i=0;i<nDim_field;i++) sigma[i]=(plasma_length[i]-plateau[i])/(2*cut[i]);
        }
    }

    else if(plasma_geometry=="polygonal"){
            ifile.extract("plasma_length", plasma_length);
            ifile.extract("vacuum_length", vacuum_length);
            ifile.extract("x_density_coor",x_density_coor);
            ifile.extract("density_rel_values_x",density_rel_values_x);
            if(x_density_coor.size()==0) ERROR("polygonal density profile not well defined");
            
        }
    
    else if(plasma_geometry=="cosine"){
        ifile.extract("plasma_length", plasma_length);
        ifile.extract("vacuum_length", vacuum_length);
        ifile.extract("mode", mode);
        ifile.extract("thetax", thetax);
        ifile.extract("ampl", ampl);
    }


     else if (plasma_geometry=="fukuda"){
        WARNING("plasma geometry: fukuda vacuum & plasma length are not used");
        ifile.extract("plasma_length", plasma_length);
        ifile.extract("vacuum_length", vacuum_length);
        vacuum_length[0] = 2.0;
        vacuum_length[1] = 0.0;
        plasma_length[0] = 112.0;
        plasma_length[1] = 40.0;
    } else if (plasma_geometry=="separated") {
        ifile.extract("plasma_length", plasma_length);
        ifile.extract("vacuum_length", vacuum_length);
    } else {
        ERROR("Unknown plasma_geometry " << plasma_geometry);
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
            WARNING("Isotropic temperature T="<< tmpSpec.temperature[0]);
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
        
        /*if ( (res_space_win_x) && (tmpSpec.bc_part_type!="supp") ) {
         WARNING( "Boundary conditions on particles don't match with moving window, modified" );
         tmpSpec.bc_part_type = "supp";
         }*/
        
        if ( !ifile.extract("ionization_model", tmpSpec.ionization_model, "species",0,n_species) )
            tmpSpec.ionization_model = "none";
        
        ifile.extract("atomic_number", tmpSpec.atomic_number, "species",0,n_species);
        
        species_param.push_back(tmpSpec);
        n_species++;
    }
    
    
    // -----------------
    // Lasers properties
    // -----------------
    n_laser=0;
    while (ifile.existGroup("laser",n_laser)) {
        
        LaserStructure tmpLaser;
        ifile.extract("a0",tmpLaser.a0 ,"laser",0,n_laser);
        
        ifile.extract("boxSide",tmpLaser.boxSide,"laser",0,n_laser);
        if ( (tmpLaser.boxSide!="west") && (tmpLaser.boxSide!="east") ) {
            ERROR("At the moment laser can enter only from West/East sides: boxSide \""
                  << tmpLaser.boxSide << "\" not defined");
        }
        
        ifile.extract("angle",tmpLaser.angle ,"laser",0,n_laser);
        ifile.extract("delta",tmpLaser.delta ,"laser",0,n_laser);
        ifile.extract("time_profile",tmpLaser.time_profile ,"laser",0,n_laser);
        ifile.extract("int_params",tmpLaser.int_params ,"laser",0,n_laser);
        ifile.extract("double_params",tmpLaser.double_params ,"laser",0,n_laser);
        ifile.extract("transv_profile",tmpLaser.transv_profile ,"laser",0,n_laser);
        ifile.extract("int_params_transv",tmpLaser.int_params_transv ,"laser",0,n_laser);
        ifile.extract("double_params_transv",tmpLaser.double_params_transv ,"laser",0,n_laser);
        
        for (unsigned int i=0; i<tmpLaser.double_params.size(); i++)
            tmpLaser.double_params[i] *= 2.0*M_PI;
        for (unsigned int i=0; i<tmpLaser.double_params_transv.size(); i++)
            tmpLaser.double_params_transv[i] *= 2.0*M_PI;
        
        laser_param.push_back(tmpLaser);
        n_laser++;
    }
    
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
    // number of time-steps
    n_time     = res_time*sim_time;
    
    // simulation time & time-step value
    sim_time  *= 2.0*M_PI;
    timestep   = 2.0*M_PI/res_time;
    
    // frozen time
    for (unsigned int i=0; i<n_species; i++) {
        species_param[i].time_frozen *= 2.0*M_PI;
    }
    
    // grid/cell properties
    n_space.resize(3);
    cell_length.resize(3);
    cell_volume=1.0;
    if (nDim_field==res_space.size() && nDim_field==sim_length.size()) {
        for (unsigned int i=0; i<nDim_field; i++) {
            n_space[i]=res_space[i]*sim_length[i];
            
            sim_length[i]*=2.0*M_PI;
            cell_length[i]=2.0*M_PI/res_space[i];
            cell_volume *= cell_length[i];
            vacuum_length[i] *= 2.0*M_PI;
            plasma_length[i] *= 2.0*M_PI;
            
            if (plasma_geometry=="trap") {
                if(slope_length.size()!=0) slope_length[i]  *= 2.0*M_PI;
                else{
                    left_slope_length[i]*= 2.0*M_PI;
                    right_slope_length[i]*= 2.0*M_PI;
                }
            }
            else if(plasma_geometry=="triangular"){
                left_slope_length[i]*= 2.0*M_PI;
                right_slope_length[i]*= 2.0*M_PI;
            }
            else if(plasma_geometry=="gaussian"){
                sigma[i]*= 2.0*M_PI;
                plateau[i]*= 2.0*M_PI;
            }
            else if(plasma_geometry=="polygonal"){
                for (i=0; i<x_density_coor.size(); i++) {
                    x_density_coor[i]*= 2.0*M_PI;
                }
            }
        }
        for (unsigned int i=nDim_field; i<3; i++) {
            n_space[i]=1;
            cell_length[i]=0.0;
        }
        
    } else {
        ERROR("This should never happen: problem with the definition of nDim_field");
    }
    
    n_space_global.resize(3, 1);	//! \todo{3 but not real size !!! Pbs in Species::Species}
    n_space.resize(3, 1);
    cell_length.resize(3, 0.);	    //! \todo{3 but not real size !!! Pbs in Species::Species}
    oversize.resize(3, 0);
    
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
        ERROR("unacceptable geometry! [" << geometry << "]");
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
    MESSAGE(1,"plasma_geometry : " << plasma_geometry);
    MESSAGE(1,"n_species       : " << n_species);
    for ( unsigned int i=0 ; i<n_species ; i++ ) {
        MESSAGE(1,"            (species_type, number of particles/cell) : ("<< species_param[i].species_type
                << ", " << species_param[i].n_part_per_cell << ") - ");
    }
    
    // Laser related parameters
    // ------------------------
    MESSAGE("Laser related parameters");
    MESSAGE(1,"n_laser        : " << n_laser);
    for ( unsigned int i=0 ; i<n_laser ; i++ ) {
        MESSAGE(2,"laser " << i << ": (boxSide, a0) : (" << laser_param[i].boxSide <<  ", " << laser_param[i].a0 <<  ")");
    }

}

