#include "PicParams.h"

#include <cmath>

#include "Tools.h"

using namespace std;

PicParams::PicParams(InputData &ifile) : restart(false), exit_after_dump(true), dump_minutes(0.0), dump_step(0) {
    //open and parse the input data file
    
	ifile.extract("dump_step", dump_step);
	ifile.extract("dump_minutes", dump_minutes);
    
	ifile.extract("exit_after_dump", exit_after_dump);
	
	ifile.extract("restart", restart);
	
	ifile.extract("res_time", res_time);
    ifile.extract("sim_time", sim_time);
    
    ifile.extract("dim", geometry);
    setDimensions();
    
    ifile.extract("interpolation_order", interpolation_order);
    if (interpolation_order!=2 && interpolation_order!=4) {
        ERROR("unacceptable order!");
    }
    if (geometry=="2d3v" && interpolation_order==4) {
        ERROR("unacceptable order for 2D simulation! (not yet implemented)");
    }
    
    ifile.extract("res_space",res_space);
    for (size_t i=0; i<res_space.size(); i++) {
        if (res_space[i] >= res_time) {
            WARNING("res_space[" << i << "] > res_time. Possible CFL problem");
        }
    }
    
    ifile.extract("sim_length",sim_length);
    if (sim_length.size()!=nDim_field) {
        ERROR("Dimension of sim_length ("<< sim_length.size() << ") != " << nDim_field << " for geometry " << geometry);
    }
    if (res_space.size()!=nDim_field) {
        ERROR("Dimension of res_space ("<< res_space.size() << ") != " << nDim_field << " for geometry " << geometry);
    }
    
    ifile.extract("wavelength_SI",wavelength_SI);
    
    ifile.extract("sim_units",sim_units);
    if (sim_units == "physicist") {
        //! \todo{change units to code units}
    }
    
    ifile.extract("plasma_geometry", plasma_geometry);
    if ( (plasma_geometry=="constant") || (plasma_geometry=="crossx") || (plasma_geometry=="crossy") ) {
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
            if (plasma_length.size()!=nDim_field || vacuum_length.size()!=nDim_field || slope_length.size()!=nDim_field) {
                ERROR("plasma_length, vacuum_length and slope_length dimension should be " << nDim_field);
            }
            
        }
        
        //not symmetric density profile
        else{
            ifile.extract("left_slope_length",left_slope_length);
            ifile.extract("right_slope_length",right_slope_length);
            if (plasma_length.size()!=nDim_field || vacuum_length.size()!=nDim_field || left_slope_length.size()!=nDim_field|| right_slope_length.size()!=nDim_field) {
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
    }else {
        ERROR("unknown plasma_geometry "<< plasma_geometry);
    }
    
    
    n_species=0;
    
    while (ifile.existGroup("species",n_species)) {
        SpeciesStructure tmpSpec;
        
        ifile.extract("species_type",tmpSpec.species_type,"species",0,n_species);
        ifile.extract("initialization_type",tmpSpec.initialization_type ,"species",0,n_species);
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
            WARNING("isotropic temperature T="<< tmpSpec.temperature[0]);
        }
        ifile.extract("dynamics_type",tmpSpec.dynamics_type ,"species",0,n_species);
        ifile.extract("time_frozen",tmpSpec.time_frozen ,"species",0,n_species);
        if (tmpSpec.time_frozen > 0 && \
            tmpSpec.initialization_type=="maxwell-juettner") {
            WARNING("For species "<< n_species << " possible conflict in maxwell-juettner initialization");
        }
        ifile.extract("radiating",tmpSpec.radiating ,"species",0,n_species);
        if (tmpSpec.dynamics_type=="rrll" && (!tmpSpec.radiating)) {
            WARNING("dynamics_type rrll forcing radiating true");
            tmpSpec.radiating=true;
        }
        ifile.extract("bc_part_type",tmpSpec.bc_part_type ,"species",0,n_species);
        
        ifile.extract("ionization_model", tmpSpec.ionization_model, "species",0,n_species);
        ifile.extract("atomic_number", tmpSpec.atomic_number, "species",0,n_species);
        
        species_param.push_back(tmpSpec);
        n_species++;
    }
    
    n_laser=0;
    while (ifile.existGroup("laser",n_laser)) {
        LaserStructure tmpLaser;
        ifile.extract("a0",tmpLaser.a0 ,"laser",0,n_laser);
        ifile.extract("angle",tmpLaser.angle ,"laser",0,n_laser);
        ifile.extract("delta",tmpLaser.delta ,"laser",0,n_laser);
        ifile.extract("time_profile",tmpLaser.time_profile ,"laser",0,n_laser);
        ifile.extract("int_params",tmpLaser.int_params ,"laser",0,n_laser);
        ifile.extract("double_params",tmpLaser.double_params ,"laser",0,n_laser);
        
        for (unsigned int i=0; i<tmpLaser.double_params.size(); i++) tmpLaser.double_params[i] *= 2.0*M_PI;
        /* DEFINITION OF THE PARAMETERS MOVED TO LASER.CPP (MG)
         if (tmpLaser.time_profile=="constant") {
         if (tmpLaser.double_params.size()<1) {
         WARNING("Laser always on");
         tmpLaser.double_params.resize(1);
         tmpLaser.double_params[0]=sim_time;
         }
         if (tmpLaser.double_params.size()>1) {
         WARNING("Too much parameters for laser "<< n_laser <<" time_profile ");
         }
         tmpLaser.double_params.resize(1);
         tmpLaser.double_params[0]*= 2.0*M_PI;
         } else {
         ERROR("Laser time_profile " << tmpLaser.time_profile << " not defined");
         }// endif laser
         */
        
        laser_param.push_back(tmpLaser);
        n_laser++;
    }

    if ( !ifile.extract("use_sort_particles", use_sort_particles) )
        use_sort_particles = true;
    if ( !ifile.extract("exchange_particles_each", exchange_particles_each) )
        exchange_particles_each = 1;
    
    compute();
    
}

/*******************************************************************************************************************
 calculate useful parameters
 ******************************************************************************************************************/
void PicParams::compute()
{
    n_time     = res_time*sim_time;
    
    sim_time  *= 2.0*M_PI;
    timestep   = 2.0*M_PI/res_time;
    
    
    for (unsigned int i=0; i<n_species; i++) {
        species_param[i].time_frozen *= 2.0*M_PI;
    }
    
    
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
        }
        for (unsigned int i=nDim_field; i<3; i++) {
            n_space[i]=1;
            cell_length[i]=0.0;
        }
    } else {
        ERROR("This should never happen!!!");
    }
    
    n_space_global.resize(3, 1);	//! \todo{3 but not real size !!! Pbs in Species::Species}
    n_space.resize(3, 1);
    cell_length.resize(3, 0.);	//! \todo{3 but not real size !!! Pbs in Species::Species}
    cell_volume = 1;
    
    oversize.resize(3, 0);
    
    
}

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

void PicParams::print()
{
    //! \todo{Display parameters at runtime}
    MESSAGE(1,"Geometry : " << geometry << "\t\t-> (nDim_particle, nDim_field) : (" << nDim_particle << ", "  << nDim_field << ")");
    MESSAGE(1,"(res_time, sim_time) : (" << res_time << ", " << sim_time << ")"
            << "\t\t-> (n_time, timestep) : (" << n_time << ", " << timestep << ")");
    
    //! \ sim_length[i]*=2.0*M_PI;
    //! \ cell_length[i]=2.0*M_PI/res_space[i];
    for ( unsigned int i=0 ; i<sim_length.size() ; i++ )
        MESSAGE(1,"dim " << i << " - (res_space, sim_length) : (" << res_space[i] << ", " << sim_length[i] << ")"
                << "\t\t-> (n_space, cell_length) : " << "(" << n_space[i] << ", " << cell_length[i] << ")");
    MESSAGE(2,"cell_volume : " << cell_volume);
    
    //! \ vacuum_length[i]*=2.0*M_PI;
    //! \ plasma_length[i]*=2.0*M_PI;
    MESSAGE(1,"plasma_geometry : " << plasma_geometry);
    for ( unsigned int i=0 ; i<plasma_length.size() ; i++ )
        MESSAGE(1,"(plasma_length, vacuum_length) : (" << plasma_length[i] << ", " << vacuum_length[i] << ")");
    MESSAGE(1,"n_species : " << n_species);
    
    MESSAGE(1,"wavelength, sim_units, n_particles : parameters not used for now");
    for ( unsigned int i=0 ; i<n_species ; i++ ) {
        MESSAGE(1,"(species_type, initialization_type, n_part_per_cell, c_part_max) : ("
                << species_param[i].species_type << ", " << species_param[i].initialization_type << ", " << species_param[i].n_part_per_cell << ", " << species_param[i].c_part_max << ") - "
                << "(mass, charge, density) : (" << species_param[i].mass << ", " <<  species_param[i].charge << ", " << species_param[i].density << ")");
        for ( unsigned int j=0 ; j<species_param[i].mean_velocity.size() ; j++ )
            MESSAGE(2,"dim " << j << " - (mean_velocity, temperature) : (" << species_param[i].mean_velocity[j] << ", " << species_param[i].temperature[j] << ")");
        MESSAGE(2," (dynamics_type, bc_part_type, time_frozen, radiating) : (" << species_param[i].dynamics_type <<  ", " << species_param[i].bc_part_type <<  ", " << species_param[i].time_frozen <<  ", " << species_param[i].radiating << ")");
    }
    
    MESSAGE(1,"n_laser : " << n_laser);
    for ( unsigned int i=0 ; i<n_laser ; i++ ) {
        MESSAGE(2,"(a0, angle, delta, time_profile) : (" << laser_param[i].a0 <<  ", " << laser_param[i].angle <<  ", " << laser_param[i].delta <<  ", " << laser_param[i].time_profile << ")");
        
        if (!laser_param[i].int_params.empty()) {
            MESSAGE(2,"int_params : ");
            for ( unsigned int j=0 ; j<laser_param[i].int_params.size() ; j++ )
                MESSAGE(3, laser_param[i].int_params[j]);
        }
        if (!laser_param[i].double_params.empty()) {
            //! \ laser_param[i].double_params[0]*= 2.0*M_PI;
            MESSAGE(2,"double_params : ");
            for ( unsigned int j=0 ; j<laser_param[i].double_params.size() ; j++ )
                MESSAGE(3,laser_param[i].double_params[j]);
        }
    }
    
    MESSAGE(1,"Interpolation_order : " <<  interpolation_order);
    
}

