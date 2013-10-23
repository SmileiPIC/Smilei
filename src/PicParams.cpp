#include "PicParams.h"
#include "Tools.h"
#include <cmath>

using namespace std;

PicParams::PicParams() {
}

void PicParams::parseInputData(InputData &ifile) {
	//open and parse the input data file
	
	DEBUGEXEC(ifile.extract("debug",debug_level));
	
	RELEASEEXEC(int debug_level=-1; ifile.extract("debug",debug_level));
	RELEASEEXEC(if (debug_level!=-1) WARNING("This is a release compilation, debug keyword is ignored"));

	ifile.extract("res_time", res_time);
	ifile.extract("sim_time", sim_time);
	
	ifile.extract("dim", geometry);
	setDimensions();

	ifile.extract("interpolation_order", interpolation_order);
	if (interpolation_order!=2) {
		ERROR("unacceptable order!");
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
	
	ifile.extract("wavelength",wavelength);	
	ifile.extract("sim_units",sim_units);
	if (sim_units == "physicist") {
		//! \todo{change units to code units}
	}
	
    ifile.extract("plasma_geometry", plasma_geometry);
	if (plasma_geometry=="constant") {
		ifile.extract("plasma_length", plasma_length);
		ifile.extract("vacuum_length", vacuum_length);
		if (plasma_length.size()!=nDim_field || vacuum_length.size()!=nDim_field) {
			ERROR("plasma_length and vacuum_length dimension should be " << nDim_field);
		}
		for (unsigned int i=0; i<nDim_field; i++) {
			if (vacuum_length[i]+plasma_length[i] > sim_length[i])
				WARNING("plasma + vacuum  dimension " << i << " > " << sim_length[i]);
		}
	} else {
		ERROR("unknown plasma_geometry "<< plasma_geometry);
	}
	
    ifile.extract("n_species", n_species);
	species_param.resize(n_species);
	
	for (unsigned int i=0; i<n_species; i++) {
		ostringstream ss;
		ss << "species " << i;
		string group=ss.str();
		
		bool found=false;
		vector<string> groups=ifile.getGroups();
		for (size_t j=0; j<groups.size(); j++) {
			if (groups[j]==group) found=true;
		}
		if (found) {			
			ifile.extract("species_type",species_param[i].species_type,group);
            ifile.extract("initialization_type",species_param[i].initialization_type ,group);
            ifile.extract("n_part_per_cell",species_param[i].n_part_per_cell,group);
            ifile.extract("c_part_max",species_param[i].c_part_max,group);
			ifile.extract("mass",species_param[i].mass ,group);
			ifile.extract("charge",species_param[i].charge ,group);
			ifile.extract("density",species_param[i].density ,group);
			ifile.extract("mean_velocity",species_param[i].mean_velocity ,group);
			if (species_param[i].mean_velocity.size()!=nDim_particle) {
				WARNING("mean_velocity size : should be " << nDim_particle);
			}
			ifile.extract("temperature",species_param[i].temperature ,group);
			if (species_param[i].temperature.size()==1) {
				species_param[i].temperature.resize(3);
				species_param[i].temperature[1]=species_param[i].temperature[2]=species_param[i].temperature[0];
				WARNING("isotropic temperature T="<< species_param[i].temperature[0]);
			}
			ifile.extract("dynamics_type",species_param[i].dynamics_type ,group);
			ifile.extract("time_frozen",species_param[i].time_frozen ,group);
			if (species_param[i].time_frozen > 0 && \
				species_param[i].initialization_type=="maxwell-juettner") {
				WARNING("For species "<<i<< " possible conflict in maxwell-juettner initialization");
			}
			ifile.extract("radiating",species_param[i].radiating ,group);
			if (species_param[i].dynamics_type=="rrll" && (!species_param[i].radiating)) {
				WARNING("dynamics_type rrll forcing radiating true");
				species_param[i].radiating=true;
			}
			ifile.extract("bc_part_type",species_param[i].bc_part_type ,group);
			
			ifile.extract("ionization", species_param[i].ionization_model, group);
			ifile.extract("atomic_number", species_param[i].atomic_number, group);
			
		} else {
			ERROR("species " << i << " not defined");
		}
	}	

	n_laser=0;
	ifile.extract("n_laser", n_laser);
	laser_param.resize(n_laser);
	
	for (unsigned int i=0; i<n_laser; i++) {
		ostringstream ss;
		ss << "laser " << i;
		string group=ss.str();
		
		bool found=false;
		vector<string> groups=ifile.getGroups();
		for (size_t j=0; j<groups.size(); j++) {
			if (groups[j]==group) found=true;
		}
		if (found) {
			ifile.extract("a0",laser_param[i].a0 ,group);
			ifile.extract("angle",laser_param[i].angle ,group);
			ifile.extract("delta",laser_param[i].delta ,group);
			ifile.extract("time_profile",laser_param[i].time_profile ,group);
			ifile.extract("int_params",laser_param[i].int_params ,group);
			ifile.extract("double_params",laser_param[i].double_params ,group);
		}
		
		if (laser_param[i].time_profile=="constant") {
			if (laser_param[i].double_params.size()<1) {
				WARNING("Laser always on");
				laser_param[i].double_params.resize(1);
				laser_param[i].double_params[0]=sim_time;
			}
			if (laser_param[i].double_params.size()>1) {
				WARNING("Too much parameters for laser "<< i <<" time_profile ");
			}
			laser_param[i].double_params.resize(1);
			laser_param[i].double_params[0]*= 2.0*M_PI;
		} else {
			ERROR("Laser time_profile " << laser_param[i].time_profile << " not defined");
		}// endif laser
	}
	
}

	/*******************************************************************************************************************
	 caclulate useful parameters
	 ******************************************************************************************************************/
void PicParams::compute()
{
	n_time=res_time*sim_time;
	//! \todo{clean this Mickael!!}
	sim_time*=2.0*M_PI;

	timestep = 2.0*M_PI/res_time;

	
	for (unsigned int i=0; i<n_species; i++) {
		species_param[i].time_frozen *= 2.0*M_PI;
	}	
	
	
	n_space.resize(3);
	cell_length.resize(3);
	cell_volume=1.0;
	if (nDim_field==res_space.size() && nDim_field==sim_length.size()) {
		for (unsigned int i=0; i<nDim_field;i++) {
			n_space[i]=res_space[i]*sim_length[i];
			//! \todo{clean this Mickael!!}
			sim_length[i]*=2.0*M_PI;
			cell_length[i]=2.0*M_PI/res_space[i];
			cell_volume *= cell_length[i];
			
			vacuum_length[i]*=2.0*M_PI;
			plasma_length[i]*=2.0*M_PI;
		}
		for (unsigned int i=nDim_field; i<3;i++) {
			n_space[i]=1;
			cell_length[i]=0.0;
		}
	} else {
		ERROR("This should never happen!!!");
	}
	
	
	
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
	cout << "\tGeometry : " << geometry
		<< "\t\t-> (nDim_particle, nDim_field) : (" << nDim_particle << ", "  << nDim_field << ")" << endl;
	cout << "\t(res_time, sim_time) : (" << res_time << ", " << sim_time << ")"
		<< "\t\t-> (n_time, timestep) : (" << n_time << ", " << timestep << ")" << endl;
	
	//! \ sim_length[i]*=2.0*M_PI;
	//! \ cell_length[i]=2.0*M_PI/res_space[i];
	for ( unsigned int i=0 ; i<sim_length.size() ; i++ )
		cout << "\tdim " << i << " - (res_space, sim_length) : (" << res_space[i] << ", " << sim_length[i] << ")"
			<< "\t\t-> (n_space, cell_length) : " << "(" << n_space[i] << ", " << cell_length[i] << ")" << endl;
	cout << "\t\tcell_volume : " << cell_volume << endl;

	//! \ vacuum_length[i]*=2.0*M_PI;
	//! \ plasma_length[i]*=2.0*M_PI;
	cout << "\tplasma_geometry : " << plasma_geometry << endl;
	for ( unsigned int i=0 ; i<plasma_length.size() ; i++ )
		cout << "\t(plasma_length, vacuum_length) : (" << plasma_length[i] << ", " << vacuum_length[i] << ")" << endl;
	cout << "\tn_species : " << n_species << endl;

	cout << "\twavelength, sim_units, n_particles : parameters not used for now" << endl;
	//! \ n_part_max : Initialised in Species::Species()
	for ( unsigned int i=0 ; i<n_species ; i++ ) {
		cout << "\t(species_type, initialization_type, n_part_per_cell, c_part_max) : ("
			<< species_param[i].species_type << ", " << species_param[i].initialization_type << ", " << species_param[i].n_part_per_cell << ", " << species_param[i].c_part_max << ") - "  
			<< "(mass, charge, density) : (" << species_param[i].mass << ", " <<  species_param[i].charge << ", " << species_param[i].density << ")" << endl;
			for ( unsigned int j=0 ; j<species_param[i].mean_velocity.size() ; j++ )
				cout << "\t\tdim " << j << " - (mean_velocity, temperature) : (" << species_param[i].mean_velocity[j] << ", " << species_param[i].temperature[j] << ")" << endl;
		cout << "\t\t (dynamics_type, bc_part_type, time_frozen, radiating) : (" << species_param[i].dynamics_type <<  ", " << species_param[i].bc_part_type <<  ", " << species_param[i].time_frozen <<  ", " << species_param[i].radiating << ")" << endl;
	}

	cout << "\tn_laser : " << n_laser << endl;
	for ( unsigned int i=0 ; i<n_laser ; i++ ) {
		cout << "\t\t(a0, angle, delta, time_profile) : (" << laser_param[i].a0 <<  ", " << laser_param[i].angle <<  ", " << laser_param[i].delta <<  ", " << laser_param[i].time_profile << ")" << endl;

		if (!laser_param[i].int_params.empty()) {
			cout << "\t\tint_params : ";
			for ( unsigned int j=0 ; j<laser_param[i].int_params.size()-1 ; j++ )
				cout << laser_param[i].int_params[j] << ", ";
			cout << laser_param[i].int_params[laser_param[i].int_params.size()-1] << endl;
		}
		if (!laser_param[i].double_params.empty()) {
			//! \ laser_param[i].double_params[0]*= 2.0*M_PI;
			cout << "\t\tdouble_params : ";
			for ( unsigned int j=0 ; j<laser_param[i].double_params.size()-1 ; j++ )
				cout << laser_param[i].double_params[j] << ", ";
			cout << laser_param[i].double_params[laser_param[i].double_params.size()-1] << endl;
		}
	}

	cout << "\tInterpolation_order : " <<  interpolation_order << endl;

}

