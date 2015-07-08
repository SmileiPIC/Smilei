#include "PicParams.h"
#include <cmath>
#include "Tools.h"
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
    
    wavelength_SI = 0.;
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
    for (unsigned int ispec = 0; ispec < ifile.nComponents("Species"); ispec++) {
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
        
        if( !ifile.extract("n_part_per_cell",tmpSpec.n_part_per_cell,"Species",ispec) ) {
            ERROR("For species #" << ispec << ", n_part_per_cell was not defined.");
        }
        
        tmpSpec.c_part_max = 1.0;// default value
        ifile.extract("c_part_max",tmpSpec.c_part_max,"Species",ispec);
        
        if( !ifile.extract("mass",tmpSpec.mass ,"Species",ispec) ) {
            ERROR("For species #" << ispec << ", mass not defined.");
        }
        
        if( !ifile.extract("charge",tmpSpec.charge ,"Species",ispec) ) {
            ERROR("For species #" << ispec << ", charge not defined.");
        }
        
        double density;
        ok = false;
        if( ifile.extract("charge_density", density ,"Species",ispec) ) {
            if( tmpSpec.charge==0. ) {
                ERROR("For species #" << ispec << ", cannot define charge_density with charge=0. Define nb_density instead.");
           }
            ok = true;
            tmpSpec.density = abs(density/((double)abs(tmpSpec.charge)));
        }
        if( ifile.extract("nb_density", density ,"Species",ispec) ) {
            if( ok ) { 
                ERROR("For species #" << ispec << ", nb_density and charge_density should not both be defined.");
            }
            ok = true;
            tmpSpec.density = density;
        }
        if( ifile.extract("density", density ,"Species",ispec) ) {
            if( ok ) { 
                ERROR("For species #" << ispec << ", choose only one of density, nb_density or charge_density.");
            }
            WARNING("For species #" << ispec << ", keyword `density` not recommended. Consider `nb_density` or `charge_density` instead.");
            ok = true;
            // Convert to number density
            if ( abs(tmpSpec.charge)==0. ) {
                tmpSpec.density = density;
                WARNING("For species #" << ispec << ", `density` assumed to be the number density.");
            } else {
                tmpSpec.density = abs(density/((double)abs(tmpSpec.charge)));
                WARNING("For species #" << ispec << ", `density` assumed to be the charge density.");
            }
        }
        if( !ok ) {
            ERROR("For species #" << ispec << ", need 'charge_density' or 'nb_density' defined.");
        }
        
        ok = ifile.extract("mean_velocity",tmpSpec.mean_velocity ,"Species",ispec);
        if (tmpSpec.mean_velocity.size()!=3) {
            if( ok ) {
                ERROR("For species #" << ispec << ", mean_velocity must be have 3 components.");
            }
            WARNING("For species #" << ispec << ", mean_velocity assumed = 0.");
            tmpSpec.mean_velocity.resize(3);
            tmpSpec.mean_velocity[0]=tmpSpec.mean_velocity[1]=tmpSpec.mean_velocity[2]=0.0;
        }
        
        ok = ifile.extract("temperature",tmpSpec.temperature ,"Species",ispec);
        if( ok ) {
            if( tmpSpec.temperature.size()==1 ) {
                tmpSpec.temperature.resize(3);
                tmpSpec.temperature[1]=tmpSpec.temperature[2]=tmpSpec.temperature[0];
                WARNING("For species #" << ispec << ", assumed isotropic temperature T = "<< tmpSpec.temperature[0]);
            } else if ( tmpSpec.temperature.size()!=3 ) {
                ERROR("For species #" << ispec << ", temperature must be have 3 components.");
            }
        } else {
            tmpSpec.temperature.resize(3);
            tmpSpec.temperature[0]=tmpSpec.temperature[1]=tmpSpec.temperature[2]=0.0;
            WARNING("For species #" << ispec << ", temperature not defined: assumed = 0.");
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
        
        tmpSpec.isTest = false; // default value
        ifile.extract("isTest",tmpSpec.isTest ,"Species",ispec);
        if (tmpSpec.ionization_model!="none" && (!tmpSpec.isTest)) {
            ERROR("Disabled for now : test & ionized");
        }
                
        // Species geometry
        // ----------------
        
        vector<double> vacuum_length;
        ok = ifile.extract("vacuum_length", vacuum_length,"Species",ispec);
        
        extractProfile(ifile, "dens"  , tmpSpec.dens_profile  , ispec, geometry, vacuum_length);
        extractProfile(ifile, "mvel_x", tmpSpec.mvel_x_profile, ispec, geometry, vacuum_length);
        extractProfile(ifile, "mvel_y", tmpSpec.mvel_y_profile, ispec, geometry, vacuum_length);
        extractProfile(ifile, "mvel_z", tmpSpec.mvel_z_profile, ispec, geometry, vacuum_length);
        extractProfile(ifile, "temp_x", tmpSpec.temp_x_profile, ispec, geometry, vacuum_length);
        extractProfile(ifile, "temp_y", tmpSpec.temp_y_profile, ispec, geometry, vacuum_length);
        extractProfile(ifile, "temp_z", tmpSpec.temp_z_profile, ispec, geometry, vacuum_length);
        
        species_param.push_back(tmpSpec);
    }
}

void PicParams::extractProfile(InputData &ifile, string prefix, ProfileStructure &P, int ispec, string geometry, vector<double> vacuum_length)
{
    
    ifile.extract(prefix+"_profile", P.profile, "Species", ispec);
    if (P.profile.empty()) {
        PyObject *mypy = ifile.extract_py(prefix+"_profile", "Species", ispec);
        if (mypy && PyCallable_Check(mypy)) {
            P.py_profile=mypy;
            P.profile="python";
        } else {
            WARNING("For species " << ispec << ", "+prefix+"_profile not defined, assumed constant.");
            P.profile = "constant";
        }
    }
    if (P.profile != "python") {
        P.vacuum_length = vacuum_length;
        ifile.extract(prefix+"_length_x", P.length_params_x, "Species", ispec);
        if ( (geometry=="2d3v") || (geometry=="3d3v") )
            ifile.extract(prefix+"_length_y", P.length_params_y, "Species", ispec);
        if ( geometry=="3d3v" )
            ifile.extract(prefix+"_length_z", P.length_params_z, "Species", ispec);
        ifile.extract(prefix+"_dbl_params", P.double_params, "Species", ispec);
        ifile.extract(prefix+"_int_params", P.int_params, "Species", ispec);
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
    timestep = conv_fac/res_time;
    sim_time = (double)(n_time) * timestep;
    
    // time after which the moving-window is turned on
    t_move_win *= conv_fac;
    
    
    // grid/cell-related parameters
    // ----------------------------
    n_space.resize(3);
    cell_length.resize(3);
    cell_volume=1.0;
    if (nDim_field==res_space.size() && nDim_field==sim_length.size()) {
        
        // compute number of cells & normalized lengths
        for (unsigned int i=0; i<nDim_field; i++) {
            cell_length[i] = conv_fac/res_space[i];
            sim_length[i] *= conv_fac;
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
        
        // --------------------------------------
        // Normalizing Species-related quantities
        // --------------------------------------
        
        SpeciesStructure * s = &(species_param[ispec]);
        
        // time during which particles are frozen
        s->time_frozen *= conv_fac;
        
        vector<ProfileStructure*> profiles;
        vector<string> prefixes;
        profiles.push_back(&(s->dens_profile  )); prefixes.push_back("dens_"  );
        profiles.push_back(&(s->mvel_x_profile)); prefixes.push_back("mvel_x_");
        profiles.push_back(&(s->mvel_y_profile)); prefixes.push_back("mvel_y_");
        profiles.push_back(&(s->mvel_z_profile)); prefixes.push_back("mvel_z_");
        profiles.push_back(&(s->temp_x_profile)); prefixes.push_back("temp_x_");
        profiles.push_back(&(s->temp_y_profile)); prefixes.push_back("temp_y_");
        profiles.push_back(&(s->temp_z_profile)); prefixes.push_back("temp_z_");
        
        for (unsigned int iprof=0; iprof<profiles.size(); iprof++) {
            
            if (profiles[iprof]->profile=="python") continue;
            
            // normalizing the vacuum lengths
            for (unsigned int i=0; i<profiles[iprof]->vacuum_length.size(); i++)
                profiles[iprof]->vacuum_length[i] *= conv_fac;
            
            // normalizing the density-related lengths
            for (unsigned int i=0; i<profiles[iprof]->length_params_x.size(); i++)
                profiles[iprof]->length_params_x[i] *= conv_fac;
            
            if ( (geometry=="2d3v") || (geometry=="3d3v") ) {
                for (unsigned int i=0; i<profiles[iprof]->length_params_y.size(); i++)
                    profiles[iprof]->length_params_y[i] *= conv_fac;
            }
            
            if ( geometry=="3d3v" ) {
                for (unsigned int i=0; i<profiles[iprof]->length_params_z.size(); i++)
                    profiles[iprof]->length_params_z[i] *= conv_fac;
            }
            
            
            // -----------------------------------------------------
            // Defining default values for species-lengths
            // (NB: here sim_length is already correctly normalized)
            // -----------------------------------------------------
            
            // defining default values for vacuum_length
            if (profiles[iprof]->vacuum_length.size()==0) {
                profiles[iprof]->vacuum_length.resize(1);
                profiles[iprof]->vacuum_length[0] = 0.0;
                //WARNING("No vacuum length defined in x-direction, automatically put to 0 for species " << ispec);
            }
            if ( (geometry=="2d3v") || (geometry=="3d3v") ) {
                if (profiles[iprof]->vacuum_length.size()<2) {
                    profiles[iprof]->vacuum_length.resize(2);
                    profiles[iprof]->vacuum_length[1] = 0.0;
                    //WARNING("No vacuum length defined in y-direction, automatically put to 0 for species " << ispec);
                }
            }
            if (geometry=="3d3v") {
                if (profiles[iprof]->vacuum_length.size()<3) {
                    profiles[iprof]->vacuum_length.resize(3);
                    profiles[iprof]->vacuum_length[2] = 0.0;
                    //WARNING("No vacuum length defined in z-direction, automatically put to 0 for species " << ispec);
                }
            }
            
            // defining default values for dens_length_{x,y,z}
            if (profiles[iprof]->length_params_x.size()==0) {
                profiles[iprof]->length_params_x.resize(1);
                profiles[iprof]->length_params_x[0] = sim_length[0] - profiles[iprof]->vacuum_length[0];
                //WARNING("No " << prefixes[iprof]<< "length_x defined, automatically put to " << profiles[iprof]->length_params_x[0]
                //        << " for species " << ispec);
            }
            if ( (geometry=="2d3v") || (geometry=="3d3v") ) {
                if (profiles[iprof]->length_params_y.size()==0) {
                    profiles[iprof]->length_params_y.resize(1);
                    profiles[iprof]->length_params_y[0] = sim_length[1] - profiles[iprof]->vacuum_length[1];
                    //WARNING("No " << prefixes[iprof]<< "length_y defined, automatically put to " << profiles[iprof]->length_params_y[0]
                    //        << " for species " << ispec);
                }
            }
            if ( geometry=="3d3v" ) {
                if (profiles[iprof]->length_params_z.size()==0) {
                    profiles[iprof]->length_params_z.resize(1);
                    profiles[iprof]->length_params_z[0] = sim_length[2] - profiles[iprof]->vacuum_length[2];
                    //WARNING("No " << prefixes[iprof]<< "length_z defined, automatically put to " << profiles[iprof]->length_params_z[0]
                    //        << " for species " << ispec);
                }
            }
            
            
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
    
    for ( unsigned int i=0 ; i<sim_length.size() ; i++ ){
        MESSAGE(1,"dimension " << i << " - (res_space, sim_length) : (" << res_space[i] << ", " << sim_length[i] << ")");
        MESSAGE(1,"            - (n_space,  cell_length) : " << "(" << n_space[i] << ", " << cell_length[i] << ")");
    }
    
    // Plasma related parameters
    // -------------------------
    MESSAGE("Plasma related parameters");
    MESSAGE(1,"number of species       : " << species_param.size());
    for ( unsigned int i=0 ; i<species_param.size() ; i++ ) {
        MESSAGE(1,"dens_profile.profile : " << species_param[i].dens_profile.profile);
        MESSAGE(1,"            (species_type, number of particles/cell) : ("<< species_param[i].species_type
                << ", " << species_param[i].n_part_per_cell << ")");
    }
    
    
}

