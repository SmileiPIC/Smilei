#include "DiagParams.h"

#include <cmath>
#include <iostream>

#include "Tools.h"

using namespace std;


inline double convertToDouble(string &s)
{
    std::istringstream i(s);
    double x;
    if (!(i >> x)) ERROR("Cannot interpret " << s << " as a number");
    return x;
}



DiagParams::DiagParams(PicParams& params, InputData &ifile) {
    
    double conv_fac = params.conv_fac; // conversion factor (see sim_units in PicParams.cpp for more details)
	
    bool ok=false;
    
    // defining default values & reading diagnostic every-parameter
    // ------------------------------------------------------------
	print_every=params.n_time/10;
    ifile.extract("print_every", print_every);
    
	fieldDump_every=0;
    ok=ifile.extract("fieldDump_every", fieldDump_every);
    if (!ok) fieldDump_every=params.global_every;
        
    avgfieldDump_every=params.res_time*10;
    ok=ifile.extract("avgfieldDump_every", avgfieldDump_every);
    if (!ok) avgfieldDump_every=params.global_every;

    //!\todo Define default behaviour : 0 or params.res_time
    //ntime_step_avg=params.res_time;
    ntime_step_avg=0;
	ifile.extract("ntime_step_avg", ntime_step_avg);
    
	particleDump_every=0;
	if (ifile.extract("particleDump_every", particleDump_every))
            WARNING("Option particleDump_every disabled");
	
	scalar_every=0;
	ok=ifile.extract("every",scalar_every,"diagnostic scalar");
    if (!ok) scalar_every=params.global_every;

    vector<double> scalar_time_range(2,0.);
    ok=ifile.extract("time_range",scalar_time_range,"diagnostic scalar");        
    if (!ok) { 
      scalar_tmin = 0.;
      scalar_tmax = params.sim_time;
    }
    else {
	scalar_tmin = scalar_time_range[0]*conv_fac;
	scalar_tmax = scalar_time_range[1]*conv_fac;
    }
	
    scalar_precision=10;
    ifile.extract("precision",scalar_precision,"diagnostic scalar");
    ifile.extract("vars",scalar_vars,"diagnostic scalar");
	
    unsigned int n_probe=0;
    while (ifile.existGroup("diagnostic probe",n_probe)) {
        probeStructure tmpStruct;
        
        tmpStruct.every=0;
        ok=ifile.extract("every",tmpStruct.every,"diagnostic probe",0,n_probe);        
        if (!ok) tmpStruct.every=params.global_every;

	vector<double> time_range(2,0.);
        ok=ifile.extract("time_range",time_range,"diagnostic probe",0,n_probe);        
        if (!ok) { 
	    tmpStruct.tmin = 0.;
	    tmpStruct.tmax = params.sim_time;
	}
	else  {
	    tmpStruct.tmin = time_range[0]*conv_fac;
	    tmpStruct.tmax = time_range[1]*conv_fac;
	}


        ifile.extract("number",tmpStruct.number,"diagnostic probe",0,n_probe);
        tmpStruct.dim=tmpStruct.number.size();
        if (tmpStruct.dim == 0) { // in 1D case you have one probe, forcing it
            tmpStruct.number.resize(1);
            tmpStruct.number[0]=1;
        }

        vector<double> pos;
        ifile.extract("pos",pos,"diagnostic probe",0,n_probe);
        for (unsigned int i=0; i<pos.size(); i++)
            pos[i] *= conv_fac;
        if (pos.size()>0) tmpStruct.pos.push_back(pos);
        
        ifile.extract("pos_first",pos,"diagnostic probe",0,n_probe);
        for (unsigned int i=0; i<pos.size(); i++)
            pos[i] *= conv_fac;
        if (pos.size()>0) tmpStruct.pos.push_back(pos);
        
        ifile.extract("pos_second",pos,"diagnostic probe",0,n_probe);
        for (unsigned int i=0; i<pos.size(); i++)
            pos[i] *= conv_fac;
        if (pos.size()>0) tmpStruct.pos.push_back(pos);
        
        ifile.extract("pos_third",pos,"diagnostic probe",0,n_probe);
        for (unsigned int i=0; i<pos.size(); i++)
            pos[i] *= conv_fac;
        if (pos.size()>0) tmpStruct.pos.push_back(pos);
        

        probeStruc.push_back(tmpStruct);

        n_probe++;
    }
    
	int n_probephase=0;
	while (ifile.existGroup("diagnostic phase",n_probephase)) {
		phaseStructure tmpPhaseStruct;
        vector<string> kind;
		ifile.extract("kind",kind,"diagnostic phase",0,n_probephase);        
        for (vector<string>::iterator it=kind.begin(); it!=kind.end();it++) {
            if (std::find(kind.begin(), it, *it) == it) {
                tmpPhaseStruct.kind.push_back(*it); 
            } else {
                WARNING("removed duplicate " << *it << " in \"diagnostic phase\" " << n_probephase);
            }
        }

        tmpPhaseStruct.every=0;
		ok=ifile.extract("every",tmpPhaseStruct.every,"diagnostic phase",0,n_probephase);
        if (!ok) {
            if (n_probephase>0) {
                tmpPhaseStruct.every=vecPhase.end()->every;
            } else {
                tmpPhaseStruct.every=params.global_every;
            }
        }

	vector<double> time_range(2,0.);
        ok=ifile.extract("time_range",time_range,"diagnostic phase",0,n_probe);        
        if (!ok) { 
	    tmpPhaseStruct.tmin = 0.;
	    tmpPhaseStruct.tmax = params.sim_time;
	}
	else {
	    tmpPhaseStruct.tmin = time_range[0]*conv_fac;
	    tmpPhaseStruct.tmax = time_range[1]*conv_fac;
	}


        ifile.extract("species",tmpPhaseStruct.species,"diagnostic phase",0,n_probephase);

        tmpPhaseStruct.deflate=0;
        ifile.extract("deflate",tmpPhaseStruct.deflate,"diagnostic phase",0,n_probephase);

		if (tmpPhaseStruct.species.size()==0) {
            WARNING("adding all species to the \"diagnostic phase\" " << n_probephase);
			for (unsigned int i=0;i<params.n_species; i++) {
				tmpPhaseStruct.species.push_back(params.species_param[i].species_type);
			}			
		}
        
		ifile.extract("pos_min",tmpPhaseStruct.pos_min,"diagnostic phase",0,n_probephase);
		ifile.extract("pos_max",tmpPhaseStruct.pos_max,"diagnostic phase",0,n_probephase);
		ifile.extract("pos_num",tmpPhaseStruct.pos_num,"diagnostic phase",0,n_probephase);
        for (unsigned int i=0; i<tmpPhaseStruct.pos_min.size(); i++) {
            tmpPhaseStruct.pos_min[i] *= conv_fac;
            tmpPhaseStruct.pos_max[i] *= conv_fac;
            if (tmpPhaseStruct.pos_min[i]==tmpPhaseStruct.pos_max[i]) {
                tmpPhaseStruct.pos_min[i] = 0.0;
                tmpPhaseStruct.pos_max[i] = params.sim_length[i];
            }
        }
        

		ifile.extract("mom_min",tmpPhaseStruct.mom_min,"diagnostic phase",0,n_probephase);
		ifile.extract("mom_max",tmpPhaseStruct.mom_max,"diagnostic phase",0,n_probephase);
		ifile.extract("mom_num",tmpPhaseStruct.mom_num,"diagnostic phase",0,n_probephase);
		
		ifile.extract("lor_min",tmpPhaseStruct.lor_min,"diagnostic phase",0,n_probephase);
		ifile.extract("lor_max",tmpPhaseStruct.lor_max,"diagnostic phase",0,n_probephase);
		ifile.extract("lor_num",tmpPhaseStruct.lor_num,"diagnostic phase",0,n_probephase);
		
		vecPhase.push_back(tmpPhaseStruct);
		n_probephase++;
	}
	
	
    // particles diagnostics start here
    // --------------------------------
    int n_diag_particles=0;
    unsigned int every, time_average, iaxis, axis_nbins;
    double axis_min, axis_max;
    string output;
    bool axis_logscale, axis_edgeinclusive;
    vector<string> species, axis;
    vector<unsigned int> species_numbers;
    DiagnosticParticlesAxis  *tmpAxis;
    vector<DiagnosticParticlesAxis*> tmpAxes;
    DiagnosticParticles * tmpDiagParticles;
    
    while (ifile.existGroup("diagnostic particles",n_diag_particles)) {
        
        // get parameter "output" that determines the quantity to sum in the output array
        output = "";
        ok = ifile.extract("output",output,"diagnostic particles",0,n_diag_particles);
        if (!ok)
            ERROR("Diagnotic Particles #" << n_diag_particles << ": parameter `output` required");
        
        // get parameter "every" which is the period (in timesteps) for getting the outputs
        every = 0;
        ok = ifile.extract("every",every,"diagnostic particles",0,n_diag_particles);
        if (!ok)
            ERROR("Diagnotic Particles #" << n_diag_particles << ": parameter `every` required");
        
        // get parameter "time_average" that determines the number of timestep to average the outputs
        time_average = 1;
        ifile.extract("time_average",time_average,"diagnostic particles",0,n_diag_particles);
        if (time_average > every)
            ERROR("Diagnotic Particles #" << n_diag_particles << ": `time_average` cannot be larger than `every`");
        if (time_average < 1) time_average=1;
        
        // get parameter "species" that determines the species to use (can be a list of species)
        species.resize(0);
        ok = ifile.extract("species",species,"diagnostic particles",0,n_diag_particles);
        if (!ok)
            ERROR("Diagnotic Particles #" << n_diag_particles << ": parameter `species` required");
        // verify that the species exist, remove duplicates and sort by number
        species_numbers = FindSpecies(species, params);
        
        // get parameter "axis" that adds one axis to the diagnostic
        //  It should contain several items:
        //      requested quantity, min value, max value ,number of bins, log (optional), edge_inclusive (optional)
        ok = true;
        iaxis = 0;
        tmpAxes.resize(0);
        while(true) { // loop in case there are several axes
        	// 1 - Find "axis" keyword and create new axis object
            axis.resize(0);
            ok = ifile.extract("axis",axis,"diagnostic particles",iaxis,n_diag_particles);
            if (!ok) break;
            if (axis.size()<4)
                ERROR("Diagnotic Particles #" << n_diag_particles << ": parameter axis needs at least 4 arguments (type, min, max, nbins)");
            tmpAxis = new DiagnosticParticlesAxis();

            // 2 - Extract axis type (e.g. 'x', 'px', etc.)
            tmpAxis->type  = axis[0];
            if (   (tmpAxis->type == "z" && params.nDim_particle <3)
                || (tmpAxis->type == "y" && params.nDim_particle <2) )
                ERROR("Diagnotic Particles #" << n_diag_particles << ": axis " << tmpAxis->type << " cannot exist in " << params.nDim_particle << "D");
            
            // 3 - Extract axis min and max
            tmpAxis->min   = convertToDouble(axis[1]);
            tmpAxis->max   = convertToDouble(axis[2]);
            
            // 4 - Extract number of bins
            tmpAxis->nbins = convertToDouble(axis[3]);
            if (tmpAxis->nbins - floor(tmpAxis->nbins) != 0.)
                ERROR("Diagnotic Particles #" << n_diag_particles << ": number of bins must be integer (not " << axis[3] << ")");

            // 5 - Check for  other keywords such as "logscale" and "edge_inclusive"
            tmpAxis->logscale = false;
            tmpAxis->edge_inclusive = false;
            for(unsigned int i=4; i<axis.size(); i++) {
                if(axis[i]=="logscale" ||  axis[i]=="log_scale" || axis[i]=="log") {
                        tmpAxis->logscale = true;
                        break;
                }
                if(axis[i]=="edges" ||  axis[i]=="edge" ||  axis[i]=="edge_inclusive" ||  axis[i]=="edges_inclusive") {
                        tmpAxis->edge_inclusive = true;
                        break;
                }
                ERROR("Diagnotic Particles #" << n_diag_particles << ": keyword `" << axis[i] << "` not understood");
            }
            // If the axis is spatial, then we need to apply the conv_fac
            if (axis[0]=="x" || axis[0]=="y" || axis[0]=="z") {
                tmpAxis->min *= conv_fac;
                tmpAxis->max *= conv_fac;
            }
            tmpAxes.push_back(tmpAxis);
            iaxis++;
        }
        if (iaxis == 0)
            ERROR("Diagnotic Particles #" << n_diag_particles << ": at least one parameter `axis` required");
        
        // create new diagnostic object
        tmpDiagParticles = new DiagnosticParticles(n_diag_particles, output, every, time_average, species_numbers, tmpAxes);
        // add this object to the list
        DiagnosticParticles::vecDiagnosticParticles.push_back(tmpDiagParticles);
        // next diagnostic
        n_diag_particles++;
    }

}


// Finds requested species in the list of existing species.
// Returns an array of the numbers of the requested species.
// Note that there might be several species that have the same "name" or "type"
//  so that we have to search for all possibilities.
vector<unsigned int> DiagParams::FindSpecies( vector<string> requested_species, PicParams& params)
{
    bool species_found;
    vector<unsigned int> result;
    unsigned int i;
    vector<string> existing_species;
    
    // Make an array of the existing species names
    existing_species.resize(0);
    for (unsigned int ispec=0 ; ispec<params.n_species ; ispec++) {
        existing_species.push_back( params.species_param[ispec].species_type );
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

