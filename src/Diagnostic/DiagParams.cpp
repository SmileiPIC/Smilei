#include "DiagParams.h"

#include <cmath>

#include "../Tools.h"

using namespace std;

DiagParams::DiagParams(PicParams& params, InputData &ifile) {
	
	print_every=params.n_time/10;
    ifile.extract("print_every", print_every);
	
	fieldDump_every=params.n_time/10;
    ifile.extract("fieldDump_every", fieldDump_every);
    
    avgfieldDump_every=params.res_time*10;
    ifile.extract("avgfieldDump_every", avgfieldDump_every);
    
    ntime_step_avg=params.res_time;
	ifile.extract("ntime_step_avg", ntime_step_avg);
    
	particleDump_every=params.n_time/10;
	if (ifile.extract("particleDump_every", particleDump_every))
            WARNING("Option particleDump_every disabled");
	
	scalar_every=0;
	ifile.extract("every",scalar_every,"diagnostic scalar");
	
    scalar_precision=10;
    ifile.extract("precision",scalar_precision,"diagnostic scalar");
    ifile.extract("vars",scalar_vars,"diagnostic scalar");
	
    unsigned int n_probe=0;
    while (ifile.existGroup("diagnostic probe",n_probe)) {
        probeStructure tmpStruct;
        
        ifile.extract("every",tmpStruct.every,"diagnostic probe",0,n_probe);
        ifile.extract("number",tmpStruct.number,"diagnostic probe",0,n_probe);
        tmpStruct.dim=tmpStruct.number.size();
        if (tmpStruct.dim == 0) { // in 1D case you have one probe, forcing it
            tmpStruct.number.resize(1);
            tmpStruct.number[0]=1;
        }
        vector<double> pos;
        ifile.extract("pos",pos,"diagnostic probe",0,n_probe);
        if (pos.size()>0) tmpStruct.pos.push_back(pos);
        ifile.extract("pos_first",pos,"diagnostic probe",0,n_probe);
        if (pos.size()>0) tmpStruct.pos.push_back(pos);
        ifile.extract("pos_second",pos,"diagnostic probe",0,n_probe);
        if (pos.size()>0) tmpStruct.pos.push_back(pos);
        ifile.extract("pos_third",pos,"diagnostic probe",0,n_probe);
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

        tmpPhaseStruct.every=1;
		ifile.extract("every",tmpPhaseStruct.every,"diagnostic phase",0,n_probephase);
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

		ifile.extract("mom_min",tmpPhaseStruct.mom_min,"diagnostic phase",0,n_probephase);
		ifile.extract("mom_max",tmpPhaseStruct.mom_max,"diagnostic phase",0,n_probephase);
		ifile.extract("mom_num",tmpPhaseStruct.mom_num,"diagnostic phase",0,n_probephase);
		
		ifile.extract("lor_min",tmpPhaseStruct.lor_min,"diagnostic phase",0,n_probephase);
		ifile.extract("lor_max",tmpPhaseStruct.lor_max,"diagnostic phase",0,n_probephase);
		ifile.extract("lor_num",tmpPhaseStruct.lor_num,"diagnostic phase",0,n_probephase);
		
		vecPhase.push_back(tmpPhaseStruct);
		n_probephase++;
	}
	
	
	
}

