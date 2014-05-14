#include "DiagParams.h"

#include <cmath>

#include "Tools.h"

using namespace std;

DiagParams::DiagParams(InputData &ifile, PicParams& params) {
	
	print_every=params.n_time/10;
    ifile.extract("print_every", print_every);
	
	fieldDump_every=params.n_time/10;
    ifile.extract("fieldDump_every", fieldDump_every);
	
	particleDump_every=params.n_time/10;
	ifile.extract("particleDump_every", particleDump_every);
	
	scalar_every=0;
	ifile.extract("every",scalar_every,"diagnostic scalar");
	

    unsigned int n_probe0d=0;
    while (ifile.existGroup("diagnostic probe0d",n_probe0d)) {
        probe0DStructure tmpStruct;
        ifile.extract("every",tmpStruct.every,"diagnostic probe0d",0,n_probe0d);
        tmpStruct.pos.resize(params.nDim_field);
        ifile.extract("pos",tmpStruct.pos,"diagnostic probe0d",0,n_probe0d);
        if (params.nDim_field != tmpStruct.pos.size()) {
            ERROR("diagnostic probe0d: different dimension of nDim and pos");
        }
        probe0DStruc.push_back(tmpStruct);
        n_probe0d++;
    }
    
    unsigned int n_probe1d=0;
    while (ifile.existGroup("diagnostic probe1d",n_probe1d)) {
        probe1DStructure tmpStruct;
        
        ifile.extract("every",tmpStruct.every,"diagnostic probe1d",0,n_probe1d);
        ifile.extract("number",tmpStruct.number,"diagnostic probe1d",0,n_probe1d);
        
        ifile.extract("pos_start",tmpStruct.posStart,"diagnostic probe1d",0,n_probe1d);
        ifile.extract("pos_end",tmpStruct.posEnd,"diagnostic probe1d",0,n_probe1d);
        
        probe1DStruc.push_back(tmpStruct);
        n_probe1d++;
    }
	
    
	int n_probephase=0;
	while (ifile.existGroup("diagnostic phase",n_probephase)) {
		phaseStructure tmpPhaseStruct;
		ifile.extract("kind",tmpPhaseStruct.kind,"diagnostic phase",0,n_probephase);
		ifile.extract("every",tmpPhaseStruct.every,"diagnostic phase",0,n_probephase);
		ifile.extract("species",tmpPhaseStruct.species,"diagnostic phase",0,n_probephase);
		if (tmpPhaseStruct.species.size()==0) {
			for (unsigned int i=0;i<params.n_species; i++) {
				tmpPhaseStruct.species.push_back(params.species_param[i].species_type);
			}			
		}
		ifile.extract("pos_min",tmpPhaseStruct.pos_min,"diagnostic phase",0,n_probephase);
		transform(tmpPhaseStruct.pos_min.begin(),tmpPhaseStruct.pos_min.end(), 
                  tmpPhaseStruct.pos_min.begin(),bind1st(multiplies<double>(),2*M_PI));
		ifile.extract("pos_max",tmpPhaseStruct.pos_max,"diagnostic phase",0,n_probephase);
		transform(tmpPhaseStruct.pos_max.begin(),tmpPhaseStruct.pos_max.end(), 
                  tmpPhaseStruct.pos_max.begin(),bind1st(multiplies<double>(),2*M_PI));
		ifile.extract("pos_num",tmpPhaseStruct.pos_num,"diagnostic phase",0,n_probephase);

		ifile.extract("mom_min",tmpPhaseStruct.mom_min,"diagnostic phase",0,n_probephase);
		ifile.extract("mom_max",tmpPhaseStruct.mom_max,"diagnostic phase",0,n_probephase);
		ifile.extract("mom_num",tmpPhaseStruct.mom_num,"diagnostic phase",0,n_probephase);
		
		ifile.extract("lor_min",tmpPhaseStruct.lor_min,"diagnostic phase",0,n_probephase);
		ifile.extract("lor_max",tmpPhaseStruct.lor_max,"diagnostic phase",0,n_probephase);
		ifile.extract("lor_num",tmpPhaseStruct.lor_num,"diagnostic phase",0,n_probephase);
		
		DEBUG(tmpPhaseStruct.mom_num.size() << " " << tmpPhaseStruct.pos_num.size());
		vecPhase.push_back(tmpPhaseStruct);
		n_probephase++;
	}
	
	
	
}

