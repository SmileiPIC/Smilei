#include "DiagnosticPhaseSpace.h"

#include <iomanip>
#include <string>
#include <iomanip>

#include "Params.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Field1D.h"
#include "Field.h"

using namespace std;

DiagnosticPhaseSpace::~DiagnosticPhaseSpace() {
}

void DiagnosticPhaseSpace::close() {
    //! check if we're on the master (the only one that opened the file)
	if (fileId != 0) {
        H5Fclose(fileId);
	}
}



DiagnosticPhaseSpace::DiagnosticPhaseSpace() : fileId(0) {}

void DiagnosticPhaseSpace::run(int timestep, std::vector<Species*>& vecSpecies) {
	//! check which diagnosticPhase to run at this timestep
	vector<DiagnosticPhase*> vecDiagPhaseActiveTimestep;	
	for (vector<DiagnosticPhase*>::const_iterator diag=vecDiagPhase.begin() ; diag != vecDiagPhase.end(); diag++) {
		if (timestep % (*diag)->every==0) vecDiagPhaseActiveTimestep.push_back(*diag);
	}
	
	if (vecDiagPhaseActiveTimestep.size()>0) {
        for (vector<Species*>::const_iterator mySpec=vecSpecies.begin(); mySpec!= vecSpecies.end(); mySpec++) {
	    if (!(*mySpec)->particles.isTestParticles) {
			
		//! check which diagnosticPhase to run for the species 
		vector<DiagnosticPhase*> vecDiagPhaseToRun;
		for (vector<DiagnosticPhase*>::const_iterator diag=vecDiagPhaseActiveTimestep.begin() ; diag != vecDiagPhaseActiveTimestep.end(); diag++) {
		    if(find((*diag)->my_species.begin(), (*diag)->my_species.end(), (*mySpec)->species_type) != (*diag)->my_species.end()) {
			vecDiagPhaseToRun.push_back(*diag);
		    }
		}
		if (vecDiagPhaseToRun.size()>0) {
                
		    //! cycle over all the particles
		    for (unsigned int ibin = 0 ; ibin < (*mySpec)->bmin.size() ; ibin++) {
			for (int iPart=(*mySpec)->bmin[ibin] ; iPart<(*mySpec)->bmax[ibin]; iPart++ ) {
			    //! fill the my_part structure
			    for(unsigned int k=0;k<my_part.pos.size();k++) {
				my_part.pos[k]=(*mySpec)->particles.position(k,iPart);
			    }
			    for(unsigned int k=0;k<3;k++) {
				my_part.mom[k]=(*mySpec)->particles.momentum(k,iPart);
			    }
			    my_part.weight=(*mySpec)->particles.weight(iPart);
			    my_part.charge=(*mySpec)->particles.charge(iPart);
			    for (vector<DiagnosticPhase*>::const_iterator diag=vecDiagPhaseToRun.begin() ; diag != vecDiagPhaseToRun.end(); diag++) {
				//! do something with each particle
				(*diag)->run(my_part);
			    }						
			}
		    }
		}
	    }
	} // END for (vector<Species*>::const_iterator mySpec
        //! and finally write the data (reduce data on 1 proc, write it and clear memory for future usage)
        for (vector<DiagnosticPhase*>::const_iterator diag=vecDiagPhaseActiveTimestep.begin() ; diag != vecDiagPhaseActiveTimestep.end(); diag++) {
            (*diag)->writeData();
        }
	}
}
