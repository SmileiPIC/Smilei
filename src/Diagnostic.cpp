#include "Diagnostic.h"
#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Species.h"

#include <string>

using namespace std;

Diagnostic::Diagnostic( PicParams* params,  DiagParams* diagparams, SmileiMPI* smpi ) {
	every = diagparams->scalar_every;

	smpi_=smpi;
	if (smpi_->isMaster()) fout.open("scalars.txt");
	
	data_.resize(2*params->n_species);
}

void Diagnostic::compute (int itime, ElectroMagn* EMfields, std::vector<Species*>& vecSpecies) {

	for (unsigned int ispec=0; ispec<vecSpecies.size(); ispec++) {
		data_[2*ispec]=(double) vecSpecies[ispec]->getNbrOfParticles();
		data_[2*ispec+1]=vecSpecies[ispec]->meanCharge();
	}
		
	
	if (every && itime % every == 0) {
		
		// collect data from nodes
		
		
		smpi_->barrier();
		if (smpi_->isMaster()) {
			fout << itime ;
			for (int i=0;i<data_.size();i++) fout << "\t" << data_[i];
			fout << std::endl;
		}
	}
}
