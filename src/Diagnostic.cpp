#include "Diagnostic.h"
#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Species.h"


Diagnostic::Diagnostic( PicParams* params,  DiagParams* diagparams, SmileiMPI* smpi ) {
	every = diagparams->scalar_every;
	fout.open("scalars.txt");
	data_.resize(params->n_species);
}

void Diagnostic::compute (int itime, ElectroMagn* EMfields, std::vector<Species*>& vecSpecies) {
	
	for (unsigned int ispec; ispec<vecSpecies.size(); ispec++) {
		data_[ispec]=vecSpecies[ispec]->meanCharge();
	}
	
	if (itime % every == 0) {
		DEBUG("here " << itime);
		fout << itime ;
		for (int i=0;i<data_.size();i++) fout << "\t" << data_[i];
		fout << std::endl ;
	}
	
}
