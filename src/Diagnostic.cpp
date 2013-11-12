#include "Diagnostic.h"
#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Species.h"

#include <string>
#include <sstream>

using namespace std;

Diagnostic::Diagnostic( PicParams* params,  DiagParams* diagparams, SmileiMPI* smpi ) {
	every = diagparams->scalar_every;

	stringstream name("");
	name << "scalars-" << smpi->getRank() << ".txt" ;

	fout.open(name.str().c_str());
	data_.resize(params->n_species);
}

void Diagnostic::compute (int itime, ElectroMagn* EMfields, std::vector<Species*>& vecSpecies) {

	for (unsigned int ispec=0; ispec<vecSpecies.size(); ispec++) {
		data_[ispec]=vecSpecies[ispec]->meanCharge();
//		DEBUG(data_[ispec]);
	}
		
	if (every && itime % every == 0) {
		fout << itime ;
		for (int i=0;i<data_.size();i++) fout << "\t" << data_[i];
		fout << std::endl;
	}
}
