#include "DiagnosticPhase1D.h"

#include <iomanip>
#include <string>
#include <iomanip>

#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Field1D.h"
#include "Field.h"

using namespace std;

DiagnosticPhase1D::~DiagnosticPhase1D() {
}

void DiagnosticPhase1D::close() {
}

DiagnosticPhase1D::DiagnosticPhase1D(PicParams* params, DiagParams* diagParams, SmileiMPI* smpi, unsigned int n) :
    smpi_(smpi) 
{

	every=diagParams->phase1D[n].every;
}

void DiagnosticPhase1D::run(int timestep, std::vector<Species*>& vecSpecies) {
	
	
	//! momentum min
	std::vector< std::vector<double> > momentum_min;
	
	//! momentum max
	std::vector< std::vector<double> > momentum_max;
	
	//! gamma min
	std::vector<double> lorentz_factor_min;
	
	//! gamma max
	std::vector<double> lorentz_factor_max;

	//! momentum min
	momentum_min.resize(vecSpecies.size());
	momentum_max.resize(vecSpecies.size());
	
	
	for (unsigned int i=0; i<vecSpecies.size(); i++) {
		momentum_min[i].resize(3);
		momentum_max[i].resize(3);
	}
	
	
	lorentz_factor_min.resize(vecSpecies.size());
	lorentz_factor_max.resize(vecSpecies.size());
	
	
	for (unsigned int j=0; j < vecSpecies.size(); j++) {
		
		// initialize momentum min/max for phase space diags
		if (vecSpecies[j]->particles.size()>0) {
			for (unsigned int i=0; i<3; i++) {
				momentum_min[j][i] = momentum_max[j][i] = vecSpecies[j]->particles.momentum(i,0);
			}
			lorentz_factor_min[j]=lorentz_factor_max[j]=vecSpecies[j]->particles.lor_fac(0);
		}
		
		for (unsigned int ibin = 0 ; ibin < vecSpecies[j]->bmin.size() ; ibin++) {
			for (int iPart=vecSpecies[j]->bmin[ibin] ; iPart<vecSpecies[j]->bmax[ibin]; iPart++ ) {
				
				for (unsigned int i=0; i<3; i++) {
					if (vecSpecies[j]->particles.momentum(i,iPart) < momentum_min[j][i]) momentum_min[j][i] = vecSpecies[j]->particles.momentum(i,iPart);
					if (vecSpecies[j]->particles.momentum(i,iPart) > momentum_max[j][i]) momentum_max[j][i] = vecSpecies[j]->particles.momentum(i,iPart);					
				}
				if (vecSpecies[j]->particles.lor_fac(iPart) < lorentz_factor_min[j]) lorentz_factor_min[j] = vecSpecies[j]->particles.lor_fac(iPart);
				if (vecSpecies[j]->particles.lor_fac(iPart) > lorentz_factor_max[j]) lorentz_factor_max[j] = vecSpecies[j]->particles.lor_fac(iPart);
				
			}
		}
	}
	
	for (unsigned int j=0; j < vecSpecies.size(); j++) {
		for (unsigned int i=0; i<3; i++) {
			double tempmin=momentum_min[j][i];
			MPI_Allreduce(&tempmin,&momentum_min[j][i],1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
			double tempmax=momentum_max[j][i];
			MPI_Allreduce(&tempmax,&momentum_max[j][i],1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
		}
		double gammamin=lorentz_factor_min[j];
		MPI_Allreduce(&gammamin,&lorentz_factor_min[j],1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
		double gammamax=lorentz_factor_max[j];
		MPI_Allreduce(&gammamax,&lorentz_factor_max[j],1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

		if (smpi_->isMaster()) {
			for (unsigned int k=0; k <momentum_min[j].size(); k++) {				
				DEBUG(timestep << " spec. " << j << " min: " << momentum_min[j][k] << " max: "<< momentum_max[j][k]);
			}
			DEBUG(timestep << " spec. " << j << " Lmin: " << lorentz_factor_min[j] << " Lmax: "<< lorentz_factor_max[j]);
		}
	}

}
