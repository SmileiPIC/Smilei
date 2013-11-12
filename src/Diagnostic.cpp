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
	
	//data_.resize(2*params->n_species);
}

void Diagnostic::compute (int itime, ElectroMagn* EMfields, std::vector<Species*>& vecSpecies) {
	vector<double> mean_values_all;
	vector<unsigned int> mean_weight_all;
	mean_values.resize(0);
	mean_weight.resize(0);
	for (unsigned int ispec=0; ispec<vecSpecies.size(); ispec++) {
		mean_values.push_back(vecSpecies[ispec]->meanCharge());
	
		mean_weight.push_back(vecSpecies[ispec]->getNbrOfParticles());
		
		//data_[2*ispec]=(double) vecSpecies[ispec]->getNbrOfParticles();
		//data_[2*ispec+1]=vecSpecies[ispec]->meanCharge();
	}
		
		num_CPUs=smpi_->getSize();
		
	if (every && itime % every == 0) {
	
	if(smpi_->isMaster())
	{
	mean_values_all.resize(mean_values.size()*(num_CPUs));
			
	mean_weight_all.resize(mean_weight.size()*(num_CPUs));

	}
		// collect data from nodes
		MPI_Gather(&mean_values[0],mean_values.size(),MPI_DOUBLE,&mean_values_all[0],mean_values.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Gather(&mean_weight[0],mean_weight.size(),MPI_UNSIGNED,&mean_weight_all[0],mean_weight.size(),MPI_UNSIGNED,0,MPI_COMM_WORLD);
		
		smpi_->barrier();
		if (smpi_->isMaster()) {

			vector<double> final_mean;
			unsigned int norm;
			final_mean.resize(mean_values.size());
			
			for(unsigned int iDiag=0; iDiag<mean_values.size();iDiag++)
				{
				final_mean[iDiag]=0.0;
				norm=0.0;
				for(unsigned int iCPU=0;iCPU<num_CPUs;iCPU++)
					{
						final_mean[iDiag] += mean_values_all[iDiag*num_CPUs] * (double)(mean_weight_all[iDiag*num_CPUs]);
						norm              += mean_weight_all[iDiag*num_CPUs];
						cout<<"final "<<final_mean[iDiag]<<" norm "<<norm<<endl;
					}
					if (norm!=0) final_mean[iDiag]/=(double)norm;
				}
			
		
			fout << itime ;
			for (int i=0;i<final_mean.size();i++) fout << "\t" << final_mean[i];
			fout << std::endl;
		}
	}
}
