#include "DiagnosticScalar.h"

#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"

#include <string>

using namespace std;

// constructor
DiagnosticScalar::DiagnosticScalar(PicParams* params, SmileiMPI* smpi) {
	smpi_=smpi;
	if (smpi_->isMaster()) {
		fout.open("scalars.txt");
		if (!fout.is_open()) ERROR("can't open scalar file");
	}
}

// wrapper of the methods
void DiagnosticScalar::run(int timestep, ElectroMagn* EMfields, vector<Species*>& vecSpecies){
	compute_gather(timestep,EMfields,vecSpecies);
	write(timestep,vecSpecies);
}


// it contains all to manage the communication of data. It is "transparent" to the user.
void DiagnosticScalar::compute_gather (int itime, ElectroMagn* EMfields, vector<Species*>& vecSpecies) {
	// 	it fills the map on each specie
	for (unsigned int ispec=0; ispec<vecSpecies.size(); ispec++) {
		vecSpecies[ispec]->computeScalar();		
	}
	
	// definition of the memory allocation vector
	vector<double> oneProc;
	for (unsigned int ispec=0; ispec<vecSpecies.size(); ispec++) {
		for (map<string, double>::iterator iter = vecSpecies[ispec]->scalars.begin(); iter != vecSpecies[ispec]->scalars.end(); iter++) {
			oneProc.push_back(iter->second);
        }
	}
	
	// 	it constructs the receiving structure on the master processor	
    vector<double> allProcs;
    if(smpi_->isMaster()){
    	allProcs.resize(smpi_->getSize()*oneProc.size());
    }
	// gathering chars of the master processor
	MPI_Gather(&oneProc[0],oneProc.size(),MPI_DOUBLE,&allProcs[0],oneProc.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	smpi_->barrier();
	
	// 	method to reconstruct the information on the master processor
	if(smpi_->isMaster()){
		unsigned int count=0;
		mpi_spec_scalars.resize(smpi_->getSize());
		for(int iCPU=0;iCPU<smpi_->getSize();iCPU++){
			mpi_spec_scalars[iCPU].resize(vecSpecies.size());
			for (unsigned int ispec=0; ispec<vecSpecies.size(); ispec++) {
				for (map<string, double>::iterator iter = vecSpecies[ispec]->scalars.begin(); iter != vecSpecies[ispec]->scalars.end(); ++iter) {
					mpi_spec_scalars[iCPU][ispec][iter->first] = allProcs[count];
					count++;
				}
			}			
		}
		if (count!=allProcs.size()) ERROR("problem here " << itime << " " << count << " " << allProcs.size());
	}
}

// Each single method of the diagnostic scalar must implemented here. it writes also on a file.
void DiagnosticScalar::write(int itime,std::vector<Species*>& vecSpecies){
	
	if(smpi_->isMaster()){
		vector<double> vecScalar;
		for(unsigned int ispec=0; ispec<vecSpecies.size();++ispec){
			double charge_tot=0;
			unsigned int part_tot=0;
			for(int iCPU=0;iCPU<smpi_->getSize();iCPU++){
				charge_tot+=mpi_spec_scalars[iCPU][ispec]["charge_tot"];
 				part_tot+=mpi_spec_scalars[iCPU][ispec]["part_number"];
			}
			if (part_tot) charge_tot/=part_tot;
			
			vecScalar.push_back(charge_tot);
			vecScalar.push_back(part_tot);
		}
		
		//		for(int iCPU=0;iCPU<smpi_->getSize();iCPU++){
		//			double totosum=0.0;
		//			for(unsigned int ispec=0; ispec<vecSpecies.size();++ispec){
		// 				totosum+=mpi_spec_scalars[iCPU][ispec]["toto"];
		//			}
		//			vecScalar.push_back(totosum);
		//		}
		
		fout << itime;
		for (unsigned int k=0;k<vecScalar.size();k++) fout << "\t" << vecScalar[k];
		fout << endl;
	}
}	
