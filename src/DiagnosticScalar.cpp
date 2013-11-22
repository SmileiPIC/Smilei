#include "DiagnosticScalar.h"

#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Species.h"

#include <string>

using namespace std;

// constructor
DiagnosticScalar::DiagnosticScalar(SmileiMPI* smpi) {
	smpi_=smpi;
	if (smpi_->isMaster()) fout.open("scalars.txt");
	
	// translation of the struct data type in a MPI_Datatype (parameters must be modified if the number of scalar diagnostics changes)	
	num_CPUs=smpi_->getSize();
	int nitems=2;
    int blocklengths[2] = {1,1};
    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_UNSIGNED};
    MPI_Aint offsets[2];
    offsets[0] = offsetof(spec_scalar_data, mean_charge);
    offsets[1] = offsetof(spec_scalar_data, part_number);
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_data_scalar);
    MPI_Type_commit(&mpi_data_scalar);
}

// wrapper of the methods
void DiagnosticScalar::run(int timestep, ElectroMagn* EMfields, vector<Species*>& vecSpecies){
	compute(timestep,EMfields,vecSpecies);
	write(timestep,vecSpecies);
}

// it contains all the methods to evaluate the scalar data
void DiagnosticScalar::compute (int itime, ElectroMagn* EMfields, vector<Species*>& vecSpecies) {
	
	// 	it fills the structure "spec_scalar_data" on each specie	
	for (unsigned int ispec=0; ispec<vecSpecies.size(); ispec++) {
		vecSpecies[ispec]->computeScalar();
	}
	
	// 	it constructs the receiving structure on the master processor	
    vector<spec_scalar_data> mpi_data_scalar_all;
    if(smpi_->isMaster()){
    	mpi_data_scalar_all.resize(num_CPUs*vecSpecies.size());
    }
	
	// it constructs a buffer of structs to be sent
    vector<spec_scalar_data> mpi_data_scalar_node(vecSpecies.size());
	
	for (unsigned int ispec=0; ispec<vecSpecies.size(); ispec++) {
		mpi_data_scalar_node[ispec]=vecSpecies[ispec]->scalar_struct();
	}
	
	// gathering of structs of the master processor
	MPI_Gather(&mpi_data_scalar_node[0],vecSpecies.size(),mpi_data_scalar,&mpi_data_scalar_all[0],vecSpecies.size(),mpi_data_scalar,0,MPI_COMM_WORLD);
	
	smpi_->barrier();
	
	// 	method to evaluate the mean charge. It is on master processor. 
	if(smpi_->isMaster()){
		vecScalar.clear();
		for(unsigned int ispec=0; ispec<vecSpecies.size();++ispec){
			double charge_tot=0;
			unsigned int part_tot=0;
			for(unsigned int iCPU=0;iCPU<num_CPUs;iCPU++){
				int k=ispec+vecSpecies.size()*iCPU;
		 		charge_tot+=mpi_data_scalar_all[k].mean_charge*mpi_data_scalar_all[k].part_number;
 				part_tot+=mpi_data_scalar_all[k].part_number;
			}
			if (part_tot) charge_tot/=part_tot;
			vecScalar.push_back(charge_tot);
		}
	}	
}

// It writes the scalar data on a file (it is a sequential method on master processor)
void DiagnosticScalar::write(int itime,std::vector<Species*>& vecSpecies){
	if(smpi_->isMaster()){
		fout << itime;
		for (unsigned int k=0;k<vecScalar.size();k++) fout << "\t" << vecScalar[k];
		fout << endl;
	}
}	
