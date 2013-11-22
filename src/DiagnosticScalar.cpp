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
}

// init scalar
void DiagnosticScalar::init (ElectroMagn* EMfields, vector<Species*>& vecSpecies) {
	typedef map<string, double>::iterator map_iter_dbl;
	typedef map<string, unsigned int>::iterator map_iter_uint;
	typedef map<string, int>::iterator map_iter_int;

	for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) {
		vecSpecies[ispec]->initScalar();
	}
	if (smpi_->isMaster()) {
		fout.open("scalars.txt");
		mpi_spec_scalar_data.resize(smpi_->getSize());
		for (int iCPU=0; iCPU < smpi_->getSize(); iCPU++) {
			mpi_spec_scalar_data[iCPU].resize(vecSpecies.size());
			for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) {
				mpi_spec_scalar_data[iCPU][ispec]=vecSpecies[ispec]->scalar_data;
			}
		}
	}
}

// wrapper of the methods
void DiagnosticScalar::run(int timestep, ElectroMagn* EMfields, vector<Species*>& vecSpecies){
	compute_gather(timestep,EMfields,vecSpecies);
	write(timestep,vecSpecies);
}


// it contains all the methods to evaluate the scalar data
void DiagnosticScalar::compute_gather (int itime, ElectroMagn* EMfields, vector<Species*>& vecSpecies) {
	
	// 	it fills the structure "spec_scalar_data" on each specie
	unsigned int totsize_struct_char=0;
	for (unsigned int ispec=0; ispec<vecSpecies.size(); ispec++) {
		vecSpecies[ispec]->computeScalar();
		totsize_struct_char+=sizeof(double)*vecSpecies[ispec]->scalar_data.map_dbl.size();
		totsize_struct_char+=sizeof(unsigned int)*vecSpecies[ispec]->scalar_data.map_uint.size();
		totsize_struct_char+=sizeof(int)*vecSpecies[ispec]->scalar_data.map_int.size();
	}
		
	typedef map<string, double>::iterator map_iter_dbl;
	typedef map<string, unsigned int>::iterator map_iter_uint;
	typedef map<string, int>::iterator map_iter_int;

	vector<char> struct_char_transl(totsize_struct_char);
	unsigned int count=0;
	for (unsigned int ispec=0; ispec<vecSpecies.size(); ispec++) {
		for (map_iter_dbl iter = vecSpecies[ispec]->scalar_data.map_dbl.begin(); 
			 iter != vecSpecies[ispec]->scalar_data.map_dbl.end(); iter++) {
			memcpy(&struct_char_transl[count],(char*)(&(iter->second)),sizeof(double));
			count+=sizeof(double);
        }
		
		for (map_iter_uint iter = vecSpecies[ispec]->scalar_data.map_uint.begin(); 
			 iter != vecSpecies[ispec]->scalar_data.map_uint.end(); iter++) {
			memcpy(&struct_char_transl[count],(char*)(&(iter->second)),sizeof(unsigned int));
			count+=sizeof(unsigned int);
        }
		for (map_iter_int iter = vecSpecies[ispec]->scalar_data.map_int.begin(); iter != vecSpecies[ispec]->scalar_data.map_int.end(); ++iter) {
			memcpy(&struct_char_transl[count],(char*)(&(iter->second)),sizeof(int));
			count+=sizeof(int);
        }
	}
	if (count!=struct_char_transl.size()) ERROR("Something wrong here " << count << " " << struct_char_transl.size());
	
	
	// 	it constructs the receiving structure on the master processor	
    vector<char> mpi_struct_char_transl;
    if(smpi_->isMaster()){
    	mpi_struct_char_transl.resize(smpi_->getSize()*totsize_struct_char);
    }
	
	// gathering chars of the master processor
	MPI_Gather(&struct_char_transl[0],struct_char_transl.size(),MPI_CHAR,&mpi_struct_char_transl[0],struct_char_transl.size(),MPI_CHAR,0,MPI_COMM_WORLD);
	
	smpi_->barrier();
	
	// 	method to evaluate the mean charge. It is on master processor. 
	if(smpi_->isMaster()){
		unsigned int count=0;
		for(int iCPU=0;iCPU<smpi_->getSize();iCPU++){
			for (unsigned int ispec=0; ispec<vecSpecies.size(); ispec++) {
				for (map_iter_dbl iter = mpi_spec_scalar_data[iCPU][ispec].map_dbl.begin(); iter != mpi_spec_scalar_data[iCPU][ispec].map_dbl.end(); ++iter) {
					iter->second = (*(double*)(&mpi_struct_char_transl[count]));
					count+=sizeof(double);
				}				
				for (map_iter_uint iter = mpi_spec_scalar_data[iCPU][ispec].map_uint.begin(); iter != mpi_spec_scalar_data[iCPU][ispec].map_uint.end(); ++iter) {
					iter->second = (*(unsigned int*)(&mpi_struct_char_transl[count]));
					count+=sizeof(unsigned int);
				}
				for (map_iter_int iter = mpi_spec_scalar_data[iCPU][ispec].map_int.begin(); iter != mpi_spec_scalar_data[iCPU][ispec].map_int.end(); ++iter) {
					iter->second = (*(int*)(&mpi_struct_char_transl[count]));
					count+=sizeof(int);
				}				
			}
		}
		if (count!=mpi_struct_char_transl.size()) ERROR("problem here " << itime << " " << count << " " << mpi_struct_char_transl.size());
		
	}
	smpi_->barrier();
}

// It writes the scalar data on a file (it is a sequential method on master processor)
void DiagnosticScalar::write(int itime,std::vector<Species*>& vecSpecies){

	if(smpi_->isMaster()){
		vector<double> vecScalar;
		for(unsigned int ispec=0; ispec<vecSpecies.size();++ispec){
			double charge_tot=0;
			unsigned int part_tot=0;
			for(int iCPU=0;iCPU<smpi_->getSize();iCPU++){
				charge_tot+=mpi_spec_scalar_data[iCPU][ispec].map_dbl["charge_tot"];
 				part_tot+=mpi_spec_scalar_data[iCPU][ispec].map_uint["part_number"];
			}
			if (part_tot) charge_tot/=part_tot;
			vecScalar.push_back(charge_tot);
			vecScalar.push_back(part_tot);
		}

		fout << itime;
		for (unsigned int k=0;k<vecScalar.size();k++) fout << "\t" << vecScalar[k];
		fout << endl;
	}
}	
