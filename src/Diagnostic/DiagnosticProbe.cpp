#include "DiagnosticProbe.h"

#include <iomanip>
#include <string>
#include <iomanip>
#include <sstream>

#include "PicParams.h"
#include "SmileiMPI.h"
#include "Patch.h"
#include "SmileiMPI_Cart1D.h"
#include "SmileiMPI_Cart2D.h"
#include "ElectroMagn.h"
#include "Field1D.h"
#include "Field2D.h"
#include "Field.h"
#include "DiagParams.h"

using namespace std;

DiagnosticProbe::DiagnosticProbe(PicParams &params, DiagParams &diagParams, Patch* patch):
cpuRank((int)patch->Hindex()),
probeSize(10), 
fileId(0) {

    every.resize(diagParams.probeStruc.size());
    probeParticles.resize(diagParams.probeStruc.size());

    probesArray.resize(diagParams.probeStruc.size());
    probesStart.resize(diagParams.probeStruc.size());

    // Loop on all Probes
    for (unsigned int np=0; np<diagParams.probeStruc.size(); np++) {
	every[np]=diagParams.probeStruc[np].every;
	unsigned int dimProbe=diagParams.probeStruc[np].dim+2;
	unsigned int ndim=params.nDim_particle;
            
	vector<unsigned int> vecNumber=diagParams.probeStruc[np].number;
            
	unsigned int totPart=1;
	for (unsigned int iDimProbe=0; iDimProbe<diagParams.probeStruc[np].dim; iDimProbe++) {
	    totPart *= vecNumber[iDimProbe];
	}

	// -----------------------------------------------------
	// Start definition of probeParticles (probes positions)
	// -----------------------------------------------------            
	probeParticles[np].initialize(totPart, ndim);
            
	vector<double> partPos(ndim*totPart,0.0);
	for(unsigned int ipart=0; ipart!=totPart; ++ipart) {
	    int found=cpuRank;
	    for(unsigned int iDim=0; iDim!=ndim; ++iDim) {
		partPos[iDim+ipart*ndim]=diagParams.probeStruc[np].pos[0][iDim];
		// the particle position is a linear combiantion of the point pos with posFirst or posSecond or posThird
		for (unsigned int iDimProbe=0; iDimProbe<diagParams.probeStruc[np].dim; iDimProbe++) {
		    partPos[iDim+ipart*ndim] += (ipart%vecNumber[iDimProbe])*(diagParams.probeStruc[np].pos[iDimProbe+1][iDim]-diagParams.probeStruc[np].pos[0][iDim])/(vecNumber[iDimProbe]-1);
		}
		probeParticles[np].position(iDim,ipart) = partPos[iDim+ipart*ndim];
                    
		//!fixme this is awful: we add one cell if we're on the upper border
		double maxToCheck=patch->getDomainLocalMax(iDim);                    
		if (ndim==1) {
		    if (patch->isEastern()) {
			maxToCheck+=params.cell_length[iDim];
		    }
		} else if (ndim==2) {
		    if ((iDim == 0 && patch->isEastern()) ||
			(iDim == 1 && patch->isNorthern())) {
			maxToCheck+=params.cell_length[iDim];
		    }                        
		} else {
		    ERROR("implement here");
		}

		if (probeParticles[np].position(iDim,ipart) < patch->getDomainLocalMin(iDim) ||
		    probeParticles[np].position(iDim,ipart) >= maxToCheck) {
		    found=-1;
		}
	    }
	}

	nProbeTot = probeParticles[np].size();
	//cout << " \t Before " << cpuRank << " nprobes : " << probeParticles[0].size() << endl;

	for ( int ipb=nProbeTot-1 ; ipb>=0 ; ipb--) {
	    if (!probeParticles[np].is_part_in_domain(ipb, patch))
		probeParticles[np].erase_particle(ipb);
	    //cout << patch->getDomainLocalMin(0) << " " << probeParticles[np].position(0, ipb) << " " << patch->getDomainLocalMax(0) << endl;
	    //cout << patch->getDomainLocalMin(1) << " " << probeParticles[np].position(1, ipb) << " " << patch->getDomainLocalMax(1) << endl;
	}
	//cout << " \t After " << cpuRank << " nprobes : " << probeParticles[0].size() << endl;
	// ---------------------------------------------------
	// End definition of probeParticles (probes positions)
	// ---------------------------------------------------

	// probesArray : np vectors x 10 vectors x probeParticles[np].size() double
	vector<unsigned int> probesArraySize(2);
	probesArraySize[0] = probeParticles[np].size();
	probesArraySize[1] = probeSize;
	probesArray[np] = new Field2D(probesArraySize);


	// ---------------------------------------------------
	// Start file split definition
	// ---------------------------------------------------
	int nPatches(1);
	for (int iDim=0;iDim<params.nDim_field;iDim++) nPatches*=params.number_of_patches[iDim];
	// probesStart
	probesStart[np] = 0;
	MPI_Status status;
	stringstream rtag("");
	rtag << cpuRank-1 << "0" << cpuRank;
	int tag(0);
	rtag >> tag;

	if (cpuRank>0) {
	    //cout << patch->Hindex() << " Recv from " << patch->getMPIRank(cpuRank-1) << " with tag " << tag << endl;
	    MPI_Recv( &(probesStart[np]), 1, MPI_INTEGER, patch->getMPIRank(cpuRank-1), tag, MPI_COMM_WORLD, &status );
	}
	    
	int probeEnd = probesStart[np]+probeParticles[np].size();
	stringstream stag("");
	stag << cpuRank << "0" << cpuRank+1;
	tag = 0;
	stag >> tag;
	if (cpuRank!=nPatches-1) {
	    //cout << patch->Hindex() << " Send to " << patch->getMPIRank(cpuRank+1) << " with tag " << tag << endl;
	    MPI_Send( &probeEnd, 1, MPI_INTEGER, patch->getMPIRank(cpuRank+1), tag, MPI_COMM_WORLD );

	}
	// ---------------------------------------------------
	// End file split definition
	// ---------------------------------------------------
    }
    //cout << " nprobes : " << probeParticles[0].size() << endl;
 

}

// Done by patch master only
void DiagnosticProbe::createFile(DiagParams &diagParams)
{
    if (!every.size())
	return;
    
    hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
    fileId = H5Fcreate( "Probes.h5", H5F_ACC_TRUNC, H5P_DEFAULT, pid);
    H5Pclose(pid);
        
    string ver(__VERSION);
    hid_t sid = H5Screate(H5S_SCALAR);
    hid_t tid = H5Tcopy(H5T_C_S1);
    H5Tset_size(tid, ver.size());
    H5Tset_strpad(tid,H5T_STR_NULLTERM);

    hid_t aid;
    aid = H5Acreate(fileId, "Version", tid, sid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(aid, tid, ver.c_str());       
    H5Aclose(aid);
    H5Sclose(sid);
    H5Tclose(tid);

    // Loop on all Probes
    for (unsigned int probe_id=0; probe_id<every.size(); probe_id++) {

	// Open group de write in	
	hid_t group_id = H5Gcreate(fileId, probeName(probe_id).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	hid_t sid;
	sid = H5Screate(H5S_SCALAR);	
	hid_t aid = H5Acreate(group_id, "every", H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(aid, H5T_NATIVE_UINT, &every[probe_id]);
	H5Aclose(aid);
	H5Sclose(sid);
            
	sid = H5Screate(H5S_SCALAR);	
	aid = H5Acreate(group_id, "dimension", H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(aid, H5T_NATIVE_UINT, &diagParams.probeStruc[probe_id].dim);
	H5Aclose(aid);
	H5Sclose(sid);

	// Close the group
	H5Gclose(group_id);

    }

}


void DiagnosticProbe::setFile(hid_t masterFileId)
{
    fileId = masterFileId;
}


void DiagnosticProbe::writePositionIn( PicParams &params, DiagParams &diagParams )
{
    // Loop on all Probes
    for (unsigned int probe_id=0; probe_id<every.size(); probe_id++) {

	// Open group de write in	
        hid_t group_id = H5Gopen(fileId, probeName(probe_id).c_str() ,H5P_DEFAULT);
	// Write in the did group
	writePositions(probe_id, params.nDim_particle, diagParams.probeStruc[probe_id].dim, group_id);
	// Close the group
	H5Gclose(group_id);

    }
    
}


void DiagnosticProbe::writePositions(int probe_id, int ndim_Particles, int probeDim, hid_t group_id ) {
    //cout << probeParticles[probe_id].size() << endl;
    //if (!probeParticles[probe_id].size()) return;

    vector<unsigned int> posArraySize(2);
    posArraySize[0] = probeParticles[probe_id].size();
    posArraySize[1] = ndim_Particles;
    Field2D* posArray = new Field2D(posArraySize);
    for ( int ipb=0 ; ipb<probeParticles[probe_id].size() ; ipb++) {
	for (int idim=0 ; idim<ndim_Particles  ; idim++ )
	    posArray->data_2D[ipb][idim] = probeParticles[probe_id].position(idim,ipb);
    }


    // memspace OK : 1 block 
    hsize_t     chunk_parts[2];
    chunk_parts[0] = probeParticles[probe_id].size();
    chunk_parts[1] = 2; 
    hid_t memspace  = H5Screate_simple(2, chunk_parts, NULL);
    // filespace :
    hsize_t dimsf[2], offset[2], stride[2], count[2];
    dimsf[0] = nProbeTot;
    dimsf[1] = 2;
    hid_t filespace = H5Screate_simple(2, dimsf, NULL);
    offset[0] = probesStart[probe_id];
    offset[1] = 0;
    stride[0] = 1;
    stride[1] = 1;
    count[0] = 1;
    count[1] = 1;
    hsize_t     block[2];
    block[0] = probeParticles[probe_id].size();
    block[1] = ndim_Particles;
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);

    //define , write_plist
    hid_t write_plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_INDEPENDENT);
    hid_t plist_id;
    hid_t dset_id;
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    //else 
    //plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    //cout << "Before create" << endl;
    if ( !H5Lexists( group_id, "positions", H5P_DEFAULT ) )
	dset_id = H5Dcreate(group_id, "positions", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    else
	dset_id = H5Dopen(group_id, "positions", H5P_DEFAULT);
    //cout << "After create" << endl;


    H5Pclose(plist_id);
    H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, write_plist, &(posArray->data_2D[0][0]) );
    H5Dclose(dset_id);
    H5Pclose( write_plist );

    H5Sclose(filespace);
    H5Sclose(memspace);



    delete posArray;

}


DiagnosticProbe::~DiagnosticProbe()
{

    for ( int np=0 ; np < probesArray.size() ; np++ )
	delete probesArray[np];

}

void DiagnosticProbe::close() {
    if (fileId>0) {
        H5Fclose(fileId);
    }
    /*cout << " I 'm here " << endl;
    H5Fflush(fileId, H5F_SCOPE_GLOBAL );
    H5Fclose(fileId);*/
}

string DiagnosticProbe::probeName(int p) {
    ostringstream prob_name("");
    prob_name << "p" << setfill('0') << setw(4) << p;
    return prob_name.str();
}


void DiagnosticProbe::run( unsigned int timestep, ElectroMagn* EMfields, Interpolator* interp )
{
    // Loop on all Probes
    for (unsigned int probe_id=0; probe_id<every.size(); probe_id++) {

	if (every[probe_id] && timestep % every[probe_id] == 0) {

	    // Open group de write in	
	    hid_t group_id = H5Gopen(fileId, probeName(probe_id).c_str() ,H5P_DEFAULT);
	    // Write in group_id
	    compute(probe_id, timestep, EMfields, interp);
	    write(probe_id, timestep, group_id);
	    // Close the group
	    H5Gclose(group_id);

	    if (fileId) H5Fflush(fileId, H5F_SCOPE_GLOBAL );

	}

    }
    
}


void DiagnosticProbe::compute(int probe_id, unsigned int timestep, ElectroMagn* EMfields, Interpolator* interp) {
    for (int iprob=0; iprob <probeParticles[probe_id].size(); iprob++) {             
	(*interp)(EMfields,probeParticles[probe_id],iprob,&Eloc_fields,&Bloc_fields,&Jloc_fields,&probesArray[probe_id]->data_2D[iprob][probeSize-1]);

	//! here we fill the probe data!!!
	probesArray[probe_id]->data_2D[iprob][0]=Eloc_fields.x;
	probesArray[probe_id]->data_2D[iprob][1]=Eloc_fields.y;
	probesArray[probe_id]->data_2D[iprob][2]=Eloc_fields.z;
	probesArray[probe_id]->data_2D[iprob][3]=Bloc_fields.x;
	probesArray[probe_id]->data_2D[iprob][4]=Bloc_fields.y;
	probesArray[probe_id]->data_2D[iprob][5]=Bloc_fields.z;
	probesArray[probe_id]->data_2D[iprob][6]=Jloc_fields.x;
	probesArray[probe_id]->data_2D[iprob][7]=Jloc_fields.y;
	probesArray[probe_id]->data_2D[iprob][8]=Jloc_fields.z;          
                
    } // End for iprob

}


void DiagnosticProbe::write(int probe_id, unsigned int timestep, hid_t group_id) {
    // memspace OK : 1 block 
    hsize_t     chunk_parts[2];
    chunk_parts[0] = probeParticles[probe_id].size();
    chunk_parts[1] = probeSize; 
    hid_t memspace  = H5Screate_simple(2, chunk_parts, NULL);
    // filespace :
    hsize_t dimsf[2], offset[2], stride[2], count[2];
    dimsf[0] = nProbeTot;
    dimsf[1] = probeSize;
    hid_t filespace = H5Screate_simple(2, dimsf, NULL);
    offset[0] = probesStart[probe_id];
    offset[1] = 0;
    stride[0] = 1;
    stride[1] = 1;
    count[0] = 1;
    count[1] = 1;
    hsize_t     block[2];
    block[0] = probeParticles[probe_id].size();
    block[1] = probeSize;
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);

    // define filespace, memspace
    hid_t write_plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_INDEPENDENT);
    hid_t did = H5Gopen2(fileId, probeName(probe_id).c_str(), H5P_DEFAULT);
    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
    ostringstream name_t;
    name_t.str("");
    name_t << "/" << probeName(probe_id).c_str() << "/" << setfill('0') << setw(10) << timestep;

    hid_t dset_id;
    htri_t status = H5Lexists( group_id, name_t.str().c_str(), H5P_DEFAULT ); 
    if (!status)
	dset_id  = H5Dcreate(group_id, name_t.str().c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    else
	dset_id = H5Dopen(group_id, name_t.str().c_str(), H5P_DEFAULT);		

    H5Pclose(plist_id);
    H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, write_plist, &(probesArray[probe_id]->data_2D[0][0]) );
    H5Dclose(dset_id);
    H5Gclose(did);
    H5Pclose( write_plist );

    H5Sclose(filespace);
    H5Sclose(memspace);

}
