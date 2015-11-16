#include "DiagnosticProbe.h"

#include <iomanip>
#include <string>
#include <iomanip>
#include <sstream>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Field1D.h"
#include "Field2D.h"
#include "Field.h"
#include "H5.h"

using namespace std;

//DiagnosticProbe::DiagnosticProbe(): fileId(0) {}
DiagnosticProbe::DiagnosticProbe(Params &params, Patch* patch):
cpuRank((int)patch->Hindex()),
probeSize(10), 
fileId(0) {

    unsigned  numProbes=PyTools::nComponents("DiagProbe");
    every.resize(numProbes);
    probeParticles.resize(numProbes);

    probesArray.resize(numProbes);
    probesStart.resize(numProbes);

    // Loop on all Probes
    for (unsigned int n_probe=0; n_probe<numProbes; n_probe++) {

	bool ok;
        ok=PyTools::extract("every",every,"DiagProbe",n_probe);        
        if (!ok) every[n_probe]=params.global_every;

	unsigned int ndim=params.nDim_particle;
            
        // Extract "number" (number of points you have in each dimension of the probe,
        // which must be smaller than the code dimensions)
        vector<unsigned int> vecNumber; 
        PyTools::extract("number",vecNumber,"DiagProbe",n_probe);

        // Dimension of the probe grid
        unsigned int dimProbe=vecNumber.size();
        if (dimProbe > params.nDim_particle) {
            ERROR("Probe #"<<n_probe<<": probe dimension is greater than simulation dimension")
        }
        
        // If there is no "number" argument provided, then it corresponds to
        // a zero-dimensional probe (one point). In this case, we say the probe
        // has actually one dimension with only one point.
        unsigned int dim=vecNumber.size();
        if (vecNumber.size() == 0) {
            vecNumber.resize(1);
            vecNumber[0]=1;
        }
            
	unsigned int totPart=1;
	for (unsigned int iDimProbe=0; iDimProbe<dim; iDimProbe++) {
	    totPart *= vecNumber[iDimProbe];
	}

	// -----------------------------------------------------
	// Start definition of probeParticles (probes positions)
	// -----------------------------------------------------            
	probeParticles[n_probe].initialize(totPart, params, ndim);
            
	vector<double> partPos(ndim*totPart,0.0);
	for(unsigned int ipart=0; ipart!=totPart; ++ipart) {
	    int found=cpuRank;
	    for(unsigned int iDim=0; iDim!=ndim; ++iDim) {
		MESSAGE ( "Attention, new definition in Diagnostic::initProbes !!!" );
		partPos[iDim+ipart*ndim]=0.;//diagParams.probeStruc[n_probe].pos[0][iDim];
		// the particle position is a linear combiantion of the point pos with posFirst or posSecond or posThird
		for (unsigned int iDimProbe=0; iDimProbe<dim; iDimProbe++) {
		    partPos[iDim+ipart*ndim] += 0.;//(ipart%vecNumber[iDimProbe])*(diagParams.probeStruc[n_probe].pos[iDimProbe+1][iDim]-diagParams.probeStruc[n_probe].pos[0][iDim])/(vecNumber[iDimProbe]-1);
		}
		probeParticles[n_probe].position(iDim,ipart) = partPos[iDim+ipart*ndim];
                    
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

		if (probeParticles[n_probe].position(iDim,ipart) < patch->getDomainLocalMin(iDim) ||
		    probeParticles[n_probe].position(iDim,ipart) >= maxToCheck) {
		    found=-1;
		}
	    }
	}

	nProbeTot = probeParticles[n_probe].size();
	//cout << " \t Before " << cpuRank << " nprobes : " << probeParticles[0].size() << endl;

	for ( int ipb=nProbeTot-1 ; ipb>=0 ; ipb--) {
	    if (!probeParticles[n_probe].is_part_in_domain(ipb, patch))
		probeParticles[n_probe].erase_particle(ipb);
	    //cout << patch->getDomainLocalMin(0) << " " << probeParticles[n_probe].position(0, ipb) << " " << patch->getDomainLocalMax(0) << endl;
	    //cout << patch->getDomainLocalMin(1) << " " << probeParticles[n_probe].position(1, ipb) << " " << patch->getDomainLocalMax(1) << endl;
	}
	//cout << " \t After " << cpuRank << " nprobes : " << probeParticles[0].size() << endl;
	// ---------------------------------------------------
	// End definition of probeParticles (probes positions)
	// ---------------------------------------------------

	// probesArray : n_probe vectors x 10 vectors x probeParticles[n_probe].size() double
	vector<unsigned int> probesArraySize(2);
	probesArraySize[0] = probeParticles[n_probe].size();
	probesArraySize[1] = probeSize;
	probesArray[n_probe] = new Field2D(probesArraySize);



    }
    //cout << " nprobes : " << probeParticles[0].size() << endl;
 

}

// Done by patch master only
void DiagnosticProbe::createFile()
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


        vector<unsigned int> vecNumber; 
	PyTools::extract("number",vecNumber,"DiagProbe",probe_id);
        // Dimension of the probe grid
        unsigned int dimProbe=vecNumber.size();

	sid = H5Screate(H5S_SCALAR);	
	aid = H5Acreate(group_id, "dimension", H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(aid, H5T_NATIVE_UINT, &dimProbe);
	H5Aclose(aid);
	H5Sclose(sid);

	// Close the group
	H5Gclose(group_id);

    }

}


void DiagnosticProbe::setFile(hid_t masterFileId, Patch* patch, Params& params)
{
    fileId = masterFileId;

    // ---------------------------------------------------
    // Start file split definition
    // ---------------------------------------------------
    int nPatches(1);
    for (int iDim=0;iDim<params.nDim_field;iDim++) nPatches*=params.number_of_patches[iDim];
    // probesStart
    unsigned  numProbes=PyTools::nComponents("DiagProbe");
    for (unsigned int np=0; np<numProbes; np++) {
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
    }

    // ---------------------------------------------------
    // End file split definition
    // ---------------------------------------------------

}

void DiagnosticProbe::setFile(hid_t masterFileId)
{
    fileId = masterFileId;  
}

void DiagnosticProbe::writePositionIn( Params &params )
{
    // Loop on all Probes
    for (unsigned int probe_id=0; probe_id<every.size(); probe_id++) {

        // Dimension of the probe grid
	vector<unsigned int> vecNumber; 
        PyTools::extract("number",vecNumber,"DiagProbe",probe_id);
        unsigned int dimProbe=vecNumber.size();

	// Open group de write in	
        hid_t group_id = H5Gopen(fileId, probeName(probe_id).c_str() ,H5P_DEFAULT);
	// Write in the did group
	writePositions(probe_id, params.nDim_particle, dimProbe, group_id);
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
    for ( unsigned int np=0 ; np < probesArray.size() ; np++ )
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


void DiagnosticProbe::run(unsigned int timestep, ElectroMagn* EMfields, Interpolator* interp) {

    double time = (double)timestep * dt;
    
    // Loop probes
    for (unsigned int np=0; np<every.size(); np++) {
        // skip if current timestep is not requested
        if ( (every[np]  && timestep % every[np] == 0) &&
             (time <= tmax[np]) && (time >= tmin[np]) ) {

            // Open the existing HDF5 group for that probe
	    hid_t group_id = H5Gopen(fileId, probeName(np).c_str() ,H5P_DEFAULT);
            //hid_t did = H5Gopen2(fileId, probeName(np).c_str(), H5P_DEFAULT);

	    // Write in group_id
	    compute(np, timestep, EMfields, interp);
	    write(np, timestep, group_id);
	    // Close the group
	    H5Gclose(group_id);

	    if (fileId) H5Fflush(fileId, H5F_SCOPE_GLOBAL );
	}
 
   }
}
            
void DiagnosticProbe::compute(int probe_id, unsigned int timestep, ElectroMagn* EMfields, Interpolator* interp) {
    
    // Loop probe ("fake") particles
    for (int iprob=0; iprob <probeParticles[probe_id].size(); iprob++) {             
	(*interp)(EMfields,probeParticles[probe_id],iprob,&Eloc_fields,&Bloc_fields,&Jloc_fields,&Rloc_fields);

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
	probesArray[probe_id]->data_2D[iprob][probeSize-1]=Rloc_fields;
        
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
    //cout << " CPU Rank " << cpuRank << " - writing at "  << probesStart[probe_id] << endl;
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

#ifdef _NEWSTYLE
    // Open the existing HDF5 group for that probe
    hid_t did = H5Gopen2(fileId, probeName(np).c_str(), H5P_DEFAULT);
    // Write the positions array into the current HDF5 group
    H5::matrix_MPI(did, name_t.str(), probesArray[np]->data_2D[0][0], nPart_total[np], nFields[np], probesStart[np], nPart_local);
    // Close the group
    H5Gclose(did);
#endif

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
