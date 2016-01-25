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

DiagnosticProbe::DiagnosticProbe(Params &params, Patch* patch):
cpuRank((int)patch->Hindex()),
probeSize(10), 
fileId(0) {

    unsigned  numProbes=PyTools::nComponents("DiagProbe");
    every.resize(numProbes);
    probeParticles.resize(numProbes);

    probesArray.resize(numProbes);
    probesStart.resize(numProbes);

    dt = params.timestep;
    every         .clear();
    tmin          .clear();
    tmax          .clear();
    probeParticles.clear();
    nPart_total   .clear();
    probesArray   .clear();
    probesStart   .clear();
    fieldname     .clear();
    fieldlocation .clear();
    nFields       .clear();

    // Loop on all Probes
    for (unsigned int n_probe=0; n_probe<numProbes; n_probe++) {

	bool ok;

	// Extract "every" (number of timesteps between each output)
	unsigned int my_every=0;
	ok=PyTools::extract("every",my_every,"DiagProbe",n_probe);
	if (!ok) my_every=params.global_every;
	every.push_back(my_every);

	// Extract "time_range" (tmin and tmax of the outputs)
	vector<double> time_range(2,0.);
	double my_tmin,my_tmax;
	ok=PyTools::extract("time_range",time_range,"DiagProbe",n_probe);
	if (!ok) {
	    my_tmin = 0.;
	    my_tmax = params.sim_time;
	} else {
	    my_tmin = time_range[0];
	    my_tmax = time_range[1];
	}
	tmin.push_back(my_tmin);
	tmax.push_back(my_tmax);

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


#ifndef _NEW_STYLE

            // Dimension of the simulation
            //unsigned int ndim=params.nDim_particle;
            
            // Extract "pos", "pos_first", "pos_second" and "pos_third"
            // (positions of the vertices of the grid)
            vector< vector<double> > allPos;
            vector<double> pos;
            
            if (PyTools::extract("pos",pos,"DiagProbe",n_probe)) {
                if (pos.size()!=ndim) {
                    ERROR("Probe #"<<n_probe<<": pos size(" << pos.size() << ") != ndim(" << ndim<< ")");
                }
                allPos.push_back(pos);
            }
            
            if (PyTools::extract("pos_first",pos,"DiagProbe",n_probe)) {
                if (pos.size()!=ndim) {
                    ERROR("Probe #"<<n_probe<<": pos_first size(" << pos.size() << ") != ndim(" << ndim<< ")");
                }
                allPos.push_back(pos);
            }
            
            if (PyTools::extract("pos_second",pos,"DiagProbe",n_probe)) {
                if (pos.size()!=ndim) {
                    ERROR("Probe #"<<n_probe<<": pos_second size(" << pos.size() << ") != ndim(" << ndim<< ")");
                }
                allPos.push_back(pos);
            }
            
            if (PyTools::extract("pos_third",pos,"DiagProbe",n_probe)) {
                if (pos.size()!=ndim) {
                    ERROR("Probe #"<<n_probe<<": pos_third size(" << pos.size() << ") != ndim(" << ndim<< ")");
                }
                allPos.push_back(pos);
            }
            
            // Extract the list of requested fields
            vector<string> fs;
            if(!PyTools::extract("fields",fs,"DiagProbe",n_probe)) {
                fs.resize(10);
                fs[0]="Ex"; fs[1]="Ey"; fs[2]="Ez";
                fs[3]="Bx"; fs[4]="By"; fs[5]="Bz";
                fs[6]="Jx"; fs[7]="Jy"; fs[8]="Jz"; fs[9]="Rho";
            }
            vector<unsigned int> locations;
            locations.resize(10);
            for( unsigned int i=0; i<10; i++) locations[i] = fs.size();
            for( unsigned int i=0; i<fs.size(); i++) {
                for( unsigned int j=0; j<i; j++) {
                    if( fs[i]==fs[j] ) {
                        ERROR("Probe #"<<n_probe<<": field "<<fs[i]<<" appears twice");
                    }
                }
                if     ( fs[i]=="Ex" ) locations[0] = i;
                else if( fs[i]=="Ey" ) locations[1] = i;
                else if( fs[i]=="Ez" ) locations[2] = i;
                else if( fs[i]=="Bx" ) locations[3] = i;
                else if( fs[i]=="By" ) locations[4] = i;
                else if( fs[i]=="Bz" ) locations[5] = i;
                else if( fs[i]=="Jx" ) locations[6] = i;
                else if( fs[i]=="Jy" ) locations[7] = i;
                else if( fs[i]=="Jz" ) locations[8] = i;
                else if( fs[i]=="Rho") locations[9] = i;
                else {
                    ERROR("Probe #"<<n_probe<<": unknown field "<<fs[i]);
                }
            }
            fieldlocation.push_back(locations);
            fieldname.push_back(fs);
            nFields.push_back(fs.size());
            
            // Calculate the total number of points in the grid
            // Each point is actually a "fake" macro-particle
            unsigned int my_nPart=1;
            for (unsigned int iDimProbe=0; iDimProbe<dimProbe; iDimProbe++) {
                my_nPart *= vecNumber[iDimProbe];
            }
            nPart_total.push_back(my_nPart);
            
            
            // Initialize the list of "fake" particles just as actual macro-particles
            Particles my_parts;
            my_parts.initialize(my_nPart, params.nDim_particle);
            
            // For each grid point, calculate its position and assign that position to the particle
            // The particle position is a linear combination of the `pos` with `pos_first` or `pos_second`, etc.
            double partPos, dx;
            vector<unsigned int> ipartND (dimProbe);
            for(unsigned int ipart=0; ipart<my_nPart; ++ipart) { // for each particle
                // first, convert the index `ipart` into N-D indexes
                unsigned int i = ipart;
                for (unsigned int iDimProbe=0; iDimProbe<dimProbe; iDimProbe++) {
                    ipartND[iDimProbe] = i%vecNumber[iDimProbe];
                    i = i/vecNumber[iDimProbe]; // integer division
                }
                // Now assign the position of the particle
                for(unsigned int iDim=0; iDim!=ndim; ++iDim) { // for each dimension of the simulation
                    partPos = allPos[0][iDim]; // position of `pos`
                    for (unsigned int iDimProbe=0; iDimProbe<dimProbe; iDimProbe++) { // for each of `pos`, `pos_first`, etc.
                        dx = (allPos[iDimProbe+1][iDim]-allPos[0][iDim])/(vecNumber[iDimProbe]-1); // distance between 2 gridpoints
                        partPos += ipartND[iDimProbe] * dx;
                    }
                    my_parts.position(iDim,ipart) = partPos;
                }
            }
            
            
            // Remove particles out of the domain
            for ( int ipb=my_nPart-1 ; ipb>=0 ; ipb--) {
                if (!my_parts.is_part_in_domain(ipb, patch))
                    my_parts.erase_particle(ipb);
            }
            probeParticles.push_back(my_parts);
            
            unsigned int nPart_local = my_parts.size(); // number of fake particles for this proc
            
            // Make the array that will contain the data
            // probesArray : 10 x nPart_tot
            vector<unsigned int> probesArraySize(2);
            probesArraySize[1] = nPart_local; // number of particles
            probesArraySize[0] = nFields[n_probe] + 1; // number of fields (Ex, Ey, etc) +1 for garbage
            Field2D *myfield = new Field2D(probesArraySize);
            probesArray.push_back(myfield);
#ifdef _TOMOVE            
            // Exchange data between MPI cpus so that they can figure out which part
            // of the grid they have to manage
            MPI_Status status;
            // Receive the location where to start from the previous node
            int my_Start = 0;
            if (patch->getRank()>0) MPI_Recv( &(my_Start), 1, MPI_INTEGER, patch->getRank()-1, 0, MPI_COMM_WORLD, &status );
            // Send the location where to end to the next node
            int probeEnd = my_Start+nPart_local;
            if (patch->getRank()!=patch->getSize()-1) MPI_Send( &probeEnd, 1, MPI_INTEGER, patch->getRank()+1, 0, MPI_COMM_WORLD );
            
            // Create group for the current probe
            ostringstream prob_name("");
            prob_name << "p" << setfill('0') << setw(4) << n_probe;
            
            // Create an array to hold the positions of local probe particles
            Field2D fieldPosProbe;
            fieldPosProbe.allocateDims(ndim,nPart_local);
            
            for (unsigned int ipb=0 ; ipb<nPart_local ; ipb++)
                for (unsigned int idim=0 ; idim<ndim  ; idim++)
                    fieldPosProbe(idim,ipb) = my_parts.position(idim,ipb);

#endif
#else




	probeParticles[n_probe].initialize(totPart, params.nDim_particle);

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

#endif



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
        
    // Write the version of the code as an attribute
    H5::attr(fileId, "Version", string(__VERSION));
    H5::attr(fileId, "CommitDate", string(__COMMITDATE));

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


	// Extract "pos", "pos_first", "pos_second" and "pos_third"
	// (positions of the vertices of the grid)
	vector< vector<double> > allPos;
	vector<double> pos;

	if (PyTools::extract("pos",pos,"DiagProbe",probe_id)) {
	    allPos.push_back(pos);
	}
            
	if (PyTools::extract("pos_first",pos,"DiagProbe",probe_id)) {
	    allPos.push_back(pos);
	}
            
	if (PyTools::extract("pos_second",pos,"DiagProbe",probe_id)) {
	    allPos.push_back(pos);
	}
            
	if (PyTools::extract("pos_third",pos,"DiagProbe",probe_id)) {
	    allPos.push_back(pos);
	}

	// Add arrays "p0", "p1", ... to the current group
	ostringstream pk;
	for (unsigned int iDimProbe=0; iDimProbe<=dimProbe; iDimProbe++) {
	    pk.str("");
	    pk << "p" << iDimProbe;
	    H5::vect(group_id, pk.str(), allPos[iDimProbe]);
	}
            
	// Add array "number" to the current group
	H5::vect(group_id, "number", vecNumber);


	vector<string> fs;
	if(!PyTools::extract("fields",fs,"DiagProbe",probe_id)) {
	    fs.resize(10);
	    fs[0]="Ex"; fs[1]="Ey"; fs[2]="Ez";
	    fs[3]="Bx"; fs[4]="By"; fs[5]="Bz";
	    fs[6]="Jx"; fs[7]="Jy"; fs[8]="Jz"; fs[9]="Rho";
	}
	vector<unsigned int> locations;
	locations.resize(10);
	for( unsigned int i=0; i<10; i++) locations[i] = fs.size();
	for( unsigned int i=0; i<fs.size(); i++) {
	    for( unsigned int j=0; j<i; j++) {
		if( fs[i]==fs[j] ) {
		    ERROR("Probe #"<<probe_id<<": field "<<fs[i]<<" appears twice");
		}
	    }
	    if     ( fs[i]=="Ex" ) locations[0] = i;
	    else if( fs[i]=="Ey" ) locations[1] = i;
	    else if( fs[i]=="Ez" ) locations[2] = i;
	    else if( fs[i]=="Bx" ) locations[3] = i;
	    else if( fs[i]=="By" ) locations[4] = i;
	    else if( fs[i]=="Bz" ) locations[5] = i;
	    else if( fs[i]=="Jx" ) locations[6] = i;
	    else if( fs[i]=="Jy" ) locations[7] = i;
	    else if( fs[i]=="Jz" ) locations[8] = i;
	    else if( fs[i]=="Rho") locations[9] = i;
	    else {
		ERROR("Probe #"<<probe_id<<": unknown field "<<fs[i]);
	    }
	}

	// Add "fields" to the current group
	ostringstream fields("");
	fields << fs[0];
	for( unsigned int i=1; i<fs.size(); i++) fields << "," << fs[i];
	H5::attr(group_id, "fields", fields.str());

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
    rsend.resize(numProbes);
    rrecv.resize(numProbes);

    for (unsigned int np=0; np<numProbes; np++) {
	probesStart[np] = 0;
	MPI_Status status;
	stringstream rtag("");
	rtag << cpuRank-1 << "0" << cpuRank;
	int tag(0);
	rtag >> tag;

	if (cpuRank>0) {
	    //cout << patch->Hindex() << " Recv from " << patch->getMPIRank(cpuRank-1) << " with tag " << tag << endl;
	    //MPI_Recv( &(probesStart[np]), 1, MPI_INTEGER, patch->getMPIRank(cpuRank-1), tag, MPI_COMM_WORLD, &status );
	    MPI_Irecv( &(probesStart[np]), 1, MPI_INTEGER, patch->getMPIRank(cpuRank-1), tag, MPI_COMM_WORLD, &(rrecv[np]) );
	}
	    
	int probeEnd = probesStart[np]+probeParticles[np].size();
	stringstream stag("");
	stag << cpuRank << "0" << cpuRank+1;
	tag = 0;
	stag >> tag;
	if (cpuRank!=nPatches-1) {
	    //cout << patch->Hindex() << " Send to " << patch->getMPIRank(cpuRank+1) << " with tag " << tag << endl;
	    //MPI_Send( &probeEnd, 1, MPI_INTEGER, patch->getMPIRank(cpuRank+1), tag, MPI_COMM_WORLD );
	    MPI_Isend( &probeEnd, 1, MPI_INTEGER, patch->getMPIRank(cpuRank+1), tag, MPI_COMM_WORLD, &(rsend[np]) );

	}
    }

    // ---------------------------------------------------
    // End file split definition
    // ---------------------------------------------------

}

void DiagnosticProbe::waitSetFile(Params& params)
{
    MPI_Status rstat;
    MPI_Status sstat;

    for (unsigned int np=0; np<rsend.size(); np++) {
	if (cpuRank>0) 
	    MPI_Wait( &(rrecv[np]), &rstat );

    int nPatches(1);
    for (int iDim=0;iDim<params.nDim_field;iDim++) nPatches*=params.number_of_patches[iDim];
    if (cpuRank!=nPatches-1)
	MPI_Wait( &(rsend[np]), &sstat );
    }
}

void DiagnosticProbe::setFile(hid_t masterFileId)
{
    fileId = masterFileId;  
}

void DiagnosticProbe::writePositionIn( Params &params )
{
    // Loop on all Probes
    for (unsigned int probe_id=0; probe_id<every.size(); probe_id++) {

#ifdef _NEWTYLE
	hid_t gid = H5::group(fileId, prob_name.str());
	// Add array "positions" into the current HDF5 group
	H5::matrix_MPI(gid, "positions", fieldPosProbe.data_2D[0][0], my_nPart, ndim, my_Start, nPart_local);
            
	probesStart.push_back(my_Start);
            
	// Add arrays "p0", "p1", ... to the current group
	ostringstream pk;
	for (unsigned int iDimProbe=0; iDimProbe<=dimProbe; iDimProbe++) {
	    pk.str("");
	    pk << "p" << iDimProbe;
	    H5::vect(gid, pk.str(), allPos[iDimProbe]);
	}
            
	// Add array "number" to the current group
	H5::vect(gid, "number", vecNumber);
            
	// Add attribute every to the current group
	H5::attr(gid, "every", my_every);
	// Add attribute "dimension" to the current group
	H5::attr(gid, "dimension", dim);
            
	// Add "fields" to the current group
	ostringstream fields("");
	fields << fs[0];
	for( unsigned int i=1; i<fs.size(); i++) fields << "," << fs[i];
	H5::attr(gid, "fields", fields.str());
            
	// Close current group
	H5Gclose(gid);
#endif

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
    chunk_parts[1] = probeParticles[probe_id].size();
    chunk_parts[0] = ndim_Particles; 
    hid_t memspace  = H5Screate_simple(2, chunk_parts, NULL);
    // filespace :
    hsize_t dimsf[2], offset[2], stride[2], count[2];
    dimsf[1] = nPart_total[probe_id];
    dimsf[0] = ndim_Particles;
    hid_t filespace = H5Screate_simple(2, dimsf, NULL);
    offset[1] = probesStart[probe_id];
    offset[0] = 0;
    stride[0] = 1;
    stride[1] = 1;
    count[0] = 1;
    count[1] = 1;
    hsize_t     block[2];
    block[1] = probeParticles[probe_id].size();
    block[0] = ndim_Particles;
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
	probesArray[probe_id]->data_2D[fieldlocation[probe_id][0]][iprob]=Eloc_fields.x;
	probesArray[probe_id]->data_2D[fieldlocation[probe_id][1]][iprob]=Eloc_fields.y;
	probesArray[probe_id]->data_2D[fieldlocation[probe_id][2]][iprob]=Eloc_fields.z;
	probesArray[probe_id]->data_2D[fieldlocation[probe_id][3]][iprob]=Bloc_fields.x;
	probesArray[probe_id]->data_2D[fieldlocation[probe_id][4]][iprob]=Bloc_fields.y;
	probesArray[probe_id]->data_2D[fieldlocation[probe_id][5]][iprob]=Bloc_fields.z;
	probesArray[probe_id]->data_2D[fieldlocation[probe_id][6]][iprob]=Jloc_fields.x;
	probesArray[probe_id]->data_2D[fieldlocation[probe_id][7]][iprob]=Jloc_fields.y;
	probesArray[probe_id]->data_2D[fieldlocation[probe_id][8]][iprob]=Jloc_fields.z;          
	probesArray[probe_id]->data_2D[fieldlocation[probe_id][9]][iprob]=Rloc_fields;
        
    } // End for iprob

}

void DiagnosticProbe::write(int probe_id, unsigned int timestep, hid_t group_id) {
    // memspace OK : 1 block 
    hsize_t     chunk_parts[2];
    chunk_parts[1] = probeParticles[probe_id].size();
    chunk_parts[0] = probeSize; 
    hid_t memspace  = H5Screate_simple(2, chunk_parts, NULL);
    // filespace :
    hsize_t dimsf[2], offset[2], stride[2], count[2];
    dimsf[1] = nPart_total[probe_id];
    dimsf[0] = probeSize;
    hid_t filespace = H5Screate_simple(2, dimsf, NULL);
    //cout << " CPU Rank " << cpuRank << " - writing at "  << probesStart[probe_id] << endl;
    offset[1] = probesStart[probe_id];
    offset[0] = 0;
    stride[0] = 1;
    stride[1] = 1;
    count[0] = 1;
    count[1] = 1;
    hsize_t     block[2];
    block[1] = probeParticles[probe_id].size();
    block[0] = probeSize;
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

    //H5::matrix_MPI(dset_id, name_t.str(), probesArray[np]->data_2D[0][0], nPart_total[np], nFields[np], probesStart[np], nPart_local);
    H5Pclose(plist_id);
    H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, write_plist, &(probesArray[probe_id]->data_2D[0][0]) );
    H5Dclose(dset_id);
    H5Gclose(did);
    H5Pclose( write_plist );

    H5Sclose(filespace);
    H5Sclose(memspace);

}
