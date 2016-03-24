
#include <string>
#include <sstream>

#include "DiagnosticTrack.h"
#include "VectorPatch.h"

using namespace std;

DiagnosticTrack::DiagnosticTrack( Params &params, SmileiMPI* smpi, Patch* patch, int diagId ) :
Diagnostic(params,smpi,patch,diagId),
nDim_particle(params.nDim_particle)
{
    probeId_ = diagId;

    // define speciesId_
    speciesId_ = 0;

    // Count track diags instanciated
    int nTrack(0);
    for ( int idiag = 0 ; idiag < patch->localDiags.size() ; idiag++ ) {
	if ( patch->localDiags[idiag]->type_ == "Track" )
	    nTrack++;
    }

    // Go to next track species
    int nTrackable(0); 
    for (int ispec = 0 ; ispec < patch->vecSpecies.size() ; ispec++ ) {
	if ( patch->vecSpecies[ispec]->particles->tracked ) {
	    nTrackable++;
	}
	if ( nTrackable > nTrack )
	    speciesId_ = ispec;
    }
    species = patch->vecSpecies[speciesId_];

    string nameSpec="Ntot_"+ patch->vecSpecies[speciesId_]->species_type;

    // Define the transfer type (collective is faster than independent)
    transfer = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio( transfer, H5FD_MPIO_INDEPENDENT);
    iter = 0;
    
    // Get the time selection from the 
    timeSelection = patch->vecSpecies[speciesId_]->particles->track_timeSelection;

    type_ = "Track";

}


DiagnosticTrack::~DiagnosticTrack()
{
}


void DiagnosticTrack::openFile( Params& params, SmileiMPI* smpi, VectorPatch& vecPatches, bool newfile )
{
    
    
    if ( newfile ) {
	// Create HDF5 file
	ostringstream hdf_filename("");
	hdf_filename << "TrackParticles_" << species->species_type  << ".h5" ;
	filename = hdf_filename.str();

	hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
	fileId_ = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, pid);
	H5Pclose(pid);

	// Define maximum size
	//hsize_t maxDimsPart[2] = {H5S_UNLIMITED, (hsize_t)nParticles};

	dims[0] = timeSelection->numberOfEvents(0, params.n_time);
	if ( !nbrParticles_ ) {
	    ERROR("DiagTrack empty or number of Particles in diag is null");
	}
	else
	    dims[1] = nbrParticles_;

	hid_t file_space = H5Screate_simple(2, dims, NULL);
        
	// Create the overall dataset in the HDF5 file with a new chunk every timestep
	hid_t plist = H5Pcreate(H5P_DATASET_CREATE);

	// Create the datasets for x, y and z
	hid_t did;
	unsigned int nPosition = species->particles->Position.size();
	for (unsigned int i=0; i<nPosition; i++) {
	    ostringstream namePos("");
	    namePos << "Position-" << i;
	    did = H5Dcreate(fileId_, namePos.str().c_str(), H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
	    H5Dclose(did);
	}
        
	// Create the datasets for px, py and pz
	unsigned int nMomentum = species->particles->Momentum.size();
	for (unsigned int i=0; i<nMomentum; i++) {
	    ostringstream nameMom("");
	    nameMom << "Momentum-" << i;
	    did = H5Dcreate(fileId_, nameMom.str().c_str(), H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
	    H5Dclose(did);
	}
        
	// Create the datasets for weight, charge and ID
	did = H5Dcreate(fileId_, "Weight", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
	H5Dclose(did);
	did = H5Dcreate(fileId_, "Charge", H5T_NATIVE_SHORT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
	H5Dclose(did);
	did = H5Dcreate(fileId_, "Id", H5T_NATIVE_UINT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
	H5Dclose(did);
        
	H5Pclose(plist);
	H5Sclose(file_space);
    
	// Create the dataset for the time array
	plist = H5Pcreate(H5P_DATASET_CREATE);
	hsize_t ntimes[1] = { dims[0] };
	file_space = H5Screate_simple(1, ntimes, NULL);
	did = H5Dcreate(fileId_, "Times", H5T_NATIVE_INT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
	H5Dclose(did);
	H5Pclose(plist);
	H5Sclose(file_space);
    

    }
    else {
	hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
	fileId_ = H5Fopen( filename.c_str(), H5F_ACC_RDWR, pid);
	H5Pclose(pid);
   }

}


void DiagnosticTrack::setFile( hid_t fid )
{
    fileId_ = fid;
}

void DiagnosticTrack::setFile( Diagnostic* diag )
{
    fileId_ = static_cast<DiagnosticTrack*>(diag)->fileId_;  
}

void DiagnosticTrack::closeFile()
{
    H5Fclose( fileId_ );
}


void DiagnosticTrack::prepare( Patch* patch, int timestep )
{
}


void DiagnosticTrack::run( Patch* patch, int timestep )
{
    int time = timestep;

    if( ! timeSelection->theTimeIsNow(time) ) return;
    
    iter ++;

    int locNbrParticles = species->getNbrOfParticles();
    
    // Create the locator (gives locations where to store particles in the file)
    vector<hsize_t> locator (locNbrParticles*2);
    for(int i=0; i<locNbrParticles; i++) {
        locator[i*2  ] = iter-1;
        locator[i*2+1] = species->particles->id(i)-1; // because particles label Id starts at 1
    }
    // Now we increase the array size for the new timestep (new chunk)
    // It is not applied to the HDF5 file yet
    dims[0] ++;

    
    // Specify the memory dataspace (the size of the local array)
    hsize_t count[2] = {1, (hsize_t)locNbrParticles};
    hid_t mem_space = H5Screate_simple(2, count, NULL);

    // For each dataspace (x, y, z, px, py, pz, weight, charge and ID), add the local
    // array to the HDF5 file -> see function appendTestParticles() in SmileiIO.h
    ostringstream namePos("");
    for (int idim=0 ; idim<nDim_particle ; idim++) {
        namePos.str("");
        namePos << "Position-" << idim;
        append( fileId_, namePos.str(), species->particles->position(idim)[0], mem_space, locNbrParticles, H5T_NATIVE_DOUBLE, locator );
    }
    ostringstream nameMom("");
    for (int idim=0 ; idim<3 ; idim++) {
        nameMom.str("");
        nameMom << "Momentum-" << idim;
        append( fileId_, nameMom.str(), species->particles->momentum(idim)[0], mem_space, locNbrParticles, H5T_NATIVE_DOUBLE, locator );
    }
    append( fileId_, "Weight", species->particles->weight()[0], mem_space, locNbrParticles, H5T_NATIVE_DOUBLE, locator );
    append( fileId_, "Charge", species->particles->charge()[0], mem_space, locNbrParticles, H5T_NATIVE_SHORT , locator );
    append( fileId_, "Id"    , species->particles->id()    [0], mem_space, locNbrParticles, H5T_NATIVE_UINT  , locator );
    
    H5Sclose( mem_space );
    
    // Write the current timestep
    locator.resize(1, iter-1);
    hsize_t onetime[1] = { 1 };
    mem_space = H5Screate_simple(1, onetime, NULL);
    append( fileId_, "Times", time, mem_space, 1, H5T_NATIVE_INT, locator );
    
    //MPI_Barrier(MPI_COMM_WORLD); // synchro to manage differently
    H5Fflush( fileId_, H5F_SCOPE_GLOBAL );

}


void DiagnosticTrack::write(int timestep)
{
}

template <class T>
void DiagnosticTrack::append( hid_t fid, string name, T & property,  hid_t  mem_space, int nParticles, hid_t type, vector<hsize_t> &locator) {
    
    // Open existing dataset
    hid_t did = H5Dopen( fid, name.c_str(), H5P_DEFAULT );
    // Increase the size of the array with the previously defined size

    //H5Dset_extent(did, dims);

    // Get the extended file space
    hid_t file_space = H5Dget_space(did);
    
    // Select locations that this proc will write
    if(nParticles>0)
        H5Sselect_elements( file_space, H5S_SELECT_SET, nParticles, &locator[0] );
    else
        H5Sselect_none(file_space);
    
    // Write
    H5Dwrite( did, type, mem_space , file_space , transfer, &property );
    
    H5Sclose(file_space);
    H5Dclose(did);
}


void DiagnosticTrack::setFileSplitting( Params& params, SmileiMPI* smpi, VectorPatch& vecPatches )
{
    // Communicate some stuff if this is a species that has to be dumped (particles have Id)
    // Need to be placed after ALL createParticles()
    if (vecPatches(0)->vecSpecies[speciesId_]->particles->tracked) {

	// Internal patches offset

	std::vector<int> localNbrParticles( vecPatches.size(), 0 );
	localNbrParticles[0] = vecPatches(0)->vecSpecies[speciesId_]->getNbrOfParticles();
	for (unsigned int ipatch=1 ; ipatch<vecPatches.size() ; ipatch++) {
	    // number of particles up to ipatch (including)
	    localNbrParticles[ipatch] += vecPatches(ipatch)->vecSpecies[speciesId_]->getNbrOfParticles() + localNbrParticles[ipatch-1];
	    vecPatches(ipatch)->vecSpecies[speciesId_]->particles->addIdOffsets(localNbrParticles[ipatch-1]);
	}
	int locNbrParticles = localNbrParticles[vecPatches.size()-1];


	// MPI offset

	//int locNbrParticles = thisSpecies->getNbrOfParticles();
	int sz(1);
	MPI_Comm_size( MPI_COMM_WORLD, &sz );
	std::vector<int> allNbrParticles(sz);
	MPI_Allgather( &locNbrParticles, 1, MPI_INTEGER, &allNbrParticles[0], 1, MPI_INTEGER, MPI_COMM_WORLD );

	int totNbrParts(0);
	for (int irk=0 ; irk<sz ; irk++) totNbrParts += allNbrParticles[irk];
	// HDF5 file open by all patch master
	this->setGlobalNbrParticles(totNbrParts);
	vecPatches(0)->localDiags[probeId_]->openFile( params, smpi, vecPatches, true );

	// Set HDF5 context for other patches
	for (unsigned int ipatch=1 ; ipatch<vecPatches.size() ; ipatch++) {
	    DiagnosticTrack* diag = static_cast<DiagnosticTrack*>( vecPatches(ipatch)->localDiags[probeId_] );
	    diag->setGlobalNbrParticles(totNbrParts);
	}

	int nParticles(0);

	nParticles =  allNbrParticles[0];
	for (int irk=1 ; irk<sz ; irk++){
	    allNbrParticles[irk] += nParticles;
	    nParticles = allNbrParticles[irk];
	}
	for (int irk=sz-1 ; irk>0 ; irk--){
	    allNbrParticles[irk] = allNbrParticles[irk-1];
	}
	allNbrParticles[0] = 0;

	int offset(0);
	MPI_Scatter(&allNbrParticles[0], 1 , MPI_INTEGER, &offset, 1, MPI_INTEGER, 0, MPI_COMM_WORLD );
            
	for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
	    vecPatches(ipatch)->vecSpecies[speciesId_]->particles->addIdOffsets(offset);

    } // End if tracked


} // End initTrackParticles

