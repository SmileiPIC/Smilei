#include <sstream>
#include <vector>

#include "DiagnosticProbes.h"

#include "VectorPatch.h"


using namespace std;

DiagnosticProbes::DiagnosticProbes( Params &params, SmileiMPI* smpi, int n_probe )
{
    probe_n = n_probe;
    nDim_particle = params.nDim_particle;
    fileId_ = 0;
    
    // Extract "every" (time selection)
    ostringstream name("");
    name << "Probe #"<<n_probe;
    timeSelection = new TimeSelection( 
        PyTools::extract_py("every", "DiagProbe", n_probe),
        name.str()
    );
    
    // Extract "number" (number of points you have in each dimension of the probe,
    // which must be smaller than the code dimensions)
    PyTools::extract("number",vecNumber,"DiagProbe",n_probe);
    
    // Dimension of the probe grid
    dimProbe=vecNumber.size();
    if (dimProbe > nDim_particle) {
        ERROR("Probe #"<<n_probe<<": probe dimension is greater than simulation dimension")
    }
    
    // If there is no "number" argument provided, then it corresponds to
    // a zero-dimensional probe (one point). In this case, we say the probe
    // has actually one dimension with only one point.
    if (dimProbe == 0) {
        vecNumber.resize(1);
        vecNumber[0]=1;
    }
    
    // Calculate the total number of points in the grid
    nPart_total=1;
    for (unsigned int iDimProbe=0; iDimProbe<dimProbe; iDimProbe++) {
        nPart_total *= vecNumber[iDimProbe];
    }
    
    // -----------------------------------------------------
    // Start definition of probeParticles (probes positions)
    // -----------------------------------------------------            
    
    // Extract "pos", "pos_first", "pos_second" and "pos_third"
    // (positions of the vertices of the grid)
    allPos.resize(0);
    vector<string> keys(4);
    keys[0] = "pos";
    keys[1] = "pos_first";
    keys[2] = "pos_second";
    keys[3] = "pos_third";
    for( int i=0; i<3; i++) {
        vector<double> pos;
        if (PyTools::extract(keys[i],pos,"DiagProbe",n_probe)) {
            if (pos.size()!=nDim_particle)
                ERROR("Probe #"<<n_probe<<": "<<keys[i]<<" size(" << pos.size() << ") != ndim(" << nDim_particle<< ")");
            allPos.push_back(pos);
        }
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
    fieldlocation = locations;
    fieldname = fs;
    nFields = fs.size();
    
    ostringstream mystream("");
    mystream << "Probes" << n_probe << ".h5";
    filename = mystream.str();
    
    MESSAGE(1, "Diagnostic created: probe #"<<n_probe);
    
    type_ = "Probes";

} // END DiagnosticProbes::DiagnosticProbes



DiagnosticProbes::~DiagnosticProbes()
{
    delete timeSelection;
}


void DiagnosticProbes::openFile( Params& params, SmileiMPI* smpi, bool newfile )
{
    if ( newfile ) {
        // Create file
        hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
        fileId_ = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, pid);
        H5Pclose(pid);
        
        // Write the version of the code as an attribute
        H5::attr(fileId_, "Version", string(__VERSION));
        H5::attr(fileId_, "CommitDate", string(__COMMITDATE));
        
        // Dimension of the probe grid
        hid_t sid = H5Screate(H5S_SCALAR);        
        hid_t aid = H5Acreate(fileId_, "dimension", H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(aid, H5T_NATIVE_UINT, &dimProbe);
        H5Aclose(aid);
        H5Sclose(sid);
        
        // Add arrays "p0", "p1", ...
        ostringstream pk;
        for (unsigned int iDimProbe=0; iDimProbe<=dimProbe; iDimProbe++) {
            pk.str("");
            pk << "p" << iDimProbe;
            H5::vect(fileId_, pk.str(), allPos[iDimProbe]);
        }
        
        // Add array "number"
        H5::vect(fileId_, "number", vecNumber);
        
        // Add "fields"
        ostringstream fields("");
        fields << fieldname[0];
        for( unsigned int i=1; i<fieldname.size(); i++) fields << "," << fieldname[i];
        H5::attr(fileId_, "fields", fields.str());

    }
    else {
        hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
        fileId_ = H5Fopen( filename.c_str(), H5F_ACC_RDWR, pid );
        H5Pclose(pid);
    }

}


void DiagnosticProbes::closeFile()
{
    if(fileId_!=0) H5Fclose(fileId_);
    fileId_ = 0;
}


bool DiagnosticProbes::prepare( int timestep )
{
    return timeSelection->theTimeIsNow(timestep);
}



void DiagnosticProbes::init(Params& params, SmileiMPI* smpi, VectorPatch& vecPatches)
{
    int nPart_MPI = 0;
    
    // 1 - Loop patches to create particles
    for (int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        
        // First loop on particles: calculate the number of points for this patch
        unsigned int nPart_patch=0, i, ipart, iDimProbe, iDim;
        double dx;
        vector<double> partPos(nDim_particle);
        vector<unsigned int> ipartND (dimProbe);
        bool is_in_domain;
        for(ipart=0; ipart<nPart_total; ipart++) { // for each particle
            // first, convert the index `ipart` into N-D indexes
            i = ipart;
            for (iDimProbe=0; iDimProbe<dimProbe; iDimProbe++) {
                ipartND[iDimProbe] = i%vecNumber[iDimProbe];
                i = i/vecNumber[iDimProbe]; // integer division
            }
            is_in_domain = true;
            // Now calculate the position of the particle
            for(iDim=0; iDim!=nDim_particle; iDim++) { // for each dimension of the simulation
                partPos[iDim] = allPos[0][iDim]; // position of `pos`
                for (iDimProbe=0; iDimProbe<dimProbe; iDimProbe++) { // for each of `pos`, `pos_first`, etc.
                    dx = (allPos[iDimProbe+1][iDim]-allPos[0][iDim])/(vecNumber[iDimProbe]-1); // distance between 2 gridpoints
                    partPos[iDim] += ipartND[iDimProbe] * dx;
                }
                // Stop if particle not in domain
                if (partPos[iDim] <  vecPatches(ipatch)->getDomainLocalMin(iDim) 
                 || partPos[iDim] >= vecPatches(ipatch)->getDomainLocalMax(iDim) ) {
                    is_in_domain = false;
                    break;
                }
            }
            if(is_in_domain) nPart_patch++;
        }
        
        // Initialize the list of "fake" particles just as actual macro-particles
        Particles * particles = &(vecPatches(ipatch)->probes[probe_n]->particles);
        particles->initialize(nPart_patch, nDim_particle);
        // Add the local offset
        vecPatches(ipatch)->probes[probe_n]->offset_in_file = nPart_MPI;
        nPart_MPI += nPart_patch;
        
        // Second loop on particles: assign the position of each particle
        // The particle position is a linear combination of the `pos` with `pos_first` or `pos_second`, etc.
        unsigned int ipart_local = 0;
        for(ipart=0; ipart<nPart_total; ipart++) { // for each particle
            if( ipart_local>= nPart_patch ) break;
            // first, convert the index `ipart` into N-D indexes
            i = ipart;
            for (iDimProbe=0; iDimProbe<dimProbe; iDimProbe++) {
                ipartND[iDimProbe] = i%vecNumber[iDimProbe];
                i = i/vecNumber[iDimProbe]; // integer division
            }
            is_in_domain = true;
            // Now calculate the position of the particle
            for(iDim=0; iDim!=nDim_particle; iDim++) { // for each dimension of the simulation
                partPos[iDim] = allPos[0][iDim]; // position of `pos`
                for (iDimProbe=0; iDimProbe<dimProbe; iDimProbe++) { // for each of `pos`, `pos_first`, etc.
                    dx = (allPos[iDimProbe+1][iDim]-allPos[0][iDim])/(vecNumber[iDimProbe]-1); // distance between 2 gridpoints
                    partPos[iDim] += ipartND[iDimProbe] * dx;
                }
                // Stop if particle not in domain
                if (partPos[iDim] <  vecPatches(ipatch)->getDomainLocalMin(iDim) 
                 || partPos[iDim] >= vecPatches(ipatch)->getDomainLocalMax(iDim) ) {
                    is_in_domain = false;
                    break;
                }
            }
            if(is_in_domain) {
                for(iDim=0; iDim!=nDim_particle; ++iDim)
                    particles->position(iDim,ipart_local) = partPos[iDim];
                ipart_local++;
            }
        }
    }
    
    // 2 - Calculate the offset of each patch to write in the file
    
    // Get the number of particles for each MPI
    int sz = smpi->getSize();
    std::vector<int> all_nPart(sz, 0);
    MPI_Allgather( &nPart_MPI, 1, MPI_INT, &all_nPart[0], 1, MPI_INT, MPI_COMM_WORLD );
    
    // Calculate the cumulative sum
    for (int irk=1 ; irk<sz ; irk++)
        all_nPart[irk] += all_nPart[irk-1];
    
    // Add the MPI offset to all patches
    if( ! smpi->isMaster() ) {
        int offset = all_nPart[smpi->getRank()-1];
        for (int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
            vecPatches(ipatch)->probes[probe_n]->offset_in_file += offset;
    }
    
    
    // 3 - Create file and write the array of the particle positions
    
    // create the file
    openFile( params, smpi, true );
    
    // Store the positions of all particles in this MPI
    vector<unsigned int> posArraySize(2);
    posArraySize[0] = nPart_MPI;
    posArraySize[1] = nDim_particle;
    Field2D* posArray = new Field2D(posArraySize);
    int ipart = 0;
    for (int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        if( ipart>=nPart_MPI ) break;
        Particles * particles = &(vecPatches(ipatch)->probes[probe_n]->particles);
        for ( int ip=0 ; ip<particles->size() ; ip++) {
            for (int idim=0 ; idim<nDim_particle  ; idim++ )
                posArray->data_2D[ipart][idim] = particles->position(idim,ip);
            ipart++;
        }
    }
    // Define size in memory
    hsize_t mem_size[2];
    mem_size[0] = nPart_MPI;
    mem_size[1] = nDim_particle; 
    hid_t memspace  = H5Screate_simple(2, mem_size, NULL);
    // Define size and location in file
    hsize_t dimsf[2], offset[2], stride[2], count[2], block[2];
    dimsf[0] = nPart_total;
    dimsf[1] = nDim_particle;
    hid_t filespace = H5Screate_simple(2, dimsf, NULL);
    if( nPart_MPI>0 ) {
        offset[0] = vecPatches(0)->probes[probe_n]->offset_in_file;
        offset[1] = 0;
        stride[0] = 1;
        stride[1] = 1;
        count[0] = 1;
        count[1] = 1;
        block[0] = nPart_MPI;
        block[1] = nDim_particle;
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);
    } else {
        H5Sselect_none(filespace);
    }
    // Define collective transfer 
    hid_t transfer = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(transfer, H5FD_MPIO_COLLECTIVE);
    // Create dataset
    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
    hid_t dset_id = H5Dcreate(fileId_, "positions", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    // Write
    H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, transfer, &(posArray->data_2D[0][0]) );
    H5Dclose(dset_id);
    H5Pclose( transfer );
    H5Sclose(filespace);
    H5Sclose(memspace);
    delete posArray;
    
    closeFile();
}



void DiagnosticProbes::run( SmileiMPI* smpi, VectorPatch& vecPatches, int timestep )
{
    
    // Leave if this timestep has already been written
    ostringstream name_t;
    name_t.str("");
    name_t << "/" << setfill('0') << setw(10) << timestep;
    if (H5Lexists( fileId_, name_t.str().c_str(), H5P_DEFAULT )) return;
    smpi->barrier();
    
    // Calculate the number of probe particles in this MPI
    int nPart_MPI = 0;
    for (int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
        nPart_MPI += vecPatches(ipatch)->probes[probe_n]->particles.size();
    
    // Make the array that will contain the data
    vector<unsigned int> probesArraySize(2);
    probesArraySize[1] = nPart_MPI; // number of particles
    probesArraySize[0] = nFields + 1; // number of fields (Ex, Ey, etc) +1 for garbage
    Field2D* probesArray = new Field2D(probesArraySize);
    
    // Loop patches to fill the array
    int iPart_MPI = 0;
    for (int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        // Loop probe ("fake") particles of current patch
        int npart = vecPatches(ipatch)->probes[probe_n]->particles.size();
        for (int ipart=0; ipart<npart; ipart++) {             
            (*(vecPatches(ipatch)->Interp)) (
                vecPatches(ipatch)->EMfields,
                vecPatches(ipatch)->probes[probe_n]->particles,
                ipart,
                &Eloc_fields, &Bloc_fields, &Jloc_fields, &Rloc_fields
            );
            
            //! here we fill the probe data!!!
            probesArray->data_2D[fieldlocation[0]][iPart_MPI]=Eloc_fields.x;
            probesArray->data_2D[fieldlocation[1]][iPart_MPI]=Eloc_fields.y;
            probesArray->data_2D[fieldlocation[2]][iPart_MPI]=Eloc_fields.z;
            probesArray->data_2D[fieldlocation[3]][iPart_MPI]=Bloc_fields.x;
            probesArray->data_2D[fieldlocation[4]][iPart_MPI]=Bloc_fields.y;
            probesArray->data_2D[fieldlocation[5]][iPart_MPI]=Bloc_fields.z;
            probesArray->data_2D[fieldlocation[6]][iPart_MPI]=Jloc_fields.x;
            probesArray->data_2D[fieldlocation[7]][iPart_MPI]=Jloc_fields.y;
            probesArray->data_2D[fieldlocation[8]][iPart_MPI]=Jloc_fields.z;          
            probesArray->data_2D[fieldlocation[9]][iPart_MPI]=Rloc_fields;
            iPart_MPI++;
        }
    }
    
    // Define size in memory
    hsize_t mem_size[2];
    mem_size[1] = nPart_MPI;
    mem_size[0] = nFields;
    hid_t memspace  = H5Screate_simple(2, mem_size, NULL);
    // Define size and location in file
    hsize_t dimsf[2], offset[2], stride[2], count[2], block[2];
    dimsf[1] = nPart_total;
    dimsf[0] = nFields;
    hid_t filespace = H5Screate_simple(2, dimsf, NULL);
    if( nPart_MPI>0 ) {
        offset[1] = vecPatches(0)->probes[probe_n]->offset_in_file;
        offset[0] = 0;
        stride[0] = 1;
        stride[1] = 1;
        count[0] = 1;
        count[1] = 1;
        block[1] = nPart_MPI;
        block[0] = nFields;
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);
    } else {
        H5Sselect_none(filespace);
    }
    // Create new dataset for this timestep
    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_alloc_time(plist_id, H5D_ALLOC_TIME_EARLY );
    hid_t dset_id  = H5Dcreate(fileId_, name_t.str().c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    // Define transfer
    hid_t transfer = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(transfer, H5FD_MPIO_COLLECTIVE);
    // Write
    H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, transfer, &(probesArray->data_2D[0][0]) );
    H5Dclose(dset_id);
    H5Pclose( transfer );
    H5Sclose(filespace);
    H5Sclose(memspace);
    
    delete probesArray;
    
    H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
}
