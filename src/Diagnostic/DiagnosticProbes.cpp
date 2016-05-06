#include <sstream>
#include <vector>

#include "DiagnosticProbes.h"

#include "VectorPatch.h"


using namespace std;

DiagnosticProbes::DiagnosticProbes( Params &params, SmileiMPI* smpi, Patch* patch, int diagId, int n_probe )
{
    diagId_ = diagId;
    probesStart = 0;
    nDim_particle = params.nDim_particle;
    
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
    
    unsigned int totPart=1;
    for (unsigned int iDimProbe=0; iDimProbe<dimProbe; iDimProbe++) {
        totPart *= vecNumber[iDimProbe];
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
    
    initParticles(params, patch);
    
    interp_ = InterpolatorFactory::create(params, patch);
    
    ostringstream mystream("");
    mystream << "Probes" << n_probe << ".h5";
    filename = mystream.str();
    
    type_ = "Probes";

} // END DiagnosticProbes::DiagnosticProbes


// Cloning constructor
DiagnosticProbes::DiagnosticProbes( DiagnosticProbes * probe, Params& params, Patch* patch )
{
    nDim_particle  = params.nDim_particle;
    diagId_        = probe->diagId_;
    fieldlocation  = probe->fieldlocation;
    fieldname      = probe->fieldname    ;
    nFields        = probe->nFields      ;
    filename       = probe->filename;
    vecNumber      = probe->vecNumber;
    allPos         = probe->allPos;
    dimProbe       = probe->dimProbe;
    
    timeSelection = new TimeSelection(probe->timeSelection);
    
    initParticles(params, patch);
    
    interp_ = InterpolatorFactory::create(params, patch);
    
    type_ = "Probes";
}



DiagnosticProbes::~DiagnosticProbes()
{
    delete interp_;
}


// During constructor, we call this method to create the probe particles
// at the right positions
void DiagnosticProbes::initParticles(Params& params, Patch * patch)
{
    // Calculate the total number of points in the grid
    // Each point is actually a "fake" macro-particle
    unsigned int my_nPart=1;
    for (unsigned int iDimProbe=0; iDimProbe<dimProbe; iDimProbe++) {
        my_nPart *= vecNumber[iDimProbe];
    }
    nPart_total = my_nPart;
    
    // Initialize the list of "fake" particles just as actual macro-particles
    probeParticles.initialize(my_nPart, nDim_particle);
    
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
        for(unsigned int iDim=0; iDim!=nDim_particle; ++iDim) { // for each dimension of the simulation
            partPos = allPos[0][iDim]; // position of `pos`
            for (unsigned int iDimProbe=0; iDimProbe<dimProbe; iDimProbe++) { // for each of `pos`, `pos_first`, etc.
                dx = (allPos[iDimProbe+1][iDim]-allPos[0][iDim])/(vecNumber[iDimProbe]-1); // distance between 2 gridpoints
                partPos += ipartND[iDimProbe] * dx;
            }
            probeParticles.position(iDim,ipart) = partPos;
        }
    }
    
    // Remove particles out of the domain
    for ( int ipb=my_nPart-1 ; ipb>=0 ; ipb--) {
        if (!probeParticles.is_part_in_domain(ipb, patch))
            probeParticles.erase_particle(ipb);
    }
    unsigned int nPart_local = probeParticles.size(); // number of fake particles for this proc
    
    // Make the array that will contain the data
    // probesArray : 10 x nPart_tot
    vector<unsigned int> probesArraySize(2);
    probesArraySize[1] = nPart_local; // number of particles
    probesArraySize[0] = nFields + 1; // number of fields (Ex, Ey, etc) +1 for garbage
    Field2D *myfield = new Field2D(probesArraySize);
    probesArray = myfield;
}


void DiagnosticProbes::openFile( Params& params, SmileiMPI* smpi, bool newfile )
{
    if ( newfile ) {
        hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
        
        
        
        fileId_ = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, pid);
        
        // Write the version of the code as an attribute
        H5::attr(fileId_, "Version", string(__VERSION));
        H5::attr(fileId_, "CommitDate", string(__COMMITDATE));
        
        //else file is created ok 
        H5Pclose(pid);
        
        hid_t sid, aid;
        
        // Dimension of the probe grid
        sid = H5Screate(H5S_SCALAR);        
        aid = H5Acreate(fileId_, "dimension", H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT);
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
    H5Fclose(fileId_);
}


bool DiagnosticProbes::prepare( int timestep )
{
    return timeSelection->theTimeIsNow(timestep);
}


void DiagnosticProbes::run( Patch* patch, int timestep )
{
    // skip if current timestep is not requested
    compute( timestep, patch->EMfields );

}


void DiagnosticProbes::write( int timestep )
{
    // memspace OK : 1 block
    hsize_t     chunk_parts[2];
    chunk_parts[1] = probeParticles.size();
    chunk_parts[0] = nFields;
    hid_t memspace  = H5Screate_simple(2, chunk_parts, NULL);
    // filespace :
    hsize_t dimsf[2], offset[2], stride[2], count[2], block[2];
    dimsf[1] = nPart_total;
    dimsf[0] = nFields;
    hid_t filespace = H5Screate_simple(2, dimsf, NULL);
    offset[1] = probesStart;
    offset[0] = 0;
    stride[0] = 1;
    stride[1] = 1;
    count[0] = 1;
    count[1] = 1;
    block[1] = probeParticles.size();
    block[0] = nFields;
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);

    // define filespace, memspace
    hid_t write_plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_INDEPENDENT);
    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
    ostringstream name_t;
    name_t.str("");
    name_t << "/" << setfill('0') << setw(10) << timestep;

    hid_t dset_id;
    if (!H5Lexists( fileId_, name_t.str().c_str(), H5P_DEFAULT ))
        dset_id  = H5Dcreate(fileId_, name_t.str().c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    else
        dset_id = H5Dopen(fileId_, name_t.str().c_str(), H5P_DEFAULT);                

    //H5::matrix_MPI(dset_id, name_t.str(), probesArray->data_2D[0][0], nPart_total, nFields, probesStart, nPart_local);
    H5Pclose(plist_id);
    if (probeParticles.size())
        H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, write_plist, &(probesArray->data_2D[0][0]) );
    H5Dclose(dset_id);
    H5Pclose( write_plist );

    H5Sclose(filespace);
    H5Sclose(memspace);
    
    //if ( fileId_ ) H5Fflush( fileId_, H5F_SCOPE_GLOBAL );

}


void DiagnosticProbes::setFileSplitting( Params& params, SmileiMPI* smpi, VectorPatch& vecPatches )
{
    // 1 - Local (MPI) offset
    
    int localNbrParticles = 0;
    for (int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        DiagnosticProbes* diag = static_cast<DiagnosticProbes*>(vecPatches(ipatch)->localDiags[diagId_]);
        diag->probesStart = localNbrParticles;
        localNbrParticles += diag->probeParticles.size();
    }
    
    // 2 - Global (MPI) offset
    
    // Get the number of particles for each MPI
    int sz = smpi->getSize();
    std::vector<int> allNbrParticles(sz, 0);
    MPI_Allgather( &localNbrParticles, 1, MPI_INT, &allNbrParticles[0], 1, MPI_INT, MPI_COMM_WORLD );
    
    // Calculate the cumulative sum
    for (int irk=1 ; irk<sz ; irk++)
        allNbrParticles[irk] += allNbrParticles[irk-1];
    
    // Calculate the starting location in the file for each patch
    int offset = 0;
    if( ! smpi->isMaster() ) offset = allNbrParticles[smpi->getRank()-1];
    for (int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        DiagnosticProbes* diag = static_cast<DiagnosticProbes*>(vecPatches(ipatch)->localDiags[diagId_]);
        diag->probesStart += offset;
    }
}


void DiagnosticProbes::init()
{
    
    // Write the array of the particle positions
    
    vector<unsigned int> posArraySize(2);
    posArraySize[0] = probeParticles.size();
    posArraySize[1] = nDim_particle;
    Field2D* posArray = new Field2D(posArraySize);
    for ( int ipb=0 ; ipb<probeParticles.size() ; ipb++) {
        for (int idim=0 ; idim<nDim_particle  ; idim++ )
            posArray->data_2D[ipb][idim] = probeParticles.position(idim,ipb);
    }
    
    // memspace : 1 block 
    hsize_t     chunk_parts[2];
    chunk_parts[0] = probeParticles.size();
    chunk_parts[1] = nDim_particle; 
    hid_t memspace  = H5Screate_simple(2, chunk_parts, NULL);
    // filespace :
    hsize_t dimsf[2], offset[2], stride[2], count[2];
    dimsf[0] = nPart_total;
    dimsf[1] = nDim_particle;
    hid_t filespace = H5Screate_simple(2, dimsf, NULL);
    offset[0] = probesStart;
    offset[1] = 0;
    stride[0] = 1;
    stride[1] = 1;
    count[0] = 1;
    count[1] = 1;
    hsize_t     block[2];
    block[0] = probeParticles.size();
    block[1] = nDim_particle;
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);
    
    hid_t write_plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_INDEPENDENT);
    hid_t plist_id;
    hid_t dset_id;
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    if ( !H5Lexists( fileId_, "positions", H5P_DEFAULT ) )
        dset_id = H5Dcreate(fileId_, "positions", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    else
        dset_id = H5Dopen(fileId_, "positions", H5P_DEFAULT);
    
    H5Pclose(plist_id);
    if (probeParticles.size())
        H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, write_plist, &(posArray->data_2D[0][0]) );
    H5Dclose(dset_id);
    H5Pclose( write_plist );
    
    H5Sclose(filespace);
    H5Sclose(memspace);
    
    delete posArray;
}


void DiagnosticProbes::compute( unsigned int timestep, ElectroMagn* EMfields )
{
    
    // Loop probe ("fake") particles
    for (int iprob=0; iprob <probeParticles.size(); iprob++) {             
        (*interp_)(EMfields,probeParticles,iprob,&Eloc_fields,&Bloc_fields,&Jloc_fields,&Rloc_fields);

        //! here we fill the probe data!!!
        probesArray->data_2D[fieldlocation[0]][iprob]=Eloc_fields.x;
        probesArray->data_2D[fieldlocation[1]][iprob]=Eloc_fields.y;
        probesArray->data_2D[fieldlocation[2]][iprob]=Eloc_fields.z;
        probesArray->data_2D[fieldlocation[3]][iprob]=Bloc_fields.x;
        probesArray->data_2D[fieldlocation[4]][iprob]=Bloc_fields.y;
        probesArray->data_2D[fieldlocation[5]][iprob]=Bloc_fields.z;
        probesArray->data_2D[fieldlocation[6]][iprob]=Jloc_fields.x;
        probesArray->data_2D[fieldlocation[7]][iprob]=Jloc_fields.y;
        probesArray->data_2D[fieldlocation[8]][iprob]=Jloc_fields.z;          
        probesArray->data_2D[fieldlocation[9]][iprob]=Rloc_fields;
        
    } // End for iprob

}
