

#include <sstream>
#include <vector>
#include <limits>

#include "DiagnosticProbes.h"

#include "VectorPatch.h"


using namespace std;





// Calculates the inverse of a square matrix, given row by row
vector<double> matrixInverse(vector<double> A) {
    unsigned int size = A.size();
    unsigned int dim  = sqrt( size );
    
    vector<double> invA (size);
    
    if( dim == 1 ) {
        invA[0] = 1./A[0];
    } else if (dim == 2) {
        double idet = 1./(A[0]*A[3]-A[1]*A[2]);
        invA[0] =  A[3]*idet;
        invA[1] = -A[1]*idet;
        invA[2] = -A[2]*idet;
        invA[3] =  A[0]*idet;
    } else if (dim == 3) {
        double idet = 1./(A[0]*A[4]*A[8]+A[1]*A[5]*A[6]+A[2]*A[3]*A[7]-A[2]*A[4]*A[6]-A[1]*A[3]*A[8]-A[0]*A[5]*A[7]);
        invA[0] = ( A[4]*A[8]-A[5]*A[7] )*idet;
        invA[1] = ( A[2]*A[7]-A[1]*A[8] )*idet;
        invA[2] = ( A[1]*A[5]-A[2]*A[4] )*idet;
        invA[3] = ( A[5]*A[6]-A[3]*A[8] )*idet;
        invA[4] = ( A[0]*A[8]-A[2]*A[6] )*idet;
        invA[5] = ( A[2]*A[4]-A[0]*A[5] )*idet;
        invA[6] = ( A[3]*A[7]-A[4]*A[6] )*idet;
        invA[7] = ( A[1]*A[6]-A[0]*A[7] )*idet;
        invA[8] = ( A[0]*A[4]-A[1]*A[3] )*idet;
    }
    
    return invA;
}

// product between a square matrix A and a vector v
vector<double> matrixTimesVector(vector<double> A, vector<double> v) {
    unsigned int size = A.size();
    unsigned int dim  = sqrt( size );
    if( dim*dim != size ) ERROR("Matrix is not squared");
    if( v.size() != dim ) ERROR("Vector has wrong size");
    
    vector<double> w(dim,0.);
    for( unsigned int i=0; i<dim; i++ )
        for( unsigned int j=0; j<dim; j++ )
                w[i] += A[i+dim*j] * v[j];
    return w;
}




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
    
    // Extract "flush_every" (time selection for flushing the file)
    flush_timeSelection = new TimeSelection( 
        PyTools::extract_py("flush_every", "DiagProbe", n_probe),
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
    
    // Extract "pos", "pos_first", "pos_second" and "pos_third"
    // (positions of the vertices of the grid)
    allPos.resize(0);
    vector<string> keys(4);
    keys[0] = "pos";
    keys[1] = "pos_first";
    keys[2] = "pos_second";
    keys[3] = "pos_third";
    for( unsigned int i=0; i<nDim_particle+1; i++) {
        vector<double> pos;
        if (PyTools::extract(keys[i],pos,"DiagProbe",n_probe)) {
            if (pos.size()!=nDim_particle)
                ERROR("Probe #"<<n_probe<<": "<<keys[i]<<" size(" << pos.size() << ") != ndim(" << nDim_particle<< ")");
            allPos.push_back(pos);
        }
    }
    
    // calculate the coordinate system (base vectors)
    axes.resize(nDim_particle*nDim_particle, 0.);
    vector<int> usedDirections (nDim_particle, 0);
    double min, max, val;
    unsigned int jmin, jmax;
    // With available axes, fill the vector, and remember which major direction it occupies
    for( unsigned int i=0; i<dimProbe; i++) {
        max = 0.;
        min = numeric_limits<double>::max();
        for( unsigned int j=0; j<nDim_particle; j++) {
            axes[j+nDim_particle*i] = allPos[i+1][j] - allPos[0][j];
            val = abs(axes[j+nDim_particle*i]);
            if( val<min ) { min=val; jmin=j; }
            if( val>max ) { max=val; jmax=j; }
        }
        usedDirections[jmax] += 3; // avoid max
        usedDirections[jmin] -= 1; // prefer min
    }
    // Then, complete the probe's coordinate system to have as many axes as the simulation dimension
    for( unsigned int i=dimProbe; i<nDim_particle; i++) {
        // find index of the most unused direction
        unsigned int unusedDirectionIndex = min_element(usedDirections.begin(), usedDirections.end()) - usedDirections.begin();
        // and use that direction as next axis
        axes[i+nDim_particle*unusedDirectionIndex] = 1.;
        for( unsigned int j=0; j<nDim_particle; j++) usedDirections[j]--;
        usedDirections[unusedDirectionIndex] += 4;
    }
    // Calculate the inverse matrix of the probe's coordinate system
    axesInverse = matrixInverse(axes);
    
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
    
} // END DiagnosticProbes::DiagnosticProbes



DiagnosticProbes::~DiagnosticProbes()
{
    delete timeSelection;
    delete flush_timeSelection;
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
        
        // Dimension of the probe grid
        H5::attr(fileId_, "dimension", dimProbe);
        
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
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        
        unsigned int i, k, ipart, iDimProbe, iDim;
        vector<double> partPos(nDim_particle);
        vector<unsigned int> ipartND (dimProbe);
        bool is_in_domain;
        
        // The first step is to reduce the area of the probe to search in this patch
        unsigned int numCorners = 1<<nDim_particle; // number of patch corners
        vector<double> point(nDim_particle, 0.);
        vector<double> mins(nDim_particle, numeric_limits<double>::max());        
        vector<double> maxs(nDim_particle, numeric_limits<double>::min());        
        // loop patch corners
        for( i=0; i<numCorners; i++ ) {
            // Get coordinates of the current corner in terms of x,y,...
            for( k=0; k<nDim_particle; k++ )
                point[k] = ( (((i>>k)&1)==0) ? vecPatches(ipatch)->getDomainLocalMin(k) : vecPatches(ipatch)->getDomainLocalMax(k) )
                          - allPos[0][k];
            // Get position of the current corner in the probe's coordinate system
            point = matrixTimesVector( axesInverse, point );
            // Store mins and maxs
            for( k=0; k<nDim_particle; k++ ) {
                if( point[k]<mins[k] ) mins[k]=point[k];
                if( point[k]>maxs[k] ) maxs[k]=point[k];
            }
        }
        // Loop directions to figure out the range of useful indices
        vector<unsigned int> minI(nDim_particle, 0), maxI(nDim_particle, 0), nI(nDim_particle, 0);
        for( i=0; i<dimProbe; i++ ) {
            if( mins[i]<0. ) mins[i]=0.;
            if( mins[i]>1. ) mins[i]=1.;
            minI[i] = (unsigned int) floor(mins[i]*((double)(vecNumber[i]-1)));
            if( maxs[i]<0. ) maxs[i]=0.;
            if( maxs[i]>1. ) maxs[i]=1.;
            maxI[i] = (unsigned int) ceil (maxs[i]*((double)(vecNumber[i]-1)));
        }
        for( i=dimProbe; i<nDim_particle; i++ )
            if( mins[i]*maxs[i] < 0 )
                maxI[i]=0;
        // Now, minI and maxI contain the min and max indexes of the probe, useful for this patch
        // Calculate total number of useful points
        unsigned int ntot = 1;
        for( i=0; i<nDim_particle; i++ ) {
            nI[i] = maxI[i]-minI[i]+1;
            ntot *= nI[i];
        }
        
        // Initialize the list of "fake" particles (points) just as actual macro-particles
        Particles * particles = &(vecPatches(ipatch)->probes[probe_n]->particles);
        particles->initialize(ntot, nDim_particle);
        
        // Loop useful probe points
        unsigned int IP, ipart_local=0;
        for( unsigned int ip=0; ip<ntot; ip++ ) {
            // Find the coordinates of this point in the global probe array
            IP = ip;
            for( i=0; i<dimProbe; i++ ) {
                point[i] = ((double)(IP % nI[i] + minI[i])) / ((double)(vecNumber[i]-1));
                IP /= nI[i];
            }
            for( i=dimProbe; i<nDim_particle; i++ )
                point[i] = 0.;
            // Compute this point's coordinates in terms of x, y, ...
            point = matrixTimesVector( axes, point );
            // Check if point is in patch
            is_in_domain = true;
            for( i=0; i<nDim_particle; i++ ) {
                point[i] += allPos[0][i];
                if (point[i] <  vecPatches(ipatch)->getDomainLocalMin(i) 
                 || point[i] >= vecPatches(ipatch)->getDomainLocalMax(i) ) {
                    is_in_domain = false;
                    break;
                }
            }
            if(is_in_domain) {
                for(iDim=0; iDim<nDim_particle; iDim++)
                    particles->position(iDim,ipart_local) = point[iDim];
                ipart_local++;
            }
        }
        
        // Resize the array with only particles in this patch
        particles->resize(ipart_local, nDim_particle);
        particles->shrink_to_fit(nDim_particle);
        
        // Add the local offset
        vecPatches(ipatch)->probes[probe_n]->offset_in_file = nPart_MPI;
        nPart_MPI += ipart_local;
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
        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
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
    unsigned int ipart = 0;
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        if( ipart>=(unsigned int)nPart_MPI ) break;
        Particles * particles = &(vecPatches(ipatch)->probes[probe_n]->particles);
        for ( unsigned int ip=0 ; ip<particles->size() ; ip++) {
            for (unsigned int idim=0 ; idim<nDim_particle  ; idim++ ){
                posArray->data_2D[ipart][idim] = particles->position(idim,ip);
            }
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
    if ( nPart_MPI>0 )
        H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, transfer, &(posArray->data_2D[0][0]) );
    else
        H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, transfer, NULL );
    H5Dclose(dset_id);
    H5Pclose( transfer );
    H5Sclose(filespace);
    H5Sclose(memspace);
    delete posArray;
    
    H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
}



void DiagnosticProbes::run( SmileiMPI* smpi, VectorPatch& vecPatches, int timestep )
{
    unsigned int nPart_MPI;
    ostringstream name_t;
    
    unsigned int nPatches( vecPatches.size() );

    // Leave if this timestep has already been written
    #pragma omp master
    {
        name_t.str("");
        name_t << "/" << setfill('0') << setw(10) << timestep;
        status = H5Lexists( fileId_, name_t.str().c_str(), H5P_DEFAULT );
    }
    #pragma omp barrier
    if (status != 0) return;
    smpi->barrier();
    
    #pragma omp master
    {
        // Calculate the number of probe particles in this MPI
        offset.resize(vecPatches.size()+1);
        offset[0] = 0;
        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
            offset[ipatch+1] = offset[ipatch] + vecPatches(ipatch)->probes[probe_n]->particles.size();
        nPart_MPI = offset.back();
        
        // Make the array that will contain the data
        vector<unsigned int> probesArraySize(2);
        probesArraySize[1] = nPart_MPI; // number of particles
        probesArraySize[0] = nFields + 1; // number of fields (Ex, Ey, etc) +1 for garbage
        probesArray = new Field2D(probesArraySize);
    }
    #pragma omp barrier
    
    // Loop patches to fill the array
    #pragma omp for schedule(static)
    for (unsigned int ipatch=0 ; ipatch<nPatches ; ipatch++) {
        // Loop probe ("fake") particles of current patch
        unsigned int iPart_MPI = offset[ipatch];
        unsigned int npart = vecPatches(ipatch)->probes[probe_n]->particles.size();

        LocalFields Eloc_fields, Bloc_fields, Jloc_fields;
        double Rloc_fields;

        for (unsigned int ipart=0; ipart<npart; ipart++) {             
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
    
    #pragma omp master
    {
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
        
        if( flush_timeSelection->theTimeIsNow(timestep) ) H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
    }
    #pragma omp barrier
}
