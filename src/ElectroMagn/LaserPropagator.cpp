#include "Tools.h"
#include "Params.h"
#include "LaserPropagator.h"
#include "Field2D.h"
#include "Field3D.h"
#include "H5.h"

#include <cmath>
#include <string>
#include <complex>
#include <algorithm>

using namespace std;

// equivalent of the numpy function fftfreq
// start or stop define a range where the array is actually returned
void fftfreq(vector<double>& freqs, unsigned int n, double d, unsigned int start, unsigned int stop)
{
    freqs.resize(stop-start);
    unsigned int N = (n-1) / 2 + 1;
    double nm = -((double)n);
    double val = 2.*M_PI/(n*d);
    unsigned int imax = min(N,stop);
    for( unsigned int i=start       ; i<imax; i++ ) freqs[i-start] = val*i;
    for( unsigned int i=max(N,start); i<stop; i++ ) freqs[i-start] = val*(nm+i);
}

// equivalent of the numpy function linspace
// start or stop define a range where the array is actually returned
vector<double> linspace(unsigned int n, double d, unsigned int start, unsigned int stop)
{
    vector<double> ret(stop-start);
    double val = d/(n-1);
    for( unsigned int i=start; i<stop; i++ ) ret[i-start] = val*i;
    return ret;
}

// equivalent of the numpy function argsort, but reverse sorting and only the largest n numbers
vector<unsigned int> partial_reverse_argsort(vector<double> x, unsigned int nmax) {
    unsigned int n = x.size();
    std::vector<unsigned int> y(n);
    for(unsigned int i=0; i<n; i++) y[i] = i;
    std::partial_sort(  y.begin(), y.begin()+nmax, y.end(),
                [&](unsigned int i1, unsigned int i2) { return x[i1] > x[i2]; } );
    y.resize(nmax);
    return y;
}


// Call a python function
PyObject* PyCall(PyObject* callable, PyObject* args, PyObject* kwargs)
{
    PyObject* a = PyObject_Call(callable, args, kwargs);
    Py_DECREF(args);
    if( kwargs ) Py_DECREF(kwargs);
    return a;
}


void LaserPropagator::init(Params* params, SmileiMPI* smpi, unsigned int side)
{

#ifdef SMILEI_USE_NUMPY
    ndim = params->nDim_field;
    _2D = ndim==2;
    MPI_size = smpi->getSize();
    MPI_rank = smpi->getRank();
    
    N     .resize( ndim );
    L     .resize( ndim );
    Nlocal.resize( ndim );
    
    // Set the grid spatial dimensions
    for( unsigned int idim=0; idim<ndim-1; idim++ ) {
        unsigned int j = (side+idim+1)%ndim;
        N[idim] = params->n_space_global[j] + 2*params->oversize[j] + 2;
        L[idim] = N[idim] * params->cell_length[j];
    }
    
    // Set the grid temporal dimension
    N[ndim-1] = params->n_time;
    L[ndim-1] = params->simulation_time;
    
    // Make the array bigger to accommodate for the parallel FFT
    N[0] =  ((int) ceil( (double)N[0] / MPI_size )) * MPI_size;
    N[1] =  ((int) ceil( (double)N[1] / MPI_size )) * MPI_size;
    
    // Calculate some array size that relates to the parallel decomposition
    Nlocal[0] = N[0] / MPI_size;
    Nlocal[1] = N[1] / MPI_size;
    
    // Arrays of coordinates (y, z, t) owning to the current processor
    local_x.resize(ndim);
    local_x[0] = linspace(N[0], L[0], MPI_rank*Nlocal[0], (MPI_rank+1)*Nlocal[0]);
    local_x[1] = linspace(N[1], L[1], 0., N[1]);
    if( ! _2D )
        local_x[2] = linspace(N[2], L[2], 0., N[2]);
    
    // Arrays of wavenumbers (ky, kz, w) owning to the current processor
    local_k.resize(ndim);
    fftfreq(local_k[0], N[0], L[0]/N[0], 0, N[0]);
    fftfreq(local_k[1], N[1], L[1]/N[1], MPI_rank*Nlocal[1], (MPI_rank+1)*Nlocal[1]);
    if( ! _2D )
        fftfreq(local_k[2], N[2], L[2]/N[2], 0, N[2]);
    
#else
    ERROR("Cannot use LaserOffset without numpy");
#endif

}

void LaserPropagator::operator() (vector<PyObject*> profiles, vector<int> profiles_n, double offset, string file, int keep_n_best_frequencies)
{
#ifdef SMILEI_USE_NUMPY
    const complex<double> i_ (0., 1.); // the imaginary number
    
    unsigned int nprofiles = profiles.size();
    
    ostringstream ro_("");
    for(unsigned int i=0;i<MPI_rank;i++) ro_<<"             ";
    
    MPI_Barrier(MPI_COMM_WORLD);
    cout << ro_.str() << "0" << endl;;
    MESSAGE(1,"nprofiles = "<<nprofiles);
    MESSAGE(1,"N = "<<N[0]<<" "<<N[1]);
    MESSAGE(1,"Nlocal = "<<Nlocal[0]<<" "<<Nlocal[1]);
    MESSAGE(1,"L = "<<L[0]<<" "<<L[1]);
    
    // Import several functions from numpy
    PyObject* numpy = PyImport_AddModule("numpy");
    PyObject* numpyfft = PyObject_GetAttrString(numpy   , "fft"  );
    PyObject* fft      = PyObject_GetAttrString(numpyfft, "fft"  );
    PyObject* ifft     = PyObject_GetAttrString(numpyfft, "ifft" );
    PyObject* fft2     = PyObject_GetAttrString(numpyfft, "fft2" );
    PyObject* ifft2    = PyObject_GetAttrString(numpyfft, "ifft2");
    Py_DECREF(numpyfft);
    
    // Import other functions (see pyprofiles.py)
    PyObject* applyProfiles    = PyObject_GetAttrString(PyImport_AddModule("__main__"), "_applyProfiles");
    PyObject* emptyComplex     = PyObject_GetAttrString(PyImport_AddModule("__main__"), "_emptyComplex");
    PyObject* reshapeTranspose = PyObject_GetAttrString(PyImport_AddModule("__main__"), "_reshapeTranspose");
    PyObject* transposeReshape = PyObject_GetAttrString(PyImport_AddModule("__main__"), "_transposeReshape");
    
    // Obtain the number associated to the complex type in numpy
    unsigned int complex_type_num=0;
    PyObject* complex128 = PyObject_GetAttrString(numpy, "complex128");
    PyObject* complex128_dtype = PyObject_CallMethod(numpy, "dtype", "O", complex128);
    PyObject* complex128_num = PyObject_GetAttrString(complex128_dtype, "num");
    PyTools::convert( complex128_num, complex_type_num );
    Py_DECREF(complex128_num);
    Py_DECREF(complex128_dtype);
    Py_DECREF(complex128);
    
    MPI_Barrier(MPI_COMM_WORLD);
    cout << ro_.str() << 1 << endl;;
    // 1- Calculate the value of the profiles at all points (y,z,t)
    // --------------------------------
    
    // Make coordinates array
    vector<PyObject*> coords(ndim);
    for( unsigned int i=0; i<ndim; i++ ) {
        npy_intp dims = local_x[i].size();
        coords[i] = PyArray_SimpleNewFromData(1, &dims, NPY_DOUBLE, (double*)(local_x[i].data()));
    }
    // Call the python function _applyProfiles
    PyObject *ret, *tuple_coords, *tuple_profiles;
    if( _2D ) tuple_coords = Py_BuildValue("(O,O)"  , coords[0], coords[1]           );
    else      tuple_coords = Py_BuildValue("(O,O,O)", coords[0], coords[1], coords[2]);
    if( nprofiles==1 ) tuple_profiles = Py_BuildValue("(O)"  , profiles[0]             );
    else               tuple_profiles = Py_BuildValue("(O,O)", profiles[0], profiles[1]);
    ret = PyCall(applyProfiles, Py_BuildValue("(O,O)", tuple_coords, tuple_profiles), NULL);
    Py_DECREF(tuple_coords);
    Py_DECREF(tuple_profiles);
    // Recover the several arrays
    vector<PyObject*> arrays( nprofiles );
    for( unsigned int i=0; i<nprofiles; i++ )
        arrays[i] = PySequence_GetItem(ret, i);
    for( unsigned int i=0; i<ndim; i++ )
        Py_DECREF(coords[i]);
    Py_DECREF(ret);
    
    MPI_Barrier(MPI_COMM_WORLD);
    cout << ro_.str() << 2 << endl;;
    // 2- Fourier transform of the fields at destination
    // --------------------------------
    
    for( unsigned int i=0; i<nprofiles; i++ ) {
        PyObject* a;
        
        // Call FFT along the last direction(s)
        if( _2D ) a = PyCall( fft , Py_BuildValue("(O)", arrays[i]), Py_BuildValue("{s:i}"    , "axis", 1   ) );
        else      a = PyCall( fft2, Py_BuildValue("(O)", arrays[i]), Py_BuildValue("{s:(i,i)}", "axes", 1, 2) );
        Py_DECREF(arrays[i]);
        
        // Change the array shape to prepare the MPI comms
        if( _2D ) arrays[i] = PyCall( reshapeTranspose, Py_BuildValue("(O, (i,i,i), (i,i,i))"    , a, Nlocal[0], MPI_size, Nlocal[1]      , 1, 0, 2   ), NULL );
        else      arrays[i] = PyCall( reshapeTranspose, Py_BuildValue("(O, (i,i,i,i), (i,i,i,i))", a, Nlocal[0], MPI_size, Nlocal[1], N[2], 1, 0, 2, 3), NULL );
        Py_DECREF(a);
        
        // Make a new empty array for MPI comms
        if( _2D ) a = PyCall( emptyComplex, Py_BuildValue("((i,i))"  , N[0], Nlocal[1]      ), NULL );
        else      a = PyCall( emptyComplex, Py_BuildValue("((i,i,i))", N[0], Nlocal[1], N[2]), NULL );
        
        // Communicate blocks to transpose the MPI decomposition
        int block_size = Nlocal[0]*Nlocal[1]* (_2D?1:N[2]);
        MPI_Alltoall(
            PyArray_GETPTR1((PyArrayObject*) arrays[i],0), block_size, MPI_DOUBLE_COMPLEX,
            PyArray_GETPTR1((PyArrayObject*) a        ,0), block_size, MPI_DOUBLE_COMPLEX,
            MPI_COMM_WORLD
        );
        Py_DECREF(arrays[i]);
        
        // Call FFT along the first direction
        arrays[i] = PyCall( fft, Py_BuildValue("(O)", a), Py_BuildValue("{s:i}", "axis", 0) );
        Py_DECREF(a);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    cout << ro_.str() << 3 << endl;;
    // 3- Select only interesting omegas
    // --------------------------------
    
    // Get the total spectrum
    vector<double> spectrum(1,0.);
    if( _2D ) {
        // Compute the spectrum locally
        vector<double> local_spectrum( Nlocal[1], 0. );
        for( unsigned int i=0; i<nprofiles; i++ ) {
            complex<double> * z = (complex<double> *) PyArray_GETPTR1((PyArrayObject*) arrays[i],0);
            for(unsigned int j=0; j<N[0]; j++) 
                for(unsigned int k=0; k<Nlocal[1]; k++) 
                    local_spectrum[k] += abs(z[j + N[0]*k]);
        }
        // In 2D, the spectrum is scattered across processors, so we gather to root
        if( MPI_rank==0 ) spectrum.resize( N[1] );
        MPI_Gather( &local_spectrum[0], Nlocal[1], MPI_DOUBLE,
                    &      spectrum[0], Nlocal[1], MPI_DOUBLE,
                    0, MPI_COMM_WORLD);
        // Remove the negative frequencies
        if( MPI_rank==0 ) spectrum.resize( N[1]/2 );
    } else {
        // Compute the spectrum locally
        unsigned int lmax = N[2]/2;
        vector<double> local_spectrum( lmax, 0. );
        for( unsigned int i=0; i<nprofiles; i++ ) {
            complex<double> * z = (complex<double> *) PyArray_GETPTR1((PyArrayObject*) arrays[i],0);
            for(unsigned int j=0; j<N[0]; j++) 
                for(unsigned int k=0; k<Nlocal[1]; k++) 
                    for(unsigned int l=0; l<lmax; l++) 
                        local_spectrum[l] += abs(z[j + N[0]*( k + Nlocal[1]*l )]);
        }
        // In 3D, each processor has the full spectrum, so we sum all contributions
        if( MPI_rank==0 ) spectrum.resize( lmax );
        MPI_Reduce( &local_spectrum[0], &spectrum[0], lmax, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    
    // Rank 0 finds the most intense points of the spectrum
    unsigned int n_omega, n_omega_local;
    vector<unsigned int> indices;
    if( MPI_rank == 0 ) {
        indices = partial_reverse_argsort(spectrum, keep_n_best_frequencies);
        sort(indices.begin(), indices.end());
        n_omega = indices.size();
    }
    
    // Broadcast the number of selected points
    MPI_Bcast(&n_omega, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    indices.resize(n_omega);
    
    // Broadcast the selected indices to all processes
    MPI_Bcast(&indices[0], n_omega, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    
    // Find the indices belonging to this process
    if( _2D ) {
        // Find the range [imin, imax[ of indices in this proc
        unsigned int imin=0;
        unsigned int proc_min= MPI_rank   *Nlocal[1];
        unsigned int proc_max=(MPI_rank+1)*Nlocal[1];
        while( imin < n_omega && indices[imin] < proc_min ) imin++;
        unsigned int imax = imin;
        while( imax < n_omega && indices[imax] < proc_max ) imax++;
        n_omega_local = imax-imin;
        cout<<ro_.str()<<imin<<" ! " <<imax<<endl;
        // Restrict to that range
        unsigned int i0=0;
        for( unsigned int i1=imin; i1<imax; i1++ ) {
            indices[i0] = indices[i1] - proc_min;
            i0++;
        }
        indices.resize(n_omega_local);
    } else {
        n_omega_local = n_omega;
    }
    
    // Make a list of the omegas
    vector<double> omega(n_omega_local);
    for(unsigned int i=0; i<n_omega_local; i++)
        omega[i] = local_k[ndim-1][indices[i]];
    
    ostringstream t("");
    t<<"\trank "<<MPI_rank<<" indices ( ";
    for( unsigned int i=0;i<n_omega_local;i++ ) 
        t<<indices[i]<<" ";
    t<<" )"<<endl;
    MPI_Barrier(MPI_COMM_WORLD);
    for( unsigned int i=0; i<MPI_size; i++) {
        if( i==MPI_rank ) cout<< t.str();
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    t.str("");
    t<<"\trank "<<MPI_rank<<" omega ( ";
    for( unsigned int i=0;i<n_omega_local;i++ ) 
        t<<omega[i]<<" ";
    t<<" )"<<endl;
    for( unsigned int i=0; i<MPI_size; i++) {
        if( i==MPI_rank ) cout<< t.str();
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    cout << ro_.str() << "4" << endl;;
    MESSAGE(1,"n_omega = "<<n_omega);
    MESSAGE(1,"n_omega_local = "<<n_omega_local);
    // 4- Multiply by the propagation term, while extracting only selected omegas
    // --------------------------------
        
    for( unsigned int i=0; i<nprofiles; i++ ) {
        if( _2D ) {
            npy_intp dims[2] = {N[0], n_omega_local};
            PyObject* a = PyArray_New(&PyArray_Type, 2, dims, complex_type_num, NULL, NULL, 128, 1, NULL);
            complex<double> * z0 = (complex<double> *) PyArray_GETPTR1((PyArrayObject*) arrays[i],0);
            complex<double> * z  = (complex<double> *) PyArray_GETPTR1((PyArrayObject*) a        ,0);
            for(unsigned int j=0; j<N[0]; j++) 
                for(unsigned int k=0; k<n_omega_local; k++) 
                    z[j + N[0]*k] = z0[j + N[0]*indices[k]]
                        * exp( i_ * offset * sqrt( max(0., pow(omega[k],2) - pow(local_k[0][j],2)) ) );
            Py_DECREF(arrays[i]);
            arrays[i] = a;
        } else {
            npy_intp dims[3] = {N[0], Nlocal[1], n_omega_local};
            PyObject* a = PyArray_New(&PyArray_Type, 3, dims, complex_type_num, NULL, NULL, 128, 1, NULL);
            complex<double> * z0 = (complex<double> *) PyArray_GETPTR1((PyArrayObject*) arrays[i],0);
            complex<double> * z  = (complex<double> *) PyArray_GETPTR1((PyArrayObject*) a        ,0);
            for(unsigned int j=0; j<N[0]; j++) 
                for(unsigned int k=0; k<Nlocal[1]; k++) 
                    for(unsigned int l=0; l<n_omega_local; l++) 
                        z[j + N[0]*( k + Nlocal[1]*l )] = z0[j + N[0]*( k + Nlocal[1]*indices[l] )]
                            * exp( i_ * offset * sqrt( max(0., pow(omega[l],2) - pow(local_k[0][j],2) - pow(local_k[1][k],2)) ) );
            Py_DECREF(arrays[i]);
            arrays[i] = a;
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    cout << ro_.str() << 5 << endl;;
    MESSAGE(1," is fortran ? "<< PyArray_ISFORTRAN((PyArrayObject*)arrays[0]));
    // 5- Fourier transform back to real space, excluding the omega axis
    // --------------------------------
    
    for( unsigned int i=0; i<nprofiles; i++ ) {
        PyObject* a;
        
        // Call FFT along the first direction
        a = PyCall( ifft, Py_BuildValue("(O)", arrays[i]), Py_BuildValue("{s:i}", "axis", 0) );
        Py_DECREF(arrays[i]);
        
        if( _2D ) {
            arrays[i] = a;
        } else {
            // Make a new empty array for MPI comms
            arrays[i] = PyCall( emptyComplex, Py_BuildValue("((i,i,i,i))", MPI_size, Nlocal[0], Nlocal[1], n_omega_local), NULL );
            
            // Communicate blocks to transpose the MPI decomposition
            int block_size = Nlocal[0]*Nlocal[1]*n_omega_local;
            MPI_Alltoall(
                PyArray_GETPTR1((PyArrayObject*) a        ,0), block_size, MPI_DOUBLE_COMPLEX,
                PyArray_GETPTR1((PyArrayObject*) arrays[i],0), block_size, MPI_DOUBLE_COMPLEX,
                MPI_COMM_WORLD
            );
            Py_DECREF(a);
            
            // Change the array shape to prepare the MPI comms
            a = PyCall( transposeReshape, Py_BuildValue("(O, (i,i,i), (i,i,i,i))", arrays[i], Nlocal[0], N[1], n_omega_local, 1, 0, 2, 3), NULL );
            Py_DECREF(arrays[i]);
            
            // Call FFT along the second direction
            arrays[i] = PyCall( ifft2, Py_BuildValue("(O)", a), Py_BuildValue("{s:(i)}", "axes", 1) );
            Py_DECREF(a);
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    cout << ro_.str() << 6 << endl;;
    // 6- Obtain the magnitude and the phase of the complex values
    // --------------------------------
    unsigned int local_size = (_2D ? N[0] : Nlocal[0]*N[1]) * n_omega_local;
    vector<vector<double> > magnitude(nprofiles), phase(nprofiles);
    double coeff_magnitude = 1./L.back(); // divide by omega increment
    
    for( unsigned int i=0; i<nprofiles; i++ ) {
        complex<double> * z = (complex<double> *) PyArray_GETPTR1((PyArrayObject*) arrays[i],0);
        magnitude[i].resize( local_size );
        phase    [i].resize( local_size );
        unsigned int i0, i1;
        if( _2D) {
            for(unsigned int j=0; j<N[0]; j++) {
                for(unsigned int k=0; k<n_omega_local; k++) {
                    // Note that we transpose for converting fortran order to C order
                    i0 = j + N[0]*k;
                    i1 = k + n_omega_local*j;
                    magnitude[i][i1] = abs(z[i0]) * coeff_magnitude;
                    phase    [i][i1] = arg(z[i0] * exp(-i_ * omega[k] * offset));
                }
            }
        } else {
            for(unsigned int j=0; j<Nlocal[0]; j++) {
                for(unsigned int k=0; k<N[1]; k++) {
                    for(unsigned int l=0; l<n_omega_local; l++) {
                        // Note that we transpose for converting fortran order to C order
                        i0 = j + Nlocal[0]*(k + N[1]*l);
                        i1 = l + n_omega_local*(k + N[1]*j);
                        magnitude[i][i1] = abs(z[i0]) * coeff_magnitude;
                        phase    [i][i1] = arg(z[i0] * exp(-i_ * omega[l] * offset));
                    }
                }
            }
        }
        Py_DECREF( arrays[i] );
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    cout << ro_.str() << 7 << endl;;
    // 7- Store all info in HDF5 file
    // --------------------------------
    
    // Set the file region where each proc will write
    hsize_t dims[ndim], start[ndim], count[ndim];
    if( _2D ) {
        unsigned int MPI_omega_offset;
        MPI_Scan(&n_omega_local, &MPI_omega_offset, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
        MPI_omega_offset -= n_omega_local;
        dims[0] = N[0]   ; start[0] = 0               ; count[0] = N[0]         ;
        dims[1] = n_omega; start[1] = MPI_omega_offset; count[1] = n_omega_local;
        cout<<MPI_rank<<" dims  " << dims[0]  << " " << dims[1] <<endl;
        cout<<MPI_rank<<" start " << start[0] << " " << start[1]<<endl;
        cout<<MPI_rank<<" count " << count[0] << " " << count[1]<<endl;
    } else {
        dims[0] = N[0]   ; start[0] = MPI_rank*Nlocal[0]; count[0] = Nlocal[0];
        dims[1] = N[1]   ; start[1] = 0                 ; count[1] = N[1]     ;
        dims[2] = n_omega; start[2] = 0                 ; count[2] = n_omega  ;
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    cout << ro_.str() << 8 << endl;;
    // Create File with parallel access
    hid_t faid = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(faid, MPI_COMM_WORLD, MPI_INFO_NULL);
    hid_t fid  = H5Fcreate( file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, faid);
    H5Pclose(faid);
    
    MPI_Barrier(MPI_COMM_WORLD);
    cout << ro_.str() << 9 << endl;;
    // Store "omega" dataset
    hsize_t len = n_omega, len_local = n_omega_local;
    hid_t filespace = H5Screate_simple(1, &len, NULL);
    hid_t memspace  = H5Screate_simple(1, &len_local, NULL);
    hid_t dcid = H5Pcreate(H5P_DATASET_CREATE);
    hid_t did = H5Dcreate( fid, "omega", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, dcid, H5P_DEFAULT);
    if( (_2D && n_omega_local > 0) || (!_2D && MPI_rank==0) ) { // Only rank 0 writes in 3D
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &start[1], NULL, &count[1], NULL);
    } else {
        H5Sselect_none(filespace);
    }
    hid_t transfer = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(transfer, H5FD_MPIO_COLLECTIVE);
    H5Dwrite( did, H5T_NATIVE_DOUBLE, memspace, filespace, transfer, &omega[0] );
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Dclose(did);
    
    MPI_Barrier(MPI_COMM_WORLD);
    cout << ro_.str() << 10 << endl;;
    // Store the magnitude and the phase
    filespace = H5Screate_simple(ndim, dims, NULL);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);
    memspace  = H5Screate_simple(ndim, count, NULL);
    for( unsigned int i=0; i<nprofiles; i++ ) {
        // Magnitude
        ostringstream name("");
        name << "magnitude" << profiles_n[i];
        did = H5Dcreate( fid, name.str().c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, dcid, H5P_DEFAULT);
        H5Dwrite( did, H5T_NATIVE_DOUBLE, memspace, filespace, transfer, &(magnitude[i][0]));
        H5Dclose(did);
        // Phase
        name.str("");
        name << "phase" << profiles_n[i];
        did = H5Dcreate( fid, name.str().c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, dcid, H5P_DEFAULT);
        H5Dwrite( did, H5T_NATIVE_DOUBLE, memspace, filespace, transfer, &(phase[i][0]));
        H5Dclose(did);
    }
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(dcid);
    H5Pclose(transfer);
    H5Fclose(fid);
    
    Py_DECREF(fft  );
    Py_DECREF(ifft );
    Py_DECREF(fft2 );
    Py_DECREF(ifft2);
    Py_DECREF(applyProfiles   );
    Py_DECREF(emptyComplex    );
    Py_DECREF(reshapeTranspose);
    Py_DECREF(transposeReshape);
    
#else
    return 0;
#endif
}

