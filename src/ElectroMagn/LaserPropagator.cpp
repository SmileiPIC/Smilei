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
void fftfreq( vector<double> &freqs, unsigned int n, double d, unsigned int start, unsigned int stop )
{
    freqs.resize( stop-start );
    unsigned int N = ( n-1 ) / 2 + 1;
    double nm = -( ( double )n );
    double val = 2.*M_PI/( n*d );
    unsigned int imax = min( N, stop );
    for( unsigned int i=start       ; i<imax; i++ ) {
        freqs[i-start] = val*i;
    }
    for( unsigned int i=max( N, start ); i<stop; i++ ) {
        freqs[i-start] = val*( nm+i );
    }
}

// equivalent of the numpy function linspace
// imin or imax define a range where the array is actually returned
vector<double> linspace( double xmin, double xmax, unsigned int n, unsigned int imin, unsigned int imax )
{
    vector<double> ret( imax-imin );
    double val = (xmax-xmin)/( n-1 );
    for( unsigned int i=imin; i<imax; i++ ) {
        ret[i-imin] = val*i + xmin;
    }
    return ret;
}

// equivalent of the numpy function argsort, but reverse sorting and only the largest n numbers
vector<unsigned int> partial_reverse_argsort( vector<double> x, unsigned int nmax )
{
    unsigned int n = x.size();
    std::vector<unsigned int> y( n );
    for( unsigned int i=0; i<n; i++ ) {
        y[i] = i;
    }
    nmax = min( nmax, n );
    std::partial_sort(
        y.begin(),
        y.begin()+nmax,
        y.end(),
        [&]( unsigned int i1, unsigned int i2 ) {
            return x[i1] > x[i2];
        }
    );
    y.resize( nmax );
    return y;
}


// Call a python function
PyObject *PyCall( PyObject *callable, PyObject *args, PyObject *kwargs )
{
    PyObject *a = PyObject_Call( callable, args, kwargs );
    Py_DECREF( args );
    if( kwargs ) {
        Py_DECREF( kwargs );
    }
    return a;
}


void LaserPropagator::init( Params *params, SmileiMPI *smpi, unsigned int side )
{

#ifdef SMILEI_USE_NUMPY
    ndim = params->nDim_field;
    _2D = ndim==2;
    MPI_size = smpi->getSize();
    MPI_rank = smpi->getRank();

    N     .resize( 3, 0 ); // third value not used in 2D
    L     .resize( ndim );
    Nlocal.resize( ndim );
    o     .resize( ndim, 0. );

    // Set the grid spatial dimensions
    for( unsigned int idim=0; idim<ndim-1; idim++ ) {
        unsigned int j = ( side+idim+1 )%ndim;
        N[idim] = params->n_space_global[j] + 2*params->oversize[j] + 2;
        L[idim] = N[idim] * params->cell_length[j];
        o[idim] = params->oversize[j] * params->cell_length[j];
    }
    ox = params->oversize[0] * params->cell_length[0];

    // Set the grid temporal dimension
    N[ndim-1] = params->n_time;
    L[ndim-1] = params->simulation_time;

    // Make the array bigger to accommodate for the parallel FFT
    double old_L = L[0];
    for( unsigned int idim=0; idim<2; idim++ ) {
        double old_N = ( double )N[idim];
        N[idim]  = ( ( int ) ceil( ( double )N[idim] / MPI_size ) ) * MPI_size;
        L[idim] *= ( ( double )N[idim]/old_N );
    }
    double O = L[0] - old_L;
    
    // Calculate some array size that relates to the parallel decomposition
    Nlocal[0] = N[0] / MPI_size;
    Nlocal[1] = N[1] / MPI_size;

    // Arrays of coordinates (y, z, t) owning to the current processor
    local_x.resize( ndim );
    local_x[0] = linspace( -o[0]-O, L[0]-o[0]-O, N[0], MPI_rank*Nlocal[0], ( MPI_rank+1 )*Nlocal[0] );
    local_x[1] = linspace( -o[1], L[1]-o[1], N[1], 0, N[1] );
    if( ! _2D ) {
        local_x[2] = linspace( 0, L[2], N[2], 0, N[2] );
    }

    // Arrays of wavenumbers (ky, kz, w) owning to the current processor
    local_k.resize( ndim );
    fftfreq( local_k[0], N[0], L[0]/N[0], 0, N[0] );
    fftfreq( local_k[1], N[1], L[1]/N[1], MPI_rank*Nlocal[1], ( MPI_rank+1 )*Nlocal[1] );
    if( ! _2D ) {
        fftfreq( local_k[2], N[2], L[2]/N[2], 0, N[2] );
    }

#else
    ERROR( "Cannot use LaserOffset without numpy" );
#endif

}

void LaserPropagator::operator()( vector<PyObject *> profiles, vector<int> profiles_n, double offset, string file, int keep_n_strongest_modes, double angle_z )
{
#ifdef SMILEI_USE_NUMPY
    //const complex<double> i_ (0., 1.); // the imaginary number

    unsigned int nprofiles = profiles.size();

    // Import several functions from numpy
    PyObject *numpy = PyImport_AddModule( "numpy" );
    PyObject *numpyfft = PyObject_GetAttrString( numpy, "fft" );
    PyObject *fft      = PyObject_GetAttrString( numpyfft, "fft" );
    PyObject *ifft     = PyObject_GetAttrString( numpyfft, "ifft" );
    PyObject *fft2     = PyObject_GetAttrString( numpyfft, "fft2" );
    Py_DECREF( numpyfft );
    int complex_type_num=0;

    // 1- Calculate the value of the profiles at all points (y,z,t)
    // --------------------------------

    // Make coordinates array
    vector<PyObject *> coords( ndim );
    for( unsigned int i=0; i<ndim; i++ ) {
        npy_intp dims = local_x[i].size();
        coords[i] = PyArray_SimpleNewFromData( 1, &dims, NPY_DOUBLE, ( double * )( local_x[i].data() ) );
    }

    // Make a "meshgrid" using the coordinates
    PyObject *meshgrid = PyObject_GetAttrString( numpy, "meshgrid" );
    PyObject *m;
    if( _2D ) {
        m = PyCall( meshgrid, Py_BuildValue( "(O,O)", coords[0], coords[1] ), Py_BuildValue( "{s:s}", "indexing", "ij" ) );
    } else {
        m = PyCall( meshgrid, Py_BuildValue( "(O,O,O)", coords[0], coords[1], coords[2] ), Py_BuildValue( "{s:s}", "indexing", "ij" ) );
    }
    PyObject *mesh = PySequence_Tuple( m );
    for( unsigned int i=0; i<ndim; i++ ) {
        Py_DECREF( coords[i] );
    }
    Py_DECREF( meshgrid );
    Py_DECREF( m );

    // Apply each profile
    vector<PyObject *> arrays( nprofiles );
    for( unsigned int i=0; i<nprofiles; i++ ) {
        // Vectorize the profile
        PyObject *profile = PyObject_CallMethod( numpy, "vectorize", "O", profiles[i] );
        // Apply to the mesh
        arrays[i] = PyObject_CallObject( profile, mesh );
        Py_DECREF( profile );
    }
    Py_DECREF( mesh );

    // 2- Fourier transform of the fields at destination
    // --------------------------------

    for( unsigned int i=0; i<nprofiles; i++ ) {
        PyObject *a;

        // Call FFT along the last direction(s)
        if( _2D ) {
            a = PyCall( fft, Py_BuildValue( "(O)", arrays[i] ), Py_BuildValue( "{s:i}", "axis", 1 ) );
        } else {
            a = PyCall( fft2, Py_BuildValue( "(O)", arrays[i] ), Py_BuildValue( "{s:(i,i)}", "axes", 1, 2 ) );
        }
        Py_DECREF( arrays[i] );

        // Change the array shape to prepare the MPI comms (reshape, transpose, C ordering)
        npy_intp D[4] = {Nlocal[0], MPI_size, Nlocal[1], N[2]};
        PyArray_Dims shape;
        shape.len=ndim+1;
        shape.ptr=D;
        arrays[i] = PyArray_Newshape( ( PyArrayObject * )a, &shape, NPY_CORDER );
        Py_DECREF( a );
        npy_intp P[4] = {1, 0, 2, 3};
        PyArray_Dims permute;
        permute.len=ndim+1;
        permute.ptr=P;
        a = PyArray_Transpose( ( PyArrayObject * )arrays[i], &permute );
        Py_DECREF( arrays[i] );
        arrays[i] = PyArray_FROM_OF( a, NPY_ARRAY_C_CONTIGUOUS );
        Py_DECREF( a );


        // Obtain the number associated to the complex type in numpy
        complex_type_num = PyArray_TYPE( ( PyArrayObject * ) arrays[i] );

        // Make a new empty array for MPI comms
        npy_intp dims[ndim];
        dims[0] = N[0];
        dims[1] = Nlocal[1];
        if( ! _2D ) {
            dims[2] = N[2];
        }
        a = PyArray_EMPTY( ndim, dims, complex_type_num, 0 );

        // Communicate blocks to transpose the MPI decomposition
        int block_size = Nlocal[0]*Nlocal[1]* ( _2D?1:N[2] );
        MPI_Alltoall(
            PyArray_GETPTR1( ( PyArrayObject * ) arrays[i], 0 ), 2*block_size, MPI_DOUBLE,
            PyArray_GETPTR1( ( PyArrayObject * ) a, 0 ), 2*block_size, MPI_DOUBLE,
            MPI_COMM_WORLD
        );
        Py_DECREF( arrays[i] );

        // Call FFT along the first direction
        arrays[i] = PyCall( fft, Py_BuildValue( "(O)", a ), Py_BuildValue( "{s:i}", "axis", 0 ) );
        Py_DECREF( a );
    }

    // 3- Select only interesting omegas
    // --------------------------------

    // Get the total spectrum
    vector<double> spectrum( 1, 0. );
    if( _2D ) {
        // Compute the spectrum locally
        vector<double> local_spectrum( Nlocal[1], 0. );
        for( unsigned int i=0; i<nprofiles; i++ ) {
            complex<double> *z = ( complex<double> * ) PyArray_GETPTR1( ( PyArrayObject * ) arrays[i], 0 );
            for( unsigned int k=0; k<Nlocal[1]; k++ )
                for( unsigned int j=0; j<N[0]; j++ ) {
                    local_spectrum[k] += abs( z[j + N[0]*k] );
                }
        }
        // In 2D, the spectrum is scattered across processors, so we gather to root
        if( MPI_rank==0 ) {
            spectrum.resize( N[1] );
        }
        MPI_Gather( &local_spectrum[0], Nlocal[1], MPI_DOUBLE,
                    &      spectrum[0], Nlocal[1], MPI_DOUBLE,
                    0, MPI_COMM_WORLD );
        // Remove the negative frequencies
        if( MPI_rank==0 ) {
            spectrum.resize( N[1]/2 );
        }
    } else {
        // Compute the spectrum locally
        unsigned int lmax = N[2]/2;
        vector<double> local_spectrum( lmax, 0. );
        for( unsigned int i=0; i<nprofiles; i++ ) {
            complex<double> *z = ( complex<double> * ) PyArray_GETPTR1( ( PyArrayObject * ) arrays[i], 0 );
            for( unsigned int l=0; l<lmax; l++ )
                for( unsigned int k=0; k<Nlocal[1]; k++ )
                    for( unsigned int j=0; j<N[0]; j++ ) {
                        local_spectrum[l] += abs( z[j + N[0]*( k + Nlocal[1]*l )] );
                    }
        }
        // In 3D, each processor has the full spectrum, so we sum all contributions
        if( MPI_rank==0 ) {
            spectrum.resize( lmax );
        }
        MPI_Reduce( &local_spectrum[0], &spectrum[0], lmax, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    }

    // Rank 0 finds the most intense points of the spectrum
    unsigned int n_omega, n_omega_local;
    vector<unsigned int> indices;
    if( MPI_rank == 0 ) {
        indices = partial_reverse_argsort( spectrum, keep_n_strongest_modes );
        sort( indices.begin(), indices.end() );
        n_omega = indices.size();
    }

    // Broadcast the number of selected points
    MPI_Bcast( &n_omega, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    indices.resize( n_omega );

    // Broadcast the selected indices to all processes
    MPI_Bcast( &indices[0], n_omega, MPI_UNSIGNED, 0, MPI_COMM_WORLD );

    MESSAGE( 2, "keeping the "<<n_omega<<" strongest temporal modes" );

    // Find the indices belonging to this process
    if( _2D ) {
        // Find the range [imin, imax[ of indices in this proc
        unsigned int imin=0;
        unsigned int proc_min= MPI_rank   *Nlocal[1];
        unsigned int proc_max=( MPI_rank+1 )*Nlocal[1];
        while( imin < n_omega && indices[imin] < proc_min ) {
            imin++;
        }
        unsigned int imax = imin;
        while( imax < n_omega && indices[imax] < proc_max ) {
            imax++;
        }
        n_omega_local = imax-imin;
        // Restrict to that range
        unsigned int i0=0;
        for( unsigned int i1=imin; i1<imax; i1++ ) {
            indices[i0] = indices[i1] - proc_min;
            i0++;
        }
        indices.resize( n_omega_local );
    } else {
        n_omega_local = n_omega;
    }

    // Make a list of the omegas
    vector<double> omega( n_omega_local );
    for( unsigned int i=0; i<n_omega_local; i++ ) {
        omega[i] = local_k[ndim-1][indices[i]];
    }

    // Make arrays of ky^2, kz^2
    vector<double> ky2( N[0] ), kz2;
    for( unsigned int i=0; i<N[0]; i++ ) {
        ky2[i] = local_k[0][i] * local_k[0][i];
    }
    if( ! _2D ) {
        kz2.resize( Nlocal[1] );
        for( unsigned int i=0; i<Nlocal[1]; i++ ) {
            kz2[i] = local_k[1][i] * local_k[1][i];
        }
    }

    // 4- Propagate and rotate in k-space, while extracting only selected omegas
    // --------------------------------

    // Define the rotation parameters
    double cz = cos( angle_z ), sz = sin( angle_z );
    double coeff_y = L[0] / 2. / M_PI;
    double kx, j_, j_max=( double )( ( N[0]-1 )/2 ), j_min=-( double )( N[0]/2 ), j_tot=N[0];
    unsigned int i1, i0;

    // Interpolate arrays on new k space
    double omega2;
    for( unsigned int i=0; i<nprofiles; i++ ) {
        if( _2D ) {
            npy_intp dims[2] = {N[0], n_omega_local};
            PyObject *a = PyArray_New( &PyArray_Type, 2, dims, complex_type_num, NULL, NULL, 128, NPY_ARRAY_F_CONTIGUOUS, NULL );
            complex<double> *z0 = ( complex<double> * ) PyArray_GETPTR1( ( PyArrayObject * ) arrays[i], 0 );
            complex<double> *z  = ( complex<double> * ) PyArray_GETPTR1( ( PyArrayObject * ) a, 0 );
            for( unsigned int k=0; k<n_omega_local; k++ ) {
                omega2 = omega[k] * omega[k];
                i1 = N[0]*k;
                i0 = N[0]*indices[k];
                for( unsigned int j=0; j<N[0]; j++ ) {
                    if( ky2[j] >= omega2 ) {
                        z[j + i1] = 0.;
                    } else {
                        kx = sqrt( omega2 - ky2[j] );
                        j_ = coeff_y * ( sz * kx + cz * local_k[0][j] ); // index of ky before rotation
                        if( j_ < j_min || j_ > j_max ) { // out of bounds of FFT
                            z[j + i1] = 0.;
                        } else {
                            if( j_ < 0. ) {
                                j_ += j_tot;    // unwrap negative frequencies
                            }
                            z[j + i1] = complex_interpolate( &z0[i0], N[0], j_, // interpolate
                                                             abs( cz - sz * local_k[0][j]/kx ), // compensate for the k-space compression
                                                             (offset+ox)* kx  // add the propagation term to the phase
                                                           );
                        }
                    }
                }
            }
            Py_DECREF( arrays[i] );
            arrays[i] = a;
        } else {
            npy_intp dims[3] = {N[0], Nlocal[1], n_omega_local};
            PyObject *a = PyArray_New( &PyArray_Type, 3, dims, complex_type_num, NULL, NULL, 128, NPY_ARRAY_F_CONTIGUOUS, NULL );
            complex<double> *z0 = ( complex<double> * ) PyArray_GETPTR1( ( PyArrayObject * ) arrays[i], 0 );
            complex<double> *z  = ( complex<double> * ) PyArray_GETPTR1( ( PyArrayObject * ) a, 0 );
            for( unsigned int l=0; l<n_omega_local; l++ ) {
                omega2 = omega[l] * omega[l];
                for( unsigned int k=0; k<Nlocal[1]; k++ ) {
                    i1 = N[0]*( k + Nlocal[1]*l );
                    i0 = N[0]*( k + Nlocal[1]*indices[l] );
                    for( unsigned int j=0; j<N[0]; j++ ) {
                        if( ky2[j] + kz2[k] >= omega2 ) {
                            z[j + i1] = 0.;
                        } else {
                            kx = sqrt( omega2 - ky2[j] - kz2[k] );
                            j_ = coeff_y * ( sz * kx + cz * local_k[0][j] ); // index of ky before rotation
                            if( j_ < j_min || j_ > j_max ) { // out of bounds of FFT
                                z[j + i1] = 0.;
                            } else {
                                if( j_ < 0. ) {
                                    j_ += j_tot;    // unwrap negative frequencies
                                }
                                z[j + i1] = complex_interpolate( &z0[i0], N[0], j_, // interpolate
                                                                 abs( cz - sz * local_k[0][j]/kx ), // compensate for the k-space compression
                                                                 (offset+ox)* kx  // add the propagation term to the phase
                                                               );
                            }
                        }
                    }
                }
            }
            Py_DECREF( arrays[i] );
            arrays[i] = a;
        }
    }

    // 5- Fourier transform back to real space, excluding the omega axis
    // --------------------------------

    for( unsigned int i=0; i<nprofiles; i++ ) {
        PyObject *a;

        // Call FFT along the first direction
        a = PyCall( ifft, Py_BuildValue( "(O)", arrays[i] ), Py_BuildValue( "{s:i}", "axis", 0 ) );
        Py_DECREF( arrays[i] );

        // Convert the array to C order
        arrays[i] = PyArray_FROM_OF( a, NPY_ARRAY_C_CONTIGUOUS );
        Py_DECREF( a );

        if( ! _2D ) {

            // Make a new empty array for MPI comms
            npy_intp dims[4];
            dims[0] = MPI_size;
            dims[1] = Nlocal[0];
            dims[2] = Nlocal[1];
            dims[3] = n_omega_local;
            a = PyArray_EMPTY( 4, dims, complex_type_num, 0 );

            // Communicate blocks to transpose the MPI decomposition
            int block_size = Nlocal[0]*Nlocal[1]*n_omega_local;
            MPI_Alltoall(
                PyArray_GETPTR1( ( PyArrayObject * ) arrays[i], 0 ), 2*block_size, MPI_DOUBLE,
                PyArray_GETPTR1( ( PyArrayObject * ) a, 0 ), 2*block_size, MPI_DOUBLE,
                MPI_COMM_WORLD
            );
            Py_DECREF( arrays[i] );

            // Change the array shape to accomodate the MPI comms
            npy_intp P[4] = {1, 0, 2, 3};
            PyArray_Dims permute;
            permute.len=4;
            permute.ptr=P;
            arrays[i] = PyArray_Transpose( ( PyArrayObject * )a, &permute );
            Py_DECREF( a );
            npy_intp D[3] = {Nlocal[0], N[1], n_omega_local};
            PyArray_Dims shape;
            shape.len=3;
            shape.ptr=D;
            a = PyArray_Newshape( ( PyArrayObject * )arrays[i], &shape, NPY_CORDER );
            Py_DECREF( arrays[i] );
            arrays[i] = PyArray_FROM_OF( a, NPY_ARRAY_C_CONTIGUOUS );
            Py_DECREF( a );

            // Call FFT along the second direction
            a = PyCall( ifft, Py_BuildValue( "(O)", arrays[i] ), Py_BuildValue( "{s:i}", "axis", 1 ) );
            Py_DECREF( arrays[i] );

            // Convert the array to C order
            arrays[i] = PyArray_FROM_OF( a, NPY_ARRAY_C_CONTIGUOUS );
            Py_DECREF( a );
        }
    }

    // 6- Obtain the magnitude and the phase of the complex values
    // --------------------------------
    unsigned int local_size = ( _2D ? N[0] : Nlocal[0]*N[1] ) * n_omega_local;
    vector<vector<double> > magnitude( nprofiles ), phase( nprofiles );

    for( unsigned int i=0; i<nprofiles; i++ ) {
        complex<double> *z = ( complex<double> * ) PyArray_GETPTR1( ( PyArrayObject * ) arrays[i], 0 );
        magnitude[i].resize( local_size );
        phase    [i].resize( local_size );
        double coeff_magnitude = 2./N[ndim-1]; // multiply by omega increment
        if( profiles_n[i]==1 ) {
            coeff_magnitude *= cz;    // multiply by cosine for By only
        }
        for( unsigned int j=0; j<local_size; j++ ) {
            magnitude[i][j] = abs( z[j] ) * coeff_magnitude;
            phase    [i][j] = arg( z[j] );
        }
        Py_DECREF( arrays[i] );
    }

    // 7- Store all info in HDF5 file
    // --------------------------------

    // Set the file region where each proc will write
    hsize_t dims[ndim], start[ndim], count[ndim];
    if( _2D ) {
        unsigned int MPI_omega_offset;
        MPI_Scan( &n_omega_local, &MPI_omega_offset, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );
        MPI_omega_offset -= n_omega_local;
        dims[0] = N[0];
        start[0] = 0;
        count[0] = N[0];
        dims[1] = n_omega;
        start[1] = MPI_omega_offset;
        count[1] = n_omega_local;
    } else {
        dims[0] = N[0];
        start[0] = MPI_rank*Nlocal[0];
        count[0] = Nlocal[0];
        dims[1] = N[1];
        start[1] = 0;
        count[1] = N[1];
        dims[2] = n_omega;
        start[2] = 0;
        count[2] = n_omega  ;
    }

    // Create File with parallel access
    hid_t faid = H5Pcreate( H5P_FILE_ACCESS );
    H5Pset_fapl_mpio( faid, MPI_COMM_WORLD, MPI_INFO_NULL );
    hid_t fid  = H5Fcreate( file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, faid );
    H5Pclose( faid );

    // Store "omega" dataset
    hsize_t len = n_omega, len_local = n_omega_local;
    hid_t filespace = H5Screate_simple( 1, &len, NULL );
    hid_t memspace  = H5Screate_simple( 1, &len_local, NULL );
    hid_t dcid = H5Pcreate( H5P_DATASET_CREATE );
    hid_t did = H5Dcreate( fid, "omega", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, dcid, H5P_DEFAULT );
    if( _2D && n_omega_local > 0 ) {
        H5Sselect_hyperslab( filespace, H5S_SELECT_SET, &start[1], NULL, &count[1], NULL );
    } else if( !_2D && MPI_rank==0 ) { // Only rank 0 writes in 3D
        H5Sselect_hyperslab( filespace, H5S_SELECT_SET, &start[2], NULL, &count[2], NULL );
    } else {
        H5Sselect_none( filespace );
        H5Sselect_none( memspace );
    }
    hid_t transfer = H5Pcreate( H5P_DATASET_XFER );
    H5Pset_dxpl_mpio( transfer, H5FD_MPIO_COLLECTIVE );
    H5Dwrite( did, H5T_NATIVE_DOUBLE, memspace, filespace, transfer, &omega[0] );
    H5Sclose( filespace );
    H5Sclose( memspace );
    H5Dclose( did );

    // Store the magnitude and the phase
    filespace = H5Screate_simple( ndim, dims, NULL );
    H5Sselect_hyperslab( filespace, H5S_SELECT_SET, start, NULL, count, NULL );
    memspace  = H5Screate_simple( ndim, count, NULL );
    for( unsigned int i=0; i<nprofiles; i++ ) {
        // Magnitude
        ostringstream name( "" );
        name << "magnitude" << profiles_n[i];
        did = H5Dcreate( fid, name.str().c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, dcid, H5P_DEFAULT );
        H5Dwrite( did, H5T_NATIVE_DOUBLE, memspace, filespace, transfer, &( magnitude[i][0] ) );
        H5Dclose( did );
        // Phase
        name.str( "" );
        name << "phase" << profiles_n[i];
        did = H5Dcreate( fid, name.str().c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, dcid, H5P_DEFAULT );
        H5Dwrite( did, H5T_NATIVE_DOUBLE, memspace, filespace, transfer, &( phase[i][0] ) );
        H5Dclose( did );
    }
    H5Sclose( filespace );
    H5Sclose( memspace );
    H5Pclose( dcid );
    H5Pclose( transfer );
    H5Fclose( fid );

    Py_DECREF( fft );
    Py_DECREF( ifft );
    Py_DECREF( fft2 );

#endif
}
