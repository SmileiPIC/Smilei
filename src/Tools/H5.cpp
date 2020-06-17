#include "H5.h"

//! Open HDF5 file + location
H5::H5( std::string file, unsigned access, bool parallel, bool _raise )
{
    parallel_ = parallel;
    
    // Analyse file string : separate file name and tree inside hdf5 file
    size_t l = file.length();
    size_t i_h5 = file.find( ".h5" ) + 3;
    size_t i_hdf5 = file.find( ".hdf5" ) + 5;
    filepath = file;
    grouppath = "";
    if( i_h5 < std::string::npos ) {
        if( i_h5 == l || ( i_h5 < l && file[i_h5] == '/' ) ) {
            filepath = file.substr( 0, i_h5 );
            grouppath = file.substr( i_h5 );
        }
    } else if( i_hdf5 < std::string::npos ) {
        if( i_hdf5 == l || ( i_hdf5 < l && file[i_hdf5] == '/' ) ) {
            filepath = file.substr( 0, i_hdf5 );
            grouppath = file.substr( i_hdf5 );
        }
    }
    
    H5E_auto2_t old_func;
    void *old_client_data;
    if( ! _raise ) {
        // Backup default error printing
        H5Eget_auto( H5E_DEFAULT, &old_func, &old_client_data );
        H5Eset_auto( H5E_DEFAULT, NULL, NULL );
    }
    
    // Open or create
    hid_t fapl = H5Pcreate( H5P_FILE_ACCESS );
    if( parallel ) {
        H5Pset_fapl_mpio( fapl, MPI_COMM_WORLD, MPI_INFO_NULL );
    }
    if( access == H5F_ACC_RDWR ) {
        fid = H5Fcreate( filepath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl );
    } else {
        fid = H5Fopen( filepath.c_str(), access, fapl );
    }
    
    // Check error
    if( H5Eget_num( H5E_DEFAULT ) > 0 ) {
        fid = -1;
        id = -1;
    } else {
        // Open group from file
        id = fid;
        if( ! grouppath.empty() ) {
            id = H5Gopen( fid, grouppath.c_str(), H5P_DEFAULT );
        }
        dxpl = H5Pcreate( H5P_DATASET_XFER );
        dcr = H5Pcreate( H5P_DATASET_CREATE );
        if( parallel ) {
            H5Pset_dxpl_mpio( dxpl, H5FD_MPIO_COLLECTIVE );
            H5Pset_alloc_time( dcr, H5D_ALLOC_TIME_EARLY );
        }
        // Check error
        if( H5Eget_num( H5E_DEFAULT ) > 0 ) {
            id = -1;
        }
    }
    
    if( ! _raise ) {
        // Restore previous error printing
        H5Eset_auto( H5E_DEFAULT, old_func, old_client_data );
    } else if( fid < 0 || id < 0 ) {
        ERROR( "Cannot open file " << filepath );
    }
}


//! Location already opened
H5::H5( hid_t ID, bool parallel ) : fid( -1 ), id( ID ), parallel_( parallel )
{
    dxpl = H5Pcreate( H5P_DATASET_XFER );
    dcr = H5Pcreate( H5P_DATASET_CREATE );
    if( parallel ) {
        H5Pset_dxpl_mpio( dxpl, H5FD_MPIO_COLLECTIVE );
        H5Pset_alloc_time( dcr, H5D_ALLOC_TIME_EARLY );
    }
}

H5::~H5()
{
    if( H5Iget_type( id ) == H5I_GROUP ) {
        H5Gclose( id );
    } else if( H5Iget_type( id ) == H5I_DATASET ) {
        H5Dclose( id );
    }
    H5Pclose( dxpl );
    H5Pclose( dcr );
    if( fid >= 0 ) {
        if( H5Fclose( fid ) < 0 ) {
            ERROR( "Can't close file " << filepath );
        }
    }
}