#include "H5.h"

//! Open HDF5 file + location
H5::H5( std::string file, unsigned access, MPI_Comm * comm, bool _raise )
{
    init( file, access, comm, _raise );
}

void H5::init( std::string file, unsigned access, MPI_Comm * comm, bool _raise )
{
    
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
    if( comm ) {
        H5Pset_fapl_mpio( fapl, *comm, MPI_INFO_NULL );
    }
    if( access == H5F_ACC_RDWR ) {
        fid_ = H5Fcreate( filepath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl );
    } else {
        fid_ = H5Fopen( filepath.c_str(), access, fapl );
    }
    H5Pclose( fapl );
    
    // Check error
    if( H5Eget_num( H5E_DEFAULT ) > 0 ) {
        fid_ = -1;
        id_ = -1;
        dxpl_ = -1;
        dcr_ = -1;
    } else {
        // Open group from file
        id_ = fid_;
        if( ! grouppath.empty() ) {
            id_ = H5Gopen( fid_, grouppath.c_str(), H5P_DEFAULT );
        }
        dxpl_ = H5Pcreate( H5P_DATASET_XFER );
        dcr_ = H5Pcreate( H5P_DATASET_CREATE );
        if( comm ) {
            H5Pset_dxpl_mpio( dxpl_, H5FD_MPIO_COLLECTIVE );
            H5Pset_alloc_time( dcr_, H5D_ALLOC_TIME_EARLY );
        }
        // Check error
        if( H5Eget_num( H5E_DEFAULT ) > 0 ) {
            id_ = -1;
        }
    }
    
    if( _raise ) {
        if( fid_ < 0 || id_ < 0 ) {
            ERROR( "Cannot open file " << filepath );
        }
    } else {
        // Restore previous error printing
        H5Eset_auto( H5E_DEFAULT, old_func, old_client_data );
    }
}


//! Location already opened
H5::H5( hid_t id, hid_t dcr, hid_t dxpl ) : fid_( -1 ), id_( id ), dcr_( dcr ), dxpl_( dxpl )
{
}

H5::~H5()
{
    if( H5Iget_type( id_ ) == H5I_GROUP ) {
        H5Gclose( id_ );
    } else if( H5Iget_type( id_ ) == H5I_DATASET ) {
        H5Dclose( id_ );
    }
    if( fid_ >= 0 ) {
        H5Pclose( dxpl_ );
        H5Pclose( dcr_ );
        herr_t err = H5Fclose( fid_ );
        if( err < 0 ) {
            H5Eprint2( H5E_DEFAULT, NULL );
            ERROR( "Can't close file " << filepath );
        }
    }
}


//! 1D
H5Space::H5Space( hsize_t size ) {
    dims_ = { size };
    global_ = size;
    sid = H5Screate_simple( 1, &size, NULL );
    if( size <= 0 ) {
        H5Sselect_none( sid );
    }
    chunk_.resize(0);
}

//! 1D
H5Space::H5Space( hsize_t size, hsize_t offset, hsize_t npoints, hsize_t chunk ) {
    dims_ = { size };
    global_ = size;
    sid = H5Screate_simple( 1, &size, NULL );
    if( size <= 0 || npoints <= 0 ) {
        H5Sselect_none( sid );
    } else {
        hsize_t count = 1;
        H5Sselect_hyperslab( sid, H5S_SELECT_SET, &offset, NULL, &count, &npoints );
    }
    if( chunk > 1 ) {
        chunk_.resize( 1, chunk );
    } else {
        chunk_.resize( 0 );
    }
}

//! ND
H5Space::H5Space( std::vector<hsize_t> size, std::vector<hsize_t> offset, std::vector<hsize_t> npoints, std::vector<hsize_t> chunk ) {
    dims_ = size;
    global_ = 1;
    for( unsigned int i=0; i<size.size(); i++ ) {
        global_ *= size[i];
    }
    sid = H5Screate_simple( size.size(), &size[0], NULL );
    if( global_ <= 0 ) {
        H5Sselect_none( sid );
    } else if( ! offset.empty() || ! npoints.empty() ) {
        unsigned int i;
        for( i=0; i<npoints.size() && npoints[i]>0; i++ ) {}
        if( i < npoints.size() ) {
            H5Sselect_none( sid );
        } else {
            std::vector<hsize_t> count( size.size(), 1 );
            H5Sselect_hyperslab( sid, H5S_SELECT_SET, &offset[0], NULL, &count[0], &npoints[0] );
        }
    }
    chunk_ = chunk;
}
