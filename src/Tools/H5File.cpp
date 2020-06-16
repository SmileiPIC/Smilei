#include "H5File.h"

//! Group already opened
H5Group::H5Group( hid_t ID ) : id( ID ), empty( false )
{
}

//! Open group from file
H5Group::H5Group( hid_t lid, std::string grouppath ) {
    id = lid;
    empty = true;
    if( ! grouppath.empty() && lid >= 0 ) {
        id = H5Gopen( lid, grouppath.c_str(), H5P_DEFAULT );
        empty = false;
    }
}

H5Group::~H5Group() {
    if( ! empty ) {
        H5Gclose( id );
    }
}

H5File::H5File( std::string file, unsigned access, bool _raise ) {
    
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
    
    // Open
    fid = H5Fopen( filepath.c_str(), access, H5P_DEFAULT );
    
    // Check error stack size
    if( H5Eget_num( H5E_DEFAULT ) > 0 ) {
        fid = -1;
    }
    
    if( ! _raise ) {
        // Restore previous error printing
        H5Eset_auto( H5E_DEFAULT, old_func, old_client_data );
    }
    
    if( _raise && fid < 0 ) {
        ERROR( "Cannot open file " << filepath );
    }
}

