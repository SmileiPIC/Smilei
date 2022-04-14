#ifndef H5_H
#define H5_H

#include <hdf5.h>
#include <string>
#include <sstream>
#include <vector>
#include "Tools.h"

#if ! H5_HAVE_PARALLEL == 1
#error "HDF5 was not built with --enable-parallel option"
#endif

class DividedString
{
public:
    DividedString( unsigned int w ) : width( std::max( w, ( unsigned int )1 ) ), numstr( 0 ), str( "" ) {};
    ~DividedString() {};
    
    void addString( std::string s )
    {
        if( s.size() <= width ) {
            str.append( s );
            std::ostringstream t( "" );
            for( unsigned int i=0; i<width-s.size(); i++ ) {
                t<<"\0";
            }
            str.append( t.str() );
        } else {
            str.append( s.substr( 0, width ) );
        }
        numstr++;
    };
    
    const char *c_str()
    {
        return str.c_str();
    };
    
    unsigned int width, numstr;
    std::string str;
};

class H5Space
{
public:
    
    //! 1D
    H5Space( hsize_t size );
    H5Space( hsize_t size, hsize_t offset, hsize_t npoints, hsize_t chunk = 0 );
    
    //! ND
    H5Space( std::vector<hsize_t> size, std::vector<hsize_t> offset = {}, std::vector<hsize_t> npoints = {}, std::vector<hsize_t> chunk = {} );
    
    ~H5Space() {
        H5Sclose( sid );
    }
    
    hid_t sid;
    std::vector<hsize_t> dims_;
    std::vector<hsize_t> chunk_;
    hsize_t global_;
    
};

class H5
{
public:
    //! Empty HDF5 object
    H5() {
        fid_ = -1;
        id_ = -1;
        dxpl_ = -1;
        dcr_ = -1;
    };
    
    //! Open HDF5 file + location
    H5( std::string file, unsigned access, MPI_Comm * comm, bool _raise );
    
    ~H5();
    
    void init( std::string file, unsigned access, MPI_Comm * comm, bool _raise );
    
    bool valid() {
        return id_ >= 0;
    }
    
    void flush() {
        H5Fflush( id_, H5F_SCOPE_GLOBAL );
    }
    
    //! Check if group exists
    bool has( std::string group_name )
    {
        return H5Lexists( id_, group_name.c_str(), H5P_DEFAULT ) > 0;
    }
    
    //! Check if attribute exists
    bool hasAttr( std::string attribute_name )
    {
        return H5Aexists( id_, attribute_name.c_str() ) > 0;
    }
    
protected:
    //! Constructor when location already opened
    H5( hid_t ID, hid_t dcr, hid_t dxpl );
    
    std::string filepath;
    std::string grouppath;
    hid_t fid_;
    hid_t id_;
    hid_t dcr_;
    hid_t dxpl_;
    
    hid_t newGroupId( std::string group_name ) {
        if( H5Lexists( id_, group_name.c_str(), H5P_DEFAULT ) > 0 ) {
            return H5Oopen( id_, group_name.c_str(), H5P_DEFAULT );
        } else {
            return H5Gcreate( id_, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
        }
    }
    
    hid_t open( std::string name ) {
        if( H5Lexists( id_, name.c_str(), H5P_DEFAULT ) > 0 ) {
            return H5Oopen( id_, name.c_str(), H5P_DEFAULT );
        } else {
            ERROR( "In HDF5 file, "<< name << " does not exist" );
        }
        return -1;
    }
};

class H5Write : public H5
{
public:
    //! Open HDF5 file + location
    H5Write( std::string file, MPI_Comm * comm = NULL, bool _raise = true )
     : H5( file, H5F_ACC_RDWR, comm, _raise ) {};
    
    //! Create group given H5Write location
    H5Write( H5Write *loc, std::string group_name )
     : H5( loc->newGroupId( group_name ), loc->dcr_, loc->dxpl_ ) {};
    
    //! Create or open (not write) a dataset given H5Write location
    H5Write( H5Write *loc, std::string name, hid_t type, H5Space *filespace )
     : H5( -1, loc->dcr_, loc->dxpl_ )
    {
        H5D_layout_t layout = H5Pget_layout( dcr_ );
        if( ! filespace->chunk_.empty() ) {
            H5Pset_chunk( dcr_, filespace->chunk_.size(), &filespace->chunk_[0] );
        }
        if( H5Lexists( loc->id_, name.c_str(), H5P_DEFAULT ) == 0 ) {
            id_  = H5Dcreate( loc->id_, name.c_str(), type, filespace->sid, H5P_DEFAULT, dcr_, H5P_DEFAULT );
        } else {
            hid_t pid = H5Pcreate( H5P_DATASET_ACCESS );
            id_ = H5Dopen( loc->id_, name.c_str(), pid );
            H5Pclose( pid );
        }
        H5Pset_layout( dcr_, layout );
    }
    
    //! Location already opened
    H5Write( hid_t id, hid_t dcr, hid_t dxpl ) : H5( id, dcr, dxpl ) {};
    
    ~H5Write() {};
    
    //! Make or open a group
    H5Write group( std::string group_name )
    {
        return H5Write( this, group_name );
    }
    
    //! Open an existing dataset
    H5Write dataset( std::string name )
    {
        return H5Write( open( name ), dcr_, dxpl_ );
    }
    
    //! Write a string as an attribute
    void attr( std::string attribute_name, std::string attribute_value )
    {
        hid_t atype = H5Tcopy( H5T_C_S1 );
        if( attribute_value.size() == 0 ) {
            H5Tset_size( atype, 1 );
        } else {
            H5Tset_size( atype, attribute_value.size() );
        }
        const char *tmp_var = attribute_value.c_str();
        attr( attribute_name, *tmp_var, atype );
        H5Tclose( atype );
    }
    
    //! Write an unsigned int as an attribute
    void attr( std::string attribute_name, unsigned int attribute_value )
    {
        attr( attribute_name, attribute_value, H5T_NATIVE_UINT );
    }
    
    //! write unsigned long int as an attribute
    void attr( std::string attribute_name, unsigned long int attribute_value )
    {
        attr( attribute_name, attribute_value, H5T_NATIVE_ULONG );
    }
    
    //! write an int as an attribute
    void attr( std::string attribute_name, int attribute_value )
    {
        attr( attribute_name, attribute_value, H5T_NATIVE_INT );
    }
    
    //! write a double as an attribute
    void attr( std::string attribute_name, double attribute_value )
    {
        attr( attribute_name, attribute_value, H5T_NATIVE_DOUBLE );
    }
    
    //! write anything as an attribute
    template<class T>
    void attr( std::string attribute_name, T &attribute_value, hid_t type )
    {
        hid_t sid = H5Screate( H5S_SCALAR );
        hid_t aid = H5Acreate( id_, attribute_name.c_str(), type, sid, H5P_DEFAULT, H5P_DEFAULT );
        H5Awrite( aid, type, &attribute_value );
        H5Sclose( sid );
        H5Aclose( aid );
    }
    
    //! write a vector<unsigned int> as an attribute
    void attr( std::string attribute_name, std::vector<unsigned int> attribute_value )
    {
        attr( attribute_name, attribute_value, H5T_NATIVE_UINT );
    }
    
    //! write an vector<double> as an attribute
    void attr( std::string attribute_name, std::vector<double> attribute_value )
    {
        attr( attribute_name, attribute_value, H5T_NATIVE_DOUBLE );
    }
    
    //! write a vector<string> as an attribute
    void attr( std::string attribute_name, std::vector<std::string> attribute_value )
    {
        std::vector<const char *> tmp_vec( attribute_value.size(), nullptr );
        for( unsigned int i=0; i<attribute_value.size(); i++ ) {
            tmp_vec[i] = attribute_value[i].c_str();
        }
        attr( attribute_name, tmp_vec );
    }
    
    //! write a vector<const char*> as an attribute
    void attr( std::string attribute_name, std::vector<const char *> attribute_value )
    {
        hid_t atype = H5Tcopy( H5T_C_S1 );
        H5Tset_size( atype, H5T_VARIABLE );
        attr( attribute_name, attribute_value, atype );
        H5Tclose( atype );
    }
    
    //! write a vector<anything> as an attribute
    template<class T>
    void attr( std::string attribute_name, std::vector<T> &attribute_value, hid_t type )
    {
        hsize_t dims = attribute_value.size();
        hid_t sid = H5Screate_simple( 1, &dims, NULL );
        hid_t aid = H5Acreate( id_, attribute_name.c_str(), type, sid, H5P_DEFAULT, H5P_DEFAULT );
        H5Awrite( aid, type, &( attribute_value[0] ) );
        H5Aclose( aid );
        H5Sclose( sid );
    }
    
    //! write a DividedString as an attribute
    void attr( std::string attribute_name, DividedString &attribute_value )
    {
        hid_t atype = H5Tcopy( H5T_C_S1 );
        H5Tset_size( atype, attribute_value.width );
        hsize_t n = attribute_value.numstr;
        hid_t sid = H5Screate_simple( 1, &n, NULL );
        hid_t aid = H5Acreate( id_, attribute_name.c_str(), atype, sid, H5P_DEFAULT, H5P_DEFAULT );
        H5Awrite( aid, atype, attribute_value.c_str() );
        H5Tclose( atype );
        H5Aclose( aid );
        H5Sclose( sid );
    }
    
    //! Write a vector<int>
    H5Write vect( std::string name, std::vector<int> v, hsize_t offset=0, hsize_t npoints=0 )
    {
        return vect( name, v[0], v.size(), H5T_NATIVE_INT, offset, npoints );
    }
    
    //! write a vector<unsigned int>
    H5Write vect(std::string name, std::vector<unsigned int> v,  hsize_t offset=0, hsize_t npoints=0 )
    {
        return vect(name, v[0], v.size(), H5T_NATIVE_UINT, offset, npoints );
    }
    
    //! write a vector<short>
    H5Write vect( std::string name, std::vector<short> v, hsize_t offset=0, hsize_t npoints=0 )
    {
        return vect( name, v[0], v.size(), H5T_NATIVE_SHORT, offset, npoints );
    }
    
    //! write a vector<doubles>
    H5Write vect( std::string name, std::vector<double> v, hsize_t offset=0, hsize_t npoints=0 )
    {
        return vect( name, v[0], v.size(), H5T_NATIVE_DOUBLE, offset, npoints );
    }
    
    //! write any vector
    template<class T>
    H5Write vect( std::string name, std::vector<T> v, hid_t type, hsize_t offset=0, hsize_t npoints=0 )
    {
        return vect( name, v[0], v.size(), type, offset, npoints );
    }
    
    //! Write a portion of a vector
    template<class T>
    H5Write vect( std::string name, T &v, int size, hid_t type, hsize_t offset=0, hsize_t npoints=0 )
    {
        // create dataspace for 1D array with good number of elements
        hsize_t dim = size;
//        if( deflate>0 ) {
//            H5Pset_chunk( dcr_, 1, &dim );
//            H5Pset_deflate( dcr_, std::min( 9, deflate ) );
//        } else {
//            H5Premove_filter( dcr_, H5Z_FILTER_DEFLATE );
//        }
        // Select portion
        if( npoints == 0 ) {
            npoints = dim - offset;
        }
        hid_t memspace = H5Screate_simple( 1, &npoints, NULL );
        hid_t filespace = H5Screate_simple( 1, &dim, NULL );
        if( offset > 0 || npoints < dim ) {
            hsize_t o = offset;
            hsize_t c = 1;
            hsize_t n = npoints;
            H5Sselect_hyperslab( filespace, H5S_SELECT_SET, &o, NULL, &c, &n );
        }
        // create dataset
        hid_t did = H5Dcreate( id_, name.c_str(), type, filespace, H5P_DEFAULT, dcr_, H5P_DEFAULT );
        // write vector in dataset
        H5Dwrite( did, type, memspace, filespace, dxpl_, &v );
        // close all
        H5Sclose( filespace );
        H5Sclose( memspace );
        return H5Write( did, dcr_, dxpl_ );
    }
    
    //! Create or open (not write) a dataset
    H5Write dataset( std::string name, hid_t type, H5Space *filespace )
    {
        return H5Write( this, name, type, filespace );
    }
    
    // Write to an open dataset
    template<class T>
    void write( T &v, hid_t type, H5Space *filespace, H5Space *memspace, bool independent = false ) {
        if( independent ) {
            if( memspace->global_ > 0 ) {
                H5FD_mpio_xfer_t xfer;
                H5Pget_dxpl_mpio( dxpl_, &xfer );
                H5Pset_dxpl_mpio( dxpl_, H5FD_MPIO_INDEPENDENT );
                H5Dwrite( id_, type, memspace->sid, filespace->sid, dxpl_, &v );
                H5Pset_dxpl_mpio( dxpl_, xfer );
            }
        } else {
            if( filespace->global_ > 0 ) {
                H5Dwrite( id_, type, memspace->sid, filespace->sid, dxpl_, &v );
            }
        }
    }
    
    // read to an open dataset
    template<class T>
    void read( T &v, hid_t type, H5Space *filespace, H5Space *memspace, bool independent = false ) {
        if( independent ) {
            if( memspace->global_ > 0 ) {
                H5FD_mpio_xfer_t xfer;
                H5Pget_dxpl_mpio( dxpl_, &xfer );
                H5Pset_dxpl_mpio( dxpl_, H5FD_MPIO_INDEPENDENT );
                H5Dread( id_, type, memspace->sid, filespace->sid, dxpl_, &v );
                H5Pset_dxpl_mpio( dxpl_, xfer );
            }
        } else {
            if( filespace->global_ > 0 ) {
                H5Dread( id_, type, memspace->sid, filespace->sid, dxpl_, &v );
            }
        }
    }
    
    //! Write a multi-dimensional array of uints
    H5Write array( std::string name, unsigned int &v, H5Space *filespace, H5Space *memspace, bool independent = false )
    {
        return array( name, v, H5T_NATIVE_UINT, filespace, memspace );
    }
    
    //! Write a multi-dimensional array of doubles
    H5Write array( std::string name, double &v, H5Space *filespace, H5Space *memspace, bool independent = false )
    {
        return array( name, v, H5T_NATIVE_DOUBLE, filespace, memspace );
    }
    
    //! Write a multi-dimensional array
    template<class T>
    H5Write array( std::string name, T &v, hid_t type, H5Space *filespace, H5Space *memspace, bool independent = false )
    {
        H5Write d = dataset( name, type, filespace );
        d.write( v, type, filespace, memspace, independent );
        return d;
    }
    
};

class H5Read : public H5
{
public:
    H5Read() : H5() {};
    
    //! Open HDF5 file + location
    H5Read( std::string file, MPI_Comm * comm = NULL, bool _raise = true )
     : H5( file, H5F_ACC_RDONLY, comm, _raise ) {};
    
    //! Location already opened
    H5Read( hid_t id, hid_t dcr, hid_t dxpl ) : H5( id, dcr, dxpl ) {};
    
    ~H5Read() {};
    
    void init( std::string file, MPI_Comm * comm = NULL, bool _raise = true )
    {
        H5::init( file, H5F_ACC_RDONLY, comm, _raise );
    }
    
    //! Open group
    H5Read group( std::string group_name )
    {
        return H5Read( open( group_name ), dcr_, dxpl_ );
    }
    
    //! Open an existing dataset
    H5Read dataset( std::string name )
    {
        return H5Read( open( name ), dcr_, dxpl_ );
    }
    
    //! Retrieve a double attribute
    void attr( std::string attribute_name, double &attribute_value )
    {
        attr( attribute_name, attribute_value, H5T_NATIVE_DOUBLE );
    }
    
    //! retrieve a unsigned int attribute
    void attr( std::string attribute_name, unsigned int &attribute_value )
    {
        attr( attribute_name, attribute_value, H5T_NATIVE_UINT );
    }
    
    //! retrieve a int attribute
    void attr( std::string attribute_name, int &attribute_value )
    {
        attr( attribute_name, attribute_value, H5T_NATIVE_INT );
    }
    
    //! retrieve a string attribute (specialized)
    void attr( std::string attribute_name, std::string &attribute_value )
    {
        if( H5Aexists( id_, attribute_name.c_str() )>0 ) {
            hid_t aid = H5Aopen_name( id_, attribute_name.c_str() );
            hid_t attr_type = H5Aget_type( aid );
            int sdim = H5Tget_size( attr_type );
            hid_t mem_type = H5Tcopy( H5T_C_S1 );
            H5Tset_size( mem_type, sdim );
            std::vector<char> tmpchar( sdim );
            // line below would crash (don't know why)
            // char* tmpchar= new char(sdim);
            if( H5Aread( aid, mem_type, &tmpchar[0] ) < 0 ) {
                WARNING( "Can't read string "<< attribute_name );
            } else {
                attribute_value = std::string( tmpchar.begin(), tmpchar.end() );
            }
            DEBUG( attribute_name << " >" << attribute_value << "< " << sdim );
            H5Tclose( mem_type );
            H5Tclose( attr_type );
            H5Aclose( aid );
        } else {
            WARNING( "Cannot find attribute " << attribute_name );
        }
    }
    
    //! retrieve anything (but string) as an attribute
    template<class T>
    void attr( std::string attribute_name, T &attribute_value, hid_t type )
    {
        if( H5Aexists( id_, attribute_name.c_str() )>0 ) {
            hid_t aid = H5Aopen( id_, attribute_name.c_str(), H5P_DEFAULT );
            H5Aread( aid, type, &( attribute_value ) );
            H5Aclose( aid );
        } else {
            WARNING( "Cannot find attribute " << attribute_name );
        }
    }
    
    //! retrieve a vector<double> as an attribute
    void attr( std::string attribute_name, std::vector<double> &attribute_value )
    {
        attr( attribute_name, attribute_value, H5T_NATIVE_DOUBLE );
    }
    
    //! retrieve a vector<anything> as an attribute
    template<class T>
    void attr( std::string attribute_name, std::vector<T> &attribute_value, hid_t type )
    {
        if( H5Aexists( id_, attribute_name.c_str() )>0 ) {
            hid_t aid = H5Aopen( id_, attribute_name.c_str(), H5P_DEFAULT );
            hid_t sid = H5Aget_space( aid );
            hssize_t npoints = H5Sget_simple_extent_npoints( sid );
            attribute_value.resize( npoints );
            H5Sclose( sid );
            H5Aread( aid, type, &( attribute_value[0] ) );
            H5Aclose( aid );
        } else {
            WARNING( "Cannot find attribute " << attribute_name );
        }
    }
    
    //! read attribute size (when vector)
    int attrSize( std::string attribute_name )
    {
        if( H5Aexists( id_, attribute_name.c_str() )>0 ) {
            hid_t aid = H5Aopen( id_, attribute_name.c_str(), H5P_DEFAULT );
            if( aid < 0 ) {
                return -1;
            }
            hid_t sid = H5Aget_space( aid );
            hssize_t npoints = H5Sget_simple_extent_npoints( sid );
            H5Sclose( sid );
            H5Aclose( aid );
            return npoints;
        } else {
            return -1;
        }
    }
    
    //! retrieve a double vector
    void vect( std::string vect_name,  std::vector<double> &v, bool resizeVect=false, hsize_t offset=0, hsize_t npoints=0 )
    {
        vect( vect_name, v, H5T_NATIVE_DOUBLE, resizeVect, offset, npoints );
    }
    
    //! retrieve an unsigned int vector
    void vect( std::string vect_name,  std::vector<unsigned int> &v, bool resizeVect=false, hsize_t offset=0, hsize_t npoints=0 )
    {
        vect( vect_name, v, H5T_NATIVE_UINT, resizeVect, offset, npoints );
    }
    
    //! retrieve a int vector
    void vect( std::string vect_name,  std::vector<int> &v, bool resizeVect=false, hsize_t offset=0, hsize_t npoints=0 )
    {
        vect( vect_name, v, H5T_NATIVE_INT, resizeVect, offset, npoints );
    }
    
    //! retrieve a short vector
    void vect( std::string vect_name,  std::vector<short> &v, bool resizeVect=false, hsize_t offset=0, hsize_t npoints=0 )
    {
        vect( vect_name, v, H5T_NATIVE_SHORT, resizeVect, offset, npoints );
    }
    
    //! template to read generic 1d vector (optionally offset and npoints)
    template<class T>
    void vect( std::string vect_name, std::vector<T> &v, hid_t type, bool resizeVect=false, hsize_t offset=0, hsize_t npoints=0 )
    {
        if( resizeVect ) {
            std::vector<hsize_t> s = shape( vect_name );
            hsize_t n = 1;
            for( unsigned int i=0; i<s.size(); i++ ) {
                n *= s[i];
            }
            if( npoints == 0 ) {
                npoints = n - offset ;
            }
            v.resize( npoints );
        }
        vect( vect_name, v[0], type, offset, npoints );
    }
    
    //! template to read generic 1d data (optionally offset and npoints)
    template<class T>
    void vect( std::string vect_name, T &v, hid_t type, hsize_t offset=0, hsize_t npoints=0 )
    {
        // Open dataset
        hid_t did = H5Dopen( id_, vect_name.c_str(), H5P_DEFAULT );
        if( did < 0 ) {
            ERROR( "Cannot read dataset " << vect_name );
        }
        if( offset != 0 || npoints != 0 ) {
            // Get shape and resize vector
            hid_t sid = H5Dget_space( did );
            int sdim = H5Sget_simple_extent_ndims( sid );
            if( sdim!=1 ) {
                ERROR( "Reading vector " << vect_name << " is not 1D but " <<sdim << "D" );
            }
            hsize_t dim = H5Sget_simple_extent_npoints( sid );
            H5Sclose( sid );
            if( npoints == 0 ) {
                npoints = dim - offset;
            }
            // Select portion
            hid_t memspace = H5Screate_simple( 1, &npoints, NULL );
            hid_t filespace = H5Screate_simple( 1, &dim, NULL );
            hsize_t o = offset;
            hsize_t c = 1;
            hsize_t n = npoints;
            H5Sselect_hyperslab( filespace, H5S_SELECT_SET, &o, NULL, &c, &n );
            // Read data
            H5Dread( did, type, memspace, filespace, dxpl_, &v );
            H5Sclose( filespace );
            H5Sclose( memspace );
        } else {
            H5Dread( did, type, H5S_ALL, H5S_ALL, dxpl_, &v );
        }
        H5Dclose( did );
    }
    
    //! Read a multi-dimensional array of uints
    void array( std::string name, unsigned int &v, H5Space *filespace, H5Space *memspace )
    {
        return array( name, v, H5T_NATIVE_UINT, filespace, memspace );
    }
    
    //! Read a multi-dimensional array of doubles
    void array( std::string name, double &v, H5Space *filespace, H5Space *memspace )
    {
        return array( name, v, H5T_NATIVE_DOUBLE, filespace, memspace );
    }
    
    //! Read a multi-dimensional array
    template<class T>
    void array( std::string name, T &v, hid_t type, H5Space *filespace, H5Space *memspace )
    {
        hid_t did = open( name );
        if( filespace->global_ > 0 ) {
            H5Dread( did, type, memspace->sid, filespace->sid, dxpl_, &v );
        }
        H5Dclose( did );
    }
    
    std::vector<hsize_t> shape( std::string name )
    {
        std::vector<hsize_t> shape( 0 );
        if( H5Lexists( id_, name.c_str(), H5P_DEFAULT ) >0 ) {
            hid_t did = H5Dopen( id_, name.c_str(), H5P_DEFAULT );
            if( did >= 0 ) {
                hid_t sid = H5Dget_space( did );
                int sdim = H5Sget_simple_extent_ndims( sid );
                shape.resize( sdim );
                H5Sget_simple_extent_dims( sid, &shape[0], NULL );
                H5Sclose( sid );
                H5Dclose( did );
            }
        }
        return shape;
    }
    
    int vectSize( std::string vect_name )
    {
        std::vector<hsize_t> s = shape( vect_name );
        if( s.size() == 0 ) {
            return -1;
        }
        if( s.size() != 1 ) {
            ERROR( "Reading vector " << vect_name << " is not 1D but " << s.size() << "D" );
        }
        return s[0];
    }
};


#endif
