#ifndef H5File_H
#define H5File_H

#include <hdf5.h>
#include <H5.h>
#include <string>
#include <sstream>
#include <vector>
#include "Tools.h"

#if ! H5_HAVE_PARALLEL == 1
#error "HDF5 was not built with --enable-parallel option"
#endif

class H5Group
{
public:
    //! Group already opened
    H5Group( hid_t ID );
    
    //! Open group from location
    H5Group( hid_t lid, std::string grouppath );
    
    ~H5Group();
    
    bool valid() {
        return id >= 0;
    }
    
protected:
    hid_t id;
    bool empty;
};

class H5GroupWrite : public H5Group
{
public:
    H5GroupWrite( hid_t id ) : H5Group( id ) {};
    H5GroupWrite( hid_t lid, std::string grouppath ) : H5Group( lid, grouppath ) {};
    
    //! Make an empty group
    H5GroupWrite group( std::string group_name )
    {
        return H5GroupWrite( H5Gcreate( id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT ) );
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
        hid_t aid = H5Acreate( id, attribute_name.c_str(), type, sid, H5P_DEFAULT, H5P_DEFAULT );
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
        hid_t aid = H5Acreate( id, attribute_name.c_str(), type, sid, H5P_DEFAULT, H5P_DEFAULT );
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
        hid_t aid = H5Acreate( id, attribute_name.c_str(), atype, sid, H5P_DEFAULT, H5P_DEFAULT );
        H5Awrite( aid, atype, attribute_value.c_str() );
        H5Tclose( atype );
        H5Aclose( aid );
        H5Sclose( sid );
    }
    
    //! Write a vector<int>
    void vect( std::string name, std::vector<int> v, int deflate=0 )
    {
        vect( name, v[0], v.size(), H5T_NATIVE_INT, deflate );
    }
    
    //! write a vector<unsigned int>
    void vect(std::string name, std::vector<unsigned int> v, int deflate=0 )
    {
        vect(name, v[0], v.size(), H5T_NATIVE_UINT, deflate );
    }
    
    //! write a vector<short>
    void vect( std::string name, std::vector<short> v, int deflate=0 )
    {
        vect( name, v[0], v.size(), H5T_NATIVE_SHORT, deflate );
    }
    
    //! write a vector<doubles>
    void vect( std::string name, std::vector<double> v, int deflate=0 )
    {
        vect( name, v[0], v.size(), H5T_NATIVE_DOUBLE, deflate );
    }
    
    //! write any vector
    template<class T>
    void vect( std::string name, std::vector<T> v, hid_t type, int deflate=0 )
    {
        vect( name, v[0], v.size(), type, deflate );
    }
    
    //! Write a portion of a vector
    template<class T>
    void vect( std::string name, T &v, int size, hid_t type, int deflate=0 )
    {
        // create dataspace for 1D array with good number of elements
        hsize_t dims = size;
        hid_t sid = H5Screate_simple( 1, &dims, NULL );
        hid_t pid = H5Pcreate( H5P_DATASET_CREATE );
        
        if( deflate>0 ) {
            H5Pset_chunk( pid, 1, &dims );
            H5Pset_deflate( pid, std::min( 9, deflate ) );
        }
        
        // create dataset
        hid_t did = H5Dcreate( id, name.c_str(), type, sid, H5P_DEFAULT, pid, H5P_DEFAULT );
        // write vector in dataset
        H5Dwrite( did, type, sid, sid, H5P_DEFAULT, &v );
        // close all
        H5Dclose( did );
        H5Pclose( pid );
        H5Sclose( sid );
    }
};

class H5GroupRead : public H5Group
{
public:
    H5GroupRead( hid_t id ) : H5Group( id ) {};
    H5GroupRead( hid_t lid, std::string grouppath ) : H5Group( lid, grouppath ) {};
    
    //! Open group
    H5GroupRead group( std::string group_name )
    {
        return H5GroupRead( H5Gopen( id, group_name.c_str(), H5P_DEFAULT ) );
    }
    
    //Check if attribute exists
    bool hasAttr( std::string attribute_name )
    {
        return H5Aexists( id, attribute_name.c_str() ) > 0;
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
        if( H5Aexists( id, attribute_name.c_str() )>0 ) {
            hid_t aid = H5Aopen_name( id, attribute_name.c_str() );
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
        if( H5Aexists( id, attribute_name.c_str() )>0 ) {
            hid_t aid = H5Aopen( id, attribute_name.c_str(), H5P_DEFAULT );
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
        if( H5Aexists( id, attribute_name.c_str() )>0 ) {
            hid_t aid = H5Aopen( id, attribute_name.c_str(), H5P_DEFAULT );
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
        if( H5Aexists( id, attribute_name.c_str() )>0 ) {
            hid_t aid = H5Aopen( id, attribute_name.c_str(), H5P_DEFAULT );
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
        // Open dataset
        hid_t did = H5Dopen( id, vect_name.c_str(), H5P_DEFAULT );
        if( did < 0 ) {
            ERROR( "Cannot read dataset " << vect_name );
        }
        // Get shape and resize vector
        hid_t sid = H5Dget_space( did );
        int sdim = H5Sget_simple_extent_ndims( sid );
        if( sdim!=1 ) {
            ERROR( "Reading vector " << vect_name << " is not 1D but " <<sdim << "D" );
        }
        hsize_t dim = H5Sget_simple_extent_npoints( sid );
        H5Sclose( sid );
        if( npoints == 0 ) {
            npoints = std::max( (hsize_t)0, dim - offset );
        }
        if( npoints != v.size() ) {
            if( resizeVect ) {
                v.resize( npoints );
            } else {
                ERROR( "Reading vector " << vect_name << " mismatch " << v.size() << " != " << npoints );
            }
        }
        // Select portion
        hid_t memspace = H5Screate_simple( 1, &npoints, NULL );
        hid_t filespace = H5Screate_simple( 1, &dim, NULL );
        if( offset > 0 || npoints < dim ) {
            hsize_t o = offset;
            hsize_t c = 1;
            hsize_t n = npoints;
            H5Sselect_hyperslab( filespace, H5S_SELECT_SET, &o, NULL, &c, &n );
        }
        // Read data
        H5Dread( did, type, memspace, filespace, H5P_DEFAULT, &v[0] );
        H5Dclose( did );
        H5Sclose( filespace );
        H5Sclose( memspace );
    }
    
    std::vector<hsize_t> shape( std::string name )
    {
        std::vector<hsize_t> shape( 0 );
        if( H5Lexists( id, name.c_str(), H5P_DEFAULT ) >0 ) {
            hid_t did = H5Dopen( id, name.c_str(), H5P_DEFAULT );
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


class H5File
{
public:
    H5File( std::string file, unsigned access, bool _raise );
    ~H5File() {
        if( fid >= 0 ) {
            H5Fclose( fid );
        }
    };
protected:
    std::string filepath;
    std::string grouppath;
    hid_t fid;
};

class H5FileWrite : public H5File, public H5GroupWrite
{
public:
    H5FileWrite( std::string file, bool _raise = true ) : H5File( file, H5F_ACC_RDWR, _raise ), H5GroupWrite( fid, grouppath ) {};
    ~H5FileWrite() {};
};

class H5FileRead : public H5File, public H5GroupRead
{
public:
    H5FileRead( std::string file, bool _raise = true ) : H5File( file, H5F_ACC_RDONLY, _raise ), H5GroupRead( fid, grouppath ) {};
    ~H5FileRead() {};
};

#endif
