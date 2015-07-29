
// -------------------
// Some HDF5 overlays
// -------------------

#ifndef H5_H
#define H5_H

#include <hdf5.h>
#include <string>

//! HDF5 help functions
class H5 {
    
    public:
    
    //! Make an empty group
    // Returns the group ID
    static hid_t group(hid_t locationId, std::string group_name) {
        return H5Gcreate(locationId, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    
    //! write a string as an attribute
    static void attr(hid_t locationId, std::string attribute_name, std::string attribute_value) {
        hid_t atype = H5Tcopy(H5T_C_S1);
        H5Tset_size(atype, attribute_value.size());
        H5Tset_strpad(atype,H5T_STR_NULLTERM);
        
        attr(locationId, attribute_name, *(attribute_value.c_str()), atype);
        
        H5Tclose(atype);
    }
    
    //! write an unsigned int as an attribute
    static void attr(hid_t locationId, std::string attribute_name, unsigned int attribute_value) {
        attr(locationId, attribute_name, attribute_value, H5T_NATIVE_UINT);}
    
    //! write an int as an attribute
    static void attr(hid_t locationId, std::string attribute_name, int attribute_value) {
        attr(locationId, attribute_name, attribute_value, H5T_NATIVE_INT);}
    
    //! write a double as an attribute
    static void attr(hid_t locationId, std::string attribute_name, double attribute_value) {
        attr(locationId, attribute_name, attribute_value, H5T_NATIVE_DOUBLE);}
    
    //! write anything as an attribute
    template<class T>
    static void attr(hid_t locationId, std::string attribute_name, T & attribute_value, hid_t type) {
        hid_t sid = H5Screate(H5S_SCALAR);
        hid_t aid = H5Acreate(locationId, attribute_name.c_str(), type, sid, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(aid, type, &attribute_value);
        H5Sclose(sid);
        H5Aclose(aid);
    }
    
    //! write a vector of unsigned ints
    //! v is the vector
    //! size is the number of elements in the vector
    static void vector(hid_t locationId, std::string name, unsigned int& v, int size) {
        vector(locationId, name, v, size, H5T_NATIVE_UINT);
    }
    
    //! write a vector of ints
    static void vector(hid_t locationId, std::string name, int& v, int size) {
        vector(locationId, name, v, size, H5T_NATIVE_INT);
    }
    
    //! write a vector of doubles
    static void vector(hid_t locationId, std::string name, double& v, int size) {
        vector(locationId, name, v, size, H5T_NATIVE_DOUBLE);
    }
    
    //! write any vector
    //! type is the h5 type (H5T_NATIVE_DOUBLE, H5T_NATIVE_INT, etc.)
    template<class T>
    static void vector(hid_t locationId, std::string name, T & v, int size, hid_t type) {
        // create dataspace for 1D array with good number of elements
        hsize_t dims = size;
        hid_t sid = H5Screate_simple(1, &dims, NULL);
        hid_t pid = H5Pcreate(H5P_DATASET_CREATE); // property list
        // create dataset 
        hid_t did = H5Dcreate(locationId, name.c_str(), type, sid, H5P_DEFAULT, pid, H5P_DEFAULT);
        // write vector in dataset
        H5Dwrite(did, type, sid, sid, H5P_DEFAULT, &v);
        // close all
        H5Dclose(did);
        H5Pclose(pid);
        H5Sclose(sid);
    }
    
    
    
    //! write a 2-D array of doubles in parallel (several MPI nodes)
    //! m is the matrix (2D array)
    //! sizex, sizey is the number of elements in both axes of the matrix
    //! offset is the x-location where the current node will start to write
    //! numel  is the x-number of elements for the current node
    static void matrix_MPI(hid_t locationId, std::string name, double& m,
                           int sizex, int sizey, int offset, int numel    ) {
        // Create a HDF5 memory space to hold the data
        hsize_t chunk_parts[2];
        chunk_parts[0] = numel;
        chunk_parts[1] = sizey;
        hid_t memspace = H5Screate_simple(2, chunk_parts, NULL);
        // Create the HDF5 filespace
        hsize_t dimsf[2];
        dimsf[1] = sizex;
        dimsf[0] = sizey;
        hid_t filespace = H5Screate_simple(2, dimsf, NULL);
        // Choose the hyperslab, which is the region where the current node will write
        hsize_t offs[2], stride[2], count[2], block[2];
        offs[1] = offset;
        offs[0] = 0;
        stride[0] = 1;
        stride[1] = 1;
        count[0] = 1;
        count[1] = 1;
        block[1] = numel;
        block[0] = sizey;
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offs, stride, count, block);
        // Open the pre-existing group and write the data inside
        hid_t write_plist = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_INDEPENDENT);
        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        hid_t dset_id  = H5Dcreate(locationId, name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
        H5Pclose(plist_id);
        H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, write_plist, &m );
        H5Dclose(dset_id);
        H5Pclose( write_plist );
        H5Sclose(filespace);
        H5Sclose(memspace);
    }
    
};

#endif
