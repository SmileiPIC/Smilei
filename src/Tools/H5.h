
// -------------------
// Some HDF5 overlays
// -------------------

#ifndef H5_H
#define H5_H

#include <hdf5.h>


class H5 {
    
    public:
    
    // write a string as an attribute
    static void attr(hid_t locationId, std::string attribute_name, std::string attribute_value) {
        hid_t sid  = H5Screate(H5S_SCALAR);
        hid_t atype = H5Tcopy(H5T_C_S1);
        H5Tset_size(atype, attribute_value.size());
        H5Tset_strpad(atype,H5T_STR_NULLTERM);
        hid_t aid = H5Acreate(locationId, attribute_name.c_str(), atype, sid, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(aid, atype, attribute_value.c_str());
        H5Aclose(aid);
        H5Sclose(sid);
        H5Tclose(atype);
    }
    
    // write an unsigned int as an attribute
    static void attr(hid_t locationId, std::string attribute_name, unsigned int attribute_value) {
        hid_t sid = H5Screate(H5S_SCALAR);
        hid_t aid = H5Acreate(locationId, attribute_name.c_str(), H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(aid, H5T_NATIVE_UINT, &attribute_value);
        H5Sclose(sid);
        H5Aclose(aid);
    }
    
    // write an int as an attribute
    static void attr(hid_t locationId, std::string attribute_name, int attribute_value) {
        hid_t sid = H5Screate(H5S_SCALAR);
        hid_t aid = H5Acreate(locationId, attribute_name.c_str(), H5T_NATIVE_INT, sid, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(aid, H5T_NATIVE_INT, &attribute_value);
        H5Sclose(sid);
        H5Aclose(aid);
    }
    
    // write a double as an attribute
    static void attr(hid_t locationId, std::string attribute_name, double attribute_value) {
        hid_t sid = H5Screate(H5S_SCALAR);
        hid_t aid = H5Acreate(locationId, attribute_name.c_str(), H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(aid, H5T_NATIVE_DOUBLE, &attribute_value);
        H5Sclose(sid);
        H5Aclose(aid);
    }
    
    // write a vector of doubles
    // v is the vector
    // size is the number of elements in the vector
    static void vector(hid_t locationId, std::string name, double& v, int size) {
        // create dataspace for 1D array with good number of elements
        hsize_t dims = size;
        hid_t sid = H5Screate_simple(1, &dims, NULL);
        hid_t pid = H5Pcreate(H5P_DATASET_CREATE); // property list
        // create dataset 
        hid_t did = H5Dcreate(locationId, name.c_str(), H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, pid, H5P_DEFAULT);
        // write vector in dataset
        H5Dwrite(did, H5T_NATIVE_DOUBLE, sid, sid, H5P_DEFAULT, &v);
        // close all
        H5Dclose(did);
        H5Pclose(pid);
        H5Sclose(sid);
    }
    
    // write a 2-D array in parallel (several MPI nodes)
    // m is the matrix (2D array)
    // sizex, sizey is the number of elements in both axes of the matrix
    // offset is the x-location where the current node will start to write
    // numel  is the x-number of elements for the current node
    static void matrix_MPI(hid_t locationId, std::string name, double& m,
                           int sizex, int sizey, int offset, int numel    ) {
        // Create a HDF5 memory space to hold the data
        hsize_t chunk_parts[2];
        chunk_parts[0] = numel;
        chunk_parts[1] = sizey;
        hid_t memspace = H5Screate_simple(2, chunk_parts, NULL);
        // Create the corresponding HDF5 filespace
        hsize_t dimsf[2], offs[2], stride[2], count[2];
        dimsf[0] = sizex;
        dimsf[1] = sizey;
        hid_t filespace = H5Screate_simple(2, dimsf, NULL);
        offs[0] = offset;
        offs[1] = 0;
        stride[0] = 1;
        stride[1] = 1;
        count[0] = 1;
        count[1] = 1;
        // Choose the hyperslab, which is the region where the current node will write
        hsize_t block[2];
        block[0] = numel;
        block[1] = sizey;
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