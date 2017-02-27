
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
        H5Tset_size(atype, H5T_VARIABLE);
        attr(locationId, attribute_name, attribute_value, atype);
        H5Tclose(atype);
    }
    
    //! write an unsigned int as an attribute
    static void attr(hid_t locationId, std::string attribute_name, unsigned int attribute_value) {
        attr(locationId, attribute_name, attribute_value, H5T_NATIVE_UINT);}
    
    //! write an uint32 as an attribute
    static void attr_uint32(hid_t locationId, std::string attribute_name, uint32_t attribute_value) {
        attr(locationId, attribute_name, attribute_value, H5T_NATIVE_UINT32);}
    
    //! write size_t as an attribute
    static void attr(hid_t locationId, std::string attribute_name, size_t attribute_value) {
        attr(locationId, attribute_name, (unsigned int) attribute_value);}
    
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
    
    
    //! write a vector<unsigned int> as an attribute
    static void attr(hid_t locationId, std::string attribute_name, std::vector<unsigned int> attribute_value) {
        attr(locationId, attribute_name, attribute_value, H5T_NATIVE_UINT);
    }
    
    //! write a vector<double> as an attribute
    static void attr(hid_t locationId, std::string attribute_name, std::vector<double> attribute_value) {
        attr(locationId, attribute_name, attribute_value, H5T_NATIVE_DOUBLE);
    }
    
    //! write a vector<string> as an attribute
    static void attr(hid_t locationId, std::string attribute_name, std::vector<std::string> attribute_value) {
        hid_t atype = H5Tcopy(H5T_C_S1);
        H5Tset_size(atype, H5T_VARIABLE);
        std::vector<char*> tmp_vec(attribute_value.size(),nullptr);
        for (unsigned int i=0; i<attribute_value.size();i++) {
            tmp_vec[i] = &attribute_value[i][0];
        }
        attr(locationId, attribute_name, tmp_vec, atype);
        H5Tclose(atype);

    }
    
    //! write a vector<anything> as an attribute
    template<class T>
    static void attr(hid_t locationId, std::string attribute_name, std::vector<T>& attribute_value, hid_t type) {
        hsize_t dims = attribute_value.size();
        hid_t sid = H5Screate_simple(1, &dims, NULL);
        hid_t aid = H5Acreate (locationId, attribute_name.c_str(), type, sid, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(aid, type, &(attribute_value[0]));
        H5Aclose(aid);
        H5Sclose(sid);
    }
    
    
    //READ ATTRIBUTES
    
    //! retrieve a double attribute
    static void getAttr(hid_t locationId, std::string attribute_name, double &attribute_value) {
        getAttr(locationId, attribute_name, attribute_value, H5T_NATIVE_DOUBLE);
    }
    
    //! retrieve a unsigned int attribute
    static void getAttr(hid_t locationId, std::string attribute_name, unsigned int &attribute_value) {
        getAttr(locationId, attribute_name, attribute_value, H5T_NATIVE_UINT);
    }
    
    //! retrieve a int attribute
    static void getAttr(hid_t locationId, std::string attribute_name, int &attribute_value) {
        getAttr(locationId, attribute_name, attribute_value, H5T_NATIVE_INT);
    }
    
    //! retrieve a string attribute
    static void getAttr(hid_t locationId, std::string attribute_name, std::string &attribute_value) {
        if (H5Aexists(locationId,attribute_name.c_str())>0) {
            hid_t attr_id = H5Aopen_name(locationId, attribute_name.c_str());
            hid_t attr_type = H5Aget_type(attr_id);
            int sdim = H5Tget_size(attr_type);
            hid_t mem_type = H5Tcopy(H5T_C_S1);
            H5Tset_size(mem_type, sdim);
            std::vector<char> tmpchar(sdim);
            // line below would crash (don't know why)
            // char* tmpchar= new char(sdim);
            if (H5Aread(attr_id, mem_type, &tmpchar[0]) < 0) {
                WARNING("Can't read string "<< attribute_name);
            } else {
                attribute_value = std::string(tmpchar.begin(),tmpchar.end());
            }
            H5Tclose(mem_type);
            H5Tclose(attr_type);
            H5Aclose(attr_id);
        } else {
            WARNING("Cannot find attribute " << attribute_name);
        }
    }
    
    //! retrieve anything as an attribute
    template<class T>
    static void getAttr(hid_t locationId, std::string attribute_name, T &attribute_value, hid_t type) {
        if (H5Aexists(locationId,attribute_name.c_str())>0) {
            hid_t aid = H5Aopen(locationId, attribute_name.c_str(), H5P_DEFAULT);
            H5Aread(aid, type, &(attribute_value));
            H5Aclose(aid);
        } else {
            WARNING("Cannot find attribute " << attribute_name);
        }
    }
    
    //! retrieve a vector<double> as an attribute
    template<class T>
    static void getAttr(hid_t locationId, std::string attribute_name, std::vector<T>& attribute_value) {
        getAttr(locationId, attribute_name, attribute_value, H5T_NATIVE_DOUBLE);
    }
    
    //! retrieve a vector<anything> as an attribute
    template<class T>
    static void getAttr(hid_t locationId, std::string attribute_name, std::vector<T>& attribute_value, hid_t type) {
        if (H5Aexists(locationId,attribute_name.c_str())>0) {
            hid_t aid = H5Aopen(locationId, attribute_name.c_str(), H5P_DEFAULT);
            hid_t sid = H5Aget_space(aid);
            hssize_t npoints = H5Sget_simple_extent_npoints( sid );
            attribute_value.resize( npoints );
            H5Sclose( sid );
            H5Aread(aid, type, &(attribute_value[0]));
            H5Aclose(aid);
        } else {
            WARNING("Cannot find attribute " << attribute_name);
        }
    }
    
    
    //! write a vector of unsigned ints
    //! v is the vector
    //! size is the number of elements in the vector
    
    //! write a vector<int>
    static void vect(hid_t locationId, std::string name, std::vector<int> v, int deflate=0) {
        vect(locationId, name, v[0], v.size(), H5T_NATIVE_INT, deflate);
    }
    
    //! write a vector<unsigned int>
    static void vect(hid_t locationId, std::string name, std::vector<unsigned int> v, int deflate=0) {
        vect(locationId, name, v[0], v.size(), H5T_NATIVE_UINT, deflate);
    }
    
    //! write a vector<short>
    static void vect(hid_t locationId, std::string name, std::vector<short> v, int deflate=0) {
        vect(locationId, name, v[0], v.size(), H5T_NATIVE_SHORT, deflate);
    }
    
    //! write a vector<doubles>
    static void vect(hid_t locationId, std::string name, std::vector<double> v, int deflate=0) {
        vect(locationId, name, v[0], v.size(), H5T_NATIVE_DOUBLE, deflate);
    }
    
    
    //! write any vector
    //! type is the h5 type (H5T_NATIVE_DOUBLE, H5T_NATIVE_INT, etc.)
    template<class T>
    static void vect(hid_t locationId, std::string name, T & v, int size, hid_t type, int deflate=0) {
        // create dataspace for 1D array with good number of elements
        hsize_t dims = size;
        hid_t sid = H5Screate_simple(1, &dims, NULL);
        hid_t pid = H5Pcreate(H5P_DATASET_CREATE); // property list
        
        if (deflate>0) {
            H5Pset_chunk(pid, 1, &dims);
            H5Pset_deflate (pid, std::min(9,deflate));
        }
        
        // create dataset
        hid_t did = H5Dcreate(locationId, name.c_str(), type, sid, H5P_DEFAULT, pid, H5P_DEFAULT);
        // write vector in dataset
        H5Dwrite(did, type, sid, sid, H5P_DEFAULT, &v);
        // close all
        H5Dclose(did);
        H5Pclose(pid);
        H5Sclose(sid);
    }
    
    
    
    //! retrieve a double vector
    static void getVect(hid_t locationId, std::string vect_name,  std::vector<double> &vect, bool resizeVect=false) {
        getVect(locationId, vect_name, vect, H5T_NATIVE_DOUBLE,resizeVect);
    }
    
    //! retrieve an unsigned int vector
    static void getVect(hid_t locationId, std::string vect_name,  std::vector<unsigned int> &vect, bool resizeVect=false) {
        getVect(locationId, vect_name, vect, H5T_NATIVE_UINT,resizeVect);
    }
    
    //! retrieve a int vector
    static void getVect(hid_t locationId, std::string vect_name,  std::vector<int> &vect, bool resizeVect=false) {
        getVect(locationId, vect_name, vect, H5T_NATIVE_INT,resizeVect);
    }
    
    //! retrieve a short vector
    static void getVect(hid_t locationId, std::string vect_name,  std::vector<short> &vect, bool resizeVect=false) {
        getVect(locationId, vect_name, vect, H5T_NATIVE_SHORT,resizeVect);
    }
    
    //! template to read generic 1d vector
    template<class T>
    static void getVect(hid_t locationId, std::string vect_name, std::vector<T> &vect, hid_t type, bool resizeVect=false) {
        hid_t did = H5Dopen(locationId, vect_name.c_str(), H5P_DEFAULT);
        hid_t sid = H5Dget_space(did);
        int sdim = H5Sget_simple_extent_ndims(sid);
        if (sdim!=1) {
            ERROR("Reading vector " << vect_name << " is not 1D but " <<sdim << "D");
        }
        hsize_t dim[1];
        H5Sget_simple_extent_dims(sid,dim,NULL);
        if (dim[0] != vect.size()) {
            if (resizeVect) {
                vect.resize(dim[0]);
            } else {
                ERROR("Reading vector " << vect_name << " mismatch " << vect.size() << " != " << dim);
            }
        }
        H5Sclose(sid);
        H5Dread(did, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vect[0]);
        H5Dclose(did);
    }
    
};

#endif
