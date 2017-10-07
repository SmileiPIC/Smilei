#ifndef DIAGNOSTICTRACK_H
#define DIAGNOSTICTRACK_H

#include "Diagnostic.h"

class Patch;
class Params;
class SmileiMPI;


class DiagnosticTrack : public Diagnostic {

public :
    //! Default constructor
    DiagnosticTrack( Params &params, SmileiMPI* smpi, VectorPatch& vecPatches, unsigned int, unsigned int, OpenPMDparams& );
    //! Default destructor
    ~DiagnosticTrack() override;
    
    void openFile( Params& params, SmileiMPI* smpi, bool newfile ) override;
    
    void closeFile() override;
    
    void init(Params& params, SmileiMPI* smpi, VectorPatch& vecPatches) override;
    
    bool prepare( int itime ) override;
    
    void run( SmileiMPI* smpi, VectorPatch& vecPatches, int itime, SimWindow* simWindow ) override;
    
    //! Get memory footprint of current diagnostic
    int getMemFootPrint() override {
        return 0;
    }
    
    //! Fills a buffer with the required particle property
    template<typename T> void fill_buffer(VectorPatch& vecPatches, unsigned int iprop, std::vector<T>& buffer);
    
    //! Write a scalar dataset with the given buffer
    template<typename T> void write_scalar( hid_t, std::string, T&, hid_t, hid_t, hid_t, hid_t, unsigned int, unsigned int );
    
    //! Write a vector component dataset with the given buffer
    template<typename T> void write_component( hid_t, std::string, T&, hid_t, hid_t, hid_t, hid_t, unsigned int, unsigned int );
    
    //! Set a given patch's particles with the required IDs
    void setIDs(Patch *);
    
    //! Set a given particles with the required IDs
    void setIDs(Particles&);
    
    //! Index of the species used
    unsigned int speciesId_;
    
    //! Last ID assigned to a particle by this MPI domain
    uint64_t latest_Id;
    
private :
    
    //! Flag to test whether IDs have been set already
    bool IDs_done;
    
    //! HDF5 objects
    hid_t data_group_id, transfer;
     
    //! Number of spatial dimensions
    unsigned int nDim_particle;
    
    //! Current particle partition among the patches own by current MPI
    std::vector<unsigned int> patch_start;
    
    //! Tells whether this diag includes a particle filter
    bool has_filter;
    
    //! Tells whether this diag includes a particle filter
    PyObject* filter;
    
    //! Selection of the filtered particles in each patch
    std::vector<std::vector<unsigned int> > patch_selection;
    
    //! Buffer for the output of double array
    std::vector<double> data_double;
    //! Buffer for the output of short array
    std::vector<short> data_short;
    //! Buffer for the output of uint64 array
    std::vector<uint64_t> data_uint64;
};



#ifdef SMILEI_USE_NUMPY
#include <numpy/arrayobject.h>

// class for exposing Particles to numpy
class ParticleData {
public:
    ParticleData(unsigned int nparticles) {
       // Use the empty python class "Particles" to store data
       PyObject* particlesClass = PyObject_GetAttrString(PyImport_AddModule("__main__"),"Particles");
       particles = PyObject_CallObject(particlesClass, NULL);
       Py_DECREF(particlesClass);
       
       dims[0] = nparticles;
    };
    ~ParticleData() {
        clear();
    };
    
    inline void resize(unsigned int nparticles) {
        dims[0] = nparticles;
    };
    
    // Expose a vector to numpy
    inline PyArrayObject* vector2numpy( std::vector<double> &vec ) {
        return (PyArrayObject*) PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, (double*)(vec.data()));
    };
    inline PyArrayObject* vector2numpy( std::vector<uint64_t> &vec ) {
        return (PyArrayObject*) PyArray_SimpleNewFromData(1, dims, NPY_UINT64, (uint64_t*)(vec.data()));
    };
    inline PyArrayObject* vector2numpy( std::vector<short> &vec ) {
        return (PyArrayObject*) PyArray_SimpleNewFromData(1, dims, NPY_SHORT, (short*)(vec.data()));
    };
    
    // Add a C++ vector as an attribute, but exposed as a numpy array
    template <typename T>
    inline void setVectorAttr( std::vector<T> &vec, std::string name ) {
        PyArrayObject* numpy_vector = vector2numpy( vec );
        PyObject_SetAttrString(particles, name.c_str(), (PyObject*)numpy_vector);
        attrs.push_back( numpy_vector );
    };
    
    // Remove python references of all attributes
    inline void clear() {
        unsigned int nattr = attrs.size();
        for( unsigned int i=0; i<nattr; i++) Py_DECREF(attrs[i]);
        attrs.resize(0);
    };
    
    inline PyObject* get() {
        return particles;
    };
    
private:
    PyObject* particles;
    std::vector<PyArrayObject*> attrs;
    npy_intp dims[1];
};


// Set a python variable to itime so that it can be accessed in the filter
inline void setIteration( int itime ) {
    PyObject* Main = PyObject_GetAttrString(PyImport_AddModule("__main__"),"Main");
    PyObject* iteration = PyInt_FromLong(itime);
    PyObject_SetAttrString(Main, "iteration", iteration);
    Py_DECREF(iteration);
    Py_DECREF(Main);
}


#endif



#endif

