
#ifndef LaserPropagator_H
#define LaserPropagator_H

#include "PyTools.h"
#include "Field2D.h"
#include "SmileiMPI.h"

#include <vector>
#include <string>
#include <cmath>

class Params;

class LaserPropagator {
public:
    LaserPropagator() {};
    ~LaserPropagator() {};
    
    void init(Params* params, SmileiMPI* smpi, unsigned int side);
    
    // Propagates the fields profiles with some offset, and writes result to file
    unsigned int operator() (std::vector<PyObject*> profiles, std::vector<int> profiles_n, double offset, std::string file);
    
protected:
    
private:
    unsigned int ndim;
    unsigned int MPI_size, MPI_rank;
    
    //! Size of the arrays in pixels
    std::vector<unsigned int> N;
    
    //! Physical size of the arrays
    std::vector<double> L;
    
    //! Array size that relates to the parallel decomposition
    std::vector<unsigned int> Nlocal;
    
    //! Arrays of coordinates (y, z, t) owning to the current processor
    std::vector<std::vector<double> > local_x;
    
    //! Arrays of wavenumbers (kx, ky, ky) owning to the current processor
    std::vector<std::vector<double> > local_k;
    
    //! True if 2D geometry, False if 3D
    bool _2D;
};



#endif
