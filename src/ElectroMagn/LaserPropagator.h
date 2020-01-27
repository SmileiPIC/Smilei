
#ifndef LaserPropagator_H
#define LaserPropagator_H

#include "PyTools.h"
#include "Field2D.h"
#include "SmileiMPI.h"

#include <vector>
#include <string>
#include <cmath>
#include <complex>

class Params;



// Interpolate a complex array z of size jmax at location j
inline std::complex<double> complex_interpolate( std::complex<double> *z, unsigned int jmax, double j, double extra_magnitude, double extra_phase )
{
    double remainder, intpart;
    remainder = std::modf( j, &intpart ); // integer and remainder parts
    unsigned int j0 = ( unsigned int ) intpart;
    unsigned int j1 = ( j0 + 1 )%jmax;
    extra_magnitude *= ( 1.-remainder )*std::abs( z[j0] ) + remainder*std::abs( z[j1] ); // interpolate the magnitude first
    double phase0 = std::arg( z[j0] );
    double phase1 = std::arg( z[j1] );
    if( phase1<phase0 ) {
        phase1 +=  6.283185307179586476925286766559;    // 2 pi
    }
    extra_phase += ( 1.-remainder )*phase0 + remainder*phase1; // interpolate the phase
    return std::polar( extra_magnitude, extra_phase );
};



class LaserPropagator
{
public:
    LaserPropagator() {};
    ~LaserPropagator() {};
    
    void init( Params *params, SmileiMPI *smpi, unsigned int side );
    
    // Propagates the fields profiles with some offset, and writes result to file
    void operator()( std::vector<PyObject *>, std::vector<int>, double, std::string, int, double );
    
protected:

private:
    unsigned int ndim;
    unsigned int MPI_size, MPI_rank;
    
    //! Size of the arrays in pixels
    std::vector<unsigned int> N;
    
    //! Physical size of the arrays
    std::vector<double> L;
    
    //! Physical oversize
    std::vector<double> o;
    
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
