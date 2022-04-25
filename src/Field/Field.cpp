#include "Field.h"

#include <iostream>
#ifdef _GPU
    #include <openacc.h>
#endif
#include "gpu.h"

void Field::put_to( double val )
{
    if( data_ ) {
        double* hostptr = data_;
#if defined( _GPU )
        void* dataptr = acc_deviceptr( hostptr );
        // Test if data exists on GPU, put_to can be used on CPU and GPU during a simulation
        #pragma acc  parallel present(hostptr[0:globalDims_]) if (dataptr!=NULL)
        #pragma acc loop gang worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    std::cout << "IsHostPointerMappedOnDevice( hostptr ) "
              << smilei::tools::HostDeviceMemoryManagment::IsHostPointerMappedOnDevice( hostptr ) << "\n";
    #pragma omp target if( smilei::tools::HostDeviceMemoryManagment::IsHostPointerMappedOnDevice( hostptr ) ) \
        defaultmap( none )                                                                                    \
            map( to                                                                                           \
                 : globalDims_, val )                                                                         \
                map( tofrom                                                                                   \
                     : hostptr [0:globalDims_] )
    #pragma omp            teams /* num_teams(xxx) thread_limit(xxx) */ // TODO(Etienne M): WG/WF tuning
    #pragma omp distribute parallel for
#endif
        for( unsigned int i=0; i<globalDims_; i++ ) {
            hostptr[i] = val;
        }
    }
}
