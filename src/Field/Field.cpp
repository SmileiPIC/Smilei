#include "Field.h"

#include "gpu.h"

void Field::put_to( double val )
{
    // If openmp or openacc are enabled in smilei, it'll return false (data_ == nullptr) or (data_ is not mapped).
    // If not (openmp or openacc) are enabled, it'll be equivalent to  data_ != nullptr.
    const bool is_hostptr_mapped_on_device = smilei::tools::gpu::HostDeviceMemoryManagment::IsHostPointerMappedOnDevice( data_ );

    if( is_hostptr_mapped_on_device ) {
#if defined( _GPU )
    // Test if data exists on GPU, put_to can be used on CPU and GPU during a simulation
    #pragma acc parallel         present( data_ [0:globalDims_] )
    #pragma acc loop gang worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target defaultmap( none ) \
        map( to                           \
             : globalDims_, val )         \
            map( tofrom                   \
                 : data_ [0:globalDims_] )
    #pragma omp            teams /* num_teams(xxx) thread_limit(xxx) */ // TODO(Etienne M): WG/WF tuning
    #pragma omp distribute parallel for
#endif
        for( unsigned int i = 0; i < globalDims_; i++ ) {
            data_[i] = val;
        }
    }
}
