#include "Field.h"

#include "gpu.h"

void Field::put_to( double val )
{
    // If openmp or openacc are enabled in smilei, it'll return false (data_ == nullptr) or (data_ is not mapped).
    // If not (openmp or openacc) are enabled, it'll be equivalent to  data_ != nullptr.
    const bool is_hostptr_mapped_on_device = smilei::tools::gpu::HostDeviceMemoryManagment::IsHostPointerMappedOnDevice( data_ );
    SMILEI_UNUSED( is_hostptr_mapped_on_device );

    if( data_ ) {
        // OpenACC needs that redundant pointeur value
        double* an_other_data_pointer = data_;
#if defined( _GPU )
    // Test if data exists on GPU, put_to can be used on CPU and GPU during a simulation
    #pragma acc parallel         present( an_other_data_pointer [0:globalDims_] ) if( is_hostptr_mapped_on_device )
    #pragma acc loop gang worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target if( is_hostptr_mapped_on_device )
    #pragma omp teams /* num_teams(xxx) thread_limit(xxx) */ // TODO(Etienne M): WG/WF tuning
    #pragma omp distribute parallel for
#endif
        for( unsigned int i = 0; i < globalDims_; i++ ) {
            an_other_data_pointer[i] = val;
        }
    }
}
