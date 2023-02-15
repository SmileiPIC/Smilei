#include "Field.h"
#include "gpu.h"

void Field::put_to( double val )
{
    SMILEI_ASSERT( data_ != nullptr );

#if defined( SMILEI_ACCELERATOR_GPU_OMP ) || defined( SMILEI_OPENACC_MODE )
    const bool is_hostptr_mapped_on_device = smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( data_ );
#endif

    // NVCC's OpenACC needs that redundant pointer value
    double* an_other_data_pointer = data_;

#if defined( SMILEI_OPENACC_MODE )
    // Test if data exists on GPU, put_to can be used on CPU and GPU during a simulation
    #pragma acc parallel present( an_other_data_pointer [0:size()] ) if( is_hostptr_mapped_on_device )
    #pragma acc loop gang worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target if( is_hostptr_mapped_on_device )
    #pragma omp teams distribute parallel for
#endif
    for( unsigned int i = 0; i < size(); i++ ) {
        an_other_data_pointer[i] = val;
    }
}

#if defined(SMILEI_ACCELERATOR_MODE)
    //! copy the field array from Host to Device
    void Field::copyFromHostToDevice()
    {
        smilei::tools::gpu::HostDeviceMemoryManagement::CopyHostToDevice( data_, number_of_points_ );
    };

    //! copy from Device to Host
    void Field::copyFromDeviceToHost()
    {
        smilei::tools::gpu::HostDeviceMemoryManagement::CopyDeviceToHost( data_, number_of_points_ );
    };

     //! Allocate only on device (without copy ore init)
    void Field::allocateOnDevice()
    {
        smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( data_, number_of_points_ );
    };

    //! allocate and copy from Device to Host
    void Field::allocateAndCopyFromHostToDevice()
    {
        smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocateAndCopyHostToDevice( data_, number_of_points_ );
    };

    //! Return if the field grid is mapped on device
    bool Field::IsOnDevice()
    {
        return smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( data_ );
    };

#endif
