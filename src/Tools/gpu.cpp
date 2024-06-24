#include "gpu.h"

#if defined( SMILEI_ACCELERATOR_GPU_OMP ) && defined( SMILEI_ACCELERATOR_GPU_OACC )
    #error "You can not enable both OpenACC and OpenMP GPU support"
#endif

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #if defined( _OPENMP )
        // NOTE: Nvidia offers poor support for OpenMP (they 
        // prefer OpenACC).
        // - NVHPC-sdk 21.3 only support OpenMP 4.0.
        // - We require OpenMP 4.5 and up (Clang has partial 
        // support for 5.0 and AMD GPU).

        // Ask for OpenMP 5.0 so the user has the feature he 
        // would expect from a up to date OpenMP implementation
        #if _OPENMP < 201811
            // 200505 2.5
            // 200805 3.0
            // 201107 3.1
            // 201307 4.0
            // 201511 4.5
            // 201811 5.0
            // 202011 5.1
            #error "OpenMP 5.0 is required"
        #endif

        #include <omp.h>
    #else
        #error "Asking for OpenMP support without enabling compiler support for OpenMP"
    #endif
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #if defined( _OPENACC )
        #include <openacc.h>
    #else
        #error "Asking for OpenACC support without enabling compiler support for OpenACC"
    #endif
#endif

namespace smilei {
    namespace tools {
        namespace gpu {
            void HostDeviceMemoryManagement::DoDeviceAllocate( const void* a_host_pointer, std::size_t a_count, std::size_t an_object_size )
            {
                const unsigned char* byte_array = static_cast<const unsigned char*>( a_host_pointer );
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target enter data map( alloc \
                                       : byte_array [0:a_count * an_object_size] )
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc enter data create( byte_array [0:a_count * an_object_size] )
#else
                SMILEI_UNUSED( a_host_pointer );
                SMILEI_UNUSED( a_count );
                SMILEI_UNUSED( an_object_size );
                SMILEI_UNUSED( byte_array );
#endif
            }

            void HostDeviceMemoryManagement::DoDeviceAllocateAndCopyHostToDevice( const void* a_host_pointer, std::size_t a_count, std::size_t an_object_size )
            {
                const unsigned char* byte_array = static_cast<const unsigned char*>( a_host_pointer );
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target enter data map( to \
                                       : byte_array [0:a_count * an_object_size] )
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc enter data copyin( byte_array [0:a_count * an_object_size] )
#else
                SMILEI_UNUSED( a_host_pointer );
                SMILEI_UNUSED( a_count );
                SMILEI_UNUSED( an_object_size );
                SMILEI_UNUSED( byte_array );
#endif
            }

            void HostDeviceMemoryManagement::DoCopyHostToDevice( const void* a_host_pointer, std::size_t a_count, std::size_t an_object_size )
            {
                const unsigned char* byte_array = static_cast<const unsigned char*>( a_host_pointer );
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target update to( byte_array [0:a_count * an_object_size] )
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc update device( byte_array [0:a_count * an_object_size] )
#else
                SMILEI_UNUSED( a_host_pointer );
                SMILEI_UNUSED( a_count );
                SMILEI_UNUSED( an_object_size );
                SMILEI_UNUSED( byte_array );
#endif
            }

            void HostDeviceMemoryManagement::DoCopyDeviceToHost( void* a_host_pointer, std::size_t a_count, std::size_t an_object_size )
            {
                unsigned char* byte_array = static_cast<unsigned char*>( a_host_pointer );
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target update from( byte_array [0:a_count * an_object_size] )
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc update host( byte_array [0:a_count * an_object_size] )
#else
                SMILEI_UNUSED( a_host_pointer );
                SMILEI_UNUSED( a_count );
                SMILEI_UNUSED( an_object_size );
                SMILEI_UNUSED( byte_array );
#endif
            }

            void HostDeviceMemoryManagement::DoCopyDeviceToHostAndDeviceFree( void* a_host_pointer, std::size_t a_count, std::size_t an_object_size )
            {
                unsigned char* byte_array = static_cast<unsigned char*>( a_host_pointer );
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target exit data map( from \
                                      : byte_array [0:a_count * an_object_size] )
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc exit data copyout( byte_array [0:a_count * an_object_size] )
#else
                SMILEI_UNUSED( a_host_pointer );
                SMILEI_UNUSED( a_count );
                SMILEI_UNUSED( an_object_size );
                SMILEI_UNUSED( byte_array );
#endif
            }

            void HostDeviceMemoryManagement::DoDeviceFree( void* a_host_pointer, std::size_t a_count, std::size_t an_object_size )
            {
                unsigned char* byte_array = static_cast<unsigned char*>( a_host_pointer );
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target exit data map( delete \
                                      : byte_array [0:a_count * an_object_size] )
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc exit data delete( byte_array [0:a_count * an_object_size] )
#else
                SMILEI_UNUSED( a_host_pointer );
                SMILEI_UNUSED( a_count );
                SMILEI_UNUSED( an_object_size );
                SMILEI_UNUSED( byte_array );
#endif
            }

            void* HostDeviceMemoryManagement::DoGetDevicePointer( const void* a_host_pointer )
            {
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
                const int device_num = ::omp_get_default_device();

                // Omp Std 5.0: A list item in a use_device_ptr clause must hold
                // the address of an object that has a corresponding list item
                // in the device data environment.
                // To be fully compliant we need to use ::omp_target_is_present

                if( ::omp_target_is_present( a_host_pointer, device_num ) == 0 ) {
                    return nullptr;
                }

                const void* a_device_pointer = nullptr;

                // NOTE: OpenMP 5.1 offers ::omp_get_mapped_ptr to the the 
                // operation below
                #pragma omp target data use_device_ptr( a_host_pointer )
                {
                    a_device_pointer = a_host_pointer;
                }

                SMILEI_ASSERT( a_device_pointer != nullptr );

                return const_cast<void*>( a_device_pointer );
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
                //return const_cast<void*>( ::acc_deviceptr( a_host_pointer ) );
                return ::acc_deviceptr( const_cast<void*>(a_host_pointer) ) ;
#else
                return const_cast<void*>( a_host_pointer );
#endif
            }

            void HostDeviceMemoryManagement::DoDeviceMemoryCopy( void* a_destination, const void* a_source, std::size_t a_count, std::size_t an_object_size )
            {
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
                const int device_num = ::omp_get_default_device();
                if( ::omp_target_memcpy( a_destination,
                                         a_source,
                                         a_count * an_object_size, 0, 0, device_num, device_num ) != 0 ) {
                    ERROR( "omp_target_memcpy failed" );
                }
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
                // It seems that the interface of ::acc_memcpy_device does not accept ptr to array of const type !
                // https://www.openacc.org/sites/default/files/inline-files/OpenACC.2.7.pdf
                // void acc_memcpy_device( d_void* dest, d_void* src, size_t bytes );
                ::acc_memcpy_device( a_destination, const_cast<void*>( a_source ), a_count * an_object_size );
#else
                std::memcpy( a_destination, a_source, a_count * an_object_size );
#endif
            }


        } // namespace gpu
    }     // namespace tools
} // namespace smilei
