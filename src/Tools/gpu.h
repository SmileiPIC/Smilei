#ifndef SMILEI_TOOLS_GPU_H
#define SMILEI_TOOLS_GPU_H

#include <cstdlib>
#include <cstring>
#include <type_traits>

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #include <omp.h>
#elif defined( _GPU )
    #include <openacc.h>
#endif

#include "Tools.h"

namespace smilei {
    namespace tools {
        namespace gpu {

            ////////////////////////////////////////////////////////////////////////////////
            // Omp/OpenACC
            ////////////////////////////////////////////////////////////////////////////////

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #define SMILEI_ACCELERATOR_DECLARE_ROUTINE     _Pragma( "omp declare target" )
    #define SMILEI_ACCELERATOR_DECLARE_ROUTINE_END _Pragma( "omp end declare target" )
#elif defined( _GPU )
    #define SMILEI_ACCELERATOR_DECLARE_ROUTINE _Pragma( "acc routine seq" )
    #define SMILEI_ACCELERATOR_DECLARE_ROUTINE_END
#else
    #define SMILEI_ACCELERATOR_DECLARE_ROUTINE
    #define SMILEI_ACCELERATOR_DECLARE_ROUTINE_END
#endif


            ////////////////////////////////////////////////////////////////////////////////
            // NonInitializingVector
            ////////////////////////////////////////////////////////////////////////////////

            /// Trivial container that does not initialize the memory after allocating.
            /// This differ from the traditionnal std::vector which does initialize the memory,
            /// leading to a significant overhead (at initialization time).
            /// This NonInitializingVector can thus better make use of the virtual memory
            /// when used in cunjunction with the openMP/OpenACC device offloading.
            ///
            /// Note:
            /// When seeking performance, more control often means more potential performance.
            /// This NonInitializingVector provides a way to automatically free the memory
            /// allocated on the device (avoid leaks) but requires the user to explicicly to
            /// the initial device allocation which, when done correctly, often source of
            /// speedup (async host<->device copies for instance).
            ///
            template <typename T,
                      bool do_device_free = false>
            class NonInitializingVector
            {
            public:
                using value_type = T;

            public:
                NonInitializingVector();

                /// Check HostAlloc() for more info
                ///
                explicit NonInitializingVector( std::size_t size );

                /// Named HostAlloc instead of just Alloc so that the user knows
                /// that it does nothing on the device!
                ///
                /// Note:
                /// Does not initialize memory, meaning, due to how the virtual
                /// memory works, that only when the memory is "touched"/set will the
                /// process' true memory usage increase. If you map to the device and never
                /// touch the host memory, the host memory usage not physically increase.
                ///
                void HostAlloc( std::size_t size );
                void Free();

                T*       begin();
                const T* cbegin() const;

                T*       end();
                const T* cend() const;

                T*       data();
                const T* data() const;

                T&
                operator[]( std::size_t index );
                const T&
                operator[]( std::size_t index ) const;

                std::size_t size() const;

                ~NonInitializingVector();

            protected:
                std::size_t           size_;
                T* /* __restrict__ */ data_;
            };


            ////////////////////////////////////////////////////////////////////////////////
            // HostDeviceMemoryManagment
            ////////////////////////////////////////////////////////////////////////////////

            /// Exploits the OpenMP/OpenACC host/device memory mapping capabilities.
            /// These function requires you to already have a pointer, this pointer will be
            /// mapped to a pointer pointing to a chunk of a given size on the device memory.
            /// This mapping is stored inside the OpenMP/OpenACC runtime.
            ///
            /// Note:
            /// - The OpenACC implementation is not complete !
            /// - The host pointer may not point to valid data for the OpenMP mapping to
            /// occur.
            /// - You may exploit the virtual memory and allocate a large chunk of memory on
            /// the host (malloc) and not use it. In this case no memory page will be marked
            /// "in use" and you host memory consumption will be zero ! You may exploit this
            /// fact using NonInitializingVector to easily produce GPU/CPU memory optimized
            /// software without worrying about the host memory consumption when offloading
            /// an Omp kernel to the GPU.
            /// If fact, you could say that malloc can be used as an excuse to get
            /// a unique value, that is, the returned pointer (as long as it's not freed).
            /// This unique value can be mapped to a valid memory chunk allocated on the GPU.
            /// - Does not support async operation. If you need that it's probably
            /// better if you do it yourself (not using HostDeviceMemoryManagment) as it can be
            /// quite tricky. HostDeviceMemoryManagment is best for allocating/copying big chunks
            /// at the start of the program.
            ///
            struct HostDeviceMemoryManagment
            {
            public:
                template <typename T>
                static void DeviceAllocate( const T* a_pointer, std::size_t a_size );
                template <typename Container>
                static void DeviceAllocate( const Container& a_vector );

                template <typename T>
                static void DeviceAllocateAndCopyHostToDevice( const T* a_pointer, std::size_t a_size );
                template <typename Container>
                static void DeviceAllocateAndCopyHostToDevice( const Container& a_vector );

                template <typename T>
                static void CopyHostToDevice( const T* a_pointer, std::size_t a_size );
                template <typename Container>
                static void CopyHostToDevice( const Container& a_vector );

                template <typename T>
                static void CopyDeviceToHost( T* a_pointer, std::size_t a_size );
                template <typename Container>
                static void CopyDeviceToHost( Container& a_vector );

                template <typename T>
                static void CopyDeviceToHostAndDeviceFree( T* a_pointer, std::size_t a_size );
                template <typename Container>
                static void CopyDeviceToHostAndDeviceFree( Container& a_vector );

                template <typename T>
                static void DeviceFree( T* a_pointer, std::size_t a_size );
                template <typename Container>
                static void DeviceFree( Container& a_vector );

                /// If OpenMP or OpenACC are enabled and if a_pointer is mapped, returns the pointer on the device.
                ///                                      else return nullptr
                /// else return a_pointer (untouched)
                ///
                /// Note:
                /// the nvidia compiler of the NVHPC 21.3 stack has a bug in ::omp_target_is_present. You can't use this 
                /// function unless you first maek the runtime "aware" (explicit mapping) of the pointer!
                ///
                /// #if defined( __NVCOMPILER )
                ///     No-op workaround to prevent from a bug in Nvidia's OpenMP implementation:
                ///     https://forums.developer.nvidia.com/t/nvc-v21-3-omp-target-is-present-crashes-the-program/215585
                /// #else
                ///
                template <typename T>
                static T* GetDevicePointer( T* a_pointer );

                template <typename T>
                static bool IsHostPointerMappedOnDevice( const T* a_pointer );

                /// Expects host pointers passed through GetDevicePointer. a_count T's are copied (dont specify the byte
                /// count only object count).
                ///
                /// ie:
                /// DeviceMemoryCopy(GetDevicePointer(a + 5), GetDevicePointer(b) + 10, 10);
                ///
                template <typename T>
                static void DeviceMemoryCopy( T* a_destination, const T* a_source, std::size_t a_count );
            };


            ////////////////////////////////////////////////////////////////////////////////
            // Macros
            ////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////
/// @def SMILEI_GPU_ASSERT_MEMORY_ON_DEVICE
///
/// Makes sure the host pointer is mapped on the device through OpenACC/OpenMP.
/// This can be used to simulate the present() clause of OpenACC in an OpenMP
/// context. There is not present() clause in OpenMP
///
/// Example usage:
///
///    #pragma omp target teams distribute parallel for
///    for(...) { ... }
///
//////////////////////////////////////
#define SMILEI_GPU_ASSERT_MEMORY_IS_ON_DEVICE( a_host_pointer ) SMILEI_ASSERT( smilei::tools::gpu::HostDeviceMemoryManagment::IsHostPointerMappedOnDevice( a_host_pointer ) )


            ////////////////////////////////////////////////////////////////////////////////
            // NonInitializingVector methods definition
            ////////////////////////////////////////////////////////////////////////////////

            template <typename T,
                      bool do_device_free>
            NonInitializingVector<T, do_device_free>::NonInitializingVector()
                : size_{}
                , data_{ nullptr }
            {
                // EMPTY
            }

            template <typename T,
                      bool do_device_free>
            NonInitializingVector<T, do_device_free>::NonInitializingVector( std::size_t size )
                : NonInitializingVector{}
            {
                HostAlloc( size );
            }

            template <typename T,
                      bool do_device_free>
            void NonInitializingVector<T, do_device_free>::HostAlloc( std::size_t size )
            {
                SMILEI_ASSERT_VERBOSE( size_ == 0 && data_ == nullptr,
                                       "NonInitializingVector::Alloc, allocation before deallocating." );

                data_ = static_cast<T*>( std::malloc( sizeof( T ) * size ) );

                SMILEI_ASSERT_VERBOSE( data_ != nullptr,
                                       "NonInitializingVector::Alloc, std::malloc() out of memory." );

                size_ = size;
            }

            template <typename T,
                      bool do_device_free>
            void NonInitializingVector<T, do_device_free>::Free()
            {
                // According to the C++ standard, if data_ == nullptr, the function does nothing
                std::free( data_ );

                if( do_device_free &&
                    // Unlike std::free, we check to avoid nullptr freeing
                    data_ != nullptr ) {
                    HostDeviceMemoryManagment::DeviceFree( *this );
                }

                data_ = nullptr;
                size_ = 0;
            }

            template <typename T,
                      bool do_device_free>
            T* NonInitializingVector<T, do_device_free>::begin()
            {
                return data_;
            }

            template <typename T,
                      bool do_device_free>
            const T* NonInitializingVector<T, do_device_free>::cbegin() const
            {
                return data_;
            }

            template <typename T,
                      bool do_device_free>
            T* NonInitializingVector<T, do_device_free>::end()
            {
                return data_ + size_;
            }

            template <typename T,
                      bool do_device_free>
            const T* NonInitializingVector<T, do_device_free>::cend() const
            {
                return data_ + size_;
            }

            template <typename T,
                      bool do_device_free>
            T* NonInitializingVector<T, do_device_free>::data()
            {
                return data_;
            }

            template <typename T,
                      bool do_device_free>
            const T* NonInitializingVector<T, do_device_free>::data() const
            {
                return data_;
            }

            template <typename T,
                      bool do_device_free>
            T& NonInitializingVector<T, do_device_free>::operator[]( std::size_t index )
            {
                return data()[index];
            }

            template <typename T,
                      bool do_device_free>
            const T&
            NonInitializingVector<T, do_device_free>::operator[]( std::size_t index ) const
            {
                return data()[index];
            }

            template <typename T,
                      bool do_device_free>
            std::size_t NonInitializingVector<T, do_device_free>::size() const
            {
                return size_;
            }

            template <typename T,
                      bool do_device_free>
            NonInitializingVector<T, do_device_free>::~NonInitializingVector()
            {
                Free();
            }


            ////////////////////////////////////////////////////////////////////////////////
            // HostDeviceMemoryManagment methods definition
            ////////////////////////////////////////////////////////////////////////////////

            template <typename T>
            void HostDeviceMemoryManagment::DeviceAllocate( const T* a_pointer, std::size_t a_size )
            {
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target enter data map( alloc \
                                       : a_pointer [0:a_size] )
#elif defined( _GPU )
    #pragma acc enter data create( a_pointer [0:a_size] )
#else
                SMILEI_UNUSED( a_pointer );
                SMILEI_UNUSED( a_size );
#endif
            }

            template <typename Container>
            void HostDeviceMemoryManagment::DeviceAllocate( const Container& a_vector )
            {
                DeviceAllocate( a_vector.data(), a_vector.size() );
            }

            template <typename T>
            void HostDeviceMemoryManagment::DeviceAllocateAndCopyHostToDevice( const T* a_pointer, std::size_t a_size )
            {
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target enter data map( to \
                                       : a_pointer [0:a_size] )
#elif defined( _GPU )
    #pragma acc enter data copyin( a_pointer [0:a_size] )
#else
                SMILEI_UNUSED( a_pointer );
                SMILEI_UNUSED( a_size );
#endif
            }

            template <typename Container>
            void HostDeviceMemoryManagment::DeviceAllocateAndCopyHostToDevice( const Container& a_vector )
            {
                DeviceAllocAndCopyHostToDevice( a_vector.data(), a_vector.size() );
            }

            template <typename T>
            void HostDeviceMemoryManagment::CopyHostToDevice( const T* a_pointer, std::size_t a_size )
            {
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target update to( a_pointer [0:a_size] )
#elif defined( _GPU )
    #pragma acc update device( a_pointer [0:a_size] )
#else
                SMILEI_UNUSED( a_pointer );
                SMILEI_UNUSED( a_size );
#endif
            }

            template <typename Container>
            void HostDeviceMemoryManagment::CopyHostToDevice( const Container& a_vector )
            {
                CopyHostToDevice( a_vector.data(), a_vector.size() );
            }

            template <typename T>
            void HostDeviceMemoryManagment::CopyDeviceToHost( T* a_pointer, std::size_t a_size )
            {
                static_assert( !std::is_const<T>::value, "" );
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target update from( a_pointer [0:a_size] )
#elif defined( _GPU )
    #pragma acc update host( a_pointer [0:a_size] )
#else
                SMILEI_UNUSED( a_pointer );
                SMILEI_UNUSED( a_size );
#endif
            }

            template <typename Container>
            void HostDeviceMemoryManagment::CopyDeviceToHost( Container& a_vector )
            {
                CopyDeviceToHost( a_vector.data(), a_vector.size() );
            }

            template <typename T>
            void HostDeviceMemoryManagment::CopyDeviceToHostAndDeviceFree( T* a_pointer, std::size_t a_size )
            {
                static_assert( !std::is_const<T>::value, "" );
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target exit data map( from \
                                      : a_pointer [0:a_size] )
#elif defined( _GPU )
    #pragma acc exit data copyout( a_pointer [0:a_size] )
#else
                SMILEI_UNUSED( a_pointer );
                SMILEI_UNUSED( a_size );
#endif
            }

            template <typename Container>
            void HostDeviceMemoryManagment::CopyDeviceToHostAndDeviceFree( Container& a_vector )
            {
                CopyDeviceToHostAndDeviceFree( a_vector.data(), a_vector.size() );
            }

            template <typename T>
            void HostDeviceMemoryManagment::DeviceFree( T* a_pointer, std::size_t a_size )
            {
                static_assert( !std::is_const<T>::value, "" );
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target exit data map( delete \
                                      : a_pointer [0:a_size] )
#elif defined( _GPU )
    #pragma acc exit data delete( a_pointer [0:a_size] )
#else
                SMILEI_UNUSED( a_pointer );
                SMILEI_UNUSED( a_size );
#endif
            }

            template <typename Container>
            void HostDeviceMemoryManagment::DeviceFree( Container& a_vector )
            {
                DeviceFree( a_vector.data(), a_vector.size() );
            }

            template <typename T>
            T* HostDeviceMemoryManagment::GetDevicePointer( T* a_host_pointer )
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

                T* a_device_pointer = nullptr;

    #pragma omp target data use_device_ptr( a_host_pointer )
                {
                    a_device_pointer = a_host_pointer;
                }

                return a_device_pointer;
#elif defined( _GPU )
                return static_cast<T*>( ::acc_deviceptr( a_host_pointer ) );
#else
                return a_host_pointer;
#endif
            }

            template <typename T>
            bool HostDeviceMemoryManagment::IsHostPointerMappedOnDevice( const T* a_host_pointer )
            {
                // We could optimize the omp version by only using ::omp_target_is_present()
                return GetDevicePointer( a_host_pointer ) != nullptr;
            }

            template <typename T>
            void HostDeviceMemoryManagment::DeviceMemoryCopy( T* a_destination, const T* a_source, std::size_t a_count )
            {
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
                const int device_num = ::omp_get_default_device();
                if( ::omp_target_memcpy( a_destination,
                                         a_source,
                                         a_count * sizeof( T ), 0, 0, device_num, device_num ) != 0 ) {
                    ERROR( "omp_target_memcpy failed" );
                }
#elif defined( _GPU )
                // It seems that the interface of ::acc_memcpy_device does not accept ptr to array of const type !
                // https://www.openacc.org/sites/default/files/inline-files/OpenACC.2.7.pdf
                // void acc_memcpy_device( d_void* dest, d_void* src, size_t bytes );
                ::acc_memcpy_device( a_destination, const_cast<T*>( a_source ), a_count * sizeof( T ) );
#else
                std::memcpy( a_destination, a_source, a_count * sizeof( T ) );
#endif
            }

        } // namespace gpu
    }     // namespace tools
} // namespace smilei

#endif
