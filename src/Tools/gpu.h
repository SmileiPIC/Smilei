#ifndef SMILEI_TOOLS_GPU_H
#define SMILEI_TOOLS_GPU_H

#include <cstdlib>
#include <cstring>
#include <type_traits>

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
    #define SMILEI_ACCELERATOR_ATOMIC _Pragma( "omp atomic update" )
    #define SMILEI_ACCELERATOR_ASYNC_POLYCY thrust::hip::par_nosync
    #define SMILEI_ACCELERATOR_DEVICE_SYNC() hipDeviceSynchronize()
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #define SMILEI_ACCELERATOR_DECLARE_ROUTINE _Pragma( "acc routine seq" )
    #define SMILEI_ACCELERATOR_DECLARE_ROUTINE_END
    #define SMILEI_ACCELERATOR_ATOMIC _Pragma( "acc atomic" )
    #define SMILEI_ACCELERATOR_ASYNC_POLYCY thrust::cuda::par_nosync
    #define SMILEI_ACCELERATOR_DEVICE_SYNC() cudaDeviceSynchronize()
#else
    #define SMILEI_ACCELERATOR_DECLARE_ROUTINE
    #define SMILEI_ACCELERATOR_DECLARE_ROUTINE_END
    #define SMILEI_ACCELERATOR_ATOMIC
#endif


            ////////////////////////////////////////////////////////////////////////////////
            // NonInitializingVector
            ////////////////////////////////////////////////////////////////////////////////

            /// Trivial container that does not initialize the memory after allocating.
            /// This differ from the traditional std::vector which does initialize the memory,
            /// leading to a significant overhead (at initialization time).
            /// This NonInitializingVector can thus better make use of the virtual memory
            /// when used in conjunction with the openMP/OpenACC device offloading.
            ///
            /// NOTE:
            /// When seeking performance, more control often means more potential performance.
            /// This NonInitializingVector provides a way to automatically free the memory
            /// allocated on the device (avoid leaks) but requires the user to explicitly to
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
                /// NOTE:
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
            // HostDeviceMemoryManagement
            ////////////////////////////////////////////////////////////////////////////////

            /// Exploits the host and device memory mapping capabilities of OpenMP/OpenACC.
            /// These functions require that you already have a pointer, this pointer will be
            /// mapped to a pointer pointing to a chunk of a given size on the device memory.
            /// This mapping is stored in the OpenMP/OpenACC runtime.
            /// You should be able to use this header whether or not you ahve GPU support
            /// enabled.
            ///
            /// Do not allocate classes using non trivial constructor/destructor !
            ///
            /// NOTE:
            /// - Do not resize in between allocations or you will confuse OpenMP and risk memory corruption!
            /// - The OpenACC implementation is not complete!
            /// - You can exploit virtual memory and allocate a large part of the memory on the
            /// the host (malloc) and not use it. The OS will allocate address space and not physical
            /// memory until you touch the page. Before touching the page, the host physical memory allocation
            /// will be zero! You can exploit this fact by using NonInitializingVector to easily produce
            /// software optimized for GPU/CPU memory without worrying about consuming host memory when offloading
            /// an kernel to the GPU. In fact, one could say that malloc can be used as an excuse to get
            /// a unique value, i.e. the returned pointer (as long as it is not freed).
            /// This unique value can be mapped to a valid chunk of memory allocated on the GPU.
            /// - Does not support asynchronous operations. If you need it, it is probably
            /// better if you do it yourself (without using HostDeviceMemoryManagement) because it can be
            /// quite tricky. HostDeviceMemoryManagement is the best solution to allocate/copy large chunks
            /// at the beginning of the program.
            /// - Everything is hidden in gpu.cpp so we dont get conflicts between GPU specific languages (HIP/Cuda)
            /// and OpenMP/OpenACC (the cray compiler can't enable both hip and OpenMP support at the same time).
            /// - The is_device_ptr() clause has support for device pointers created outside of OpenMP but the behavior
            /// is implementation defined. For us, it behaves as expected.
            ///
            struct HostDeviceMemoryManagement
            {
            public:
                template <typename T>
                static void DeviceAllocate( const T* a_host_pointer, std::size_t a_count );
                template <typename Container>
                static void DeviceAllocate( const Container& a_vector );

                template <typename T>
                static void DeviceAllocateAndCopyHostToDevice( const T* a_host_pointer, std::size_t a_count );
                template <typename Container>
                static void DeviceAllocateAndCopyHostToDevice( const Container& a_vector );

                template <typename T>
                static void CopyHostToDevice( const T* a_host_pointer, std::size_t a_count );
                template <typename Container>
                static void CopyHostToDevice( const Container& a_vector );

                template <typename T>
                static void CopyDeviceToHost( T* a_host_pointer, std::size_t a_count );
                template <typename Container>
                static void CopyDeviceToHost( Container& a_vector );

                template <typename T>
                static void CopyDeviceToHostAndDeviceFree( T* a_host_pointer, std::size_t a_count );
                template <typename Container>
                static void CopyDeviceToHostAndDeviceFree( Container& a_vector );

                template <typename T>
                static void DeviceFree( T* a_host_pointer, std::size_t a_count );
                template <typename Container>
                static void DeviceFree( Container& a_vector );

                /// If OpenMP or OpenACC are enabled and if a_host_pointer is mapped, returns the pointer on the device.
                ///                                      else return nullptr
                /// else return a_host_pointer (untouched)
                ///
                /// NOTE:
                /// the nvidia compiler of the NVHPC 21.3 stack has a bug in ::omp_target_is_present. You can't use this
                /// function unless you first make the runtime "aware" (explicit mapping) of the pointer!
                ///
                /// #if defined( __NVCOMPILER )
                ///     No-op workaround to prevent from a bug in Nvidia's OpenMP implementation:
                ///     https://forums.developer.nvidia.com/t/nvc-v21-3-omp-target-is-present-crashes-the-program/215585
                /// #else
                ///
                template <typename T>
                static T* GetDevicePointer( T* a_host_pointer );

                /// Smilei's code does a lot of runtime checking to know if we are using GPU or CPU data.
                /// Sometimes, we just want to get the GPU pointer if it exist, or the host pointer if no GPU equivalent
                /// exists. ie: MPI_Isend in Patch::initSumField() There, we dont know which buffer we deal with under
                /// the contract that the data will be exchanged will be from the GPU if it exists there, or else, on
                /// from the CPU.
                ///
                template <typename T>
                static T* GetDeviceOrHostPointer( T* a_host_pointer );

                template <typename T>
                static bool IsHostPointerMappedOnDevice( const T* a_host_pointer );

                /// Expects host pointers passed through GetDevicePointer. a_count T's are copied (dont specify the byte
                /// count only object count).
                ///
                /// ie:
                /// DeviceMemoryCopy(GetDevicePointer(a + 5), GetDevicePointer(b) + 10, 10);
                ///
                template <typename T>
                static void DeviceMemoryCopy( T* a_destination, const T* a_source, std::size_t a_count );

            protected:
                static void  DoDeviceAllocate( const void* a_host_pointer, std::size_t a_count, std::size_t an_object_size );
                static void  DoDeviceAllocateAndCopyHostToDevice( const void* a_host_pointer, std::size_t a_count, std::size_t an_object_size );
                static void  DoCopyHostToDevice( const void* a_host_pointer, std::size_t a_count, std::size_t an_object_size );
                static void  DoCopyDeviceToHost( void* a_host_pointer, std::size_t a_count, std::size_t an_object_size );
                static void  DoCopyDeviceToHostAndDeviceFree( void* a_host_pointer, std::size_t a_count, std::size_t an_object_size );
                static void  DoDeviceFree( void* a_host_pointer, std::size_t a_count, std::size_t an_object_size );
                static void* DoGetDevicePointer( const void* a_host_pointer );
                static void  DoDeviceMemoryCopy( void* a_destination, const void* a_source, std::size_t a_count, std::size_t an_object_size );
            };


            ////////////////////////////////////////////////////////////////////////////////
            // Macros
            ////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////
/// @def SMILEI_GPU_ASSERT_MEMORY_ON_DEVICE
/// @def SMILEI_GPU_ASSERT_MEMORY_NOT_ON_DEVICE
///
/// Makes sure the host pointer is mapped on or is not mapped on the device
/// through OpenACC/OpenMP. This can be used to simulate the present() clause of
///  OpenACC in an OpenMP context. There is no present() clause in OpenMP.
///
/// Example usage:
///     smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( a_host_pointer, 10 );
///     SMILEI_GPU_ASSERT_MEMORY_IS_ON_DEVICE( a_host_pointer );  // Succeed
///     SMILEI_GPU_ASSERT_MEMORY_NOT_ON_DEVICE( a_host_pointer ); // Fails
///     smilei::tools::gpu::HostDeviceMemoryManagement::DeviceFree( a_host_pointer, 10 );
///     SMILEI_GPU_ASSERT_MEMORY_IS_ON_DEVICE( a_host_pointer );  // Fails
///     SMILEI_GPU_ASSERT_MEMORY_NOT_ON_DEVICE( a_host_pointer ); // Succeed
///
//////////////////////////////////////
#define SMILEI_GPU_ASSERT_MEMORY_IS_ON_DEVICE( a_host_pointer ) SMILEI_ASSERT( smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( a_host_pointer ) )

#define SMILEI_GPU_ASSERT_MEMORY_NOT_ON_DEVICE( a_host_pointer ) SMILEI_ASSERT( !smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( a_host_pointer ) )

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
                    HostDeviceMemoryManagement::DeviceFree( *this );
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
            // HostDeviceMemoryManagement methods definition
            ////////////////////////////////////////////////////////////////////////////////

            template <typename T>
            void HostDeviceMemoryManagement::DeviceAllocate( const T* a_host_pointer, std::size_t a_count )
            {
                static_assert( std::is_pod<T>::value, "" );
                DoDeviceAllocate( a_host_pointer, a_count, sizeof( T ) );
            }

            template <typename Container>
            void HostDeviceMemoryManagement::DeviceAllocate( const Container& a_vector )
            {
                // Dont use capacity(), because that does not play nice with 
                // functions like DeviceAllocateAndCopyHostToDevice (copy
                // uninitialized values).
                DeviceAllocate( a_vector.data(), a_vector.size() );
            }

            template <typename T>
            void HostDeviceMemoryManagement::DeviceAllocateAndCopyHostToDevice( const T* a_host_pointer, std::size_t a_count )
            {
                static_assert( std::is_pod<T>::value, "" );
                DoDeviceAllocateAndCopyHostToDevice( a_host_pointer, a_count, sizeof( T ) );
            }

            template <typename Container>
            void HostDeviceMemoryManagement::DeviceAllocateAndCopyHostToDevice( const Container& a_vector )
            {
                DeviceAllocateAndCopyHostToDevice( a_vector.data(), a_vector.size() );
            }

            template <typename T>
            void HostDeviceMemoryManagement::CopyHostToDevice( const T* a_host_pointer, std::size_t a_count )
            {
                static_assert( std::is_pod<T>::value, "" );
                DoCopyHostToDevice( a_host_pointer, a_count, sizeof( T ) );
            }

            template <typename Container>
            void HostDeviceMemoryManagement::CopyHostToDevice( const Container& a_vector )
            {
                CopyHostToDevice( a_vector.data(), a_vector.size() );
            }

            template <typename T>
            void HostDeviceMemoryManagement::CopyDeviceToHost( T* a_host_pointer, std::size_t a_count )
            {
                static_assert( !std::is_const<T>::value, "" );
                static_assert( std::is_pod<T>::value, "" );
                DoCopyDeviceToHost( a_host_pointer, a_count, sizeof( T ) );
            }

            template <typename Container>
            void HostDeviceMemoryManagement::CopyDeviceToHost( Container& a_vector )
            {
                CopyDeviceToHost( a_vector.data(), a_vector.size() );
            }

            template <typename T>
            void HostDeviceMemoryManagement::CopyDeviceToHostAndDeviceFree( T* a_host_pointer, std::size_t a_count )
            {
                static_assert( !std::is_const<T>::value, "" );
                static_assert( std::is_pod<T>::value, "" );
                DoCopyDeviceToHostAndDeviceFree( a_host_pointer, a_count, sizeof( T ) );
            }

            template <typename Container>
            void HostDeviceMemoryManagement::CopyDeviceToHostAndDeviceFree( Container& a_vector )
            {
                CopyDeviceToHostAndDeviceFree( a_vector.data(), a_vector.size() );
            }

            template <typename T>
            void HostDeviceMemoryManagement::DeviceFree( T* a_host_pointer, std::size_t a_count )
            {
                static_assert( !std::is_const<T>::value, "" );
                static_assert( std::is_pod<T>::value, "" );
                DoDeviceFree( a_host_pointer, a_count, sizeof( T ) );
            }

            template <typename Container>
            void HostDeviceMemoryManagement::DeviceFree( Container& a_vector )
            {
                DeviceFree( a_vector.data(), a_vector.size() );
            }

            template <typename T>
            T* HostDeviceMemoryManagement::GetDevicePointer( T* a_host_pointer )
            {
                return static_cast<T*>( DoGetDevicePointer( static_cast<const void*>( a_host_pointer ) ) );
            }

            template <typename T>
            T* HostDeviceMemoryManagement::GetDeviceOrHostPointer( T* a_host_pointer )
            {
                T* const a_device_pointer = GetDevicePointer( a_host_pointer );
                return a_device_pointer == nullptr ?
                           a_host_pointer : // Not mapped to the GPU
                           a_device_pointer;
            }

            template <typename T>
            bool HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( const T* a_host_pointer )
            {
                // We could optimize the omp version by only using ::omp_target_is_present()
                return GetDevicePointer( a_host_pointer ) != nullptr;
            }

            template <typename T>
            void HostDeviceMemoryManagement::DeviceMemoryCopy( T* a_destination, const T* a_source, std::size_t a_count )
            {
                DoDeviceMemoryCopy( a_destination, a_source, a_count, sizeof( T ) );
            }

        } // namespace gpu
    }     // namespace tools
} // namespace smilei

#endif
