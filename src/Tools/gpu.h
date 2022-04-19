#ifndef SMILEI_TOOLS_GPU_H
#define SMILEI_TOOLS_GPU_H

#include "Tools.h"

namespace smilei {
    namespace tools {

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

        template <typename T,
                  bool do_device_alloc = false>
        class NonInitializingVector
        {
        public:
            using value_type = T;

        public:
            NonInitializingVector();

            explicit NonInitializingVector( std::size_t size );

            void Alloc( std::size_t size );
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
        ///
        struct HostDeviceMemoryManagment
        {
        public:
            template <typename T>
            static void DeviceAlloc( const T* a_pointer, std::size_t a_size );
            template <typename Container>
            static void DeviceAlloc( const Container& a_vector );

            template <typename T>
            static void DeviceAllocAndCopyHostToDevice( const T* a_pointer, std::size_t a_size );
            template <typename Container>
            static void DeviceAllocAndCopyHostToDevice( const Container& a_vector );

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
        };


        ////////////////////////////////////////////////////////////////////////////////
        // NonInitializingVector methods definition
        ////////////////////////////////////////////////////////////////////////////////

        template <typename T,
                  bool do_device_alloc>
        NonInitializingVector<T, do_device_alloc>::NonInitializingVector()
            : size_{}
            , data_{ nullptr }
        {
            // EMPTY
        }

        template <typename T,
                  bool do_device_alloc>
        NonInitializingVector<T, do_device_alloc>::NonInitializingVector( std::size_t size )
            : NonInitializingVector{}
        {
            Alloc( size );
        }

        template <typename T,
                  bool do_device_alloc>
        void NonInitializingVector<T, do_device_alloc>::Alloc( std::size_t size )
        {
            if( size_ != 0 || data_ != nullptr ) {
                ERROR( "NonInitializingVector::Alloc, allocation before dealloc" );
            }

            data_ = static_cast<T*>( std::malloc( sizeof( T ) * size ) );
            if( data_ == nullptr ) {
                ERROR( "NonInitializingVector::Alloc, std::malloc() out of memory." );
            }

            size_ = size;

            if( do_device_alloc ) {
                HostDeviceMemoryManagment::DeviceAlloc( *this );
            }
        }

        template <typename T,
                  bool do_device_alloc>
        void NonInitializingVector<T, do_device_alloc>::Free()
        {
            std::free( data_ );

            if( do_device_alloc ) {
                HostDeviceMemoryManagment::DeviceFree( *this );
            }

            data_ = nullptr;
            size_ = 0;
        }

        template <typename T,
                  bool do_device_alloc>
        T* NonInitializingVector<T, do_device_alloc>::begin()
        {
            return data_;
        }

        template <typename T,
                  bool do_device_alloc>
        const T* NonInitializingVector<T, do_device_alloc>::cbegin() const
        {
            return data_;
        }

        template <typename T,
                  bool do_device_alloc>
        T* NonInitializingVector<T, do_device_alloc>::end()
        {
            return data_ + size_;
        }

        template <typename T,
                  bool do_device_alloc>
        const T* NonInitializingVector<T, do_device_alloc>::cend() const
        {
            return data_ + size_;
        }

        template <typename T,
                  bool do_device_alloc>
        T* NonInitializingVector<T, do_device_alloc>::data()
        {
            return data_;
        }

        template <typename T,
                  bool do_device_alloc>
        const T* NonInitializingVector<T, do_device_alloc>::data() const
        {
            return data_;
        }

        template <typename T,
                  bool do_device_alloc>
        T& NonInitializingVector<T, do_device_alloc>::operator[]( std::size_t index )
        {
            return data()[index];
        }

        template <typename T,
                  bool do_device_alloc>
        const T&
        NonInitializingVector<T, do_device_alloc>::operator[]( std::size_t index ) const
        {
            return data()[index];
        }

        template <typename T,
                  bool do_device_alloc>
        std::size_t NonInitializingVector<T, do_device_alloc>::size() const
        {
            return size_;
        }

        template <typename T,
                  bool do_device_alloc>
        NonInitializingVector<T, do_device_alloc>::~NonInitializingVector()
        {
            Free();
        }


        ////////////////////////////////////////////////////////////////////////////////
        // HostDeviceMemoryManagment methods definition
        ////////////////////////////////////////////////////////////////////////////////

        template <typename T>
        void HostDeviceMemoryManagment::DeviceAlloc( const T* a_pointer, std::size_t a_size )
        {
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target enter data map( alloc \
                                       : a_pointer [0:a_size] )
// #elif defined( _GPU )
#else
            SMILEI_UNUSED( a_pointer );
            SMILEI_UNUSED( a_size );
#endif
        }

        template <typename Container>
        void HostDeviceMemoryManagment::DeviceAlloc( const Container& a_vector )
        {
            DeviceAlloc( a_vector.data(), a_vector.size() );
        }

        template <typename T>
        void HostDeviceMemoryManagment::DeviceAllocAndCopyHostToDevice( const T* a_pointer, std::size_t a_size )
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
        void HostDeviceMemoryManagment::DeviceAllocAndCopyHostToDevice( const Container& a_vector )
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
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target exit data map( from \
                                      : a_pointer [0:a_size] )
#elif defined( _GPU )
    #pragma acc exit data delete( a_pointer [0:a_size] )
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
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target exit data map( delete \
                                      : a_pointer [0:a_size] )
// #elif defined( _GPU )
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

    } // namespace tools
} // namespace smilei

#endif
