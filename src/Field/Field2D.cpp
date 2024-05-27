#include "Field2D.h"

#include <cstring>
#include <iostream>
#include <vector>

#include "Params.h"
#include "Patch.h"
#include "SmileiMPI.h"
#include "gpu.h"

// ---------------------------------------------------------------------------------------------------------------------
// Creators for Field2D
// ---------------------------------------------------------------------------------------------------------------------

// with no input argument
Field2D::Field2D() : Field()
{
    sendFields_.resize(4,NULL);
    recvFields_.resize(4,NULL);
}

// with the dimensions as input argument
Field2D::Field2D( std::vector<unsigned int> dims ) : Field( dims )
{
    allocateDims( dims );
    sendFields_.resize(4,NULL);
    recvFields_.resize(4,NULL);
}

// with the dimensions and output (dump) file name as input argument
Field2D::Field2D( std::vector<unsigned int> dims, std::string name_in ) : Field( dims, name_in )
{
    allocateDims( dims );
    sendFields_.resize(4,NULL);
    recvFields_.resize(4,NULL);
}

// with the dimensions as input argument
Field2D::Field2D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal ) : Field( dims, mainDim, isPrimal )
{
    allocateDims( dims, mainDim, isPrimal );
    sendFields_.resize(4,NULL);
    recvFields_.resize(4,NULL);
}

// with the dimensions and output (dump) file name as input argument
Field2D::Field2D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, std::string name_in ) : Field( dims, mainDim, isPrimal, name_in )
{
    allocateDims( dims, mainDim, isPrimal );
    sendFields_.resize(4,NULL);
    recvFields_.resize(4,NULL);
}

// without allocating
Field2D::Field2D( std::string name_in, std::vector<unsigned int> dims ) : Field( dims, name_in )
{
    dims_ = dims;
    number_of_points_ = dims_[0]*dims_[1];
    sendFields_.resize(4,NULL);
    recvFields_.resize(4,NULL);
}



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Field2D
// ---------------------------------------------------------------------------------------------------------------------
Field2D::~Field2D()
{
    for (int iside=0 ; iside<(int)(sendFields_.size()) ; iside++ ) {
        if ( sendFields_[iside] != NULL ) {

#if defined ( SMILEI_ACCELERATOR_GPU )
            if ( sendFields_[iside]->isOnDevice() )
            {
                sendFields_[iside]->deleteOnDevice();
            }
            
            if ( recvFields_[iside]->isOnDevice() )
            {
                recvFields_[iside]->deleteOnDevice();
            }
#endif

            delete sendFields_[iside];
            sendFields_[iside] = NULL;
            delete recvFields_[iside];
            recvFields_[iside] = NULL;
        }
    }
    if( data_!=NULL ) {
        delete [] data_;
        delete [] data_2D;
    }
}


// ---------------------------------------------------------------------------------------------------------------------
// Method used for allocating the dimension of a Field2D
// ---------------------------------------------------------------------------------------------------------------------
void Field2D::allocateDims()
{
    if( dims_.size()!=2 ) {
        ERROR( "Alloc error must be 2 : " << dims_.size() );
    }
    if( data_!=NULL ) {
        delete [] data_;
    }
    
    isDual_.resize( dims_.size(), 0 );
    
    data_ = new double[dims_[0]*dims_[1]];
    //! \todo{check row major order!!! (JD)}
    
    data_2D = new double *[dims_[0]];

    // Try not to use data_2D, it breaks the homogeneity of the code
    for( unsigned int i = 0; i < dims_[0]; i++ ) {
        data_2D[i] = data_ + i * dims_[1];
    }
    
    number_of_points_ = dims_[0]*dims_[1];
    
    Field::put_to(0.0);
}

void Field2D::deallocateDataAndSetTo( Field* f )
{
    delete [] data_;
    data_ = NULL;
    delete [] data_2D;
    data_2D = NULL;

    data_   = f->data_;
    data_2D = (static_cast<Field2D *>(f))->data_2D;
    
}

void Field2D::allocateDims( unsigned int dims1, unsigned int dims2 )
{
    std::vector<unsigned int> dims( 2 );
    dims[0]=dims1;
    dims[1]=dims2;
    allocateDims( dims );
}

// ---------------------------------------------------------------------------------------------------------------------
// Method used for allocating the dimension of a Field2D
// ---------------------------------------------------------------------------------------------------------------------
void Field2D::allocateDims( unsigned int mainDim, bool isPrimal )
{
    if( dims_.size()!=2 ) {
        ERROR( "Alloc error must be 2 : " << dims_.size() );
    }
    if( data_ ) {
        delete [] data_;
    }
    
    // isPrimal define if mainDim is Primal or Dual
    isDual_.resize( dims_.size(), 0 );
    for( unsigned int j=0 ; j<dims_.size() ; j++ ) {
        if( ( j==mainDim ) && ( !isPrimal ) ) {
            isDual_[j] = 1;
        } else if( ( j!=mainDim ) && ( isPrimal ) ) {
            isDual_[j] = 1;
        }
    }
    
    for( unsigned int j=0 ; j<dims_.size() ; j++ ) {
        dims_[j] += isDual_[j];
    }
    
    data_ = new double[dims_[0]*dims_[1]];
    //! \todo{check row major order!!! (JD)}

    data_2D = new double *[dims_[0]];

    // Try not to use data_2D, it breaks the homogeneity of the code
    for( unsigned int i = 0; i < dims_[0]; i++ ) {
        data_2D[i] = data_ + i * dims_[1];
    }

    number_of_points_ = dims_[0]*dims_[1];

    Field::put_to(0.0);
}



// ---------------------------------------------------------------------------------------------------------------------
// Method to shift field in space
// ---------------------------------------------------------------------------------------------------------------------
void Field2D::shift_x( unsigned int delta )
{
    memmove( &( data_2D[0][0] ), &( data_2D[delta][0] ), ( dims_[1]*dims_[0]-delta*dims_[1] )*sizeof( double ) );
    memset( &( data_2D[dims_[0]-delta][0] ), 0, delta*dims_[1]*sizeof( double ) );
    
}

double Field2D::norm2( unsigned int istart[3][2], unsigned int bufsize[3][2] )
{
    double nrj( 0. );
    
    int idxlocalstart[2];
    int idxlocalend[2];
    for( int i=0 ; i<2 ; i++ ) {
        idxlocalstart[i] = istart[i][isDual_[i]];
        idxlocalend[i]   = istart[i][isDual_[i]]+bufsize[i][isDual_[i]];
    }

    for( int i=idxlocalstart[0] ; i<idxlocalend[0] ; i++ ) {
        #pragma omp simd reduction(+:nrj)
        for( int j=idxlocalstart[1] ; j<idxlocalend[1] ; j++ ) {
            nrj += data_2D[i][j]*data_2D[i][j];
        }
    }
    
    return nrj;
}

//! Perform the norm2 on Device
#if defined(SMILEI_ACCELERATOR_GPU)
double Field2D::norm2OnDevice( unsigned int istart[3][2], unsigned int bufsize[3][2] )
{

    //std::cout << name << std::endl;

    double nrj( 0. );
    
    int idxlocalstart[2];
    int idxlocalend[2];
    for( int i=0 ; i<2 ; i++ ) {
        idxlocalstart[i] = istart[i][isDual_[i]];
        idxlocalend[i]   = istart[i][isDual_[i]]+bufsize[i][isDual_[i]];
    }
    
    const unsigned iystart = idxlocalstart[1];
    const unsigned iyend = idxlocalend[1];
    const unsigned ny = dims_[1];

    const double *const __restrict__ field = data();

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target teams distribute parallel for collapse(2) \
		      map(tofrom: nrj)  \
                      map(to: ny, idxlocalstart[0], idxlocalstart[1], iystart, iyend) \
		      /* is_device_ptr( data_ )*/ \
		      reduction(+:nrj) 
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel present(field) //deviceptr( data_ )
    #pragma acc loop gang worker vector collapse(2) reduction(+:nrj)
#endif

    for( unsigned int i=idxlocalstart[0] ; i< idxlocalend[0] ; i++) {
        for( unsigned int j=iystart ; j< iyend ; j++ ) {
            nrj += field[i*ny + j]*field[i*ny + j];
            //nrj += 1;
        }
    }

    //std::cout << nrj << std::endl;

    return nrj;
}
#endif

void Field2D::put( Field *outField, Params &params, Patch *thisPatch, Patch *outPatch )
{
    Field2D *out2D = static_cast<Field2D *>( outField );
    
    std::vector<unsigned int> dual =  this->isDual_;
    
    int iout = thisPatch->Pcoordinates[0]*params.patch_size_[0] - ( outPatch->getCellStartingGlobalIndex(0) + params.region_oversize[0] ) ;
    int jout = thisPatch->Pcoordinates[1]*params.patch_size_[1] - ( outPatch->getCellStartingGlobalIndex(1) + params.region_oversize[1] ) ;
    
    //for ( unsigned int i = params.oversize[0] ; i < this->dims_[0]-params.oversize[0] ; i++ ) {
    //    for ( unsigned int j = params.oversize[1] ; j < this->dims_[1]-params.oversize[1] ; j++ ) {
    for( unsigned int i = 0 ; i < params.patch_size_[0]+1+dual[0]+2*params.oversize[0] ; i++ ) {
        for( unsigned int j = 0 ; j < params.patch_size_[1]+1+dual[1]+2*params.oversize[1] ; j++ ) {
            ( *out2D )( iout+i+params.region_oversize[0]-params.oversize[0], jout+j+params.region_oversize[1]-params.oversize[1] ) = ( *this )( i, j );
        }
    }
    
}


void Field2D::add( Field *outField, Params &params, Patch *thisPatch, Patch *outPatch )
{
    Field2D *out2D = static_cast<Field2D *>( outField );
    
    std::vector<unsigned int> dual =  this->isDual_;
    
    int iout = thisPatch->Pcoordinates[0]*params.patch_size_[0] - ( outPatch->getCellStartingGlobalIndex(0) + params.region_oversize[0] ) ;
    int jout = thisPatch->Pcoordinates[1]*params.patch_size_[1] - ( outPatch->getCellStartingGlobalIndex(1) + params.region_oversize[1] ) ;
    
    //for ( unsigned int i = params.oversize[0] ; i < this->dims_[0]-params.oversize[0] ; i++ ) {
    //    for ( unsigned int j = params.oversize[1] ; j < this->dims_[1]-params.oversize[1] ; j++ ) {
    for( unsigned int i = 0 ; i < params.patch_size_[0]+1+dual[0]+2*params.oversize[0] ; i++ ) {
        for( unsigned int j = 0 ; j < params.patch_size_[1]+1+dual[1]+2*params.oversize[1] ; j++ ) {
            ( *out2D )( iout+i+params.region_oversize[0]-params.oversize[0], jout+j+params.region_oversize[1]-params.oversize[1] ) += ( *this )( i, j );
        }
    }
    
}

void Field2D::get( Field *inField, Params &params, Patch *inPatch, Patch *thisPatch )
{
    Field2D *in2D  = static_cast<Field2D *>( inField );
    
    std::vector<unsigned int> dual =  in2D->isDual_;
    
    int iin = thisPatch->Pcoordinates[0]*params.patch_size_[0] - ( inPatch->getCellStartingGlobalIndex(0) + params.region_oversize[0] );
    int jin = thisPatch->Pcoordinates[1]*params.patch_size_[1] - ( inPatch->getCellStartingGlobalIndex(1) + params.region_oversize[1] );
    
    //for ( unsigned int i = params.oversize[0] ; i < out2D->dims_[0]-params.oversize[0] ; i++ ) {
    //    for ( unsigned int j = params.oversize[1] ; j < out2D->dims_[1]-params.oversize[1] ; j++ ) {
    for( unsigned int i = 0 ; i < params.patch_size_[0]+1+dual[0]+2*params.oversize[0] ; i++ ) {
        for( unsigned int j = 0 ; j < params.patch_size_[1]+1+dual[1]+2*params.oversize[1] ; j++ ) {
            ( *this )( i, j ) = ( *in2D )( iin+i+params.region_oversize[0]-params.oversize[0], jin+j+params.region_oversize[1]-params.oversize[1] );
            //( *out2D )( i, j ) = in2D->hindex;
        }
    }
    
}

void Field2D::create_sub_fields( int iDim, int iNeighbor, int ghost_size )
{
    std::vector<unsigned int> size = dims_;
    size[iDim] = ghost_size;
    if ( sendFields_[iDim*2+iNeighbor] == NULL ) {
        // Allocate only the first time we call this function
        sendFields_[iDim*2+iNeighbor] = new Field2D(size);
        recvFields_[iDim*2+iNeighbor] = new Field2D(size);

#if defined( SMILEI_ACCELERATOR_GPU ) 
       if( ( name[0] == 'B' ) || ( name[0] == 'J' || name[0] == 'R' ) ) {
           sendFields_[iDim * 2 + iNeighbor]->allocateAndCopyFromHostToDevice();
           recvFields_[iDim * 2 + iNeighbor]->allocateAndCopyFromHostToDevice();
       }
#endif
    } 
    else if ( ghost_size != (int)(sendFields_[iDim*2+iNeighbor]->dims_[iDim]) ) {
#if defined( SMILEI_ACCELERATOR_GPU_OACC ) || defined( SMILEI_ACCELERATOR_GPU_OMP )
        ERROR( "To Do GPU : envelope" );
#endif
        delete sendFields_[iDim*2+iNeighbor];
        sendFields_[iDim*2+iNeighbor] = new Field2D(size);
        delete recvFields_[iDim*2+iNeighbor];
        recvFields_[iDim*2+iNeighbor] = new Field2D(size);
    }
}

void Field2D::extract_fields_exch( int iDim, int iNeighbor, int ghost_size )
{
    std::vector<unsigned int> size = dims_;
    size[iDim] = ghost_size;

    std::vector<int> idx( 2, 0 );
    idx[iDim] = 1;
    int istart = iNeighbor * ( dims_[iDim]- ( 2*ghost_size+1+isDual_[iDim] ) ) + ( 1-iNeighbor ) * ( ghost_size + 1 + isDual_[iDim] );
    int ix = idx[0]*istart;
    int iy = idx[1]*istart;

    unsigned int NX = size[0];
    unsigned int NY = size[1];

    int dimY = dims_[1];

    double *__restrict__ sub         = sendFields_[iDim * 2 + iNeighbor]->data_;
    const double *__restrict__ field = data_;

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    // At initialization, this data is NOT on the GPU
    const bool should_manipulate_gpu_memory = name[0] == 'B' &&
                                              smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( sub );
    SMILEI_ASSERT( smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( field ) ==
                   smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( sub ) );
    const unsigned field_first = ix * dimY + iy;
    const unsigned field_last  = ( ix + NX - 1 ) * dimY + iy + NY;

    #pragma omp target if( should_manipulate_gpu_memory )
    #pragma omp teams distribute parallel for collapse( 2 )
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    const int subSize = sendFields_[iDim*2+iNeighbor]->size();
    const int fSize = number_of_points_;
    bool fieldName( (name.substr(0,1) == "B") );
    #pragma acc parallel present( field[0:fSize], sub[0:subSize] ) if (fieldName)
    #pragma acc loop gang
#endif
    for( unsigned int i=0; i<NX; i++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
	#pragma acc loop worker vector
#endif
        for( unsigned int j=0; j<NY; j++ ) {
            sub[i*NY+j] = field[ (ix+i)*dimY+(iy+j) ];
        }
    }
}

void Field2D::inject_fields_exch ( int iDim, int iNeighbor, int ghost_size )
{
    std::vector<unsigned int> size = dims_;
    size[iDim] = ghost_size;

    std::vector<int> idx( 2, 0 );
    idx[iDim] = 1;
    int istart = ( ( iNeighbor+1 )%2 ) * ( dims_[iDim] - 1- ( ghost_size-1 ) ) + ( 1-( iNeighbor+1 )%2 ) * ( 0 )  ;
    int ix = idx[0]*istart;
    int iy = idx[1]*istart;

    unsigned int NX = size[0];
    unsigned int NY = size[1];

    int dimY = dims_[1];

    const double *__restrict__ sub = recvFields_[iDim * 2 + ( iNeighbor + 1 ) % 2]->data_;
    double *__restrict__ field     = data_;

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    // At initialization, this data is NOT on the GPU
    const bool should_manipulate_gpu_memory = name[0] == 'B' &&
                                              smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( sub );
    // SMILEI_ASSERT( /*iteration == 0 && */!smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( field ) );
    const unsigned field_first = ix * dimY + iy;
    const unsigned field_last  = ( ix + NX - 1 ) * dimY + iy + NY;

    #pragma omp target if( should_manipulate_gpu_memory ) \
        map( tofrom                                       \
             : field [field_first:field_last - field_first] )
    #pragma omp teams distribute parallel for collapse( 2 )
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    int subSize = recvFields_[iDim*2+(iNeighbor+1)%2]->size();
    const int fSize = number_of_points_;
    bool fieldName( name.substr(0,1) == "B" );
    #pragma acc parallel present( field[0:fSize], sub[0:subSize] ) if (fieldName)
    #pragma acc loop gang
#endif
    for( unsigned int i=0; i<NX; i++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
	#pragma acc loop worker vector
#endif
        for( unsigned int j=0; j<NY; j++ ) {
            field[ (ix+i)*dimY+(iy+j) ] = sub[i*NY+j];
        }
    }
}

void Field2D::extract_fields_sum ( int iDim, int iNeighbor, int ghost_size )
{
    std::vector<unsigned int> size = dims_;
    size[iDim] = 2*ghost_size+1+isDual_[iDim];

    std::vector<int> idx( 2, 0 );
    idx[iDim] = 1;
    int istart = iNeighbor * ( dims_[iDim]- ( 2*ghost_size+1+isDual_[iDim] ) ) + ( 1-iNeighbor ) * 0;
    int ix = idx[0]*istart;
    int iy = idx[1]*istart;

    unsigned int NX = size[0];
    unsigned int NY = size[1];

    int dimY = dims_[1];

    double *__restrict__ sub         = sendFields_[iDim * 2 + iNeighbor]->data_;
    const double *__restrict__ field = data_;

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    // At initialization, this data is NOT on the GPU
    const bool should_manipulate_gpu_memory = (name[0] == 'J' || name[0] == 'R') &&
                                              smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( sub );
    // SMILEI_ASSERT( iteration == 0 && !smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( field ) );
    const unsigned field_first = ix * dimY + iy;
    const unsigned field_last  = ( ix + NX - 1 ) * dimY + iy + NY;

    #pragma omp target if( should_manipulate_gpu_memory ) \
        map( to                                           \
             : field [field_first:field_last - field_first] )
    #pragma omp teams distribute parallel for collapse( 2 )
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    const int subSize = sendFields_[iDim*2+iNeighbor]->size();
    const int fSize = number_of_points_;
    bool fieldName( ((name.substr(0,1) == "J") || (name.substr(0,1) == "R") ) && smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( sub ));
    #pragma acc parallel copy(field[0:fSize]) present(  sub[0:subSize] ) if (fieldName)
    //#pragma acc parallel present( field[0:fSize], sub[0:subSize] ) if (fieldName)
    #pragma acc loop gang
#endif
    for( unsigned int i=0; i<NX; i++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
	#pragma acc loop worker vector
#endif
        for( unsigned int j=0; j<NY; j++ ) {
            sub[i*NY+j] = field[ (ix+i)*dimY+(iy+j) ];
        }
    }
}

void Field2D::inject_fields_sum  ( int iDim, int iNeighbor, int ghost_size )
{
    std::vector<unsigned int> size = dims_;
    size[iDim] = 2*ghost_size+1+isDual_[iDim];

    std::vector<int> idx( 2, 0 );
    idx[iDim] = 1;
    int istart = ( ( iNeighbor+1 )%2 ) * ( dims_[iDim] - ( 2*ghost_size+1+isDual_[iDim] ) ) + ( 1-( iNeighbor+1 )%2 ) * ( 0 )  ;
    int ix = idx[0]*istart;
    int iy = idx[1]*istart;

    unsigned int NX = size[0];
    unsigned int NY = size[1];

    int dimY = dims_[1];

    const double *__restrict__ sub = recvFields_[iDim * 2 + ( iNeighbor + 1 ) % 2]->data_;
    double *__restrict__ field     = data_;

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    // At initialization, this data is NOT on the GPU
    const bool should_manipulate_gpu_memory = (name[0] == 'J' || name[0] == 'R') &&
                                              smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( sub );
    // SMILEI_ASSERT( iteration == 0 && !smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( field ) );
    const unsigned field_first = ix * dimY + iy;
    const unsigned field_last  = ( ix + NX - 1 ) * dimY + iy + NY;

    #pragma omp target if( should_manipulate_gpu_memory ) \
        map( tofrom                                       \
             : field [field_first:field_last - field_first] )
    #pragma omp teams distribute parallel for collapse( 2 )
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    int subSize = recvFields_[iDim*2+(iNeighbor+1)%2]->size();
    int fSize = number_of_points_;
    bool fieldName( name.substr(0,1) == "J" || name.substr(0,1) == "R");
    #pragma acc parallel copy(field[0:fSize]) present(  sub[0:subSize] ) if (fieldName)
    //#pragma acc parallel present( field[0:fSize], sub[0:subSize] ) if (fieldName)
    #pragma acc loop gang
#endif
    for( unsigned int i=0; i<NX; i++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
	#pragma acc loop worker vector
#endif
        for( unsigned int j=0; j<NY; j++ ) {
            field[ (ix+i)*dimY+(iy+j) ] += sub[i*NY+j];
        }
    }
}
