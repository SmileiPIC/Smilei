#include "Field1D.h"

#include <iostream>
#include <vector>
#include <cstring>

#include "Params.h"
#include "SmileiMPI.h"
#include "Patch.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creators for Field1D
// ---------------------------------------------------------------------------------------------------------------------

// with no input argument
Field1D::Field1D() : Field()
{
    sendFields_.resize(2,NULL);
    recvFields_.resize(2,NULL);
}

// with the dimensions as input argument
Field1D::Field1D( vector<unsigned int> dims ) : Field( dims )
{
    allocateDims( dims );
    sendFields_.resize(2,NULL);
    recvFields_.resize(2,NULL);
}

// with the dimensions and output (dump) file name as input argument
Field1D::Field1D( vector<unsigned int> dims, string name_in ) : Field( dims, name_in )
{
    allocateDims( dims );
    sendFields_.resize(2,NULL);
    recvFields_.resize(2,NULL);
}

// with the dimensions as input argument
Field1D::Field1D( vector<unsigned int> dims, unsigned int mainDim, bool isPrimal ) : Field( dims, mainDim, isPrimal )
{
    allocateDims( dims, mainDim, isPrimal );
    sendFields_.resize(2,NULL);
    recvFields_.resize(2,NULL);
}

// with the dimensions and output (dump) file name as input argument
Field1D::Field1D( vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, string name_in ) : Field( dims, mainDim, isPrimal, name_in )
{
    allocateDims( dims, mainDim, isPrimal );
    sendFields_.resize(2,NULL);
    recvFields_.resize(2,NULL);
}

// without allocating
Field1D::Field1D( string name_in, vector<unsigned int> dims ) : Field( dims, name_in )
{
    dims_ = dims;
    number_of_points_ = dims_[0];
    sendFields_.resize(2,NULL);
    recvFields_.resize(2,NULL);
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Field1D
// ---------------------------------------------------------------------------------------------------------------------
Field1D::~Field1D()
{
    for( unsigned int iside=0 ; iside<sendFields_.size() ; iside++ ) {
        if ( sendFields_[iside] != NULL ) {
            delete sendFields_[iside];
            sendFields_[iside] = NULL;
            delete recvFields_[iside];
            recvFields_[iside] = NULL;
        }
    }
    if( data_!=NULL ) {
        delete [] data_;
    }
}


// ---------------------------------------------------------------------------------------------------------------------
// Method used for allocating the dimension of a Field1D
// ---------------------------------------------------------------------------------------------------------------------
void Field1D::allocateDims()
{
    // for Field1D only
    if( dims_.size()!=1 ) {
        ERROR( "Alloc error must be 1 : " << dims_.size() );
    }
    
    isDual_.resize( dims_.size(), 0 );
    
    data_ = new double[ dims_[0] ];
    //! \todo{change to memset (JD)}
    for( unsigned int i=0; i<dims_[0]; i++ ) {
        data_[i]=0.0;
    }
    
    number_of_points_ = dims_[0];
    
}

void Field1D::deallocateDataAndSetTo( Field* f )
{
    delete [] data_;
    data_=NULL;

    data_ = f->data_;
}


void Field1D::allocateDims( unsigned int dims1 )
{
    vector<unsigned int> dims( 1 );
    dims[0]=dims1;
    allocateDims( dims );
}



// ---------------------------------------------------------------------------------------------------------------------
// Method used for allocating the dimension of a Field1D
// ---------------------------------------------------------------------------------------------------------------------
void Field1D::allocateDims( unsigned int mainDim, bool isPrimal )
{
    // for Field1D only
    if( dims_.size()!=1 ) {
        ERROR( "Alloc error must be 1 : " << dims_.size() );
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
    
    data_ = new double[ dims_[0] ];
    //! \todo{change to memset (JD)}
    for( unsigned int i=0; i<dims_[0]; i++ ) {
        data_[i]=0.0;
    }
    
    number_of_points_ = dims_[0];
    
}


// ---------------------------------------------------------------------------------------------------------------------
// Method to shift field in space
// ---------------------------------------------------------------------------------------------------------------------
void Field1D::shift_x( unsigned int delta )
{
    memmove( &( data_[0] ), &( data_[delta] ), ( dims_[0]-delta )*sizeof( double ) );
    //memset ( &(data_[dims_[0]-delta]), 0, delta*sizeof(double));
    for( int i=dims_[0]-delta; i<( int )dims_[0]; i++ ) {
        data_[i] = 0.;
    }
    
}

double Field1D::norm2( unsigned int istart[3][2], unsigned int bufsize[3][2] )
{

    double nrj( 0. );
    
    int idxlocalstart[1];
    int idxlocalend[1];
    for( int i=0 ; i<1 ; i++ ) {
        idxlocalstart[i] = istart[i][isDual_[i]];
        idxlocalend[i]   = istart[i][isDual_[i]]+bufsize[i][isDual_[i]];
    }
    for( int i=idxlocalstart[0] ; i<idxlocalend[0] ; i++ ) {
        nrj += data_[i]*data_[i];
    }
    
    return nrj;
}
//! Perform the norm2 on Device
#if defined(SMILEI_ACCELERATOR_GPU)
double Field1D::norm2OnDevice( unsigned int istart[3][2], unsigned int bufsize[3][2] )
{

    double nrj( 0. );
    
    int idxlocalstart[1];
    int idxlocalend[1];
    idxlocalstart[0] = istart[0][isDual_[0]];
    idxlocalend[0]   = istart[0][isDual_[0]]+bufsize[0][isDual_[0]];

    const double *const __restrict__ field = data();

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target teams distribute parallel for\
		      map(tofrom: nrj)  \
                      map(to: idxlocalstart[0]) \
		      /* is_device_ptr( data_ )*/ \
		      reduction(+:nrj) 
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel present(field) //deviceptr( data_ )
    #pragma acc loop gang worker vector reduction(+:nrj)
#endif

    for( unsigned int i=idxlocalstart[0] ; i< idxlocalend[0] ; i++) {
            nrj += field[i]*field[i];
    }

    return nrj;

}
#endif

void Field1D::put( Field *outField, Params &params, Patch *thisPatch, Patch *outPatch )
{
    Field1D *out1D = static_cast<Field1D *>( outField );
    
    std::vector<unsigned int> dual =  this->isDual_;
    
    int iout = thisPatch->Pcoordinates[0]*params.patch_size_[0] - ( outPatch->getCellStartingGlobalIndex(0) + params.oversize[0] ) ;
    
    for( unsigned int i = 0 ; i < this->dims_[0] ; i++ ) {
        ( *out1D )( iout+i ) = ( *this )( i );
    }
    
}

void Field1D::add( Field *outField, Params &params, Patch *thisPatch, Patch *outPatch )
{
    Field1D *out1D = static_cast<Field1D *>( outField );
    
    std::vector<unsigned int> dual =  this->isDual_;
    
    int iout = thisPatch->Pcoordinates[0]*params.patch_size_[0] - ( outPatch->getCellStartingGlobalIndex(0) + params.oversize[0] ) ;
    
    for( unsigned int i = 0 ; i < this->dims_[0] ; i++ ) {
        ( *out1D )( iout+i ) += ( *this )( i );
    }
    
}

void Field1D::get( Field *inField, Params &params, Patch *inPatch, Patch *thisPatch )
{
    Field1D *in1D  = static_cast<Field1D *>( inField );
    
    std::vector<unsigned int> dual =  in1D->isDual_;
    
    int iin = thisPatch->Pcoordinates[0]*params.patch_size_[0] - ( inPatch->getCellStartingGlobalIndex(0) + params.oversize[0] );
    
    for( unsigned int i = 0 ; i < this->dims_[0] ; i++ ) {
        ( *this )( i ) = ( *in1D )( iin+i );
    }
    
}

void Field1D::create_sub_fields  ( int iDim, int iNeighbor, int ghost_size )
{
    std::vector<unsigned int> size = dims_;
    size[iDim] = ghost_size;
    if ( sendFields_[iDim*2+iNeighbor] == NULL ) {
        sendFields_[iDim*2+iNeighbor] = new Field1D(size);
        recvFields_[iDim*2+iNeighbor] = new Field1D(size);
#if defined( SMILEI_ACCELERATOR_GPU ) 
       if( ( name[0] == 'B' ) || ( name[0] == 'J' || name[0] == 'R' ) ) {
           sendFields_[iDim * 2 + iNeighbor]->allocateAndCopyFromHostToDevice();
           recvFields_[iDim * 2 + iNeighbor]->allocateAndCopyFromHostToDevice();
       }
#endif
    }
    else if( ghost_size != (int) sendFields_[iDim*2+iNeighbor]->dims_[iDim] ) {
#if defined( SMILEI_ACCELERATOR_GPU )
        ERROR( "To Do GPU : envelope" );
#endif
        delete sendFields_[iDim*2+iNeighbor];
        sendFields_[iDim*2+iNeighbor] = new Field1D(size);
        delete recvFields_[iDim*2+iNeighbor];
        recvFields_[iDim*2+iNeighbor] = new Field1D(size);
    }
}
void Field1D::extract_fields_exch( int iDim, int iNeighbor, int ghost_size )
{
    std::vector<unsigned int> size = dims_;
    size[iDim] = ghost_size;

    vector<int> idx( 1, 0 );
    idx[iDim] = 1;
    int istart = iNeighbor * ( dims_[iDim]- ( 2*ghost_size+1+isDual_[iDim] ) ) + ( 1-iNeighbor ) * ( ghost_size + 1 + isDual_[iDim] );
    int ix = idx[0]*istart;

    unsigned int NX = size[0];

    double *__restrict__ sub = sendFields_[iDim*2+iNeighbor]->data_;
    const double*__restrict__ field = data_;

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    // At initialization, this data is NOT on the GPU
    const bool should_manipulate_gpu_memory = name[0] == 'B' &&
                                              smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( sub );
    SMILEI_ASSERT( smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( field ) ==
                   smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( sub ) );
    const unsigned field_first = ix;
    const unsigned field_last  = ix + NX - 1;
    #pragma omp target if( should_manipulate_gpu_memory )
    #pragma omp teams distribute parallel for
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    const int subSize = sendFields_[iDim*2+iNeighbor]->size();
    const int fSize = number_of_points_;
    bool fieldName( (name.substr(0,1) == "B") );
    #pragma acc parallel present( field[0:fSize], sub[0:subSize] ) if (fieldName)
    #pragma acc loop gang worker vector
#endif
    for( unsigned int i=0; i<NX; i++ ) {
        sub[i] = field[ (ix+i) ];
    }
}
void Field1D::inject_fields_exch ( int iDim, int iNeighbor, int ghost_size )
{
    std::vector<unsigned int> size = dims_;
    size[iDim] = ghost_size;

    vector<int> idx( 1, 0 );
    idx[iDim] = 1;
    int istart = ( ( iNeighbor+1 )%2 ) * ( dims_[iDim] - 1- ( ghost_size-1 ) ) + ( 1-( iNeighbor+1 )%2 ) * ( 0 )  ;
    int ix = idx[0]*istart;

    unsigned int NX = size[0];

    const double *__restrict__ sub = recvFields_[iDim*2+(iNeighbor+1)%2]->data_;
    double *__restrict__ field = data_;

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    // At initialization, this data is NOT on the GPU
    const bool should_manipulate_gpu_memory = name[0] == 'B' &&
                                              smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( sub );
    const unsigned field_first = ix;
    const unsigned field_last  = ix + NX - 1;
    #pragma omp target if( should_manipulate_gpu_memory ) \
        map( tofrom : field [field_first:field_last - field_first] )
    #pragma omp teams distribute parallel for
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    int subSize = recvFields_[iDim*2+(iNeighbor+1)%2]->size();
    const int fSize = number_of_points_;
    bool fieldName( name.substr(0,1) == "B" );
    #pragma acc parallel present( field[0:fSize], sub[0:subSize] ) if (fieldName)
    #pragma acc loop gang worker vector
#endif
    for( unsigned int i=0; i<NX; i++ ) {
        field[ (ix+i) ] = sub[i];
    }
}

void Field1D::extract_fields_sum ( int iDim, int iNeighbor, int ghost_size )
{
    std::vector<unsigned int> size = dims_;
    size[iDim] = 2*ghost_size+1+isDual_[iDim];

    vector<int> idx( 1, 0 );
    idx[iDim] = 1;
    int istart = iNeighbor * ( dims_[iDim]- ( 2*ghost_size+1+isDual_[iDim] ) ) + ( 1-iNeighbor ) * 0;
    int ix = idx[0]*istart;

    unsigned int NX = size[0];

    double *__restrict__ sub         = sendFields_[iDim*2+iNeighbor]->data_;
    const double *__restrict__ field = data_;

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    // At initialization, this data is NOT on the GPU
    const bool should_manipulate_gpu_memory = (name[0] == 'J' || name[0] == 'R') &&
                                              smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( sub );
    const unsigned field_first = ix;
    const unsigned field_last  = ix + NX - 1;
    #pragma omp target if( should_manipulate_gpu_memory ) \
        map( to : field [field_first:field_last - field_first] )
    #pragma omp teams distribute parallel for 
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    const int subSize = sendFields_[iDim*2+iNeighbor]->size();
    const int fSize = number_of_points_;
    bool fieldName( ((name.substr(0,1) == "J") || (name.substr(0,1) == "R") ) && smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( sub ));
    #pragma acc parallel copy(field[0:fSize]) present(  sub[0:subSize] ) if (fieldName)
    #pragma acc loop gang worker vector
#endif
    for( unsigned int i=0; i<NX; i++ ) {
        sub[i] = field[ (ix+i) ];
    }
}
void Field1D::inject_fields_sum  ( int iDim, int iNeighbor, int ghost_size )
{
    std::vector<unsigned int> size = dims_;
    size[iDim] = 2*ghost_size+1+isDual_[iDim];

    vector<int> idx( 1, 0 );
    idx[iDim] = 1;
    int istart = ( ( iNeighbor+1 )%2 ) * ( dims_[iDim] - ( 2*ghost_size+1+isDual_[iDim] ) ) + ( 1-( iNeighbor+1 )%2 ) * ( 0 )  ;
    int ix = idx[0]*istart;

    unsigned int NX = size[0];

    const double *__restrict__ sub = recvFields_[iDim*2+(iNeighbor+1)%2]->data_;
    double *__restrict__ field     = data_;

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    // At initialization, this data is NOT on the GPU
    const bool should_manipulate_gpu_memory = (name[0] == 'J' || name[0] == 'R') &&
                                              smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( sub );
    const unsigned field_first = ix;
    const unsigned field_last  = ix + NX - 1;
    #pragma omp target if( should_manipulate_gpu_memory ) \
        map( tofrom : field [field_first:field_last - field_first] )
    #pragma omp teams distribute parallel for
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    int subSize = recvFields_[iDim*2+(iNeighbor+1)%2]->size();
    int fSize = number_of_points_;
    bool fieldName( name.substr(0,1) == "J" || name.substr(0,1) == "R");
    #pragma acc parallel copy(field[0:fSize]) present( sub[0:subSize] ) if (fieldName)
    #pragma acc loop gang worker vector
#endif
    for( unsigned int i=0; i<NX; i++ ) {
        field[ (ix+i) ] += sub[i];
    }
}

