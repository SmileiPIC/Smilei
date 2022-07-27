// ----------------------------------------------------------------------------
//! \file Table.cpp
//
//! \brief Table class function implementation
//
// ----------------------------------------------------------------------------

#include "Table.h"

// -----------------------------------------------------------------------------
// Constructor for Table
// -----------------------------------------------------------------------------
Table::Table()
{
}

// -----------------------------------------------------------------------------
// Destructor for Table
// -----------------------------------------------------------------------------
Table::~Table()
{
    if (data_) {
        delete [] data_;
        data_ = nullptr;
    }
}

// -----------------------------------------------------------------------------
//! Allocation of data_ based on proposed size
// -----------------------------------------------------------------------------
void Table::allocate () {
    
    if (size_ == 0) {
        ERROR("Problem during table allocation, size is equal to 0");
    }
    
    data_ = new double [size_];
}

// -----------------------------------------------------------------------------
//! Bcast data_ and metadata to all MPI processes
// -----------------------------------------------------------------------------
void Table::bcast (SmileiMPI *smpi) {
    
    if (smpi->isMaster()) {
        if (size_ == 0) {
            ERROR("Problem during table bcast, size is equal to 0");
        }
        
        if (!data_) {
            ERROR("Problem during table bcast, data not allocated");
        }
    }
    
    // Position for MPI pack and unack
    int position;
    // buffer size for MPI pack and unpack
    int buf_size;

    // -------------------------------------------------------
    // Bcast of all the parameters
    // We pack everything in a buffer
    // --------------------------------------------------------

    // buffer size
    if( smpi->isMaster() ) {
        MPI_Pack_size( 1, MPI_INT, smpi->world(), &position );
        buf_size = position;
        MPI_Pack_size( 2, MPI_DOUBLE, smpi->world(), &position );
        buf_size += position;
        MPI_Pack_size( size_, MPI_DOUBLE, smpi->world(),
                       &position );
        buf_size += position;
    }

    //MESSAGE( 2,"Buffer size: " << buf_size );

    // Exchange buf_size with all ranks
    MPI_Bcast( &buf_size, 1, MPI_INT, 0, smpi->world() );

    // Packet that will contain all parameters
    char *buffer = new char[buf_size];

    // Proc 0 packs
    if( smpi->isMaster() ) {
        position = 0;
        MPI_Pack( &size_,
                  1, MPI_INT, buffer, buf_size, &position, smpi->world() );
        MPI_Pack( &min_,
                  1, MPI_DOUBLE, buffer, buf_size, &position, smpi->world() );
        MPI_Pack( &max_,
                  1, MPI_DOUBLE, buffer, buf_size, &position, smpi->world() );

        MPI_Pack( &data_[0], size_,
                  MPI_DOUBLE, buffer, buf_size, &position, smpi->world() );

    }

    // Bcast all parameters
    MPI_Bcast( &buffer[0], buf_size, MPI_PACKED, 0, smpi->world() );

    // Other ranks unpack
    if( !smpi->isMaster() ) {
        position = 0;
        MPI_Unpack( buffer, buf_size, &position,
                    &size_, 1, MPI_INT, smpi->world() );
        MPI_Unpack( buffer, buf_size, &position,
                    &min_, 1, MPI_DOUBLE, smpi->world() );
        MPI_Unpack( buffer, buf_size, &position,
                    &max_, 1, MPI_DOUBLE, smpi->world() );

        // Resize table before unpacking values
        allocate( );

        MPI_Unpack( buffer, buf_size, &position, &data_[0],
                    size_, MPI_DOUBLE, smpi->world() );

    }

    delete[] buffer;

    if (min_ != 0) {
        log10_min_ = std::log10( min_ );
    } else {
        log10_min_ = 0;
    }

    // Computation of the delta
    delta_ = ( std::log10( max_ )
                      - log10_min_ )/( size_-1 );

    // Inverse delta
    inv_delta_ = 1./delta_;
}

// -----------------------------------------------------------------------------
//! get value using linear interpolation
// -----------------------------------------------------------------------------
double Table::get (double x) {

    // Position in the niel_.table
    double d = ( std::log10( x )-log10_min_ )*inv_delta_;
    int index = int( std::floor( d ) );

    // distance for interpolation
    d = d - std::floor( d );

    // Linear interpolation
    return data_[index]*( 1.-d ) + data_[index+1]*( d );
}

// -----------------------------------------------------------------------------
//! get value using linear interpolation
// -----------------------------------------------------------------------------
void Table::set (std::vector<double> & input_data) {

    if (input_data.size() != size_) {
        ERROR("Impossible to initialize data in Table " << name_ << " because sizes do not match.")
    }

    for (unsigned int i = 0 ; i < size_ ; i++) {
        data_[i] = input_data[i];
    }

}
