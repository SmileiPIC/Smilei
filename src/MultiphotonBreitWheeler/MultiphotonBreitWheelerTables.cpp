// ----------------------------------------------------------------------------
//! \file MultiphotonBreitWheelerTables.cpp
//
//! \brief This file is for the methods of the class
//! MultiphotonBreitWheelerTables.
//! This class contains the methods and tools to generate and manage
//! the physical tables for the multiphoton Breit-wheeler process.
//
//! \details
//! The implementation is adapted from the thesis results of M. Lobet
//! See http://www.theses.fr/2015BORD0361
// ----------------------------------------------------------------------------
 
#include "MultiphotonBreitWheelerTables.h"
#include "MultiphotonBreitWheelerTablesDefault.h"
#include "H5.h"

// -----------------------------------------------------------------------------
// INITILIZATION AND DESTRUCTION
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Constructor for MultiphotonBreitWheelerTables
// -----------------------------------------------------------------------------
MultiphotonBreitWheelerTables::MultiphotonBreitWheelerTables()
{
}

// -----------------------------------------------------------------------------
// Destructor for MultiphotonBreitWheelerTables
// -----------------------------------------------------------------------------
MultiphotonBreitWheelerTables::~MultiphotonBreitWheelerTables()
{
}

// -----------------------------------------------------------------------------
//
//! Initialization of the parameters for the multiphoton Breit-Wheeler pair
//! creation
//
//! \param params Object Params for the parameters from the input script
// -----------------------------------------------------------------------------
void MultiphotonBreitWheelerTables::initialization( Params &params, SmileiMPI *smpi )
{
    if( params.has_multiphoton_Breit_Wheeler_ ) {
        TITLE( "Initializing multiphoton Breit-Wheeler" )

        // Preliminary checks
        if( params.reference_angular_frequency_SI <= 0. ) {
            ERROR( "The parameter `reference_angular_frequency_SI` needs "
                   << "to be defined and positive to compute radiation losses" );
        }

    }

    // If the namelist for Nonlinear Inverse Compton Scattering exists
    // We read the properties
    if( PyTools::nComponents( "MultiphotonBreitWheeler" ) ) {
        // Path to the databases
        PyTools::extract( "table_path", table_path_, "MultiphotonBreitWheeler"  );
    }

    // Computation of some parameters
    if( params.has_multiphoton_Breit_Wheeler_ ) {
        // Computation of the normalized Compton wavelength
        normalized_Compton_wavelength_ = params.red_planck_cst*params.reference_angular_frequency_SI
                                         / ( params.electron_mass*params.c_vacuum_*params.c_vacuum_ );
        // Computation of the factor factor_dNBW_dt_
        factor_dNBW_dt_ = params.fine_struct_cst/( normalized_Compton_wavelength_ * M_PI*std::sqrt( 3.0 ) );
    }

    // Messages and checks
    if( params.has_multiphoton_Breit_Wheeler_ ) {
        if (table_path_.size() > 0) {
            MESSAGE( 1,"Reading of the external database, path: " << table_path_ );
            readTables( params, smpi );
        } else {
            MESSAGE(1,"Default tables (stored in the code) are used:");
            MultiphotonBreitWheelerTablesDefault::setDefault( T_, xi_ );
        }

        MESSAGE( "" )
        MESSAGE( 1,"--- Table `integration_dt_dchi`:" );
        MESSAGE( 2,"Dimension quantum parameter: "
                 << T_.size_ );
        MESSAGE( 2,"Minimum photon quantum parameter chi: "
                 << T_.min_ );
        MESSAGE( 2,"Maximum photon quantum parameter chi: "
                 << T_.max_ );

        MESSAGE( "" )
        MESSAGE( 1,"--- Table `min_particle_chi_for_xi` and `xi`:" );
        MESSAGE( 2,"Dimension photon chi: " << xi_.dim_size_[0] );
        MESSAGE( 2,"Dimension particle chi: " << xi_.dim_size_[1] );
        MESSAGE( 2,"Minimum photon chi: " << xi_.min_ );
        MESSAGE( 2,"Maximum photon chi: " << xi_.max_ );

    }

}


// -----------------------------------------------------------------------------
//! Computation of the electron and positron quantum parameters for
//! the multiphoton Breit-Wheeler pair creation
//
//! \param photon_chi photon quantum parameter
// -----------------------------------------------------------------------------
void MultiphotonBreitWheelerTables::computePairQuantumParameter( const double photon_chi, 
                                                                 double * pair_chi,
                                                                 const double xip )
{

    // -----------------------------------------------------------
    // Computation of the index associated to the given photon_chi
    // -----------------------------------------------------------

    int ichiph;
    const double logchiph = std::log10( photon_chi );

    // Lower boundary of the table
    if( photon_chi < xi_.min_ ) {
        ichiph = 0;
    }
    // Upper boundary of the table
    else if( photon_chi >= xi_.max_ ) {
        ichiph = xi_.dim_size_[0]-1;
    }
    // Inside the table
    else {
        // Use floor so that photon_chi corresponding to ichiph is <= given photon_chi
        ichiph = int( std::floor( ( logchiph-xi_.log10_min_ )*( xi_.inv_delta_ ) ) );
    }

    // ---------------------------------------
    // Search of the index ichiph for photon_chi
    // ---------------------------------------

    // The array uses the symmetric properties of the T fonction,
    // Cases xip > or <= 0.5 are treated seperatly

    // If xip > 0.5, the electron will bring more energy than the positron
    const double xipp = xip > 0.5 ? 1.0 - xip : xip;    

    // index for particle chi
    int ichipa;

    // check boundaries
    // Lower bound
    if( xipp < xi_.data_[ichiph*xi_.dim_size_[1]] ) {
        ichipa = 0;
    }
    // Upper bound
    else if( xipp >= xi_.data_[( ichiph+1 )*xi_.dim_size_[1]-1] ) {
        ichipa = xi_.dim_size_[1]-2;
    } else {
        // Search for the corresponding index ichipa for xip
        ichipa = userFunctions::searchValuesInMonotonicArray(
                     &xi_.data_[ichiph*xi_.dim_size_[1]], xipp, xi_.dim_size_[1] );
    }

    // Delta for the particle_chi dimension
    const double delta_chipa = ( std::log10( 0.5*photon_chi )-xi_.axis1_min_[ichiph] )
                  * xi_.inv_dim_size_minus_one_[1];

    const int ixip = ichiph*xi_.dim_size_[1]+ ichipa;

    const double log10_chipam = ichipa*delta_chipa + xi_.axis1_min_[ichiph];
    const double log10_chipap = log10_chipam + delta_chipa;

    const double d = ( xipp - xi_.data_[ixip] ) / ( xi_.data_[ixip+1] - xi_.data_[ixip] );

    // If xip > 0.5, the electron will bring more energy than the positron
    if( xip > 0.5 ) {

        // Positron quantum parameter
        pair_chi[1] = std::pow( 10., log10_chipam*( 1.0-d ) + log10_chipap*( d ) );

        // Electron quantum parameter
        pair_chi[0] = photon_chi - pair_chi[1];
    }
    // If xip <= 0.5, the positron will bring more energy than the electron
    else {
        // Electron quantum parameter
        pair_chi[0] = std::pow( 10., log10_chipam*( 1.0-d ) + log10_chipap*( d ) );

        // Positron quantum parameter
        pair_chi[1] = photon_chi - pair_chi[0];
    }
}

// -----------------------------------------------------------------------------
//! Computation of the production rate of pairs per photon
//! \param photon_chi photon quantum parameter
//! \param photon_gamma photon normalized energy
// -----------------------------------------------------------------------------
double MultiphotonBreitWheelerTables::computeBreitWheelerPairProductionRate( 
    const double photon_chi,
    const double photon_gamma )
{
    // final value to return
    double dNBWdt;

    // Log of the photon quantum parameter particle_chi
    const double logchiph = std::log10( photon_chi );

    // Lower index for interpolation in the table integfochi
    int ichiph = int( std::floor( ( logchiph-T_.log10_min_ )
                         *T_.inv_delta_ ) );

    // If photon_chi is below the lower bound of the table
    // An asymptotic approximation is used
    if( ichiph < 0 ) {
        ichiph = 0;
        // 0.2296 * sqrt(3) * pi [MG/correction by Antony]
        dNBWdt = 1.2493450020845291*std::exp( -8.0/( 3.0*photon_chi ) ) * photon_chi*photon_chi;
    }
    // If photon_chi is above the upper bound of the table
    // An asymptotic approximation is used
    else if( (unsigned int) ichiph >= T_.size_-1 ) {
        ichiph = T_.size_-2;
        dNBWdt = 2.067731275227008 * cbrt(photon_chi*photon_chi*photon_chi*photon_chi*photon_chi);
    } else {
        // Upper and lower values for linear interpolation
        const double logchiphm = ichiph*T_.delta_ + T_.log10_min_;
        const double logchiphp = logchiphm + T_.delta_;

        // Interpolation
        dNBWdt = ( T_.data_[ichiph+1]*std::fabs( logchiph-logchiphm ) +
                   T_.data_[ichiph]*std::fabs( logchiphp - logchiph ) )*T_.inv_delta_;
    }
    return factor_dNBW_dt_*dNBWdt/(photon_chi*photon_gamma);
}


// -----------------------------------------------------------------------------
// TABLE READING
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//! Read the external table T
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
void MultiphotonBreitWheelerTables::readTableT( SmileiMPI *smpi )
{
    std::string file = table_path_ + "/multiphoton_Breit_Wheeler_tables.h5";
    if( Tools::fileExists( file ) ) {

        if( smpi->isMaster() ) {

            H5Read f( file );

            unsigned int size;

            // First, we read attributes
            H5Read d = f.dataset( "integration_dt_dchi" );
            d.attr( "size_photon_chi", size );
            d.attr( "min_photon_chi", T_.min_ );
            d.attr( "max_photon_chi", T_.max_ );
            
            // Initialize table parameters
            T_.set_size(&size);

            // Allocate data
            T_.allocate();
            
            // Fill data
            f.vect( "integration_dt_dchi", T_.data_[0], H5T_NATIVE_DOUBLE, 0, T_.size_ );

            // Resize and read table
            // T_.table_.resize( T_.size_photon_chi_ );
            // f.vect( "integration_dt_dchi", T_.table_ );

        }
    }
    else
    {
        ERROR("The table T for the nonlinear Breit-Wheeler pair process could not be read from the provided path: `"
              << table_path_<<"`. Please check that the path is correct.")
    }

    // Bcast the table to all MPI ranks
    // MultiphotonBreitWheelerTables::bcastTableT( smpi );
    T_.bcast(smpi);

}

// -----------------------------------------------------------------------------
//! Read the external table xip_chipamin and xip
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
void MultiphotonBreitWheelerTables::readTableXi( SmileiMPI *smpi )
{
    std::string file = table_path_ + "/multiphoton_Breit_Wheeler_tables.h5";
    if( Tools::fileExists( file ) ) {

        if( smpi->isMaster() ) {

            H5Read f( file );

            // First, we read attributes
            H5Read xi = f.dataset( "xi" );
            
            unsigned int dim_size[2];
            
            xi.attr( "size_photon_chi", dim_size[0] );
            xi.attr( "size_particle_chi", dim_size[1] );
            xi.attr( "min_photon_chi", xi_.min_ );
            xi.attr( "max_photon_chi", xi_.max_ );
            
            xi_.set_size(&dim_size[0]);

            // Allocate and read arrays
            // xi_.min_particle_chi_.resize( xi_.size_photon_chi_ );
            // xi_.table_.resize( xi_.size_particle_chi_*xi_.size_photon_chi_ );
            // f.vect( "min_particle_chi_for_xi", xi_.min_particle_chi_ );
            // f.vect( "xi", xi_.table_ );

            xi_.allocate();
            f.vect( "min_particle_chi_for_xi", xi_.axis1_min_[0], H5T_NATIVE_DOUBLE );
            f.vect( "xi", xi_.data_[0], H5T_NATIVE_DOUBLE );

        }
    }
    else
    {
        ERROR("The tables chipamin and xip for the nonlinear Breit-Wheeler pair"
              << " process could not be read from the provided path: `"
              << table_path_<<"`. Please check that the path is correct.")
    }

    // Bcast the table to all MPI ranks
    xi_.bcast( smpi );

}

// -----------------------------------------------------------------------------
//! Read all tables
//
//! \param params list of simulation parameters
//! \param smpi MPI parameters
// -----------------------------------------------------------------------------
void MultiphotonBreitWheelerTables::readTables( Params &params, SmileiMPI *smpi )
{
    // These tables are loaded only if if one species has Monte-Carlo Compton radiation
    // And if the h values are not computed from a numerical fit
    if( params.has_multiphoton_Breit_Wheeler_) {
        MultiphotonBreitWheelerTables::readTableT( smpi );
        MultiphotonBreitWheelerTables::readTableXi( smpi );
    }
}

// -----------------------------------------------------------------------------
// TABLE COMMUNICATIONS
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//! Bcast of the external table T
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
// void MultiphotonBreitWheelerTables::bcastTableT( SmileiMPI *smpi )
// {
//     // Position for MPI pack and unack
//     int position;
//     // buffer size for MPI pack and unpack
//     int buf_size;
// 
//     // -------------------------------------------------------
//     // Bcast of all the parameters
//     // We pack everything in a buffer
//     // --------------------------------------------------------
// 
//     // buffer size
//     if( smpi->getRank() == 0 ) {
//         MPI_Pack_size( 1, MPI_INT, smpi->world(), &position );
//         buf_size = position;
//         MPI_Pack_size( 2, MPI_DOUBLE, smpi->world(), &position );
//         buf_size += position;
//         MPI_Pack_size( T_.size_, MPI_DOUBLE, smpi->world(),
//                        &position );
//         buf_size += position;
//     }
// 
//     //MESSAGE( 2,"Buffer size: " << buf_size );
// 
//     // Exchange buf_size with all ranks
//     MPI_Bcast( &buf_size, 1, MPI_INT, 0, smpi->world() );
// 
//     // Packet that will contain all parameters
//     char *buffer = new char[buf_size];
// 
//     // Proc 0 packs
//     if( smpi->getRank() == 0 ) {
//         position = 0;
//         MPI_Pack( &T_.size_photon_chi_,
//                   1, MPI_INT, buffer, buf_size, &position, smpi->world() );
//         MPI_Pack( &T_.min_photon_chi_,
//                   1, MPI_DOUBLE, buffer, buf_size, &position, smpi->world() );
//         MPI_Pack( &T_.max_photon_chi_,
//                   1, MPI_DOUBLE, buffer, buf_size, &position, smpi->world() );
// 
//         MPI_Pack( &T_.table_[0], T_.size_photon_chi_,
//                   MPI_DOUBLE, buffer, buf_size, &position, smpi->world() );
// 
//     }
// 
//     // Bcast all parameters
//     MPI_Bcast( &buffer[0], buf_size, MPI_PACKED, 0, smpi->world() );
// 
//     // Other ranks unpack
//     if( smpi->getRank() != 0 ) {
//         position = 0;
//         MPI_Unpack( buffer, buf_size, &position,
//                     &T_.size_photon_chi_, 1, MPI_INT, smpi->world() );
//         MPI_Unpack( buffer, buf_size, &position,
//                     &T_.min_photon_chi_, 1, MPI_DOUBLE, smpi->world() );
//         MPI_Unpack( buffer, buf_size, &position,
//                     &T_.max_photon_chi_, 1, MPI_DOUBLE, smpi->world() );
// 
//         // Resize table before unpacking values
//         T_.table_.resize( T_.size_photon_chi_ );
// 
//         MPI_Unpack( buffer, buf_size, &position, &T_.table_[0],
//                     T_.size_photon_chi_, MPI_DOUBLE, smpi->world() );
// 
//     }
// 
//     delete[] buffer;
// 
//     T_.log10_min_photon_chi_ = log10( T_.min_photon_chi_ );
// 
//     // Computation of the delta
//     T_.photon_chi_delta_ = ( log10( T_.max_photon_chi_ )
//                       - T_.log10_min_photon_chi_ )/( T_.size_photon_chi_-1 );
// 
//     // Inverse delta
//     T_.photon_chi_inv_delta_ = 1.0/T_.photon_chi_delta_;
// }

// -----------------------------------------------------------------------------
//! Bcast of the external table xip_chipamin and xip
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
// void MultiphotonBreitWheelerTables::bcastTableXi( SmileiMPI *smpi )
// {
//     // Position for MPI pack and unack
//     int position = 0;
//     // buffer size for MPI pack and unpack
//     int buf_size = 0;
// 
//     // -------------------------------------------
//     // Bcast of all the parameters
//     // We pack everything in a buffer
//     // -------------------------------------------
// 
//     // Compute the buffer size
//     if( smpi->getRank() == 0 ) {
//         MPI_Pack_size( 2, MPI_INT, smpi->world(), &position );
//         buf_size = position;
//         MPI_Pack_size( 2, MPI_DOUBLE, smpi->world(), &position );
//         buf_size += position;
//         MPI_Pack_size( xi_.size_photon_chi_, MPI_DOUBLE, smpi->world(),
//                        &position );
//         buf_size += position;
//         MPI_Pack_size( xi_.size_photon_chi_*xi_.size_particle_chi_, MPI_DOUBLE,
//                        smpi->world(), &position );
//         buf_size += position;
//     }
// 
//     //MESSAGE( 2,"Buffer size for MPI exchange: " << buf_size );
// 
//     // Exchange buf_size with all ranks
//     MPI_Bcast( &buf_size, 1, MPI_INT, 0, smpi->world() );
// 
//     // Packet that will contain all parameters
//     char *buffer = new char[buf_size];
// 
//     // Proc 0 packs
//     if( smpi->getRank() == 0 ) {
//         position = 0;
//         MPI_Pack( &xi_.size_photon_chi_,
//                   1, MPI_INT, buffer, buf_size, &position, smpi->world() );
//         MPI_Pack( &xi_.size_particle_chi_,
//                   1, MPI_INT, buffer, buf_size, &position, smpi->world() );
//         MPI_Pack( &xi_.min_photon_chi_,
//                   1, MPI_DOUBLE, buffer, buf_size, &position, smpi->world() );
//         MPI_Pack( &xi_.max_photon_chi_,
//                   1, MPI_DOUBLE, buffer, buf_size, &position, smpi->world() );
// 
//         MPI_Pack( &xi_.min_particle_chi_[0], xi_.size_photon_chi_,
//                   MPI_DOUBLE, buffer, buf_size, &position, smpi->world() );
// 
//         MPI_Pack( &xi_.table_[0], xi_.size_particle_chi_*xi_.size_photon_chi_,
//                   MPI_DOUBLE, buffer, buf_size, &position, smpi->world() );
//     }
// 
//     // Bcast all parameters
//     MPI_Bcast( &buffer[0], buf_size, MPI_PACKED, 0, smpi->world() );
// 
//     // Other ranks unpack
//     if( smpi->getRank() != 0 ) {
//         position = 0;
//         MPI_Unpack( buffer, buf_size, &position,
//                     &xi_.size_photon_chi_, 1, MPI_INT, smpi->world() );
//         MPI_Unpack( buffer, buf_size, &position,
//                     &xi_.size_particle_chi_, 1, MPI_INT, smpi->world() );
//         MPI_Unpack( buffer, buf_size, &position,
//                     &xi_.min_photon_chi_, 1, MPI_DOUBLE, smpi->world() );
//         MPI_Unpack( buffer, buf_size, &position,
//                     &xi_.max_photon_chi_, 1, MPI_DOUBLE, smpi->world() );
// 
//         // Resize tables before unpacking values
//         xi_.min_particle_chi_.resize( xi_.size_photon_chi_ );
//         xi_.table_.resize( xi_.size_particle_chi_*xi_.size_photon_chi_ );
// 
//         MPI_Unpack( buffer, buf_size, &position, &xi_.min_particle_chi_[0],
//                     xi_.size_photon_chi_, MPI_DOUBLE, smpi->world() );
// 
//         MPI_Unpack( buffer, buf_size, &position, &xi_.table_[0],
//                     xi_.size_particle_chi_*xi_.size_photon_chi_, MPI_DOUBLE, smpi->world() );
//     }
// 
//     delete[] buffer;
// 
//     // Log10 of xi_.min_photon_chi_ for efficiency
//     xi_.log10_min_photon_chi_ = log10( xi_.min_photon_chi_ );
// 
//     // Computation of the delta
//     xi_.photon_chi_delta_ = ( log10( xi_.max_photon_chi_ )
//                         - xi_.log10_min_photon_chi_ )/( xi_.size_photon_chi_-1 );
// 
//     // Inverse of delta
//     xi_.photon_chi_inv_delta_ = 1./xi_.photon_chi_delta_;
// 
//     // Inverse particle_chi discetization (regularly used)
//     xi_.inv_size_particle_chi_minus_one_ = 1./( xi_.size_particle_chi_ - 1. );
// 
// }
