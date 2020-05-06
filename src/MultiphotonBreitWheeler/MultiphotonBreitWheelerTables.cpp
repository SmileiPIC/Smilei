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

// -----------------------------------------------------------------------------
// INITILIZATION AND DESTRUCTION
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Constructor for MutliphotonBreitWheelerTables
// -----------------------------------------------------------------------------
MultiphotonBreitWheelerTables::MultiphotonBreitWheelerTables()
{
    
    setDefault();
    
}

// -----------------------------------------------------------------------------
// Destructor for MutliphotonBreitWheelerTables
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
    if( params.hasMultiphotonBreitWheeler ) {
        TITLE( "Initializing mutliphoton Breit-Wheeler" )

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
    if( params.hasMultiphotonBreitWheeler ) {
        // Computation of the normalized Compton wavelength
        normalized_Compton_wavelength_ = params.red_planck_cst*params.reference_angular_frequency_SI
                                         / ( params.electron_mass*params.c_vacuum*params.c_vacuum );

        // Computation of the factor factor_dNBW_dt_
        factor_dNBW_dt_ = params.fine_struct_cst/( normalized_Compton_wavelength_ * M_PI*std::sqrt( 3.0 ) );
    }

    // Messages and checks
    if( params.hasMultiphotonBreitWheeler ) {
        if (table_path_.size() > 0) {
            
            MESSAGE( 1,"table path: " << table_path_ );
            MESSAGE( "" )
            
            readTables( params, smpi );
            
        } else {
            MESSAGE(1,"Default tables (stored in the code) are used:");
        }
        
        MESSAGE( 1,"--- Table `integration_dt_dchi`:" );
        MESSAGE( 2,"Reading of the external database" );
        MESSAGE( 2,"Dimension quantum parameter: "
                 << T_.size_photon_chi_ );
        MESSAGE( 2,"Minimum photon quantum parameter chi: "
                 << T_.min_photon_chi_ );
        MESSAGE( 2,"Maximum photon quantum parameter chi: "
                 << T_.max_photon_chi_ );
                 
        MESSAGE( "" )
        
        MESSAGE( 1,"--- Table `min_particle_chi_for_xi` and `xi`:" );
        MESSAGE( 2,"Reading of the external database" );
        MESSAGE( 2,"Dimension photon chi: " << xi_.size_photon_chi_ );
        MESSAGE( 2,"Dimension particle chi: " << xi_.size_particle_chi_ );
        MESSAGE( 2,"Minimum photon chi: " << xi_.min_photon_chi_ );
        MESSAGE( 2,"Maximum photon chi: " << xi_.max_photon_chi_ );
        
    }

}

// -----------------------------------------------------------------------------
// PHYSICAL COMPUTATION
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//! Computation of the production rate of pairs per photon
//! \param photon_chi photon quantum parameter
//! \param gamma photon normalized energy
// -----------------------------------------------------------------------------
double MultiphotonBreitWheelerTables::computeBreitWheelerPairProductionRate( double photon_chi, double gamma )
{
    // ________________________________________
    // Parameters

    // Log of the photon quantum parameter particle_chi
    double logchiph;
    double logchiphm;
    double logchiphp;
    // Index
    int ichiph;
    // final value
    double dNBWdt;

    // ________________________________________
    // Computation

    logchiph = log10( photon_chi );

    // Lower index for interpolation in the table integfochi
    ichiph = int( floor( ( logchiph-T_.log10_min_photon_chi_ )
                         *T_.photon_chi_inv_delta_ ) );

    // If photon_chi is below the lower bound of the table
    // An asymptotic approximation is used
    if( ichiph < 0 ) {
        ichiph = 0;
        // 0.2296 * sqrt(3) * pi [MG/correction by Antony]
        dNBWdt = 1.2493450020845291*exp( -8.0/( 3.0*photon_chi ) ) * photon_chi*photon_chi;
    }
    // If photon_chi is above the upper bound of the table
    // An asymptotic approximation is used
    else if( ichiph >= T_.size_photon_chi_-1 ) {
        ichiph = T_.size_photon_chi_-2;
        dNBWdt = 2.067731275227008*pow( photon_chi, 5.0/3.0 );
    } else {
        // Upper and lower values for linear interpolation
        logchiphm = ichiph*T_.photon_chi_delta_ + T_.log10_min_photon_chi_;
        logchiphp = logchiphm + T_.photon_chi_delta_;

        // Interpolation
        dNBWdt = ( T_.table_[ichiph+1]*fabs( logchiph-logchiphm ) +
                   T_.table_[ichiph]*fabs( logchiphp - logchiph ) )*T_.photon_chi_inv_delta_;
    }
    return factor_dNBW_dt_*dNBWdt/(photon_chi*gamma);
}

// -----------------------------------------------------------------------------
//! Computation of the electron and positron quantum parameters for
//! the multiphoton Breit-Wheeler pair creation
//
//! \param photon_chi photon quantum parameter
// -----------------------------------------------------------------------------
double *MultiphotonBreitWheelerTables::computePairQuantumParameter( double photon_chi, Random * rand )
{
    // Parameters
    double *chi = new double[2];
    double logchiph;
    double log10_chipam, log10_chipap;
    double d;
    double delta_chipa;
    double xip, xipp;
    int ichiph;
    int ichipa;
    int ixip;

    // -----------------------------------------------------------
    // Computation of the index associated to the given photon_chi
    // -----------------------------------------------------------

    logchiph = log10( photon_chi );

    // Lower boundary of the table
    if( photon_chi < xi_.min_photon_chi_ ) {
        ichiph = 0;
    }
    // Upper boundary of the table
    else if( photon_chi >= xi_.max_photon_chi_ ) {
        ichiph = xi_.size_photon_chi_-1;
    }
    // Inside the table
    else {
        // Use floor so that photon_chi corresponding to ichiph is <= given photon_chi
        ichiph = int( floor( ( logchiph-xi_.log10_min_photon_chi_ )*( xi_.photon_chi_inv_delta_ ) ) );
    }

    // ---------------------------------------
    // Search of the index ichiph for photon_chi
    // ---------------------------------------

    // First, we compute a random xip in [0,1[
    // xip = Rand::uniform();
    xip = rand->uniform();

    // The array uses the symmetric properties of the T fonction,
    // Cases xip > or <= 0.5 are treated seperatly

    // If xip > 0.5, the electron will bring more energy than the positron
    if( xip > 0.5 ) {
        xipp = 1.-xip;
    } else {
        xipp = xip;
    }

    // check boundaries
    // Lower bound
    if( xipp < xi_.table_[ichiph*xi_.size_particle_chi_] ) {
        ichipa = 0;
    }
    // Upper bound
    else if( xipp >= xi_.table_[( ichiph+1 )*xi_.size_particle_chi_-1] ) {
        ichipa = xi_.size_particle_chi_-2;
    } else {
        // Search for the corresponding index ichipa for xip
        ichipa = userFunctions::searchValuesInMonotonicArray(
                     &xi_.table_[ichiph*xi_.size_particle_chi_], xipp, xi_.size_particle_chi_ );
    }

    // Delta for the particle_chi dimension
    delta_chipa = ( log10( 0.5*photon_chi )-xi_.min_particle_chi_[ichiph] )
                  * xi_.inv_size_particle_chi_minus_one_;

    ixip = ichiph*xi_.size_particle_chi_ + ichipa;

    log10_chipam = ichipa*delta_chipa + xi_.min_particle_chi_[ichiph];
    log10_chipap = log10_chipam + delta_chipa;

    d = ( xipp - xi_.table_[ixip] ) / ( xi_.table_[ixip+1] - xi_.table_[ixip] );

    // If xip > 0.5, the electron will bring more energy than the positron
    if( xip > 0.5 ) {

        // Positron quantum parameter
        chi[1] = pow( 10, log10_chipam*( 1.0-d ) + log10_chipap*( d ) );

        // Electron quantum parameter
        chi[0] = photon_chi - chi[1];
    }
    // If xip <= 0.5, the positron will bring more energy than the electron
    else {
        // Electron quantum parameter
        chi[0] = pow( 10, log10_chipam*( 1.0-d ) + log10_chipap*( d ) );

        // Positron quantum parameter
        chi[1] = photon_chi - chi[0];
    }

    return chi;
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

    if( Tools::fileExists( table_path_ + "/multiphoton_Breit_Wheeler_tables.h5" ) ) {

        if( smpi->getRank()==0 ) {

            hid_t       fileId;
            hid_t       dataset_id;
            std::string buffer;

            buffer = table_path_ + "/multiphoton_Breit_Wheeler_tables.h5";

            fileId = H5Fopen( buffer.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );

            dataset_id = H5Dopen2( fileId, "integration_dt_dchi", H5P_DEFAULT );

            // If this dataset exists, we read it
            if( dataset_id > 0 ) {

                // First, we read attributes
                H5::getAttr( dataset_id, "size_photon_chi", T_.size_photon_chi_ );
                H5::getAttr( dataset_id, "min_photon_chi", T_.min_photon_chi_ );
                H5::getAttr( dataset_id, "max_photon_chi", T_.max_photon_chi_ );

                // Resize of the array integfochi_table before reading
                T_.table_.resize( T_.size_photon_chi_ );

                // then the dataset
                H5Dread( dataset_id,
                         H5T_NATIVE_DOUBLE, H5S_ALL,
                         H5S_ALL, H5P_DEFAULT,
                         &T_.table_[0] );

                H5Dclose( dataset_id );
                H5Fclose( fileId );
            }
            else {
                ERROR("Table T does not exist in the specified file `multiphoton_Breit_Wheeler_tables.h5`")
            }
        }
    }
    else
    {
        ERROR("The table T for the nonlinear Breit-Wheeler pair process could not be read from the provided path: `"
              << table_path_<<"`. Please check that the path is correct.")
    }

    // Bcast the table to all MPI ranks
    MultiphotonBreitWheelerTables::bcastTableT( smpi );

}

// -----------------------------------------------------------------------------
//! Read the external table xip_chipamin and xip
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
void MultiphotonBreitWheelerTables::readTableXi( SmileiMPI *smpi )
{
    if( Tools::fileExists( table_path_ + "/multiphoton_Breit_Wheeler_tables.h5" ) ) {
        if( smpi->getRank()==0 ) {

            hid_t       fileId;
            hid_t       dataset_id_chipamin;
            hid_t       dataset_id_xip;
            std::string buffer;

            buffer = table_path_ + "/multiphoton_Breit_Wheeler_tables.h5";

            fileId = H5Fopen( buffer.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );

            dataset_id_chipamin = H5Dopen2( fileId, "min_particle_chi_for_xi", H5P_DEFAULT );
            dataset_id_xip = H5Dopen2( fileId, "xi", H5P_DEFAULT );

            // If this dataset exists, we read it
            if( dataset_id_chipamin > 0 && dataset_id_xip > 0 ) {

                // First, we read attributes
                H5::getAttr( dataset_id_xip, "size_photon_chi", xi_.size_photon_chi_ );
                H5::getAttr( dataset_id_xip, "size_particle_chi", xi_.size_particle_chi_ );
                H5::getAttr( dataset_id_xip, "min_photon_chi", xi_.min_photon_chi_ );
                H5::getAttr( dataset_id_xip, "max_photon_chi", xi_.max_photon_chi_ );

                // Allocation of the array xip
                xi_.min_particle_chi_.resize( xi_.size_photon_chi_ );
                xi_.table_.resize( xi_.size_particle_chi_*xi_.size_photon_chi_ );

                // then the dataset for chipamin
                H5Dread( dataset_id_chipamin,
                         H5T_NATIVE_DOUBLE, H5S_ALL,
                         H5S_ALL, H5P_DEFAULT,
                         &xi_.min_particle_chi_[0] );

                // then the dataset for xip
                H5Dread( dataset_id_xip,
                         H5T_NATIVE_DOUBLE, H5S_ALL,
                         H5S_ALL, H5P_DEFAULT,
                         &xi_.table_[0] );

                H5Dclose( dataset_id_xip );
                H5Dclose( dataset_id_chipamin );
                H5Fclose( fileId );
            }
            else {
                ERROR("Table xip or xip_chipamin does not exist in the"
                <<" specified file `multiphoton_Breit_Wheeler_tables.h5`")
            }
        }

    }
    else
    {
        ERROR("The tables chipamin and xip for the nonlinear Breit-Wheeler pair"
              << " process could not be read from the provided path: `"
              << table_path_<<"`. Please check that the path is correct.")
    }

    // Bcast the table to all MPI ranks
    MultiphotonBreitWheelerTables::bcastTableXi( smpi );

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
    if( params.hasMultiphotonBreitWheeler) {
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
void MultiphotonBreitWheelerTables::bcastTableT( SmileiMPI *smpi )
{
    // Position for MPI pack and unack
    int position;
    // buffer size for MPI pack and unpack
    int buf_size;

    // -------------------------------------------------------
    // Bcast of all the parameters
    // We pack everything in a buffer
    // --------------------------------------------------------

    // buffer size
    if( smpi->getRank() == 0 ) {
        MPI_Pack_size( 1, MPI_INT, smpi->getGlobalComm(), &position );
        buf_size = position;
        MPI_Pack_size( 2, MPI_DOUBLE, smpi->getGlobalComm(), &position );
        buf_size += position;
        MPI_Pack_size( T_.size_photon_chi_, MPI_DOUBLE, smpi->getGlobalComm(),
                       &position );
        buf_size += position;
    }

    //MESSAGE( 2,"Buffer size: " << buf_size );

    // Exchange buf_size with all ranks
    MPI_Bcast( &buf_size, 1, MPI_INT, 0, smpi->getGlobalComm() );

    // Packet that will contain all parameters
    char *buffer = new char[buf_size];

    // Proc 0 packs
    if( smpi->getRank() == 0 ) {
        position = 0;
        MPI_Pack( &T_.size_photon_chi_,
                  1, MPI_INT, buffer, buf_size, &position, smpi->getGlobalComm() );
        MPI_Pack( &T_.min_photon_chi_,
                  1, MPI_DOUBLE, buffer, buf_size, &position, smpi->getGlobalComm() );
        MPI_Pack( &T_.max_photon_chi_,
                  1, MPI_DOUBLE, buffer, buf_size, &position, smpi->getGlobalComm() );

        MPI_Pack( &T_.table_[0], T_.size_photon_chi_,
                  MPI_DOUBLE, buffer, buf_size, &position, smpi->getGlobalComm() );

    }

    // Bcast all parameters
    MPI_Bcast( &buffer[0], buf_size, MPI_PACKED, 0, smpi->getGlobalComm() );

    // Other ranks unpack
    if( smpi->getRank() != 0 ) {
        position = 0;
        MPI_Unpack( buffer, buf_size, &position,
                    &T_.size_photon_chi_, 1, MPI_INT, smpi->getGlobalComm() );
        MPI_Unpack( buffer, buf_size, &position,
                    &T_.min_photon_chi_, 1, MPI_DOUBLE, smpi->getGlobalComm() );
        MPI_Unpack( buffer, buf_size, &position,
                    &T_.max_photon_chi_, 1, MPI_DOUBLE, smpi->getGlobalComm() );

        // Resize table before unpacking values
        T_.table_.resize( T_.size_photon_chi_ );

        MPI_Unpack( buffer, buf_size, &position, &T_.table_[0],
                    T_.size_photon_chi_, MPI_DOUBLE, smpi->getGlobalComm() );

    }

    delete[] buffer;

    T_.log10_min_photon_chi_ = log10( T_.min_photon_chi_ );

    // Computation of the delta
    T_.photon_chi_delta_ = ( log10( T_.max_photon_chi_ )
                      - T_.log10_min_photon_chi_ )/( T_.size_photon_chi_-1 );

    // Inverse delta
    T_.photon_chi_inv_delta_ = 1.0/T_.photon_chi_delta_;
}

// -----------------------------------------------------------------------------
//! Bcast of the external table xip_chipamin and xip
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
void MultiphotonBreitWheelerTables::bcastTableXi( SmileiMPI *smpi )
{
    // Position for MPI pack and unack
    int position = 0;
    // buffer size for MPI pack and unpack
    int buf_size = 0;

    // -------------------------------------------
    // Bcast of all the parameters
    // We pack everything in a buffer
    // -------------------------------------------

    // Compute the buffer size
    if( smpi->getRank() == 0 ) {
        MPI_Pack_size( 2, MPI_INT, smpi->getGlobalComm(), &position );
        buf_size = position;
        MPI_Pack_size( 2, MPI_DOUBLE, smpi->getGlobalComm(), &position );
        buf_size += position;
        MPI_Pack_size( xi_.size_photon_chi_, MPI_DOUBLE, smpi->getGlobalComm(),
                       &position );
        buf_size += position;
        MPI_Pack_size( xi_.size_photon_chi_*xi_.size_particle_chi_, MPI_DOUBLE,
                       smpi->getGlobalComm(), &position );
        buf_size += position;
    }

    //MESSAGE( 2,"Buffer size for MPI exchange: " << buf_size );

    // Exchange buf_size with all ranks
    MPI_Bcast( &buf_size, 1, MPI_INT, 0, smpi->getGlobalComm() );

    // Packet that will contain all parameters
    char *buffer = new char[buf_size];

    // Proc 0 packs
    if( smpi->getRank() == 0 ) {
        position = 0;
        MPI_Pack( &xi_.size_photon_chi_,
                  1, MPI_INT, buffer, buf_size, &position, smpi->getGlobalComm() );
        MPI_Pack( &xi_.size_particle_chi_,
                  1, MPI_INT, buffer, buf_size, &position, smpi->getGlobalComm() );
        MPI_Pack( &xi_.min_photon_chi_,
                  1, MPI_DOUBLE, buffer, buf_size, &position, smpi->getGlobalComm() );
        MPI_Pack( &xi_.max_photon_chi_,
                  1, MPI_DOUBLE, buffer, buf_size, &position, smpi->getGlobalComm() );

        MPI_Pack( &xi_.min_particle_chi_[0], xi_.size_photon_chi_,
                  MPI_DOUBLE, buffer, buf_size, &position, smpi->getGlobalComm() );

        MPI_Pack( &xi_.table_[0], xi_.size_particle_chi_*xi_.size_photon_chi_,
                  MPI_DOUBLE, buffer, buf_size, &position, smpi->getGlobalComm() );
    }

    // Bcast all parameters
    MPI_Bcast( &buffer[0], buf_size, MPI_PACKED, 0, smpi->getGlobalComm() );

    // Other ranks unpack
    if( smpi->getRank() != 0 ) {
        position = 0;
        MPI_Unpack( buffer, buf_size, &position,
                    &xi_.size_photon_chi_, 1, MPI_INT, smpi->getGlobalComm() );
        MPI_Unpack( buffer, buf_size, &position,
                    &xi_.size_particle_chi_, 1, MPI_INT, smpi->getGlobalComm() );
        MPI_Unpack( buffer, buf_size, &position,
                    &xi_.min_photon_chi_, 1, MPI_DOUBLE, smpi->getGlobalComm() );
        MPI_Unpack( buffer, buf_size, &position,
                    &xi_.max_photon_chi_, 1, MPI_DOUBLE, smpi->getGlobalComm() );

        // Resize tables before unpacking values
        xi_.min_particle_chi_.resize( xi_.size_photon_chi_ );
        xi_.table_.resize( xi_.size_particle_chi_*xi_.size_photon_chi_ );

        MPI_Unpack( buffer, buf_size, &position, &xi_.min_particle_chi_[0],
                    xi_.size_photon_chi_, MPI_DOUBLE, smpi->getGlobalComm() );

        MPI_Unpack( buffer, buf_size, &position, &xi_.table_[0],
                    xi_.size_particle_chi_*xi_.size_photon_chi_, MPI_DOUBLE, smpi->getGlobalComm() );
    }

    delete[] buffer;

    // Log10 of xi_.min_photon_chi_ for efficiency
    xi_.log10_min_photon_chi_ = log10( xi_.min_photon_chi_ );

    // Computation of the delta
    xi_.photon_chi_delta_ = ( log10( xi_.max_photon_chi_ )
                        - xi_.log10_min_photon_chi_ )/( xi_.size_photon_chi_-1 );

    // Inverse of delta
    xi_.photon_chi_inv_delta_ = 1./xi_.photon_chi_delta_;

    // Inverse particle_chi discetization (regularly used)
    xi_.inv_size_particle_chi_minus_one_ = 1./( xi_.size_particle_chi_ - 1. );

}
