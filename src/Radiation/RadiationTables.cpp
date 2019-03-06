// ----------------------------------------------------------------------------
//! \file RadiationTables.cpp
//
//! \brief This class contains the tables and the functions to generate them
//! for the Nonlinear Inverse Compton Scattering
//
//! \details The implementation is adapted from the thesis results of M. Lobet
//! See http://www.theses.fr/2015BORD0361
// ----------------------------------------------------------------------------

#include "RadiationTables.h"

#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <cmath>

#include "userFunctions.h"

// -----------------------------------------------------------------------------
// INITILIZATION AND DESTRUCTION
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Constructor for RadiationTables
// -----------------------------------------------------------------------------
RadiationTables::RadiationTables()
{
    h_table.resize( 0 );
    integfochi_table.resize( 0 );
    xip_chiphmin_table.resize( 0 );
    xip_table.resize( 0 );
    
    h_computed = false;
    integfochi_computed = false;
    xip_computed = false;
    
    minimum_chi_continuous_ = 1e-3;
    minimum_chi_discontinuous_ = 1e-2;
}

// -----------------------------------------------------------------------------
// Destructor for RadiationTables
// -----------------------------------------------------------------------------
RadiationTables::~RadiationTables()
{
}

// -----------------------------------------------------------------------------
//
//! Initialization of the parameters for the nonlinear inverse Compton
//! scattering
//
//! \param params Object Params for the parameters from the input script
// -----------------------------------------------------------------------------
void RadiationTables::initializeParameters( Params &params )
{

    if( params.hasMCRadiation ||
            params.hasLLRadiation ||
            params.hasNielRadiation ) {
        TITLE( "Initializing radiation reaction" )
        
        // Preliminary checks
        if( params.reference_angular_frequency_SI <= 0. )
            ERROR( "The parameter `reference_angular_frequency_SI` needs "
                   << "to be defined and positive to compute radiation losses" );
                   
    }
    
    if( params.hasMCRadiation ) {
        MESSAGE( "        The Monte-Carlo Compton radiation module"
                 << " is requested by some species.\n" );
    }
    if( params.hasNielRadiation ) {
        MESSAGE( "        The synchrotron-like stochastic radiation module"
                 << " of Niel et al. is requested by some species.\n" );
    }
    
    // If the namelist for Nonlinear Inverse Compton Scattering exists
    // We read the properties
    if( PyTools::nComponents( "RadiationReaction" ) != 0 ) {
    
        // If stochastic radiation loss is requested
        if( params.hasNielRadiation ) {
            // Extraction of the parameter from the input file
            PyTools::extract( "h_chipa_min", h_chipa_min, "RadiationReaction" );
            PyTools::extract( "h_chipa_max", h_chipa_max, "RadiationReaction" );
            PyTools::extract( "h_dim", h_dim, "RadiationReaction" );
            PyTools::extract( "h_computation_method", h_computation_method, "RadiationReaction" );
            
            h_log10_chipa_min = log10( h_chipa_min );
        }
        
        // If Monte-Carlo radiation loss is requested
        if( params.hasMCRadiation ) {
        
            // Extraction of the parameter from the input file
            PyTools::extract( "integfochi_chipa_min", integfochi_chipa_min, "RadiationReaction" );
            PyTools::extract( "integfochi_chipa_max", integfochi_chipa_max, "RadiationReaction" );
            PyTools::extract( "integfochi_dim", integfochi_dim, "RadiationReaction" );
            
            PyTools::extract( "xip_chipa_min", xip_chipa_min, "RadiationReaction" );
            PyTools::extract( "xip_chipa_max", xip_chipa_max, "RadiationReaction" );
            PyTools::extract( "xip_power", xip_power, "RadiationReaction" );
            PyTools::extract( "xip_threshold", xip_threshold, "RadiationReaction" );
            PyTools::extract( "xip_chipa_dim", xip_chipa_dim, "RadiationReaction" );
            PyTools::extract( "xip_chiph_dim", xip_chiph_dim, "RadiationReaction" );
            
            // Discontinuous minimum threshold
            PyTools::extract( "minimum_chi_discontinuous",
                              minimum_chi_discontinuous_, "RadiationReaction" );
                              
            // Additional regularly used parameters
            xip_log10_chipa_min = log10( xip_chipa_min );
            integfochi_log10_chipa_min = log10( integfochi_chipa_min );
            xip_inv_chiph_dim_minus_one = 1./( xip_chiph_dim - 1. );
        }
        
        // With any radiation model
        if( params.hasNielRadiation || params.hasMCRadiation ) {
            // Format of the tables
            PyTools::extract( "output_format", output_format_, "RadiationReaction" );
            
            // Path to the databases
            PyTools::extract( "table_path", table_path, "RadiationReaction" );
            
            // Radiation threshold on the quantum parameter particle_chi
            PyTools::extract( "minimum_chi_continuous",
                              minimum_chi_continuous_, "RadiationReaction" );
                              
        }
    }
    
    // Computation of some parameters
    if( params.hasMCRadiation ||
            params.hasLLRadiation ||
            params.hasNielRadiation ) {
            
        // Computation of the normalized Compton wavelength
        normalized_Compton_wavelength_ = params.red_planck_cst*params.reference_angular_frequency_SI
                                         / ( params.electron_mass*params.c_vacuum*params.c_vacuum );
                                         
        // Computation of the factor factor_dNphdt
        factor_dNphdt = sqrt( 3. )*params.fine_struct_cst/( 2.*M_PI*normalized_Compton_wavelength_ );
        
        // Computation of the factor for the classical radiated power
        factor_classical_radiated_power_ = 2.*params.fine_struct_cst/( 3.*normalized_Compton_wavelength_ );
        
        MESSAGE( "        Factor classical radiated power: " << factor_classical_radiated_power_ )
        
    }
    
    // Messages...
    // Computation of some parameters
    if( params.hasMCRadiation ||
            params.hasLLRadiation ||
            params.hasNielRadiation ) {
        MESSAGE( "        Minimum quantum parameter for continuous radiation: "
                 << std::setprecision( 5 ) << minimum_chi_continuous_ );
    }
    if( params.hasMCRadiation ) {
        MESSAGE( "        Minimum quantum parameter for discontinuous radiation: "
                 << std::setprecision( 5 ) << minimum_chi_discontinuous_ );
    }
    if( params.hasMCRadiation ||
            params.hasNielRadiation ) {
        MESSAGE( "        Table path: " << table_path );
    }
    if( params.hasNielRadiation ) {
        if( h_computation_method == "table" ||
                h_computation_method == "fit5"  ||
                h_computation_method == "fit10" ||
                h_computation_method == "ridgers" ) {
            MESSAGE( "        Niel h function computation method: " << h_computation_method )
        } else {
            ERROR( " The parameter `h_computation_method` must be `table`, `fit5`, `fit10` or `ridgers`." )
        }
    }
    
    MESSAGE( "" )
    
    // Some additional checks
    if( params.hasMCRadiation ) {
        if( integfochi_chipa_min >= integfochi_chipa_max ) {
            ERROR( "integfochi_chipa_min (" << integfochi_chipa_min
                   << ") >= integfochi_chipa_max (" << integfochi_chipa_max << ")" )
        }
        if( xip_chipa_min >= xip_chipa_max ) {
            ERROR( "xip_chipa_min (" << xip_chipa_min
                   << ") >= xip_chipa_max (" << xip_chipa_max << ")" )
        }
    }
    if( params.hasNielRadiation ) {
        if( h_chipa_min >= h_chipa_max ) {
            ERROR( "h_chipa_min (" << h_chipa_min
                   << ") >= h_chipa_max (" << h_chipa_max << ")" )
        }
    }
}

// -----------------------------------------------------------------------------
// TABLE COMPUTATION
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//! Computation of the table h that is a discetization of the h function
//! in the stochastic model of Niel.
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
void RadiationTables::compute_h_table( SmileiMPI *smpi )
{

    // Parameters
    int rank; // Rank number
    // timers
    double t0, t1;
    // tabel_exists
    bool table_exists;
    
    // Get the MPI rank
    rank = smpi->getRank();
    
    MESSAGE( "        --- h table:" );
    
    // Initial timer
    t0 = MPI_Wtime();
    
    // If external tables are available, we read them
    table_exists = RadiationTables::read_h_table( smpi );
    
    // Else we compute them
    if( !table_exists ) {
    
        // Temporary particle chi value
        double particle_chi;
        // For percentages
        double pct = 0.;
        double dpct = 10.;
        // table load repartition
        int *imin_table;
        int *length_table;
        // Local array
        double *buffer;
        int nb_ranks; // Number of ranks
        //int err;  // error MPI
        
        // checks
        if( minimum_chi_continuous_ < h_chipa_min ) {
            ERROR( "Parameter `minimum_chi_continuous_` is below `h_chipa_min`,"
                   << "the lower bound of the h table should be equal or below"
                   << "the radiation threshold on chi." )
        }
        
        // Get the number of ranks
        nb_ranks = smpi->getSize();
        
        // Allocation of the array h_table
        h_table.resize( h_dim );
        
        // Allocation of the table for load repartition
        imin_table = new int[nb_ranks];
        length_table = new int[nb_ranks];
        
        // Computation of the delta
        h_chipa_delta = ( log10( h_chipa_max )
                          - h_log10_chipa_min )/( h_dim-1 );
                          
        // Inverse delta
        h_chipa_inv_delta = 1./h_chipa_delta;
        
        // Load repartition
        userFunctions::distribute_load_1d_table( nb_ranks,
                h_dim,
                imin_table,
                length_table );
                
        // Allocation of the local buffer
        buffer = new double [length_table[rank]];
        
        MESSAGE( std::setprecision( 5 ) <<std::setprecision( 5 ) <<"            Dimension quantum parameter: "
                 << h_dim );
        MESSAGE( std::setprecision( 5 ) <<"            Minimum particle quantum parameter chi: "
                 << h_chipa_min );
        MESSAGE( std::setprecision( 5 ) <<"            Maximum particle quantum parameter chi: "
                 << h_chipa_max );
        MESSAGE( "            MPI repartition:" );
        // Print repartition
        if( rank==0 ) {
            for( int i =0 ; i < nb_ranks ; i++ ) {
                MESSAGE( "            Rank: " << i
                         << " imin: " << imin_table[i]
                         << " length: " << length_table[i] );
            }
        }
        
        MESSAGE( "            Computation:" );
        dpct = std::max( dpct, 100./length_table[rank] );
        // Loop over the table values
        for( int i = 0 ; i < length_table[rank] ; i++ ) {
            particle_chi = pow( 10., ( imin_table[rank] + i )*h_chipa_delta
                                + h_log10_chipa_min );
                                
            // 100 iterations is in theory sufficient to get the convergence
            // at 1e-15
            buffer[i] = RadiationTables::computeHNiel( particle_chi, 400, 1e-15 );
            
            if( 100.*i >= length_table[rank]*pct ) {
                pct += dpct;
                MESSAGE( "            " << i + 1<< "/" << length_table[rank]
                         << " - " << ( int )( std::round( pct ) )
                         << "%" );
            }
        }
        
        // Communication of the data
        MPI_Allgatherv( &buffer[0], length_table[rank], MPI_DOUBLE,
                        &h_table[0], &length_table[0], &imin_table[0],
                        MPI_DOUBLE, smpi->getGlobalComm() );
                        
        // flag computed at true
        h_computed = true;
        
        // Free memory
        delete buffer;
        delete imin_table;
        delete length_table;
        
    }
    
    // Final timer
    t1 = MPI_Wtime();
    MESSAGE( "        done in " << ( t1 - t0 ) << "s" );
}

// -----------------------------------------------------------------------------
//! Computation of the table values of integfochi
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
void RadiationTables::compute_integfochi_table( SmileiMPI *smpi )
{

    // Parameters
    int rank; // Rank number
    // timers
    double t0, t1;
    // tabel_exists
    bool table_exists;
    
    // Get the MPI rank
    rank = smpi->getRank();
    
    MESSAGE( "        --- Integration F/particle_chi table:" );
    
    // Initial timer
    t0 = MPI_Wtime();
    
    // If external tables are available, we read them
    table_exists = RadiationTables::read_integfochi_table( smpi );
    
    // Else we compute them
    if( !table_exists ) {
    
        // Temporary particle chi value
        double particle_chi;
        // For percentages
        double pct = 0.;
        double dpct = 10.;
        // table load repartition
        int *imin_table;
        int *length_table;
        // Local array
        double *buffer;
        int nb_ranks; // Number of ranks
        //int err;  // error MPI
        
        // Get the number of ranks
        nb_ranks = smpi->getSize();
        
        // Allocation of the array integfochi_table
        integfochi_table.resize( integfochi_dim );
        
        // Allocation of the table for load repartition
        imin_table = new int[nb_ranks];
        length_table = new int[nb_ranks];
        
        // Computation of the delta
        integfochi_chipa_delta = ( log10( integfochi_chipa_max )
                                   - integfochi_log10_chipa_min )/( integfochi_dim-1 );
                                   
        // Inverse delta
        integfochi_chipa_inv_delta = 1./integfochi_chipa_delta;
        
        // Load repartition
        userFunctions::distribute_load_1d_table( nb_ranks,
                integfochi_dim,
                imin_table,
                length_table );
                
        // Allocation of the local buffer
        buffer = new double [length_table[rank]];
        
        
        MESSAGE( "            MPI repartition:" );
        // Print repartition
        if( rank==0 ) {
            for( int i =0 ; i < nb_ranks ; i++ ) {
                MESSAGE( "            Rank: " << i
                         << " imin: " << imin_table[i]
                         << " length: " << length_table[i] );
            }
        }
        
        MESSAGE( "            Computation:" );
        dpct = std::max( dpct, 100./length_table[rank] );
        // Loop over the table values
        for( int i = 0 ; i < length_table[rank] ; i++ ) {
            particle_chi = pow( 10., ( imin_table[rank] + i )*integfochi_chipa_delta
                                + integfochi_log10_chipa_min );
                                
            buffer[i] = RadiationTables::integrateSynchrotronEmissivity( particle_chi,
                        0.98e-40*particle_chi, 0.98*particle_chi, 400, 1e-15 );
                        
            //std::cout << rank << " " << buffer[i] << std::endl;
            
            if( 100.*i >= length_table[rank]*pct ) {
                pct += dpct;
                MESSAGE( "            " << i + 1<< "/" << length_table[rank] << " - " << ( int )( std::round( pct ) ) << "%" );
            }
        }
        
        // Communication of the data
        MPI_Allgatherv( &buffer[0], length_table[rank], MPI_DOUBLE,
                        &integfochi_table[0], &length_table[0], &imin_table[0],
                        MPI_DOUBLE, smpi->getGlobalComm() );
                        
        // flag computed at true
        integfochi_computed = true;
        
        // Free memory
        delete buffer;
        delete imin_table;
        delete length_table;
        
    }
    
    // Final timer
    t1 = MPI_Wtime();
    MESSAGE( "        done in " << ( t1 - t0 ) << "s" );
}

// -----------------------------------------------------------------------------
//! Computation of the minimum photon quantum parameter chiphmin
//! for the xip array and computation of the xip array.
//
//! \details Under the minimum photon_chi value, photon energy is
//! considered negligible.
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
void RadiationTables::compute_xip_table( SmileiMPI *smpi )
{

    // Parameters
    int rank; // Rank number
    // timers
    double t0, t1;
    // Flag table exists
    bool table_exists;
    
    // Get the MPI rank
    rank = smpi->getRank();
    
    t0 = MPI_Wtime();
    
    MESSAGE( "        --- Table chiphmin and xip:" );
    
    // If external tables are available, we read them
    table_exists = RadiationTables::read_xip_table( smpi );
    
    // Else we compute them
    if( !table_exists ) {
        // Parameters:
        double particle_chi; // Temporary particle chi value
        double photon_chi; // Temporary photon chi value
        double chiph_delta; // Temporary delta for photon_chi
        double logchiphmin; // Temporary log10 of photon chi value
        double xip;   // Temporary xip
        double numerator;
        double denominator;
        // For percentages
        double pct = 0.;
        double dpct = 10.;
        // table load repartition
        int *imin_table;
        int *length_table;
        // Local array
        double *buffer;
        //int err;  // error MPI
        int nb_ranks; // Number of ranks
        // Iterator
        int  k;
        
        // Get the number of ranks
        nb_ranks = smpi->getSize();
        
        // Allocation of the array xip_chiphmin_table
        xip_chiphmin_table.resize( xip_chipa_dim );
        
        // Allocation of the array xip_table
        xip_table.resize( xip_chipa_dim*xip_chiph_dim );
        
        // Allocation of the table for load repartition
        imin_table = new int[nb_ranks];
        length_table = new int[nb_ranks];
        
        // Computation of the delta
        xip_chipa_delta = ( log10( xip_chipa_max )
                            - xip_log10_chipa_min )/( xip_chipa_dim-1 );
                            
        // Inverse of delta
        xip_chipa_inv_delta = 1./xip_chipa_delta;
        
        // Load repartition
        userFunctions::distribute_load_1d_table( nb_ranks,
                xip_chipa_dim,
                imin_table,
                length_table );
                
        // Allocation of the local buffer
        buffer = new double [length_table[rank]];
        
        MESSAGE( "            MPI repartition:" );
        // Print repartition
        if( rank==0 ) {
            for( int i =0 ; i < nb_ranks ; i++ ) {
                MESSAGE( "            Rank: " << i
                         << " imin: "   << imin_table[i]
                         << " length: " << length_table[i] );
            }
        }
        
        // 1. - Computation of xip_chiphmin_table
        MESSAGE( "            Computation of log10(chiphmin):" );
        dpct = std::max( 10., 100./length_table[rank] );
        
        // Loop for chiphmin
        for( int ichipa = 0 ; ichipa < length_table[rank] ; ichipa++ ) {
        
            xip = 1;
            logchiphmin = ( imin_table[rank] + ichipa )*xip_chipa_delta
                          + xip_log10_chipa_min;
            particle_chi = pow( 10., logchiphmin );
            
            // Denominator of xip
            denominator = RadiationTables::integrateSynchrotronEmissivity( particle_chi,
                          0.99e-40*particle_chi, 0.99*particle_chi, 200, 1e-13 );
                          
            k = 0;
            while( k < xip_power ) {
                logchiphmin -= pow( 0.1, k );
                photon_chi = pow( 10., logchiphmin );
                numerator = RadiationTables::integrateSynchrotronEmissivity( particle_chi,
                            0.99e-40*photon_chi, 0.99*photon_chi, 200, 1e-13 );
                            
                if( numerator == 0 ) {
                    xip = 0;
                } else {
                    xip = numerator/denominator;
                }
                
                if( xip < xip_threshold ) {
                    logchiphmin += pow( 0.1, k );
                    k += 1;
                }
            }
            buffer[ichipa] = logchiphmin;
            
            // display percentage
            if( 100.*ichipa >= length_table[rank]*pct ) {
                pct += dpct;
                MESSAGE( "            " << ichipa + 1 << "/" << length_table[rank]
                         << " - " << ( int )( std::round( pct ) ) << "%" );
            }
        }
        
        // Communication of the xip_chiphmin table
        MPI_Allgatherv( &buffer[0], length_table[rank], MPI_DOUBLE,
                        &xip_chiphmin_table[0], &length_table[0], &imin_table[0],
                        MPI_DOUBLE, smpi->getGlobalComm() );
                        
        // 2. - Computation of the xip table
        MESSAGE( "            Computation of xip:" );
        
        // Allocation of the local buffer
        buffer = new double [length_table[rank]*xip_chiph_dim];
        
        // Loop for xip in the particle_chi dimension
        pct = 0;
        dpct = std::max( 10., 100./( length_table[rank]*xip_chiph_dim ) );
        for( int ichipa = 0 ; ichipa < length_table[rank] ; ichipa++ ) {
        
            particle_chi = ( imin_table[rank] + ichipa )*xip_chipa_delta
                           + xip_log10_chipa_min;
                           
            chiph_delta = ( particle_chi - xip_chiphmin_table[imin_table[rank] + ichipa] )
                          / ( xip_chiph_dim - 1 );
                          
            particle_chi = pow( 10., particle_chi );
            
            // Denominator of xip
            denominator = RadiationTables::integrateSynchrotronEmissivity( particle_chi,
                          1e-40*particle_chi, particle_chi, 300, 1e-15 );
                          
            // Loop in the photon_chi dimension
            for( int ichiph = 0 ; ichiph < xip_chiph_dim ; ichiph ++ ) {
                // Local photon_chi value
                photon_chi = pow( 10., ichiph*chiph_delta +
                                  xip_chiphmin_table[imin_table[rank] + ichipa] );
                                  
                // Numerator of xip
                numerator = RadiationTables::integrateSynchrotronEmissivity( particle_chi,
                            0.99e-40*photon_chi, 0.99*photon_chi, 300, 1e-15 );
                            
                // Update local buffer value
                buffer[ichipa*xip_chiph_dim + ichiph] = std::min( 1., numerator / denominator );
                
                // If buffer == 1, end of the loop with 1
                if( buffer[ichipa*xip_chiph_dim + ichiph] == 1 ) {
                    for( int i = ichiph+1 ; i < xip_chiph_dim ; i ++ ) {
                        buffer[ichipa*xip_chiph_dim + i] = 1.;
                    }
                    ichiph = xip_chiph_dim;
                }
                
                // Artificial monotony
                if( ( ichiph > 0 ) &&
                        ( buffer[ichipa*xip_chiph_dim + ichiph] <
                          buffer[ichipa*xip_chiph_dim + ichiph -1] ) ) {
                    buffer[ichipa*xip_chiph_dim + ichiph] =
                        buffer[ichipa*xip_chiph_dim + ichiph -1];
                }
                
                // display percentage
                if( 100.*( ichipa*xip_chiph_dim+ichiph )
                        >= length_table[rank]*xip_chiph_dim*pct ) {
                    pct += dpct;
                    MESSAGE( "            " << ichipa*xip_chiph_dim+ichiph + 1
                             << "/"
                             << length_table[rank]*xip_chiph_dim
                             << " - " << ( int )( std::round( pct ) )
                             << "%" );
                }
                
            }
        }
        
        // Update length_table and imin_table
        for( int i = 0 ; i < nb_ranks ; i++ ) {
            length_table[i] *= xip_chiph_dim;
            imin_table[i] *= xip_chiph_dim;
        }
        
        // Communication of the xip table
        MPI_Allgatherv( &buffer[0], length_table[rank], MPI_DOUBLE,
                        &xip_table[0], &length_table[0], &imin_table[0],
                        MPI_DOUBLE, smpi->getGlobalComm() );
                        
        // flag computed at true
        xip_computed = true;
        
        // clean temporary arrays
        delete buffer;
        delete length_table;
        delete imin_table;
        
    }
    
    t1 = MPI_Wtime();
    MESSAGE( "        done in " << ( t1 - t0 ) << "s" );
    
}

// -----------------------------------------------------------------------------
//! Output the computed tables so that thay can be read at the next run.
//
//! \param params list of simulation parameters
//! \param smpi MPI parameters
// -----------------------------------------------------------------------------
void RadiationTables::compute_tables( Params &params, SmileiMPI *smpi )
{
    // These tables are loaded only if if one species has Monte-Carlo Compton radiation
    // And if the h values are not computed from a numerical fit
    if( params.hasNielRadiation && this->h_computation_method == "table" ) {
        RadiationTables::compute_h_table( smpi );
    }
    if( params.hasMCRadiation ) {
        RadiationTables::compute_integfochi_table( smpi );
        RadiationTables::compute_xip_table( smpi );
    }
}

// -----------------------------------------------------------------------------
// TABLE OUTPUTS
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//! Ouput in a file of the table values of h for the Niel radiation model
//
// -----------------------------------------------------------------------------
void RadiationTables::output_h_table()
{

    if( output_format_ == "ascii" ) {
        std::ofstream file;
        file.open( table_path + "/tab_h.dat" );
        
        if( file.is_open() ) {
        
            file.precision( 12 );
            
            file << "Stochastic synchrotron-like radiation model of Niel \n";
            
            file << "Dimension particle_chi - particle_chi min - particle_chi max \n";
            
            file << h_dim;
            file << " "
                 << h_chipa_min << " "
                 << h_chipa_max << "\n";;
                 
            // Loop over the table values
            for( int i = 0 ; i < h_dim ; i++ ) {
                file <<  h_table[i] << "\n";
            }
            
            file.close();
        }
    } else if( output_format_ == "binary" ) {
        std::ofstream file;
        file.open( table_path + "/tab_h.bin", std::ios::binary );
        
        if( file.is_open() ) {
        
            file.write( ( char * )&h_dim, sizeof( h_dim ) );
            file.write( ( char * )&h_chipa_min, sizeof( double ) );
            file.write( ( char * )&h_chipa_max, sizeof( double ) );
            
            // Loop over the table values
            for( int i = 0 ; i < h_dim ; i++ ) {
                file.write( ( char * )&h_table[i], sizeof( double ) );
            }
            
            file.close();
        }
    }
    // HDF5
    // The table is written as a dataset
    else if( output_format_ == "hdf5" ) {
    
        hid_t       fileId;
        hid_t       datasetId;
        hid_t       dataspaceId;
        hsize_t     dims;
        std::string buffer;
        
        buffer = table_path + "/radiation_tables.h5";
        
        // We first check whether the file already exists
        // If yes, we simply open the file
        if( Tools::file_exists( buffer ) ) {
            fileId = H5Fopen( buffer.c_str(),
                              H5F_ACC_RDWR,
                              H5P_DEFAULT );
                              
        }
        // Else, we create the file
        else {
            fileId  = H5Fcreate( buffer.c_str(),
                                 H5F_ACC_TRUNC,
                                 H5P_DEFAULT,
                                 H5P_DEFAULT );
        }
        
        // Create the data space for the dataset.
        dims = h_dim;
        dataspaceId = H5Screate_simple( 1, &dims, NULL );
        
        // Creation of the dataset
        datasetId = H5Dcreate( fileId,
                               "h",
                               H5T_NATIVE_DOUBLE,
                               dataspaceId,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                               
        // Fill the dataset
        H5Dwrite( datasetId, H5T_NATIVE_DOUBLE,
                  H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  &h_table[0] );
                  
        // Create attributes
        H5::attr( datasetId, "chipa_min", h_chipa_min );
        H5::attr( datasetId, "chipa_max", h_chipa_max );
        H5::attr( datasetId, "chipa_dim", h_dim );
        
        // Close everything
        H5Dclose( datasetId );
        H5Sclose( dataspaceId );
        H5Fclose( fileId );
        
    } else {
        MESSAGE( "The table output format " << output_format_
                 << " is not recognized" );
    }
}

// -----------------------------------------------------------------------------
//! Ouput in a file of the table values of integfochi
//
// -----------------------------------------------------------------------------
void RadiationTables::output_integfochi_table()
{

    if( output_format_ == "ascii" ) {
        std::ofstream file;
        file.open( table_path + "/tab_integfochi.dat" );
        
        if( file.is_open() ) {
        
            file.precision( 12 );
            
            file << "Table Integration F(CHI)/CHI for Nonlinear Compton Scattering \n";
            
            file << "Dimension particle_chi - particle_chi min - particle_chi max \n";
            
            file << integfochi_dim ;
            file << " "
                 << integfochi_chipa_min << " "
                 << integfochi_chipa_max << "\n";;
                 
            // Loop over the table values
            for( int i = 0 ; i < integfochi_dim ; i++ ) {
                file <<  integfochi_table[i] << "\n";
            }
            
            file.close();
        }
    } else if( output_format_ == "binary" ) {
        std::ofstream file;
        file.open( table_path + "/tab_integfochi.bin", std::ios::binary );
        
        if( file.is_open() ) {
        
            double temp0, temp1;
            
            temp0 = integfochi_chipa_min;
            temp1 = integfochi_chipa_max;
            
            file.write( ( char * )&integfochi_dim, sizeof( integfochi_dim ) );
            file.write( ( char * )&temp0, sizeof( double ) );
            file.write( ( char * )&temp1, sizeof( double ) );
            
            // Loop over the table values
            for( int i = 0 ; i < integfochi_dim ; i++ ) {
                file.write( ( char * )&integfochi_table[i], sizeof( double ) );
            }
            
            file.close();
        }
    }
    // HDF5
    // The table is written as a dataset
    else if( output_format_ == "hdf5" ) {
    
        hid_t       fileId;
        hid_t       datasetId;
        hid_t       dataspaceId;
        hsize_t     dims;
        std::string buffer;
        
        buffer = table_path + "/radiation_tables.h5";
        
        // We first check whether the file already exists
        // If yes, we simply open the file
        if( Tools::file_exists( buffer ) ) {
            fileId = H5Fopen( buffer.c_str(),
                              H5F_ACC_RDWR,
                              H5P_DEFAULT );
                              
        }
        // Else, we create the file
        else {
            fileId  = H5Fcreate( buffer.c_str(),
                                 H5F_ACC_TRUNC,
                                 H5P_DEFAULT,
                                 H5P_DEFAULT );
        }
        
        // Create the data space for the dataset.
        dims = integfochi_dim;
        dataspaceId = H5Screate_simple( 1, &dims, NULL );
        
        // Creation of the dataset
        datasetId = H5Dcreate( fileId,
                               "integfochi",
                               H5T_NATIVE_DOUBLE,
                               dataspaceId,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                               
        // Fill the dataset
        H5Dwrite( datasetId, H5T_NATIVE_DOUBLE,
                  H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  &integfochi_table[0] );
                  
        // Create attributes
        H5::attr( datasetId, "chipa_min", integfochi_chipa_min );
        H5::attr( datasetId, "chipa_max", integfochi_chipa_max );
        H5::attr( datasetId, "chipa_dim", integfochi_dim );
        
        // Close everything
        H5Dclose( datasetId );
        H5Sclose( dataspaceId );
        H5Fclose( fileId );
        
    } else {
        MESSAGE( "The table output format " << output_format_
                 << " is not recognized" );
    }
}

// -----------------------------------------------------------------------------
//! File output of xip_chiphmin_table and xip_table
//
// -----------------------------------------------------------------------------
void RadiationTables::output_xip_table()
{

    if( output_format_ == "ascii" ) {
        std::ofstream file;
        file.open( table_path + "/tab_xip.dat" );
        
        if( file.is_open() ) {
        
            file.precision( 12 );
            
            file << "Table xip_chiphmin and xip for Nonlinear Compton Scattering \n";
            
            file << "Dimension particle_chi - Dimension photon_chi - particle_chi min - particle_chi max \n";
            
            file << xip_chipa_dim << " "
                 << xip_chiph_dim << " "
                 << xip_chipa_min << " "
                 << xip_chipa_max << "\n";;
                 
            // Loop over the xip values
            for( int ichipa = 0 ; ichipa < xip_chipa_dim ; ichipa++ ) {
                for( int ichiph = 0 ; ichiph < xip_chiph_dim ; ichiph++ ) {
                    file <<  xip_table[ichipa*xip_chiph_dim+ichiph] << " ";
                }
                file << "\n";
            }
            
            file.close();
        }
    } else if( output_format_ == "binary" ) {
        std::ofstream file;
        file.open( table_path + "/tab_xip.bin", std::ios::binary );
        
        if( file.is_open() ) {
        
            double temp0, temp1;
            
            temp0 = xip_chipa_min;
            temp1 = xip_chipa_max;
            
            file.write( ( char * )&xip_chipa_dim, sizeof( int ) );
            file.write( ( char * )&xip_chiph_dim, sizeof( int ) );
            file.write( ( char * )&temp0, sizeof( double ) );
            file.write( ( char * )&temp1, sizeof( double ) );
            
            // Write all table values of xip_chimin_table
            file.write( ( char * )&xip_chiphmin_table[0], sizeof( double )*xip_chipa_dim );
            
            // Write all table values of xip_table
            file.write( ( char * )&xip_table[0], sizeof( double )*xip_chipa_dim*xip_chiph_dim );
            
            file.close();
        }
    }
    // HDF5
    // The table is written as a dataset
    else if( output_format_ == "hdf5" ) {
    
        hid_t       fileId;
        hid_t       datasetId;
        hid_t       dataspaceId;
        hsize_t     dims[2];
        std::string buffer;
        
        buffer = table_path + "/radiation_tables.h5";
        
        // We first check whether the file already exists
        // If yes, we simply open the file
        if( Tools::file_exists( buffer ) ) {
            fileId = H5Fopen( buffer.c_str(),
                              H5F_ACC_RDWR,
                              H5P_DEFAULT );
                              
        }
        // Else, we create the file
        else {
            fileId  = H5Fcreate( buffer.c_str(),
                                 H5F_ACC_TRUNC,
                                 H5P_DEFAULT,
                                 H5P_DEFAULT );
        }
        
        // Creation of the datasat chiphmin
        dims[0] = xip_chipa_dim;
        dataspaceId = H5Screate_simple( 1, dims, NULL );
        datasetId = H5Dcreate( fileId,
                               "xip_chiphmin",
                               H5T_NATIVE_DOUBLE,
                               dataspaceId,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                               
        // Fill the dataset chiphmin
        H5Dwrite( datasetId, H5T_NATIVE_DOUBLE,
                  H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  &xip_chiphmin_table[0] );
                  
        // Attribute creation
        H5::attr( datasetId, "chipa_min", xip_chipa_min );
        H5::attr( datasetId, "chipa_max", xip_chipa_max );
        H5::attr( datasetId, "chipa_dim", xip_chipa_dim );
        H5::attr( datasetId, "power", xip_power );
        H5::attr( datasetId, "threshold", xip_threshold );
        
        // Creation of the datasat chiphmin xip
        dims[0] = xip_chipa_dim;
        dims[1] = xip_chiph_dim;
        dataspaceId = H5Screate_simple( 2, dims, NULL );
        datasetId = H5Dcreate( fileId,
                               "xip",
                               H5T_NATIVE_DOUBLE,
                               dataspaceId,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                               
        // Fill the dataset
        H5Dwrite( datasetId, H5T_NATIVE_DOUBLE,
                  H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  &xip_table[0] );
                  
        // Attribute creation
        H5::attr( datasetId, "chipa_min", xip_chipa_min );
        H5::attr( datasetId, "chipa_max", xip_chipa_max );
        H5::attr( datasetId, "chipa_dim", xip_chipa_dim );
        H5::attr( datasetId, "chiph_dim", xip_chiph_dim );
        H5::attr( datasetId, "power", xip_power );
        H5::attr( datasetId, "threshold", xip_threshold );
        
        // Close everything
        H5Dclose( datasetId );
        H5Sclose( dataspaceId );
        H5Fclose( fileId );
    } else {
        MESSAGE( "The output format " << output_format_ << " is not recognized" );
    }
}

// -----------------------------------------------------------------------------
//! Output the computed tables so that thay can be read at the next run.
//! Table output by the master MPI rank
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
void RadiationTables::output_tables( SmileiMPI *smpi )
{
    // Sequential output
    if( smpi->isMaster() ) {
        // If tables have been computed, they are output on the disk
        // to be used for the next run
        if( h_computed ) {
            RadiationTables::output_h_table();
        }
        if( integfochi_computed ) {
            RadiationTables::output_integfochi_table();
        }
        if( xip_computed ) {
            RadiationTables::output_xip_table();
        }
    }
}

// -----------------------------------------------------------------------------
// PHYSICAL COMPUTATION
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//! Computation of the photon quantum parameter photon_chi for emission
//! ramdomly and using the tables xip and chiphmin
//
//! \param particle_chi particle quantum parameter
// -----------------------------------------------------------------------------
double RadiationTables::computeRandomPhotonChi( double particle_chi )
{
    // Log10 of particle_chi
    double logchipa;
    double photon_chi;
    double chiph_xip_delta;
    // Random xip
    double xip;
    int ichipa;
    int ichiph;
    // For the interpolation
    double log10_chiphm;
    double log10_chiphp;
    double d;
    int ixip;
    
    logchipa = log10( particle_chi );
    
    // ---------------------------------------
    // index of particle_chi in xip_table
    // ---------------------------------------
    // Use floor so that particle_chi corresponding to ichipa is <= given particle_chi
    ichipa = int( floor( ( logchipa-xip_log10_chipa_min )*( xip_chipa_inv_delta ) ) );
    
    // Checking that ichipa is in the range of the tables
    // Else we use the values at the boundaries
    if( ichipa < 0 ) {
        ichipa = 0;
    } else if( ichipa > xip_chipa_dim-1 ) {
        ichipa = xip_chipa_dim-1;
    }
    
    // ---------------------------------------
    // Search of the index ichiph for photon_chi
    // ---------------------------------------
    
    // First, we compute a random xip in [0,1[
    xip = Rand::uniform();
    
    // If the randomly computed xip if below the first one of the row,
    // we take the first one which corresponds to the minimal photon photon_chi
    if( xip <= xip_table[ichipa*xip_chiph_dim] ) {
        ichiph = 0;
        xip = xip_table[ichipa*xip_chiph_dim];
    }
    // Above the last xip of the row, the last one corresponds
    // to the maximal photon photon_chi
    else if( xip > xip_table[( ichipa+1 )*xip_chiph_dim-2] ) {
        ichiph = xip_chiph_dim-2;
        xip = xip_table[( ichipa+1 )*xip_chiph_dim-1];
        // If nearest point: ichiph = xip_chiph_dim-1
    } else {
        // Search for the corresponding index ichiph for xip
        ichiph = userFunctions::search_elem_in_array(
                     &xip_table[ichipa*xip_chiph_dim], xip, xip_chiph_dim );
    }
    
    // Corresponding particle_chi for ichipa
    logchipa = ichipa*xip_chipa_delta+xip_log10_chipa_min;
    
    // Delta for the corresponding particle_chi
    chiph_xip_delta = ( logchipa - xip_chiphmin_table[ichipa] )
                      *xip_inv_chiph_dim_minus_one;
                      
    // --------------------------------------------------------------------
    // Compute photon_chi
    // This method is slow but more accurate than taking the nearest point
    // --------------------------------------------------------------------
    
    ixip = ichipa*xip_chiph_dim + ichiph;
    
    // Computation of the final photon_chi by interpolation
    if( xip_table[ixip+1] - xip_table[ixip] > 1e-15 ) {
        log10_chiphm = ichiph*chiph_xip_delta
                       + xip_chiphmin_table[ichipa];
        log10_chiphp = log10_chiphm + chiph_xip_delta;
        
        d = ( xip - xip_table[ixip] ) / ( xip_table[ixip+1] - xip_table[ixip] );
        
        // Chiph after linear interpolation in the logarithmic scale
        photon_chi = pow( 10., log10_chiphm*( 1.0-d ) + log10_chiphp*( d ) );
    } else
        // For integration reasons, we can have xip_table[ixip+1] = xip_table[ixip]
        // In this case, no interpolation
    {
        photon_chi = pow( 10., ichiph*chiph_xip_delta
                          + xip_chiphmin_table[ichipa] );
    }
    
    // ------------------------------------------------------------
    // Compute photon_chi
    // Fastest method using the nearest point but less accurate
    // ------------------------------------------------------------
    
    /*photon_chi = pow(10.,ichiph*chiph_xip_delta
           + xip_chiphmin_table[ichipa]);*/
    
    
    return photon_chi;
}

// ---------------------------------------------------------------------------------------------------------------------
//! Computation of the Cross Section dNph/dt which is also
//! the number of photons generated per time unit.
//
//! \param particle_chi particle quantum parameter
//! \param particle_gamma particle gamma factor
// ---------------------------------------------------------------------------------------------------------------------
double RadiationTables::computePhotonProductionYield( double particle_chi, double particle_gamma )
{

    // Log of the particle quantum parameter particle_chi
    double logchipa;
    double logchipam;
    double logchipap;
    // Index
    int ichipa;
    // final value
    double dNphdt;
    
    logchipa = log10( particle_chi );
    
    // Lower index for interpolation in the table integfochi
    ichipa = int( floor( ( logchipa-integfochi_log10_chipa_min )
                         *integfochi_chipa_inv_delta ) );
                         
    // If we are not in the table...
    if( ichipa < 0 ) {
        ichipa = 0;
        dNphdt = integfochi_table[ichipa];
    } else if( ichipa >= integfochi_dim-1 ) {
        ichipa = integfochi_dim-2;
        dNphdt = integfochi_table[ichipa];
    } else {
        // Upper and lower values for linear interpolation
        logchipam = ichipa*integfochi_chipa_delta + integfochi_log10_chipa_min;
        logchipap = logchipam + integfochi_chipa_delta;
        
        // Interpolation
        dNphdt = ( integfochi_table[ichipa+1]*fabs( logchipa-logchipam ) +
                   integfochi_table[ichipa]*fabs( logchipap - logchipa ) )*integfochi_chipa_inv_delta;
    }
    
    return factor_dNphdt*dNphdt*particle_chi/particle_gamma;
    
}

// ---------------------------------------------------------------------------------------------------------------------
//! \brief
//! Compute integration of F/chi between
//! using Gauss-Legendre for a given particle_chi value
//
//
//! \param particle_chi particle (electron for instance) quantum parameter
//! \param min_photon_chi Minimal integration value (photon quantum parameter)
//! \param max_photon_chi Maximal integration value (photon quantum parameter)
//! \param nb_iterations number of points for integration
//! \param eps integration accuracy
// ---------------------------------------------------------------------------------------------------------------------
double RadiationTables::integrateSynchrotronEmissivity( double particle_chi,
        double min_photon_chi,
        double max_photon_chi,
        int nb_iterations,
        double eps )
{

    //std::cout << "RadiationTables::integrateSynchrotronEmissivity" << std::endl;
    
    // Arrays for Gauss-Legendre integration
    double *gauleg_x = new double[nb_iterations];
    double *gauleg_w = new double[nb_iterations];
    // Photon quantum parameter
    double photon_chi;
    // Integration result
    double integ;
    // Synchrotron emissivity
    double sync_emi;
    // Iterator
    int i;
    
    // gauss Legendre coefficients
    userFunctions::gauss_legendre_coef( log10( min_photon_chi ), log10( max_photon_chi ),
                                        gauleg_x, gauleg_w, nb_iterations, eps );
                                        
    // Integration loop
    integ = 0;
    #pragma omp parallel for reduction(+:integ) private(photon_chi,sync_emi) shared(particle_chi,gauleg_w,gauleg_x)
    for( i=0 ; i< nb_iterations ; i++ ) {
        photon_chi = pow( 10., gauleg_x[i] );
        sync_emi = RadiationTables::computeRitusSynchrotronEmissivity( particle_chi, photon_chi, 200, 1e-15 );
        integ += gauleg_w[i]*sync_emi*log( 10 );
    }
    
    return integ;
    
}

// ---------------------------------------------------------------------------------------------------------------------
//! Computation of the synchrotron emissivity following the formulae of Ritus
//
//! \param particle_chi particle quantum parameter
//! \param photon_chi photon quantum parameter
//! \param nb_iterations number of iterations for the Gauss-Legendre integration
//! \param eps epsilon for the modified bessel function
// ---------------------------------------------------------------------------------------------------------------------
double RadiationTables::computeRitusSynchrotronEmissivity( double particle_chi,
        double photon_chi, int nb_iterations, double eps )
{

    //std::cout << "RadiationTables::computeRitusSynchrotronEmissivity" << std::endl;
    
    // The photon quantum parameter should be below the electron one
    if( particle_chi > photon_chi ) {
        // Arrays for Gauss-Legendre integration
        double *gauleg_x = new double[nb_iterations];
        double *gauleg_w = new double[nb_iterations];
        // Values for Bessel results
        double I, dI;
        double K, dK;
        // Parts of the formulae
        double part1, part2;
        // Iterator
        int i;
        
        double y = photon_chi/( 3.*particle_chi*( particle_chi-photon_chi ) );
        
        // Computation of Part. 1
        // Call the modified Bessel function to get K
        userFunctions::modified_bessel_IK( 2./3., 2*y, I, dI, K, dK, 50000, eps, false );
        
        part1 = ( 2. + 3.*photon_chi*y )*( K );
        
        // Computation of Part. 2
        // Using Gauss Legendre integration
        
        userFunctions::gauss_legendre_coef( log10( 2*y ), log10( y )+50., gauleg_x,
                                            gauleg_w, nb_iterations, eps );
                                            
        part2 = 0;
        for( i=0 ; i< nb_iterations ; i++ ) {
            y = pow( 10., gauleg_x[i] );
            userFunctions::modified_bessel_IK( 1./3., y, I, dI, K, dK, 50000, eps, false );
            part2 += gauleg_w[i]*K*y*log( 10. );
        }
        
        // Factor for final result
        y = 2*photon_chi/( 3*pow( particle_chi, 2. ) );
        
        return ( part1 - part2 )*y;
        
        
    } else if( particle_chi == photon_chi ) {
        return 0;
    } else {
        ERROR( "In computeRitusSynchrotronEmissivity: particle_chi " << particle_chi
               << " < photon_chi "
               << photon_chi );
        return -1.;
    }
}

// -----------------------------------------------------------------------------
//! Return the value of the function h(particle_chi) of Niel et al.
//! Use an integration of Gauss-Legendre
//
//! \param particle_chi particle quantum parameter
//! \param nb_iterations number of iterations for the Gauss-Legendre integration
//! \param eps epsilon for the modified bessel function
// -----------------------------------------------------------------------------
double RadiationTables::computeHNiel( double particle_chi,
                                      int nb_iterations, double eps )
{
    // Arrays for Gauss-Legendre integration
    double *gauleg_x = new double[nb_iterations];
    double *gauleg_w = new double[nb_iterations];
    double nu;
    double I, dI;
    double K23, K53, dK;
    double h = 0.;
    int    i;
    
    // Gauss-Legendre coefs between log10(1E-20)
    // and log10(50.) = 1.6989700043360187
    userFunctions::gauss_legendre_coef( -20., log10( 50. ), gauleg_x,
                                        gauleg_w, nb_iterations, eps );
                                        
    for( i=0 ; i< nb_iterations ; i++ ) {
    
        nu = pow( 10., gauleg_x[i] );
        
        userFunctions::modified_bessel_IK( 5./3., nu, I, dI, K53, dK, 50000, eps, false );
        
        userFunctions::modified_bessel_IK( 2./3., nu, I, dI, K23, dK, 50000, eps, false );
        
        h += gauleg_w[i]*log( 10. )*nu
             * ( ( 2.*pow( particle_chi*nu, 3 ) / pow( 2. + 3.*nu*particle_chi, 3 )*K53 )
                 + ( 54.*pow( particle_chi, 5 )*pow( nu, 4 ) / pow( 2. + 3.*nu*particle_chi, 5 )*K23 ) );
                 
    }
    
    return 9.*sqrt( 3. )/( 4.*M_PI )*h;
    
}

// -----------------------------------------------------------------------------
//! Return the value of the function h(particle_chi) of Niel et al.
//! from the computed table h_table
//! \param particle_chi particle quantum parameter
// -----------------------------------------------------------------------------
double RadiationTables::getHNielFromTable( double particle_chi )
{
    int ichipa;
    double d;
    
    // Position in the h_table
    d = ( log10( particle_chi )-h_log10_chipa_min )*h_chipa_inv_delta;
    ichipa = int( floor( d ) );
    
    // distance for interpolation
    d = d - floor( d );
    
    // Linear interpolation
    return h_table[ichipa]*( 1.-d ) + h_table[ichipa+1]*( d );
}

// -----------------------------------------------------------------------------
//! Return the stochastic diffusive component of the pusher
//! of Niel et al.
//! \param gamma particle Lorentz factor
//! \param particle_chi particle quantum parameter
// -----------------------------------------------------------------------------
double RadiationTables::getNielStochasticTerm( double gamma,
        double particle_chi,
        double sqrtdt )
{
    // Get the value of h for the corresponding particle_chi
    double h, r;
    
    h = RadiationTables::getHNielFromTable( particle_chi );
    
    // Pick a random number in the normal distribution of standard
    // deviation sqrt(dt) (variance dt)
    r = Rand::normal( sqrtdt );
    
    /*std::random_device device;
    std::mt19937 gen(device());
    std::normal_distribution<double> normal_distribution(0., sqrt(dt));
    r = normal_distribution(gen);*/
    
    return sqrt( factor_classical_radiated_power_*gamma*h )*r;
}

// -----------------------------------------------------------------------------
// TABLE READING
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//! Read the external table h
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
bool RadiationTables::read_h_table( SmileiMPI *smpi )
{
    // Flag database available
    bool table_exists = false;
    
    // Test if an external table exists, if yes we read the table...
    // Binary table
    if( Tools::file_exists( table_path + "/tab_h.bin" ) ) {
    
        table_exists = true;
        
        if( smpi->getRank()==0 ) {
        
            // Reading of the table file
            std::ifstream file;
            file.open( table_path + "/tab_h.bin", std::ios::binary );
            
            if( file.is_open() ) {
            
                // Read the header
                file.read( ( char * )&h_dim,
                           sizeof( h_dim ) );
                file.read( ( char * )&h_chipa_min,
                           sizeof( h_chipa_min ) );
                file.read( ( char * )&h_chipa_max,
                           sizeof( h_chipa_max ) );
                           
                // Resize of the array integfochi_table before reading
                h_table.resize( h_dim );
                
                // Read the table values
                file.read( ( char * )&h_table[0], sizeof( double )*h_dim );
                
                file.close();
            }
            
        }
        
    }
    // HDF5 format
    else if( Tools::file_exists( table_path + "/radiation_tables.h5" ) ) {
        hid_t       fileId;
        hid_t       datasetId;
        std::string buffer;
        
        if( smpi->getRank()==0 ) {
        
            buffer = table_path + "/radiation_tables.h5";
            
            fileId = H5Fopen( buffer.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
            
            datasetId = H5Dopen2( fileId, "h", H5P_DEFAULT );
            
            // If this dataset exists, we read it
            if( datasetId > 0 ) {
            
                table_exists = true;
                
                // First, we read attributes
                H5::getAttr( datasetId, "chipa_dim", h_dim );
                H5::getAttr( datasetId, "chipa_min", h_chipa_min );
                H5::getAttr( datasetId, "chipa_max", h_chipa_max );
                
                // Resize of the array integfochi_table before reading
                h_table.resize( h_dim );
                
                // then the dataset
                H5Dread( datasetId,
                         H5T_NATIVE_DOUBLE, H5S_ALL,
                         H5S_ALL, H5P_DEFAULT,
                         &h_table[0] );
                         
                H5Dclose( datasetId );
                H5Fclose( fileId );
            }
            // Else, we will have to compute it
            else {
                table_exists = false;
            }
        }
        
        // Bcast table_exists
        int TE = table_exists;
        MPI_Bcast( &TE, 1, MPI_INT, 0, smpi->getGlobalComm() );
        table_exists = TE;
        
    }
    
    // If the table exists, they have been read...
    if( table_exists ) {
    
        // checks
        if( minimum_chi_continuous_ < h_chipa_min ) {
            ERROR( "Parameter `minimum_chi_continuous_` is below `h_chipa_min`,"
                   << "the lower bound of the h table should be equal or below"
                   << "the radiation threshold on chi." )
        }
        
        MESSAGE( "            Reading of the external database" );
        MESSAGE( "            Dimension quantum parameter: "
                 << h_dim );
        MESSAGE( "            Minimum particle quantum parameter chi: "
                 << h_chipa_min );
        MESSAGE( "            Maximum particle quantum parameter chi: "
                 << h_chipa_max );
                 
        // Bcast the table to all MPI ranks
        RadiationTables::bcast_h_table( smpi );
    }
    
    return table_exists;
    
}

// -----------------------------------------------------------------------------
//! Read the external table integfochi
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
bool RadiationTables::read_integfochi_table( SmileiMPI *smpi )
{
    // Flag database available
    bool table_exists = false;
    
    // Test if an external table exists, if yes we read the table...
    // Binary table
    if( Tools::file_exists( table_path + "/tab_integfochi.bin" ) ) {
    
        table_exists = true;
        
        if( smpi->getRank()==0 ) {
        
            // Reading of the table file
            std::ifstream file;
            file.open( table_path + "/tab_integfochi.bin", std::ios::binary );
            
            if( file.is_open() ) {
            
                // Read the header
                file.read( ( char * )&integfochi_dim, sizeof( integfochi_dim ) );
                file.read( ( char * )&integfochi_chipa_min,
                           sizeof( integfochi_chipa_min ) );
                file.read( ( char * )&integfochi_chipa_max,
                           sizeof( integfochi_chipa_max ) );
                           
                // Resize of the array integfochi_table before reading
                integfochi_table.resize( integfochi_dim );
                
                // Read the table values
                file.read( ( char * )&integfochi_table[0], sizeof( double )*integfochi_dim );
                
                file.close();
            }
            
        }
        
    }
    // HDF5 format
    else if( Tools::file_exists( table_path + "/radiation_tables.h5" ) ) {
        hid_t       fileId;
        hid_t       datasetId;
        std::string buffer;
        
        if( smpi->getRank()==0 ) {
        
            buffer = table_path + "/radiation_tables.h5";
            
            fileId = H5Fopen( buffer.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
            
            datasetId = H5Dopen2( fileId, "integfochi", H5P_DEFAULT );
            
            // If this dataset exists, we read it
            if( datasetId > 0 ) {
            
                table_exists = true;
                
                // First, we read attributes
                H5::getAttr( datasetId, "chipa_dim", integfochi_dim );
                H5::getAttr( datasetId, "chipa_min", integfochi_chipa_min );
                H5::getAttr( datasetId, "chipa_max", integfochi_chipa_max );
                
                // Resize of the array integfochi_table before reading
                integfochi_table.resize( integfochi_dim );
                
                // then the dataset
                H5Dread( datasetId,
                         H5T_NATIVE_DOUBLE, H5S_ALL,
                         H5S_ALL, H5P_DEFAULT,
                         &integfochi_table[0] );
                         
                H5Dclose( datasetId );
                H5Fclose( fileId );
            }
            // Else, we will have to compute it
            else {
                table_exists = false;
            }
        }
        
        // Bcast table_exists
        int TE = table_exists;
        MPI_Bcast( &TE, 1, MPI_INT, 0, smpi->getGlobalComm() );
        table_exists = TE;
        
    }
    
    // If the table exists, they have been read...
    if( table_exists ) {
    
        MESSAGE( "            Reading of the external database" );
        MESSAGE( "            Dimension quantum parameter: " << integfochi_dim );
        MESSAGE( "            Minimum particle quantum parameter chi: " << integfochi_chipa_min );
        MESSAGE( "            Maximum particle quantum parameter chi: " << integfochi_chipa_max );
        
        // Bcast the table to all MPI ranks
        RadiationTables::bcast_integfochi_table( smpi );
    }
    
    return table_exists;
    
}

// -----------------------------------------------------------------------------
//! Read the external table xip_chiphmin and xip
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
bool RadiationTables::read_xip_table( SmileiMPI *smpi )
{
    // Flag database available
    bool table_exists = false;
    
    // Test if an external table exists, we read the table...
    if( Tools::file_exists( table_path + "/tab_xip.bin" ) ) {
    
        table_exists = true;
        
        if( smpi->getRank()==0 ) {
        
            // Reading of the table file
            std::ifstream file;
            file.open( table_path + "/tab_xip.bin", std::ios::binary );
            
            if( file.is_open() ) {
            
                // Read the header
                file.read( ( char * )&xip_chipa_dim, sizeof( xip_chipa_dim ) );
                file.read( ( char * )&xip_chiph_dim, sizeof( xip_chiph_dim ) );
                file.read( ( char * )&xip_chipa_min, sizeof( xip_chipa_min ) );
                file.read( ( char * )&xip_chipa_max, sizeof( xip_chipa_max ) );
                
                // Allocation of the array xip
                xip_chiphmin_table.resize( xip_chipa_dim );
                xip_table.resize( xip_chipa_dim*xip_chiph_dim );
                
                // Read the table values
                file.read( ( char * )&xip_chiphmin_table[0], sizeof( double )*xip_chipa_dim );
                
                // Read the table values
                file.read( ( char * )&xip_table[0], sizeof( double )*xip_chipa_dim*xip_chiph_dim );
                
                file.close();
            }
            
        }
    }
    // HDF5 format
    else if( Tools::file_exists( table_path + "/radiation_tables.h5" ) ) {
        if( smpi->getRank()==0 ) {
        
            hid_t       fileId;
            hid_t       datasetId_chiphmin;
            hid_t       datasetId_xip;
            std::string buffer;
            
            buffer = table_path + "/radiation_tables.h5";
            
            fileId = H5Fopen( buffer.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
            
            datasetId_chiphmin = H5Dopen2( fileId, "xip_chiphmin", H5P_DEFAULT );
            datasetId_xip = H5Dopen2( fileId, "xip", H5P_DEFAULT );
            
            // If this dataset exists, we read it
            if( datasetId_chiphmin > 0 && datasetId_xip > 0 ) {
            
                table_exists = true;
                
                // First, we read attributes
                H5::getAttr( datasetId_xip, "chipa_dim", xip_chipa_dim );
                H5::getAttr( datasetId_xip, "chiph_dim", xip_chiph_dim );
                H5::getAttr( datasetId_xip, "chipa_min", xip_chipa_min );
                H5::getAttr( datasetId_xip, "chipa_max", xip_chipa_max );
                
                // Allocation of the array xip
                xip_chiphmin_table.resize( xip_chipa_dim );
                xip_table.resize( xip_chipa_dim*xip_chiph_dim );
                
                // then the dataset for chiphmin
                H5Dread( datasetId_chiphmin,
                         H5T_NATIVE_DOUBLE, H5S_ALL,
                         H5S_ALL, H5P_DEFAULT,
                         &xip_chiphmin_table[0] );
                         
                // then the dataset for xip
                H5Dread( datasetId_xip,
                         H5T_NATIVE_DOUBLE, H5S_ALL,
                         H5S_ALL, H5P_DEFAULT,
                         &xip_table[0] );
                         
                H5Dclose( datasetId_xip );
                H5Dclose( datasetId_chiphmin );
                H5Fclose( fileId );
            }
            // Else, we will have to compute it
            else {
                table_exists = false;
            }
        }
        
        // Bcast table_exists
        MPI_Bcast( &table_exists, 1, MPI_INT, 0, smpi->getGlobalComm() );
        
    }
    
    // If the table exists, they have been read...
    if( table_exists ) {
    
        MESSAGE( "            Reading of the external database" );
        MESSAGE( "            Dimension particle chi: " << xip_chipa_dim );
        MESSAGE( "            Dimension photon chi: " << xip_chiph_dim );
        MESSAGE( "            Minimum particle chi: " << xip_chipa_min );
        MESSAGE( "            Maximum particle chi: " << xip_chipa_max );
        
        // Bcast the table to all MPI ranks
        RadiationTables::bcast_xip_table( smpi );
    }
    
    return table_exists;
}

// -----------------------------------------------------------------------------
// TABLE COMMUNICATIONS
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//! Bcast of the external table h for the Niel radiation model
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
void RadiationTables::bcast_h_table( SmileiMPI *smpi )
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
        MPI_Pack_size( h_dim, MPI_DOUBLE, smpi->getGlobalComm(),
                       &position );
        buf_size += position;
    }
    
    MESSAGE( "            Buffer size: " << buf_size );
    
    // Exchange buf_size with all ranks
    MPI_Bcast( &buf_size, 1, MPI_INT, 0, smpi->getGlobalComm() );
    
    // Packet that will contain all parameters
    char *buffer = new char[buf_size];
    
    // Proc 0 packs
    if( smpi->getRank() == 0 ) {
        position = 0;
        MPI_Pack( &h_dim,
                  1, MPI_INT, buffer, buf_size, &position, smpi->getGlobalComm() );
        MPI_Pack( &h_chipa_min,
                  1, MPI_DOUBLE, buffer, buf_size, &position, smpi->getGlobalComm() );
        MPI_Pack( &h_chipa_max,
                  1, MPI_DOUBLE, buffer, buf_size, &position, smpi->getGlobalComm() );
                  
        MPI_Pack( &h_table[0], h_dim,
                  MPI_DOUBLE, buffer, buf_size, &position, smpi->getGlobalComm() );
                  
    }
    
    // Bcast all parameters
    MPI_Bcast( &buffer[0], buf_size, MPI_PACKED, 0, smpi->getGlobalComm() );
    
    // Other ranks unpack
    if( smpi->getRank() != 0 ) {
        position = 0;
        MPI_Unpack( buffer, buf_size, &position,
                    &h_dim, 1, MPI_INT, smpi->getGlobalComm() );
        MPI_Unpack( buffer, buf_size, &position,
                    &h_chipa_min, 1, MPI_DOUBLE, smpi->getGlobalComm() );
        MPI_Unpack( buffer, buf_size, &position,
                    &h_chipa_max, 1, MPI_DOUBLE, smpi->getGlobalComm() );
                    
        // Resize table before unpacking values
        h_table.resize( h_dim );
        
        MPI_Unpack( buffer, buf_size, &position, &h_table[0],
                    h_dim, MPI_DOUBLE, smpi->getGlobalComm() );
                    
    }
    
    h_log10_chipa_min = log10( h_chipa_min );
    
    // Computation of the delta
    h_chipa_delta = ( log10( h_chipa_max )
                      - h_log10_chipa_min )/( h_dim-1 );
                      
    // Inverse delta
    h_chipa_inv_delta = 1./h_chipa_delta;
}

// -----------------------------------------------------------------------------
//! Bcast of the external table integfochi
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
void RadiationTables::bcast_integfochi_table( SmileiMPI *smpi )
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
        MPI_Pack_size( integfochi_dim, MPI_DOUBLE, smpi->getGlobalComm(),
                       &position );
        buf_size += position;
    }
    
    MESSAGE( "            Buffer size: " << buf_size );
    
    // Exchange buf_size with all ranks
    MPI_Bcast( &buf_size, 1, MPI_INT, 0, smpi->getGlobalComm() );
    
    // Packet that will contain all parameters
    char *buffer = new char[buf_size];
    
    // Proc 0 packs
    if( smpi->getRank() == 0 ) {
        position = 0;
        MPI_Pack( &integfochi_dim,
                  1, MPI_INT, buffer, buf_size, &position, smpi->getGlobalComm() );
        MPI_Pack( &integfochi_chipa_min,
                  1, MPI_DOUBLE, buffer, buf_size, &position, smpi->getGlobalComm() );
        MPI_Pack( &integfochi_chipa_max,
                  1, MPI_DOUBLE, buffer, buf_size, &position, smpi->getGlobalComm() );
                  
        MPI_Pack( &integfochi_table[0], integfochi_dim,
                  MPI_DOUBLE, buffer, buf_size, &position, smpi->getGlobalComm() );
                  
    }
    
    // Bcast all parameters
    MPI_Bcast( &buffer[0], buf_size, MPI_PACKED, 0, smpi->getGlobalComm() );
    
    // Other ranks unpack
    if( smpi->getRank() != 0 ) {
        position = 0;
        MPI_Unpack( buffer, buf_size, &position,
                    &integfochi_dim, 1, MPI_INT, smpi->getGlobalComm() );
        MPI_Unpack( buffer, buf_size, &position,
                    &integfochi_chipa_min, 1, MPI_DOUBLE, smpi->getGlobalComm() );
        MPI_Unpack( buffer, buf_size, &position,
                    &integfochi_chipa_max, 1, MPI_DOUBLE, smpi->getGlobalComm() );
                    
        // Resize table before unpacking values
        integfochi_table.resize( integfochi_dim );
        
        MPI_Unpack( buffer, buf_size, &position, &integfochi_table[0],
                    integfochi_dim, MPI_DOUBLE, smpi->getGlobalComm() );
                    
    }
    
    integfochi_log10_chipa_min = log10( integfochi_chipa_min );
    
    // Computation of the delta
    integfochi_chipa_delta = ( log10( integfochi_chipa_max )
                               - integfochi_log10_chipa_min )/( integfochi_dim-1 );
                               
    // Inverse delta
    integfochi_chipa_inv_delta = 1./integfochi_chipa_delta;
}

// -----------------------------------------------------------------------------
//! Bcast of the external table xip_chiphmin and xip
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
void RadiationTables::bcast_xip_table( SmileiMPI *smpi )
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
        MPI_Pack_size( xip_chipa_dim, MPI_DOUBLE, smpi->getGlobalComm(),
                       &position );
        buf_size += position;
        MPI_Pack_size( xip_chipa_dim*xip_chiph_dim, MPI_DOUBLE,
                       smpi->getGlobalComm(), &position );
        buf_size += position;
    }
    
    MESSAGE( "            Buffer size for MPI exchange: " << buf_size );
    
    // Exchange buf_size with all ranks
    MPI_Bcast( &buf_size, 1, MPI_INT, 0, smpi->getGlobalComm() );
    
    // Packet that will contain all parameters
    char *buffer = new char[buf_size];
    
    // Proc 0 packs
    if( smpi->getRank() == 0 ) {
        position = 0;
        MPI_Pack( &xip_chipa_dim,
                  1, MPI_INT, buffer, buf_size, &position, smpi->getGlobalComm() );
        MPI_Pack( &xip_chiph_dim,
                  1, MPI_INT, buffer, buf_size, &position, smpi->getGlobalComm() );
        MPI_Pack( &xip_chipa_min,
                  1, MPI_DOUBLE, buffer, buf_size, &position, smpi->getGlobalComm() );
        MPI_Pack( &xip_chipa_max,
                  1, MPI_DOUBLE, buffer, buf_size, &position, smpi->getGlobalComm() );
                  
        MPI_Pack( &xip_chiphmin_table[0], xip_chipa_dim,
                  MPI_DOUBLE, buffer, buf_size, &position, smpi->getGlobalComm() );
                  
        MPI_Pack( &xip_table[0], xip_chipa_dim*xip_chiph_dim,
                  MPI_DOUBLE, buffer, buf_size, &position, smpi->getGlobalComm() );
    }
    
    // Bcast all parameters
    MPI_Bcast( &buffer[0], buf_size, MPI_PACKED, 0, smpi->getGlobalComm() );
    
    // Other ranks unpack
    if( smpi->getRank() != 0 ) {
        position = 0;
        MPI_Unpack( buffer, buf_size, &position,
                    &xip_chipa_dim, 1, MPI_INT, smpi->getGlobalComm() );
        MPI_Unpack( buffer, buf_size, &position,
                    &xip_chiph_dim, 1, MPI_INT, smpi->getGlobalComm() );
        MPI_Unpack( buffer, buf_size, &position,
                    &xip_chipa_min, 1, MPI_DOUBLE, smpi->getGlobalComm() );
        MPI_Unpack( buffer, buf_size, &position,
                    &xip_chipa_max, 1, MPI_DOUBLE, smpi->getGlobalComm() );
                    
        // Resize tables before unpacking values
        xip_chiphmin_table.resize( xip_chipa_dim );
        xip_table.resize( xip_chipa_dim*xip_chiph_dim );
        
        MPI_Unpack( buffer, buf_size, &position, &xip_chiphmin_table[0],
                    xip_chipa_dim, MPI_DOUBLE, smpi->getGlobalComm() );
                    
        MPI_Unpack( buffer, buf_size, &position, &xip_table[0],
                    xip_chipa_dim*xip_chiph_dim, MPI_DOUBLE, smpi->getGlobalComm() );
    }
    
    // Log10 of xip_chipa_min for efficiency
    xip_log10_chipa_min = log10( xip_chipa_min );
    
    // Computation of the delta
    xip_chipa_delta = ( log10( xip_chipa_max )
                        - xip_log10_chipa_min )/( xip_chipa_dim-1 );
                        
    // Inverse of delta
    xip_chipa_inv_delta = 1./xip_chipa_delta;
    
    // Inverse photon_chi discetization (regularly used)
    xip_inv_chiph_dim_minus_one = 1./( xip_chiph_dim - 1. );
    
}
