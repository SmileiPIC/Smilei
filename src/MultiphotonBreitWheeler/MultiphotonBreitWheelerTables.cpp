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

#include <iomanip>

#include "MultiphotonBreitWheelerTables.h"

// -----------------------------------------------------------------------------
// INITILIZATION AND DESTRUCTION
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Constructor for MutliphotonBreitWheelerTables
// -----------------------------------------------------------------------------
MultiphotonBreitWheelerTables::MultiphotonBreitWheelerTables()
{
    // T table
    T_table.resize( 0 );
    T_computed = false;
    xip_computed = false;
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

        // Flag to activate the table computation
        PyTools::extract( "compute_table", compute_table_, "MultiphotonBreitWheeler" );

        // If we compute the table, we read the table parameters
        if (compute_table_)
        {
            // Extraction of the parameter from the input file
            PyTools::extract( "T_chiph_min", T_chiph_min, "MultiphotonBreitWheeler" );
            PyTools::extract( "T_chiph_max", T_chiph_max, "MultiphotonBreitWheeler" );
            PyTools::extract( "T_dim", T_dim, "MultiphotonBreitWheeler" );

            T_log10_chiph_min = log10( T_chiph_min );

            // xip table
            PyTools::extract( "xip_chiph_min", xip_chiph_min, "MultiphotonBreitWheeler" );
            PyTools::extract( "xip_chiph_max", xip_chiph_max, "MultiphotonBreitWheeler" );
            PyTools::extract( "xip_chiph_dim", xip_chiph_dim, "MultiphotonBreitWheeler" );
            PyTools::extract( "xip_chipa_dim", xip_chipa_dim, "MultiphotonBreitWheeler" );
            PyTools::extract( "xip_power", xip_power, "MultiphotonBreitWheeler" );
            PyTools::extract( "xip_threshold", xip_threshold, "MultiphotonBreitWheeler" );

            if( T_chiph_min >= T_chiph_max ) {
                ERROR( "T_chiph_min (" << T_chiph_min
                       << ") >= T_chiph_max (" << T_chiph_max << ")" );
            }
            if( T_dim < 1 ) {
                ERROR( "T_dim is too low" );
            }
            if( xip_chiph_min >= xip_chiph_max ) {
                ERROR( "xip_chiph_min (" << xip_chiph_min
                       << ") >= xip_chiph_max (" << xip_chiph_max << ")" );
            }
            if( xip_chiph_dim < 1 ) {
                ERROR( "In MultiphotonBreitWheelerTables, xip_chiph_dim must be above 0" );
            }
            if( xip_chipa_dim < 1 ) {
                ERROR( "In MultiphotonBreitWheelerTables, xip_chipa_dim must be above 0" );
            }

            // Additional regularly used parameters
            xip_log10_chiph_min = log10( xip_chiph_min );
            xip_inv_chipa_dim_minus_one = 1./( xip_chipa_dim - 1. );

            // Format of the tables
            PyTools::extract( "output_format", output_format_, "MultiphotonBreitWheeler" );
        }

        // Path to the databases
        PyTools::extract( "table_path", table_path_, "MultiphotonBreitWheeler" );
    }

    // Computation of some parameters
    if( params.hasMultiphotonBreitWheeler ) {
        // Computation of the normalized Compton wavelength
        normalized_Compton_wavelength_ = params.red_planck_cst*params.reference_angular_frequency_SI
                                         / ( params.electron_mass*params.c_vacuum*params.c_vacuum );

        // Computation of the factor factor_dNBWdt
        factor_dNBWdt = params.fine_struct_cst/( normalized_Compton_wavelength_ );
    }

    // Messages and checks
    if( params.hasMultiphotonBreitWheeler ) {
        MESSAGE( 1,"table path: " << table_path_ );
    }

    MESSAGE( "" )

    // If the user does not request to compute the tables, then we read them
    if (!compute_table_) {
        readTables( params, smpi );
    }
    // Else, the tables are computed and store in the provided path `table_path_`
    else {
        computeTables( params, smpi );
        outputTables( smpi );
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
double MultiphotonBreitWheelerTables::compute_dNBWdt( double photon_chi, double gamma )
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
    ichiph = int( floor( ( logchiph-T_log10_chiph_min )
                         *T_chiph_inv_delta ) );

    // If photon_chi is below the lower bound of the table
    // An asymptotic approximation is used
    if( ichiph < 0 ) {
        ichiph = 0;
        dNBWdt = 0.2296*exp( -8./( 3.*photon_chi ) );
    }
    // If photon_chi is above the upper bound of the table
    // An asymptotic approximation is used
    else if( ichiph >= T_dim-1 ) {
        ichiph = T_dim-2;
        dNBWdt = 0.38*pow( photon_chi, -1./3. );
    } else {
        // Upper and lower values for linear interpolation
        logchiphm = ichiph*T_chiph_delta + T_log10_chiph_min;
        logchiphp = logchiphm + T_chiph_delta;

        // Interpolation
        dNBWdt = ( T_table[ichiph+1]*fabs( logchiph-logchiphm ) +
                   T_table[ichiph]*fabs( logchiphp - logchiph ) )*T_chiph_inv_delta;
    }
    return factor_dNBWdt*dNBWdt*photon_chi/gamma;
}

// -----------------------------------------------------------------------------
//! Computation of the value T(photon_chi) using the approximated
//! formula of Erber
//! \param photon_chi photon quantum parameter
//! \param nb_iterations number of iteration for the Bessel evaluation
//! \param eps epsilon for the Bessel evaluation
// -----------------------------------------------------------------------------
double MultiphotonBreitWheelerTables::compute_Erber_T( double photon_chi, int nb_iterations,
        double eps )
{
    // Values for Bessel results
    //double I,dI;
    double K;

    //userFunctions::modifiedBesselIK(1./3.,4./(3.*photon_chi),I,dI,K,dK,nb_iterations,eps);
    K = userFunctions::modifiedBesselK( 1./3., 4./( 3.*photon_chi ), nb_iterations, eps, false );

    return 0.16*K*K/photon_chi;
}

// -----------------------------------------------------------------------------
//! Computation of the value of the integration of  dT(photon_chi)/dhicpa
//! using the formula of Ritus.
//! \param photon_chi photon quantum parameter
//! \param particle_chi particle quantum parameter for integration (=0.5*photon_chi for full integration)
//! \param nb_iterations number of iteration for the Gauss-Legendre integration
//! \param eps epsilon for the Bessel evaluation
// -----------------------------------------------------------------------------
double MultiphotonBreitWheelerTables::compute_integration_Ritus_dTdchi( double photon_chi,
        double particle_chi,
        int nb_iterations,
        double eps )
{

    // Parameters:
    // Arrays for Gauss-Legendre integration
    double *gauleg_x = new double[nb_iterations];
    double *gauleg_w = new double[nb_iterations];
    int i;
    double u;
    double T;

    // gauss Legendre coefficients
    userFunctions::gaussLegendreCoef( log10( particle_chi )-50., log10( particle_chi ),
                                        gauleg_x, gauleg_w, nb_iterations, eps );

    T = 0;
    for( i=0 ; i< nb_iterations ; i++ ) {
        u = pow( 10., gauleg_x[i] );
        T += gauleg_w[i]*compute_Ritus_dTdchi( photon_chi, u, nb_iterations, eps )*u*log( 10. );
    }

    delete[] gauleg_x;
    delete[] gauleg_w;

    return T;

}

// -----------------------------------------------------------------------------
//! Computation of the electron and positron quantum parameters for
//! the multiphoton Breit-Wheeler pair creation
//
//! \param photon_chi photon quantum parameter
// -----------------------------------------------------------------------------
double *MultiphotonBreitWheelerTables::compute_pair_chi( double photon_chi )
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
    if( photon_chi < xip_chiph_min ) {
        ichiph = 0;
    }
    // Upper boundary of the table
    else if( photon_chi >= xip_chiph_max ) {
        ichiph = xip_chiph_dim-1;
    }
    // Inside the table
    else {
        // Use floor so that photon_chi corresponding to ichiph is <= given photon_chi
        ichiph = int( floor( ( logchiph-xip_log10_chiph_min )*( xip_chiph_inv_delta ) ) );
    }

    // ---------------------------------------
    // Search of the index ichiph for photon_chi
    // ---------------------------------------

    // First, we compute a random xip in [0,1[
    xip = Rand::uniform();

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
    if( xipp < xip_table[ichiph*xip_chipa_dim] ) {
        ichipa = 0;
    }
    // Upper bound
    else if( xipp >= xip_table[( ichiph+1 )*xip_chipa_dim-1] ) {
        ichipa = xip_chipa_dim-2;
    } else {
        // Search for the corresponding index ichipa for xip
        ichipa = userFunctions::searchValuesInMonotonicArray(
                     &xip_table[ichiph*xip_chipa_dim], xipp, xip_chipa_dim );
    }

    // Delta for the particle_chi dimension
    delta_chipa = ( log10( 0.5*photon_chi )-xip_chipamin_table[ichiph] )
                  * xip_inv_chipa_dim_minus_one;

    ixip = ichiph*xip_chipa_dim + ichipa;

    log10_chipam = ichipa*delta_chipa + xip_chipamin_table[ichiph];
    log10_chipap = log10_chipam + delta_chipa;

    d = ( xipp - xip_table[ixip] ) / ( xip_table[ixip+1] - xip_table[ixip] );

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
//! Computation of the value dT/dchipa(photon_chi) using the formula of Ritus
//! \param photon_chi photon quantum parameter
//! \param particle_chi particle quantum parameter
//! \param nb_iterations number of iteration for the Gauss-Legendre integration
//! \param eps epsilon for the Bessel evaluation
// -----------------------------------------------------------------------------
double MultiphotonBreitWheelerTables::compute_Ritus_dTdchi( double photon_chi,
        double particle_chi, int nb_iterations, double eps )
{

    // Parameters:
    // Arrays for Gauss-Legendre integration
    double *gauleg_x = new double[nb_iterations];
    double *gauleg_w = new double[nb_iterations];
    double y, u;
    double p1, p2;
    int i;
    //double I,dI;
    //double K,dK;

    y = ( photon_chi/( 3.*particle_chi*( photon_chi-particle_chi ) ) );

    //userFunctions::modifiedBesselIK(2./3.,2.*y,I,dI,K,dK,500,1e-16);
    //p1 = (2. - 3.*photon_chi*y)*K;

    p1 = ( 2. - 3.*photon_chi*y )*userFunctions::modifiedBesselK( 2./3., 2*y, 500, 1e-16, false );

    // gauss Legendre coefficients
    userFunctions::gaussLegendreCoef( log10( 2*y ), log10( 2*y )+50.,
                                        gauleg_x, gauleg_w, nb_iterations, eps );

    // Integration loop
    p2 = 0;
    //#pragma omp parallel for reduction(+:p2) private(u) shared(gauleg_w,gauleg_x)
    for( i=0 ; i< nb_iterations ; i++ ) {
        u = pow( 10., gauleg_x[i] );
        p2 += gauleg_w[i]*userFunctions::modifiedBesselK( 1./3., u, 500, 1e-16, false )*u*log( 10. );
        //userFunctions::modifiedBesselIK(1./3.,u,I,dI,K,dK,500,1e-16);
        //p2 += gauleg_w[i]*K*u*log(10.);
    }

    delete[] gauleg_x;
    delete[] gauleg_w;

    return ( p2 - p1 )/( M_PI*sqrt( 3. )*pow( photon_chi, 2 ) );
}

// -----------------------------------------------------------------------------
// TABLE COMPUTATION
// -----------------------------------------------------------------------------

//! Computation of the table T_table that is a discetization of the T function
//! for the multiphoton Breit-Wheeler process
//! \param smpi Object of class SmileiMPI containing MPI properties
void MultiphotonBreitWheelerTables::computeTtable( SmileiMPI *smpi )
{
    // Parameters
    int rank; // Rank number
    // timers
    double t0, t1;

    // Get the MPI rank
    rank = smpi->getRank();

    MESSAGE( 1,"--- T table:" );

    // Initial timer
    t0 = MPI_Wtime();

    // Temporary photon chi value
    double photon_chi;
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

    // Allocation of the array T_table
    T_table.resize( T_dim );

    // Allocation of the table for load repartition
    imin_table = new int[nb_ranks];
    length_table = new int[nb_ranks];

    // Computation of the delta
    T_chiph_delta = ( log10( T_chiph_max )
                      - T_log10_chiph_min )/( T_dim-1 );

    // Inverse delta
    T_chiph_inv_delta = 1./T_chiph_delta;

    // Load repartition
    userFunctions::distributeArray( nb_ranks,
            T_dim,
            imin_table,
            length_table );

    // Allocation of the local buffer
    buffer = new double [length_table[rank]];


    MESSAGE( 2,"MPI repartition:" );
    // Print repartition
    if( rank==0 ) {
        for( int i =0 ; i < nb_ranks ; i++ ) {
            MESSAGE( 2,"Rank: " << i
                     << " imin: " << imin_table[i]
                     << " length: " << length_table[i] );
        }
    }

    MESSAGE( 2,"Computation:" );
    dpct = std::max( dpct, 100./length_table[rank] );
    // Loop over the table values
    for( int i = 0 ; i < length_table[rank] ; i++ ) {
        photon_chi = pow( 10., ( imin_table[rank] + i )*T_chiph_delta
                          + T_log10_chiph_min );

        buffer[i] = 2.*MultiphotonBreitWheelerTables::compute_integration_Ritus_dTdchi( photon_chi, 0.5*photon_chi, 200, 1e-15 );
        //buffer[i] = MultiphotonBreitWheelerTables::compute_Erber_T(photon_chi,500000,1e-10);

        if( 100.*i >= length_table[rank]*pct ) {
            pct += dpct;
            MESSAGE( 3,i + 1<< "/" << length_table[rank]
                     << " - " << ( int )( std::round( pct ) )
                     << "%" );
        }
    }

    // Communication of the data
    MPI_Allgatherv( &buffer[0], length_table[rank], MPI_DOUBLE,
                    &T_table[0], &length_table[0], &imin_table[0],
                    MPI_DOUBLE, smpi->getGlobalComm() );

    // flag computed at true
    T_computed = true;

    // Free memory
    delete[] length_table;
    delete[] buffer;
    delete[] imin_table;


    // Final timer
    t1 = MPI_Wtime();
    MESSAGE( 2,"Done in " << ( t1 - t0 ) << "s" );
}

// -----------------------------------------------------------------------------
//! Computation of the minimum particle quantum parameter chipamin
//! for the photon xip array and computation of the photon xip array.
//
//! \details Under the minimum particle_chi value, the particle kinetic energy is
//! considered negligible. All energy goes to the other.
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
void MultiphotonBreitWheelerTables::computeXipTable( SmileiMPI *smpi )
{

    // Parameters
    int rank; // Rank number
    // timers
    double t0, t1;

    // Get the MPI rank
    rank = smpi->getRank();

    t0 = MPI_Wtime();

    MESSAGE( 1,"--- Table chipamin and xip:" );

    // Parameters:
    double particle_chi; // Temporary particle chi value
    double photon_chi; // Temporary photon chi value
    double chipa_delta; // Temporary delta for photon_chi
    double logchipamin; // Temporary log10 of photon chi value
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
    xip_chipamin_table.resize( xip_chiph_dim );

    // Allocation of the array xip_table
    xip_table.resize( xip_chipa_dim*xip_chiph_dim );

    // Allocation of the table for load repartition
    imin_table = new int[nb_ranks];
    length_table = new int[nb_ranks];

    // Computation of the delta
    xip_chiph_delta = ( log10( xip_chiph_max )
                        - xip_log10_chiph_min )/( xip_chiph_dim-1 );

    // Inverse of delta
    xip_chiph_inv_delta = 1./xip_chiph_delta;

    // Load repartition
    userFunctions::distributeArray( nb_ranks,
            xip_chiph_dim,
            imin_table,
            length_table );

    // Allocation of the local buffer
    buffer = new double [length_table[rank]];

    MESSAGE( 2,"MPI repartition:" );
    // Print repartition
    if( rank==0 ) {
        for( int i =0 ; i < nb_ranks ; i++ ) {
            MESSAGE( 3,"Rank: " << i
                     << " imin: "   << imin_table[i]
                     << " length: " << length_table[i] );
        }
    }

    // _______________________________________________________________
    // 1. - Computation of xip_chipamin_table

    MESSAGE( 2,"Computation of log10(chipamin):" );
    dpct = std::max( 10., 100./length_table[rank] );

    // Loop for chiphmin
    for( int ichiph = 0 ; ichiph < length_table[rank] ; ichiph++ ) {

        xip = 1;
        logchipamin = ( imin_table[rank] + ichiph )*xip_chiph_delta
                      + xip_log10_chiph_min;

        photon_chi = pow( 10., logchipamin );

        logchipamin = log10( 0.5*photon_chi );

        // Denominator of xip
        denominator = MultiphotonBreitWheelerTables::compute_integration_Ritus_dTdchi( photon_chi,
                      0.5*photon_chi, 200, 1e-15 );

        k = 0;
        while( k < xip_power ) {
            logchipamin -= pow( 0.1, k );
            particle_chi = pow( 10., logchipamin );
            numerator = MultiphotonBreitWheelerTables::compute_integration_Ritus_dTdchi( photon_chi,
                        particle_chi, 200, 1e-15 );

            if( numerator == 0 || denominator == 0 ) {
                xip = 0;
            } else {
                xip = numerator/( 2.*denominator );
            }

            if( xip < xip_threshold ) {
                logchipamin += pow( 0.1, k );
                k += 1;
            }
        }
        buffer[ichiph] = logchipamin;

        // display percentage
        if( 100.*ichiph >= length_table[rank]*pct ) {
            pct += dpct;
            MESSAGE( 3,ichiph + 1 << "/" << length_table[rank]
                     << " - " << ( int )( std::round( pct ) ) << "%" );
        }
    }

    // Communication of the xip_chipamin table
    MPI_Allgatherv( &buffer[0], length_table[rank], MPI_DOUBLE,
                    &xip_chipamin_table[0], &length_table[0], &imin_table[0],
                    MPI_DOUBLE, smpi->getGlobalComm() );

    // _______________________________________________________________
    // 2. - Computation of the xip table
    MESSAGE( 2,"Computation of xip:" );

    // Allocation of the local buffer
    buffer = new double [length_table[rank]*xip_chipa_dim];

    // Loop for xip in the photon_chi dimension
    dpct = std::max( 10., 100./( length_table[rank]*xip_chipa_dim ) );
    pct = dpct;
    for( int ichiph = 0 ; ichiph < length_table[rank] ; ichiph++ ) {

        // log10(photon_chi) for the current ichiph
        photon_chi = ( imin_table[rank] + ichiph )*xip_chiph_delta
                     + xip_log10_chiph_min;

        // photon_chi for the current ichiph
        photon_chi = pow( 10., photon_chi );

        // Delta on particle_chi
        chipa_delta = ( log10( 0.5*photon_chi ) - xip_chipamin_table[imin_table[rank] + ichiph] )
                      / ( xip_chipa_dim - 1 );

        // Denominator of xip
        denominator = MultiphotonBreitWheelerTables::compute_integration_Ritus_dTdchi( photon_chi,
                      0.5*photon_chi, 200, 1e-15 );

        // Loop in the particle_chi dimension
        for( int ichipa = 0 ; ichipa < xip_chipa_dim ; ichipa ++ ) {
            // Local particle_chi value
            particle_chi = pow( 10., ichipa*chipa_delta +
                                xip_chipamin_table[imin_table[rank] + ichiph] );

            // Numerator of xip
            numerator = MultiphotonBreitWheelerTables::compute_integration_Ritus_dTdchi( photon_chi,
                        particle_chi, 200, 1e-15 );

            // Update local buffer value
            // buffer[ichiph*xip_chipa_dim + ichipa] = std::min(1.,numerator / 2*denominator);
            buffer[ichiph*xip_chipa_dim + ichipa] = numerator / ( 2*denominator );

            // If buffer == 1, end of the loop with 1
            /*if (buffer[ichiph*xip_chipa_dim + ichipa] == 1)
            {
                for (int i = ichipa+1 ; i < xip_chipa_dim ; i ++)
                {
                    buffer[ichiph*xip_chipa_dim + ichipa] = 1.;
                }
                ichipa = xip_chipa_dim;
            }*/

            // Artificial monotony
            /*if ((ichipa > 0) &&
                (buffer[ichiph*xip_chipa_dim + ichipa] <
                 buffer[ichiph*xip_chipa_dim + ichipa -1]))
            {
                buffer[ichiph*xip_chipa_dim + ichipa] =
                buffer[ichiph*xip_chipa_dim + ichipa -1];
            }*/

            // display percentage
            if( 100.*( ichiph*xip_chipa_dim+ichipa )
                >= length_table[rank]*xip_chipa_dim*pct ) {
                MESSAGE( 3,ichiph*xip_chipa_dim+ichipa
                         << "/"
                         << length_table[rank]*xip_chipa_dim
                         << " - " << ( int )( std::round( pct ) )
                         << "%" );
                pct += dpct;
            }

        }
    }

    // Update length_table and imin_table
    for( int i = 0 ; i < nb_ranks ; i++ ) {
        length_table[i] *= xip_chipa_dim;
        imin_table[i] *= xip_chipa_dim;
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

    t1 = MPI_Wtime();
    MESSAGE( 2,"Done in " << ( t1 - t0 ) << "s" );

}

// -----------------------------------------------------------------------------
//! Output the computed tables so that thay can be read at the next run.
//
//! \param params list of simulation parameters
//! \param smpi MPI parameters
// -----------------------------------------------------------------------------
void MultiphotonBreitWheelerTables::computeTables( Params &params, SmileiMPI *smpi )
{
    // These tables are loaded only if if one species has Monte-Carlo Compton radiation
    if( params.hasMultiphotonBreitWheeler ) {
        MultiphotonBreitWheelerTables::computeTtable( smpi );
        MultiphotonBreitWheelerTables::computeXipTable( smpi );
    }
}

// -----------------------------------------------------------------------------
// TABLE OUTPUTS
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//! Ouput in a file of the table values of T for the
//! mutliphoton Breit-Wheeler process
// -----------------------------------------------------------------------------
void MultiphotonBreitWheelerTables::outputTableT()
{

    if( output_format_ == "ascii" ) {
        std::ofstream file;
        file.open( table_path_ + "/tab_T.dat" );

        if( file.is_open() ) {

            file.precision( 12 );

            file << "Multiphoton Breit-Wheeler T table\n";

            file << "Dimension photon_chi - photon_chi min - photon_chi max \n";

            file << T_dim;
            file << " "
                 << T_chiph_min << " "
                 << T_chiph_max << "\n";;

            // Loop over the table values
            for( int i = 0 ; i < T_dim ; i++ ) {
                file <<  T_table[i] << "\n";
            }

            file.close();
        }
    } else if( output_format_ == "binary" ) {
        std::ofstream file;
        file.open( table_path_ + "/tab_T.bin", std::ios::binary );

        if( file.is_open() ) {

            file.write( ( char * )&T_dim, sizeof( T_dim ) );
            file.write( ( char * )&T_chiph_min, sizeof( double ) );
            file.write( ( char * )&T_chiph_max, sizeof( double ) );

            // Loop over the table values
            for( int i = 0 ; i < T_dim ; i++ ) {
                file.write( ( char * )&T_table[i], sizeof( double ) );
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

        buffer = table_path_ + "/multiphoton_Breit_Wheeler_tables.h5";

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
        dims = T_dim;
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
                  &T_table[0] );

        // Create attributes
        H5::attr( datasetId, "chiph_min", T_chiph_min );
        H5::attr( datasetId, "chiph_max", T_chiph_max );
        H5::attr( datasetId, "chiph_dim", T_dim );

        // Close everything
        H5Dclose( datasetId );
        H5Sclose( dataspaceId );
        H5Fclose( fileId );

    } else {
        ERROR( "The table output format " << output_format_
                 << " is not recognized" );
    }
}

// -----------------------------------------------------------------------------
//! File output of xip_chipamin_table and xip_table
//
// -----------------------------------------------------------------------------
void MultiphotonBreitWheelerTables::outputXipTable()
{

    if( output_format_ == "ascii" ) {
        std::ofstream file;
        file.open( table_path_ + "/tab_mBW_xip.dat" );

        if( file.is_open() ) {

            file.precision( 12 );

            file << "Table xip_chipamin and xip for multiphoton Breit-Wheeler \n";

            file << "Dimension photon_chi - Dimension particle_chi - photon_chi min - photon_chi max \n";

            file << xip_chiph_dim << " "
                 << xip_chipa_dim << " "
                 << xip_chiph_min << " "
                 << xip_chiph_max << "\n";;

            // Loop over the xip values
            for( int ichiph = 0 ; ichiph < xip_chiph_dim ; ichiph++ ) {
                for( int ichipa = 0 ; ichipa < xip_chipa_dim ; ichipa++ ) {
                    file <<  xip_table[ichiph*xip_chipa_dim+ichipa] << " ";
                }
                file << "\n";
            }

            file.close();
        }
    } else if( output_format_ == "binary" ) {
        std::ofstream file;
        file.open( table_path_ + "/tab_mBW_xip.bin", std::ios::binary );

        if( file.is_open() ) {

            double temp0, temp1;

            temp0 = xip_chiph_min;
            temp1 = xip_chiph_max;

            file.write( ( char * )&xip_chiph_dim, sizeof( int ) );
            file.write( ( char * )&xip_chipa_dim, sizeof( int ) );
            file.write( ( char * )&temp0, sizeof( double ) );
            file.write( ( char * )&temp1, sizeof( double ) );

            // Write all table values of xip_chimin_table
            file.write( ( char * )&xip_chipamin_table[0], sizeof( double )*xip_chiph_dim );

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

        buffer = table_path_ + "/multiphoton_Breit_Wheeler_tables.h5";

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
        dims[0] = xip_chiph_dim;
        dataspaceId = H5Screate_simple( 1, dims, NULL );
        datasetId = H5Dcreate( fileId,
                               "xip_chipamin",
                               H5T_NATIVE_DOUBLE,
                               dataspaceId,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

        // Fill the dataset chiphmin
        H5Dwrite( datasetId, H5T_NATIVE_DOUBLE,
                  H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  &xip_chipamin_table[0] );

        // Attribute creation
        H5::attr( datasetId, "chiph_min", xip_chiph_min );
        H5::attr( datasetId, "chiph_max", xip_chiph_max );
        H5::attr( datasetId, "chiph_dim", xip_chiph_dim );
        H5::attr( datasetId, "power", xip_power );
        H5::attr( datasetId, "threshold", xip_threshold );

        // Creation of the datasat chiphmin xip
        dims[0] = xip_chiph_dim;
        dims[1] = xip_chipa_dim;
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
        H5::attr( datasetId, "chiph_min", xip_chiph_min );
        H5::attr( datasetId, "chiph_max", xip_chiph_max );
        H5::attr( datasetId, "chiph_dim", xip_chiph_dim );
        H5::attr( datasetId, "chipa_dim", xip_chipa_dim );
        H5::attr( datasetId, "power", xip_power );
        H5::attr( datasetId, "threshold", xip_threshold );

        // Close everything
        H5Dclose( datasetId );
        H5Sclose( dataspaceId );
        H5Fclose( fileId );
    } else {
        ERROR( "The output format " << output_format_ << " is not recognized" );
    }
}

// -----------------------------------------------------------------------------
//! Output the computed tables so that thay can be read at the next run.
//! Table output by the master MPI rank
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
void MultiphotonBreitWheelerTables::outputTables( SmileiMPI *smpi )
{
    // Sequential output
    if( smpi->isMaster() ) {
        // If tables have been computed, they are output on the disk
        // to be used for the next run
        if( T_computed ) {
            outputTableT();
        }
        if( xip_computed ) {
            outputXipTable();
        }
    }
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

    // Test if an external table exists, if yes we read the table...
    // Binary table
    if( Tools::file_exists( table_path_ + "/tab_T.bin" ) ) {

        if( smpi->getRank()==0 ) {

            // Reading of the table file
            std::ifstream file;
            file.open( table_path_ + "/tab_T.bin", std::ios::binary );

            if( file.is_open() ) {

                // Read the header
                file.read( ( char * )&T_dim,
                           sizeof( T_dim ) );
                file.read( ( char * )&T_chiph_min,
                           sizeof( T_chiph_min ) );
                file.read( ( char * )&T_chiph_max,
                           sizeof( T_chiph_max ) );

                // Resize of the array integfochi_table before reading
                T_table.resize( T_dim );

                // Read the table values
                file.read( ( char * )&T_table[0], sizeof( double )*T_dim );

                file.close();
            }

        }

    }
    // HDF5 format
    else if( Tools::file_exists( table_path_ + "/multiphoton_Breit_Wheeler_tables.h5" ) ) {
        hid_t       fileId;
        hid_t       datasetId;
        std::string buffer;

        if( smpi->getRank()==0 ) {

            buffer = table_path_ + "/multiphoton_Breit_Wheeler_tables.h5";

            fileId = H5Fopen( buffer.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );

            datasetId = H5Dopen2( fileId, "h", H5P_DEFAULT );

            // If this dataset exists, we read it
            if( datasetId > 0 ) {

                // First, we read attributes
                H5::getAttr( datasetId, "chiph_dim", T_dim );
                H5::getAttr( datasetId, "chiph_min", T_chiph_min );
                H5::getAttr( datasetId, "chiph_max", T_chiph_max );

                // Resize of the array integfochi_table before reading
                T_table.resize( T_dim );

                // then the dataset
                H5Dread( datasetId,
                         H5T_NATIVE_DOUBLE, H5S_ALL,
                         H5S_ALL, H5P_DEFAULT,
                         &T_table[0] );

                H5Dclose( datasetId );
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

    MESSAGE( 1,"--- Table T:" );
    MESSAGE( 2,"Reading of the external database" );
    MESSAGE( 2,"Dimension quantum parameter: "
             << T_dim );
    MESSAGE( 2,"Minimum photon quantum parameter chi: "
             << T_chiph_min );
    MESSAGE( 2,"Maximum photon quantum parameter chi: "
             << T_chiph_max );

    // Bcast the table to all MPI ranks
    MultiphotonBreitWheelerTables::bcastTableT( smpi );

}

// -----------------------------------------------------------------------------
//! Read the external table xip_chipamin and xip
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
void MultiphotonBreitWheelerTables::readTableXip( SmileiMPI *smpi )
{

    // Test if an external table exists, we read the table...
    if( Tools::file_exists( table_path_ + "/tab_mBW_xip.bin" ) ) {

        if( smpi->getRank()==0 ) {

            // Reading of the table file
            std::ifstream file;
            file.open( table_path_ + "/tab_mBW_xip.bin", std::ios::binary );

            if( file.is_open() ) {

                // Read the header
                file.read( ( char * )&xip_chiph_dim, sizeof( xip_chiph_dim ) );
                file.read( ( char * )&xip_chipa_dim, sizeof( xip_chipa_dim ) );
                file.read( ( char * )&xip_chiph_min, sizeof( xip_chiph_min ) );
                file.read( ( char * )&xip_chiph_max, sizeof( xip_chiph_max ) );

                // Allocation of the array xip
                xip_chipamin_table.resize( xip_chiph_dim );
                xip_table.resize( xip_chiph_dim*xip_chipa_dim );

                // Read the table values
                file.read( ( char * )&xip_chipamin_table[0], sizeof( double )*xip_chiph_dim );

                // Read the table values
                file.read( ( char * )&xip_table[0], sizeof( double )*xip_chipa_dim*xip_chiph_dim );

                file.close();
            }

        }
    }
    // HDF5 format
    else if( Tools::file_exists( table_path_ + "/multiphoton_Breit_Wheeler_tables.h5" ) ) {
        if( smpi->getRank()==0 ) {

            hid_t       fileId;
            hid_t       datasetId_chipamin;
            hid_t       datasetId_xip;
            std::string buffer;

            buffer = table_path_ + "/multiphoton_Breit_Wheeler_tables.h5";

            fileId = H5Fopen( buffer.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );

            datasetId_chipamin = H5Dopen2( fileId, "xip_chipamin", H5P_DEFAULT );
            datasetId_xip = H5Dopen2( fileId, "xip", H5P_DEFAULT );

            // If this dataset exists, we read it
            if( datasetId_chipamin > 0 && datasetId_xip > 0 ) {

                // First, we read attributes
                H5::getAttr( datasetId_xip, "chiph_dim", xip_chiph_dim );
                H5::getAttr( datasetId_xip, "chipa_dim", xip_chipa_dim );
                H5::getAttr( datasetId_xip, "chiph_min", xip_chiph_min );
                H5::getAttr( datasetId_xip, "chiph_max", xip_chiph_max );

                // Allocation of the array xip
                xip_chipamin_table.resize( xip_chiph_dim );
                xip_table.resize( xip_chipa_dim*xip_chiph_dim );

                // then the dataset for chipamin
                H5Dread( datasetId_chipamin,
                         H5T_NATIVE_DOUBLE, H5S_ALL,
                         H5S_ALL, H5P_DEFAULT,
                         &xip_chipamin_table[0] );

                // then the dataset for xip
                H5Dread( datasetId_xip,
                         H5T_NATIVE_DOUBLE, H5S_ALL,
                         H5S_ALL, H5P_DEFAULT,
                         &xip_table[0] );

                H5Dclose( datasetId_xip );
                H5Dclose( datasetId_chipamin );
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

    MESSAGE( 1,"--- Table chipamin and xip:" );
    MESSAGE( 2,"Reading of the external database" );
    MESSAGE( 2,"Dimension photon chi: " << xip_chiph_dim );
    MESSAGE( 2,"Dimension particle chi: " << xip_chipa_dim );
    MESSAGE( 2,"Minimum photon chi: " << xip_chiph_min );
    MESSAGE( 2,"Maximum photon chi: " << xip_chiph_max );

    // Bcast the table to all MPI ranks
    MultiphotonBreitWheelerTables::bcastTableXip( smpi );

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
        MultiphotonBreitWheelerTables::readTableXip( smpi );
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
        MPI_Pack_size( T_dim, MPI_DOUBLE, smpi->getGlobalComm(),
                       &position );
        buf_size += position;
    }

    MESSAGE( 2,"Buffer size: " << buf_size );

    // Exchange buf_size with all ranks
    MPI_Bcast( &buf_size, 1, MPI_INT, 0, smpi->getGlobalComm() );

    // Packet that will contain all parameters
    char *buffer = new char[buf_size];

    // Proc 0 packs
    if( smpi->getRank() == 0 ) {
        position = 0;
        MPI_Pack( &T_dim,
                  1, MPI_INT, buffer, buf_size, &position, smpi->getGlobalComm() );
        MPI_Pack( &T_chiph_min,
                  1, MPI_DOUBLE, buffer, buf_size, &position, smpi->getGlobalComm() );
        MPI_Pack( &T_chiph_max,
                  1, MPI_DOUBLE, buffer, buf_size, &position, smpi->getGlobalComm() );

        MPI_Pack( &T_table[0], T_dim,
                  MPI_DOUBLE, buffer, buf_size, &position, smpi->getGlobalComm() );

    }

    // Bcast all parameters
    MPI_Bcast( &buffer[0], buf_size, MPI_PACKED, 0, smpi->getGlobalComm() );

    // Other ranks unpack
    if( smpi->getRank() != 0 ) {
        position = 0;
        MPI_Unpack( buffer, buf_size, &position,
                    &T_dim, 1, MPI_INT, smpi->getGlobalComm() );
        MPI_Unpack( buffer, buf_size, &position,
                    &T_chiph_min, 1, MPI_DOUBLE, smpi->getGlobalComm() );
        MPI_Unpack( buffer, buf_size, &position,
                    &T_chiph_max, 1, MPI_DOUBLE, smpi->getGlobalComm() );

        // Resize table before unpacking values
        T_table.resize( T_dim );

        MPI_Unpack( buffer, buf_size, &position, &T_table[0],
                    T_dim, MPI_DOUBLE, smpi->getGlobalComm() );

    }

    delete[] buffer;

    T_log10_chiph_min = log10( T_chiph_min );

    // Computation of the delta
    T_chiph_delta = ( log10( T_chiph_max )
                      - T_log10_chiph_min )/( T_dim-1 );

    // Inverse delta
    T_chiph_inv_delta = 1./T_chiph_delta;
}

// -----------------------------------------------------------------------------
//! Bcast of the external table xip_chipamin and xip
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
void MultiphotonBreitWheelerTables::bcastTableXip( SmileiMPI *smpi )
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
        MPI_Pack_size( xip_chiph_dim, MPI_DOUBLE, smpi->getGlobalComm(),
                       &position );
        buf_size += position;
        MPI_Pack_size( xip_chiph_dim*xip_chipa_dim, MPI_DOUBLE,
                       smpi->getGlobalComm(), &position );
        buf_size += position;
    }

    MESSAGE( 2,"Buffer size for MPI exchange: " << buf_size );

    // Exchange buf_size with all ranks
    MPI_Bcast( &buf_size, 1, MPI_INT, 0, smpi->getGlobalComm() );

    // Packet that will contain all parameters
    char *buffer = new char[buf_size];

    // Proc 0 packs
    if( smpi->getRank() == 0 ) {
        position = 0;
        MPI_Pack( &xip_chiph_dim,
                  1, MPI_INT, buffer, buf_size, &position, smpi->getGlobalComm() );
        MPI_Pack( &xip_chipa_dim,
                  1, MPI_INT, buffer, buf_size, &position, smpi->getGlobalComm() );
        MPI_Pack( &xip_chiph_min,
                  1, MPI_DOUBLE, buffer, buf_size, &position, smpi->getGlobalComm() );
        MPI_Pack( &xip_chiph_max,
                  1, MPI_DOUBLE, buffer, buf_size, &position, smpi->getGlobalComm() );

        MPI_Pack( &xip_chipamin_table[0], xip_chiph_dim,
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
                    &xip_chiph_dim, 1, MPI_INT, smpi->getGlobalComm() );
        MPI_Unpack( buffer, buf_size, &position,
                    &xip_chipa_dim, 1, MPI_INT, smpi->getGlobalComm() );
        MPI_Unpack( buffer, buf_size, &position,
                    &xip_chiph_min, 1, MPI_DOUBLE, smpi->getGlobalComm() );
        MPI_Unpack( buffer, buf_size, &position,
                    &xip_chiph_max, 1, MPI_DOUBLE, smpi->getGlobalComm() );

        // Resize tables before unpacking values
        xip_chipamin_table.resize( xip_chiph_dim );
        xip_table.resize( xip_chipa_dim*xip_chiph_dim );

        MPI_Unpack( buffer, buf_size, &position, &xip_chipamin_table[0],
                    xip_chiph_dim, MPI_DOUBLE, smpi->getGlobalComm() );

        MPI_Unpack( buffer, buf_size, &position, &xip_table[0],
                    xip_chipa_dim*xip_chiph_dim, MPI_DOUBLE, smpi->getGlobalComm() );
    }

    delete[] buffer;

    // Log10 of xip_chiph_min for efficiency
    xip_log10_chiph_min = log10( xip_chiph_min );

    // Computation of the delta
    xip_chiph_delta = ( log10( xip_chiph_max )
                        - xip_log10_chiph_min )/( xip_chiph_dim-1 );

    // Inverse of delta
    xip_chiph_inv_delta = 1./xip_chiph_delta;

    // Inverse particle_chi discetization (regularly used)
    xip_inv_chipa_dim_minus_one = 1./( xip_chipa_dim - 1. );

}
