// ----------------------------------------------------------------------------
//! \file NonlinearInverseComptonScattering.h
//
//! \brief This class contains methods to generate the nonlinear inverse Compton scattering tables
//
// ----------------------------------------------------------------------------

#include "NonlinearInverseComptonScattering.h"

//! Creation of the tables
void NonlinearComptonScattering::createTables(int argc, std::string * arguments)
{
    
    // _______________________________________________________________________
    // MPI
    
    int rank;
    int number_of_ranks;
    
    MPI_Comm_size( MPI_COMM_WORLD, &number_of_ranks );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    if( rank==0 ) {
        std::cout << "\n You have selected the creation of tables"
                  << " for the nonlinear inverse Compton scattering." << std::endl;
    }
    
    // _______________________________________________________________________
    // Parameters
    
    double t0;
    double t1;
    
    double min_particle_chi;
    double max_particle_chi;
    
    int size_particle_chi;
    int size_photon_chi;
    
    double particle_chi;
    double photon_chi;
    
    double log10_photon_chi;
    
    double log10_min_particle_chi;
    double log10_max_particle_chi;
    
    double delta_particle_chi;
    double delta_photon_chi;
    double inverse_delta_particle_chi;
    
    double xi_power;
    double xi_threshold;
    
    double xi;
    
    int i_particle_chi;
    int i_photon_chi;
    
    int number_of_draws;
    
    // Vector for 1d tables
    std::vector <double> table_1d;
    std::vector <double> table_2d;
    
    bool verbose;
    
    // Parameter default initialization
    size_particle_chi      = 128;
    size_photon_chi        = 128;
    min_particle_chi       = 1e-3;
    max_particle_chi       = 1e3;
    xi_power               = 4;
    xi_threshold           = 1e-3;
    number_of_draws        = 0;
    verbose                = false;
    
    std::string help_message;
    help_message =  "\n Help page specific to the nonlinear inverse Compton Scattering:\n";
    help_message += "\n";
    help_message += " List of available commands:\n";
    help_message += " -h, --help                       print a help message and exit.\n";
    help_message += " -s, --size       int int         respective size of the particle and photon chi axis. (default 128 128)\n";
    help_message += " -b, --boundaries double double   min and max of the particle chi axis. (default 1e-3 1e3)\n";
    help_message += " -e, --error      int             compute error due to discretization and use the provided int as a number of draws. (default 0)\n";
    help_message += " -t, --threshold  double          Minimum targeted value of xi in the computation the minimum particle quantum parameter. (default 1e-3)\n";
    help_message += " -p, --power      int             Maximum decrease in order of magnitude for the search for the minimum particle quantum parameter. (default 4)\n";
    help_message += " -v, --verbose                    Dump the tables\n";
    
    // Read from command line
    int i_arg = 2;
    while(i_arg < argc) {
        if (arguments[i_arg] == "-s" || arguments[i_arg] == "--size") {
            size_particle_chi = std::stoi(arguments[i_arg+1]);
            size_photon_chi = std::stoi(arguments[i_arg+2]);
            i_arg+=3;
        } else if (arguments[i_arg] == "-b" || arguments[i_arg] == "--boundaries") {
            min_particle_chi = std::stod(arguments[i_arg+1]);
            max_particle_chi = std::stod(arguments[i_arg+2]);
            i_arg+=3;
        } else if (arguments[i_arg] == "-e" || arguments[i_arg] == "--error") {
            number_of_draws = std::stoi(arguments[i_arg+1]);
            i_arg+=2;
        } else if (arguments[i_arg] == "-t" || arguments[i_arg] == "--threshold") {
            xi_threshold = std::stod(arguments[i_arg+1]);
            i_arg+=2;
        } else if (arguments[i_arg] == "-p" || arguments[i_arg] == "--power") {
            xi_power = std::stod(arguments[i_arg+1]);
            i_arg+=2;
        } else if (arguments[i_arg] == "-v" || arguments[i_arg] == "--verbose") {
            verbose = true;
            i_arg+=1;
        } else if (arguments[i_arg] == "-h" || arguments[i_arg] == "--help") {
            if (rank == 0) {
                std::cout << help_message << std::endl;
            }
            //MPI_Abort(MPI_COMM_WORLD,error);
            exit(0);
            i_arg+=2;
        } else {
            ERROR("Keywork " << arguments[i_arg] << " not recognized");
        }
    }
    
    log10_min_particle_chi = std::log10(min_particle_chi);
    log10_max_particle_chi = std::log10(max_particle_chi);
    
    delta_particle_chi = (log10_max_particle_chi - log10_min_particle_chi) / (size_particle_chi - 1);
    
    inverse_delta_particle_chi = 1.0 / delta_particle_chi;
    
    table_1d.resize(size_particle_chi);
    table_2d.resize(size_photon_chi*size_particle_chi);
    
    if( rank==0 ) {
        std::cout << " Size particle chi axis: " << size_particle_chi << "\n"
                  << " Size photon chi axis: " << size_photon_chi << "\n"
                  << " Min particle chi axis: " << min_particle_chi << "\n"
                  << " Max particle chi axis: " << max_particle_chi << "\n"
                  << " Power: " << xi_power << "\n"
                  << " Threshold: " << xi_threshold
                  << std::endl;
    }
    
    // _______________________________________________________________________
    // Load repartition
    
    // Load distribution parameters
    int *rank_first_index = new int[number_of_ranks];
    int *rank_indexes = new int[number_of_ranks];
    
    Tools::distributeArray( number_of_ranks,
        size_particle_chi,
        rank_first_index,
        rank_indexes );
        
    // Print repartition
    if( rank==0 ) {
        std::cout << std::endl;
        std::cout << " MPI load distribution:" << std::endl;
        for( int i =0 ; i < number_of_ranks ; i++ ) {
            std::cout << " - Rank " << i
                      << " - 1st index: " << rank_first_index[i]
                      << " - length: " << rank_indexes[i]
                      << std::endl;
        }
    }
    
    // _______________________________________________________________________
    // First management of the output
    
    hid_t       fileId;
    std::string path;
    
    path = "./radiation_tables.h5";

    if (rank == 0) {
        remove(path.c_str());
    }
    
    // _______________________________________________________________________
    // Computation of integrale of F/chi between 0 and chi_particles
    
    // Allocation of the local buffer
    double * buffer = new double [rank_indexes[rank]];
    
    if( rank==0 ) {
        std::cout << std::endl;
        std::cout << " Computation of integrale of F/chi between 0 and chi_particles"
                  << std::endl;
    }
    
    double delta_percentage = std::max( delta_percentage, 100.0/rank_indexes[rank] );
    double percentage = 0;
    
    t0 = MPI_Wtime();

    // Computation of integrale of F/chi between 0 and chi_particles
    for( i_particle_chi = 0 ; i_particle_chi < rank_indexes[rank] ; i_particle_chi++ ) {
        particle_chi = std::pow( 10., ( rank_first_index[rank] + i_particle_chi )* delta_particle_chi
                             + log10_min_particle_chi );

        buffer[i_particle_chi] = NonlinearComptonScattering::integrateSynchrotronEmissivity( particle_chi,
                    1e-40*particle_chi, particle_chi, 400, 1e-15 );

        t1 = MPI_Wtime();

        if( rank==0 ) {
            if( 100.0*(i_particle_chi+1) >= rank_indexes[rank]*percentage ) {
                percentage += delta_percentage;
                std::cout << " - " << i_particle_chi + 1<< "/" << rank_indexes[rank]
                                    << " - " << ( int )( std::round( percentage ) ) << "%"
                                    << " - " << t1 - t0 << " s"
                                    << std::endl;
            }
        }
    }
    
    // All data gathered
    MPI_Allgatherv( &buffer[0], rank_indexes[rank], MPI_DOUBLE,
                    &table_1d[0], &rank_indexes[0], &rank_first_index[0],
                    MPI_DOUBLE, MPI_COMM_WORLD );

    t1 = MPI_Wtime();
    if (rank==0) {
        std::cout << " Total time: " << t1 - t0 << " s" << std::endl;
    }
    
    // Error evaluation
    if (number_of_draws > 0) {
        std::default_random_engine generator;
        std::uniform_real_distribution<double> distribution(min_particle_chi,max_particle_chi);
        double value;
        double interpolated_value;
        double error;
        double local_max_error = 0;
        double max_error = 0;
        double distance;
        
        if (rank==0) std::cout << " Error computation: " << std::endl;
        
        for(int i = 0 ; i < number_of_draws; i++) {
            particle_chi = distribution(generator);
            value = NonlinearComptonScattering::integrateSynchrotronEmissivity( particle_chi,
                        1e-40*particle_chi, particle_chi, 400, 1e-15 );
            i_particle_chi = int((log10(particle_chi) - log10_min_particle_chi) * inverse_delta_particle_chi);
            distance = std::abs(log10(particle_chi) - (i_particle_chi*delta_particle_chi + log10_min_particle_chi)) * inverse_delta_particle_chi;
            interpolated_value = table_1d[i_particle_chi]*(1 - distance) + table_1d[i_particle_chi+1]*distance;
            error = std::abs(value - interpolated_value)/value;
            local_max_error = std::max(local_max_error,error);
        }
    
        MPI_Reduce(&local_max_error,&max_error,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
        
        if (rank==0) {
            std::cout << " - Maximal relative error: " << max_error << std::endl;
        }
    }
    
    // Output of the tables
    if (verbose && rank == 0) {
        std:: cout << "\n table integfochi: " << std::setprecision(std::numeric_limits<double>::digits10 + 1) <<std::endl;
        for( i_particle_chi = 0 ; i_particle_chi < size_particle_chi  ; i_particle_chi += 8 ) {
            std::cout << " ";
            for (int i = i_particle_chi ; i< std::min(i_particle_chi+8,size_particle_chi) ; i++) {
                std::cout << table_1d[i] << ", ";
            }
            std::cout << std::endl;
        }
    }
    
    // Output of the vector
    if (rank==0) {
        fileId  = H5Fcreate( path.c_str(),
                             H5F_ACC_TRUNC,
                             H5P_DEFAULT,
                             H5P_DEFAULT );
        
        std::string vect_name("integfochi");
        
        H5::vect( fileId, vect_name, table_1d, 0 );
        
        std::string attr_name("min_particle_chi");
        H5::attr( fileId, vect_name, attr_name, min_particle_chi);
        
        attr_name = "max_particle_chi";
        H5::attr( fileId, vect_name, attr_name, max_particle_chi);
        
        attr_name = "size_particle_chi";
        H5::attr( fileId, vect_name, attr_name, size_particle_chi);

        H5Fclose( fileId );
    }
    
    // _______________________________________________________________________
    // Computation of the table H for the model of Niel et al.
    
    if( rank==0 ) {
        std::cout << std::endl;
        std::cout << " Computation of the table H for Niel at al."
                  << std::endl;
    }
    
    percentage = 0;
    t0 = MPI_Wtime();
    
    for( i_particle_chi = 0 ; i_particle_chi < rank_indexes[rank] ; i_particle_chi ++ ) {
        
        particle_chi = std::pow( 10.0, ( rank_first_index[rank] + i_particle_chi )*delta_particle_chi
                            + log10_min_particle_chi );

        // 100 iterations is in theory sufficient to get the convergence
        // at 1e-15
        buffer[i_particle_chi] = NonlinearComptonScattering::computeHNiel( particle_chi, 400, 1e-15 );

        t1 = MPI_Wtime();
        
        if( rank==0 ) {
            if( 100.0*(i_particle_chi+1) >= rank_indexes[rank]*percentage ) {
                percentage += delta_percentage;
                std::cout << " - " << i_particle_chi + 1<< "/" << rank_indexes[rank]
                                  << " - " << ( int )( std::round( percentage ) ) << "%"
                                  << " - " << t1 - t0 << " s"
                                  << std::endl;
            }
        }
    }
    
    t1 = MPI_Wtime();
    if (rank==0) {
        std::cout << " Total time: " << t1 - t0 << " s" << std::endl;
    }

    // All data gathered
    MPI_Allgatherv( &buffer[0], rank_indexes[rank], MPI_DOUBLE,
                    &table_1d[0], &rank_indexes[0], &rank_first_index[0],
                    MPI_DOUBLE, MPI_COMM_WORLD );
    
    // Output of the tables
    if (verbose && rank == 0) {
        std:: cout << "\n table h for Niel: " << std::setprecision(std::numeric_limits<double>::digits10 + 1) <<std::endl;
        for( i_particle_chi = 0 ; i_particle_chi < size_particle_chi  ; i_particle_chi += 8 ) {
            std::cout << " ";
            for (int i = i_particle_chi ; i< std::min(i_particle_chi+8,size_particle_chi) ; i++) {
                std::cout << table_1d[i] << ", ";
            }
            std::cout << std::endl;
        }
    }
    
    // Output of the vector
    if (rank==0) {
        fileId = H5Fopen( path.c_str(),
                          H5F_ACC_RDWR,
                          H5P_DEFAULT );
        
        std::string vect_name("h");
        
        H5::vect( fileId, vect_name, table_1d, 0 );
        
        std::string attr_name("min_particle_chi");
        H5::attr( fileId, vect_name, attr_name, min_particle_chi);
        
        attr_name = "max_particle_chi";
        H5::attr( fileId, vect_name, attr_name, max_particle_chi);
        
        attr_name = "size_particle_chi";
        H5::attr( fileId, vect_name, attr_name, size_particle_chi);

        H5Fclose( fileId );
    }

    // _______________________________________________________________________
    // Computation of the minimal photon chi value for xi
    
    if( rank==0 ) {
        std::cout << std::endl;
        std::cout << " Computation of the minimum photon chi value for the xi table"
                  << std::endl;
    }
    
    double denominator;
    double numerator;
    int k;
    
    percentage = 0;
    t0 = MPI_Wtime();
        
    for( int i_particle_chi = 0 ; i_particle_chi < rank_indexes[rank] ; i_particle_chi++ ) {

        xi = 1;
        log10_photon_chi = ( rank_first_index[rank] + i_particle_chi )*delta_particle_chi
                      + log10_min_particle_chi;
        particle_chi = pow( 10.0, log10_photon_chi );
        
        // Denominator of xi
        denominator = NonlinearComptonScattering::integrateSynchrotronEmissivity( particle_chi,
                      0.99e-40*particle_chi, particle_chi, 200, 1e-13 );
        
        k = 0;
        while( k < xi_power ) {
            log10_photon_chi -= pow( 0.1, k );
            photon_chi = pow( 10.0, log10_photon_chi );
            numerator = NonlinearComptonScattering::integrateSynchrotronEmissivity( particle_chi,
                      0.99e-40*photon_chi, photon_chi, 200, 1e-13 );

            if( numerator == 0 ) {
                xi = 0;
            } else {
                xi = numerator/denominator;
            }

            if( xi < xi_threshold ) {
                log10_photon_chi += pow( 0.1, k );
                k += 1;
            }
        }
        buffer[i_particle_chi] = log10_photon_chi;
        
        t1 = MPI_Wtime();
        
        if( rank==0 ) {
            if( 100.0*(i_particle_chi+1) >= rank_indexes[rank]*percentage ) {
                percentage += delta_percentage;
                std::cout << " - " << i_particle_chi + 1<< "/" << rank_indexes[rank]
                                  << " - " << ( int )( std::round( percentage ) ) << "%"
                                  << " - " << t1 - t0 << " s"
                                  << std::endl;
            }
        }
        
    }

    t1 = MPI_Wtime();
    if (rank==0) {
        std::cout << " Total time: " << t1 - t0 << " s" << std::endl;
    }

    // All data gathered
    MPI_Allgatherv( &buffer[0], rank_indexes[rank], MPI_DOUBLE,
                    &table_1d[0], &rank_indexes[0], &rank_first_index[0],
                    MPI_DOUBLE, MPI_COMM_WORLD );

    // Output of the tables
    if (verbose && rank == 0) {
        std:: cout << "\n table min_photon_chi_for_xi: " << std::setprecision(std::numeric_limits<double>::digits10 + 1) <<std::endl;
        for( i_particle_chi = 0 ; i_particle_chi < size_particle_chi  ; i_particle_chi += 8 ) {
            std::cout << " ";
            for (int i = i_particle_chi ; i< std::min(i_particle_chi+8,size_particle_chi) ; i++) {
                std::cout << table_1d[i] << ", ";
            }
            std::cout << std::endl;
        }
    }

    // Output
    if (rank==0) {
        
        fileId = H5Fopen( path.c_str(),
                          H5F_ACC_RDWR,
                          H5P_DEFAULT );
        
        std::string vect_name("min_photon_chi_for_xi");
        
        H5::vect( fileId, vect_name, table_1d, 0 );
        
        std::string attr_name("min_particle_chi");
        H5::attr( fileId, vect_name, attr_name, min_particle_chi);
        
        attr_name = "max_particle_chi";
        H5::attr( fileId, vect_name, attr_name, max_particle_chi);
        
        attr_name = "size_particle_chi";
        H5::attr( fileId, vect_name, attr_name, size_particle_chi);

        attr_name = "power";
        H5::attr( fileId, vect_name, attr_name, xi_power);

        attr_name = "threshold";
        H5::attr( fileId, vect_name, attr_name, xi_threshold);

        H5Fclose( fileId );
    }

    // _______________________________________________________________________
    // Computation of xi
    
    if( rank==0 ) {
        std::cout << std::endl;
        std::cout << " Computation of the xi table"
                  << std::endl;
    }
    
    // Allocation of the local buffer
    buffer = new double [rank_indexes[rank]*size_photon_chi];

    percentage = 0;
    t0 = MPI_Wtime();

    for( int i_particle_chi = 0 ; i_particle_chi < rank_indexes[rank] ; i_particle_chi++ ) {
        
        particle_chi = ( rank_first_index[rank] + i_particle_chi )*delta_particle_chi
               + log10_min_particle_chi;

        delta_photon_chi = ( particle_chi - table_1d[rank_first_index[rank] + i_particle_chi] )
              / ( size_photon_chi - 1 );
        
        particle_chi = std::pow( 10.0, particle_chi );
        
        // Denominator of xi
        denominator = NonlinearComptonScattering::integrateSynchrotronEmissivity( particle_chi,
                    1e-40*particle_chi, particle_chi, 300, 1e-15 );
        
        bool unity_not_reached = true;
        
        // Loop over the photon chi axis
        for( i_photon_chi = 0 ; i_photon_chi < size_photon_chi ; i_photon_chi ++ ) {
            
            // Local photon_chi value
            photon_chi = std::pow( 10.0, i_photon_chi*delta_photon_chi +
                              table_1d[rank_first_index[rank] + i_particle_chi] );
                              
            // Numerator of xi
            numerator = NonlinearComptonScattering::integrateSynchrotronEmissivity( particle_chi,
                      1e-40*photon_chi, photon_chi, 300, 1e-15 );
                              
            // Update local buffer value
            buffer[i_particle_chi*size_photon_chi + i_photon_chi] = std::min( 1.0, numerator / denominator );
                      
            // Check monotony
            if ( i_photon_chi > 0 ) {
                        
                if ( buffer[i_particle_chi*size_photon_chi + i_photon_chi] <
                buffer[i_particle_chi*size_photon_chi + i_photon_chi -1] ) {
                    std::cout << "   > Monotony issue at " << i_particle_chi*size_photon_chi + i_photon_chi
                              << " " << buffer[i_particle_chi*size_photon_chi + i_photon_chi-1]
                              << " " << buffer[i_particle_chi*size_photon_chi + i_photon_chi]
                              << "." << std::endl;
                }
                        
                // buffer[i_particle_chi*size_photon_chi + i_photon_chi] =
                // buffer[i_particle_chi*size_photon_chi + i_photon_chi -1];
            }
            // First unity value
            if ( i_photon_chi < size_photon_chi - 1 ) {
              if ( buffer[i_particle_chi*size_photon_chi + i_photon_chi] >= 1 and unity_not_reached) {
                std::cerr << "   > Unity reached at " << i_photon_chi << std::endl;
                unity_not_reached = false;
              }
            }
                      
        }
        
        t1 = MPI_Wtime();
        
        if( rank==0 ) {
            if( 100.0*(i_particle_chi+1) >= rank_indexes[rank]*percentage ) {
                percentage += delta_percentage;
                std::cout << " - " << i_particle_chi + 1<< "/" << rank_indexes[rank]
                                  << " - " << ( int )( std::round( percentage ) ) << "%"
                                  << " - " << t1 - t0 << " s"
                                  << std::endl;
            }
        }
        
    }

    t1 = MPI_Wtime();
    if (rank==0) {
        std::cout << " Total time: " << t1 - t0 << " s" << std::endl;
    }

    // Update rank_indexes and rank_first_index
    for( int i = 0 ; i < number_of_ranks ; i++ ) {
        rank_indexes[i] *= size_photon_chi;
        rank_first_index[i] *= size_photon_chi;
    }

    // Communication of the xip table
    MPI_Allgatherv( &buffer[0], rank_indexes[rank], MPI_DOUBLE,
                    &table_2d[0], &rank_indexes[0], &rank_first_index[0],
                    MPI_DOUBLE, MPI_COMM_WORLD );

    // Output of the tables
    if (verbose && rank == 0) {
        //std::string message = "";
        std:: cout << "\n table xi: " << std::setprecision(std::numeric_limits<double>::digits10 + 1) <<std::endl;
        for( i_particle_chi = 0 ; i_particle_chi < size_particle_chi  ; i_particle_chi++ ) {
            std::cout << " ";
            for( i_photon_chi = 0 ; i_photon_chi < size_photon_chi ; i_photon_chi ++ ) {
                //message += std::to_string(table_2d[i_photon_chi*size_particle_chi + i_particle_chi]) + ", ";
                std::cout << table_2d[i_particle_chi*size_photon_chi + i_photon_chi] << ", ";
            }
            //message += "\n";
            std::cout << std::endl;
        }
    }

    // Output
    if (rank==0) {
        
        fileId = H5Fopen( path.c_str(),
                          H5F_ACC_RDWR,
                          H5P_DEFAULT );
        
        std::string vect_name("xi");
        
        int size[2];
        size[0] = size_particle_chi;
        size[1] = size_photon_chi;
        H5::H5Vector2D( fileId, vect_name, &size[0], table_2d);
        //H5::vect( fileId, vect_name, table_2d, 0 );
        
        std::string attr_name("min_particle_chi");
        H5::attr( fileId, vect_name, attr_name, min_particle_chi);
        
        attr_name = "max_particle_chi";
        H5::attr( fileId, vect_name, attr_name, max_particle_chi);
        
        attr_name = "size_particle_chi";
        H5::attr( fileId, vect_name, attr_name, size_particle_chi);

        attr_name = "size_photon_chi";
        H5::attr( fileId, vect_name, attr_name, size_photon_chi);

        H5Fclose( fileId );
    }

    // _______________________________________________________________________
    // Cleaning

    // Free memory
    delete buffer;
    delete rank_first_index;
    delete rank_indexes;
    
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
//! \param discretization number of points for integration
//! \param eps integration accuracy
// ---------------------------------------------------------------------------------------------------------------------
double NonlinearComptonScattering::integrateSynchrotronEmissivity( double particle_chi,
        double min_photon_chi,
        double max_photon_chi,
        int discretization,
        double eps )
{

    //std::cout << "RadiationTables::integrateSynchrotronEmissivity" << std::endl;

    // Arrays for Gauss-Legendre integration
    double *roots = new double[discretization];
    double *weights = new double[discretization];
    // Photon quantum parameter
    double photon_chi;
    // Integration result
    double integ;
    // Synchrotron emissivity
    double sync_emi;
    // Iterator
    int i;

    // gauss Legendre coefficients
    Tools::GaussLegendreQuadrature( log10( min_photon_chi ), log10( max_photon_chi ),
                                        roots, weights, discretization, eps );

    // Integration loop
    integ = 0;
    // #pragma omp parallel for reduction(+:integ) private(photon_chi,sync_emi) shared(particle_chi,roots,roots)
    for( i=0 ; i< discretization ; i++ ) {
        photon_chi = std::pow( 10., roots[i] );
        sync_emi = NonlinearComptonScattering::computeRitusSynchrotronEmissivity( particle_chi, photon_chi, 200, 1e-15 );
        integ += weights[i]*sync_emi*std::log( 10 );
    }

    delete[] roots;
    delete[] weights;

    return integ;

}

// ---------------------------------------------------------------------------------------------------------------------
//! Computation of the synchrotron emissivity following the formulae of Ritus
//
//! \param particle_chi particle quantum parameter
//! \param photon_chi photon quantum parameter
//! \param discretization number of iterations for the Gauss-Legendre integration
//! \param eps epsilon for the modified bessel function
// ---------------------------------------------------------------------------------------------------------------------
double NonlinearComptonScattering::computeRitusSynchrotronEmissivity( double particle_chi,
        double photon_chi, int discretization, double eps )
{

    //std::cout << "RadiationTables::computeRitusSynchrotronEmissivity" << std::endl;

    // The photon quantum parameter should be below the electron one
    if( particle_chi > photon_chi ) {
        // Arrays for Gauss-Legendre integration
        double *roots = new double[discretization];
        double *weights = new double[discretization];
        // Values for Bessel results
        double K;
        // Parts of the formulae
        double part1, part2;
        // Iterator
        int i;

        double y = photon_chi/( 3.0*particle_chi*( particle_chi-photon_chi ) );

        // Computation of Part. 1
        // Call the modified Bessel function to get K
        K = Tools::BesselK( 2.0/3.0, 2*y);

        part1 = ( 2.0 + 3.0*photon_chi*y )*( K );

        // Computation of Part. 2
        // Using Gauss Legendre integration

        Tools::GaussLegendreQuadrature( std::log10( 2*y ), std::log10( y )+5, roots,
                                            weights, discretization, eps );

        part2 = 0;
        for( i=0 ; i< discretization ; i++ ) {
            y = std::pow( 10., roots[i] );
            K = Tools::BesselK( 1./3., y);
            part2 += weights[i]*K*y*std::log( 10. );
        }

        // Factor for final result
        y = 2*photon_chi/( 3*std::pow( particle_chi, 2. ) );

        delete[] roots;
        delete[] weights;

        return ( part1 - part2 )*y;


    } else if( particle_chi == photon_chi ) {
        return 0;
    } else {
        std::cerr << "In computeRitusSynchrotronEmissivity: particle_chi " << particle_chi
               << " < photon_chi "
               << photon_chi << std::endl;
        return -1.;
    }
}

// -----------------------------------------------------------------------------
//! Return the value of the function h(particle_chi) of Niel et al.
//! Use an integration of Gauss-Legendre
//
//! \param particle_chi particle quantum parameter
//! \param discretization number of iterations for the Gauss-Legendre integration
//! \param eps epsilon for the modified bessel function
// -----------------------------------------------------------------------------
double NonlinearComptonScattering::computeHNiel( double particle_chi,
                                      int discretization, double eps )
{
    // Arrays for Gauss-Legendre integration
    double *roots = new double[discretization];
    double *weights = new double[discretization];
    double nu;
    double K23, K53;
    double h = 0.;
    double c53 = 5.0/3.0;
    double c23 = 2.0/3.0;
    int    i;

    // Gauss-Legendre coefs between log10(1E-20)
    // and log10(50.) = 1.6989700043360187
    Tools::GaussLegendreQuadrature( -20., std::log10( 50. ), roots,
                                        weights, discretization, eps );

    for( i=0 ; i< discretization ; i++ ) {

        nu = std::pow( 10., roots[i] );

        K53 = Tools::BesselK( c53, nu);

        K23 = Tools::BesselK( c23, nu);

        h += weights[i]*std::log( 10. )*nu
             * ( ( 2.0*std::pow( particle_chi*nu, 3 ) / std::pow( 2.0 + 3.0*nu*particle_chi, 3 )*K53 )
                 + ( 54.0*std::pow( particle_chi, 5 )*std::pow( nu, 4 ) / std::pow( 2.0 + 3.0*nu*particle_chi, 5 )*K23 ) );

    }

    delete[] roots;
    delete[] weights;

    return 9.0*std::sqrt( 3.0 )/( 4.0*M_PI )*h;

}
