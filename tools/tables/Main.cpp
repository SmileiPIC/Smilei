// ---------------------------------------------------------------------------------------------------------------------
//                                                   MAIN CODE
// ---------------------------------------------------------------------------------------------------------------------

#include <mpi.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <stdio.h>
#include <limits>
#include <iomanip>
#include "Tools.h"
#include "NonlinearInverseComptonScattering.h"
#include "MultiphotonBreitWheeler.h"
#include "H5.h"

int main( int argc, char *argv[] )
{
    // MPI parameters
    int mpi_provided;
    int rank;
    int number_of_ranks;
    // int error;

    // Initialization of MPI
    MPI_Init_thread( 0, 0, MPI_THREAD_MULTIPLE, &mpi_provided );
    MPI_Comm_size( MPI_COMM_WORLD, &number_of_ranks );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    std::string table_mode;
    std::string help_message;
    help_message =  "\n This tool generates Smilei tables for different physical mechanisms.\n";
    help_message += " The first keyword to provide represents the table mode you want ot generate:\n";
    help_message += " - 'nonlinear inverse Compton scattering' or 'nics'\n";
    help_message += " - 'multiphoton Breit Wheeler' or 'mbw'\n";
    help_message += "\n";
    help_message += " List of available commands:\n";
    help_message += " -h    print a help message and exit.\n";
    
    // _______________________________________________________________________
    // Read from command line
    
    // We pack arguments in a list of string
    std::string * arguments = new std::string[argc];
    for (int i = 0 ; i < argc ; i++) {
        arguments[i] = argv[i];
    }
    
    if( rank==0 ) {
        std::cout << " _______________________________________________________________________ \n\n"
                  << " Smilei Tables \n"
                  << " _______________________________________________________________________ " << std::endl;
    }
    
    if (argc < 2) {
        ERROR("Please, specify which table mode to generate.\n\n" << help_message );
    } else {
        if (arguments[1] == "nics" || arguments[1] == "nonlinear inverse Compton scattering") {
            table_mode = "nics";
        }
        else if (arguments[1] == "mbw" || arguments[1] == "multiphoton Breit Wheeler") {
            table_mode = "mbw";
        }
        else if (arguments[1] == "-h") {
            if (rank ==0 ) {
                std::cout << help_message << std::endl;
            }
            exit(0);
        }
        else {
            ERROR(" Unknown argument " << arguments[1]);
        }
    }
    
    if (table_mode == "nics") {
        NonlinearComptonScattering::createTables(argc, arguments);
    } else if (table_mode == "mbw") {
        MultiphotonBreitWheeler::createTables(argc, arguments);
    }

    // std::cout << std::setprecision(10) << Tools::asymptoticBesselK(2.0/3.0, 2.) << std::endl;
    // std::cout << std::cyl_bessel_k( 2.0/3.0, 2.) << std::endl;

    if (rank==0) {
        std::cout << "\n All tables generated." << std::endl;
    }

    // Finalization of MPI
    MPI_Finalize();
}
