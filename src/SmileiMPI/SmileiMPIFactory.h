#ifndef SMILEIMPIFACTORY_H
#define SMILEIMPIFACTORY_H

#include <iostream>
#include <string>

#include "SmileiMPI.h"
#include "SmileiMPI_test.h"

class SmileiMPIFactory {
public:

    static SmileiMPI* create(int* argc, char*** argv) {
        // If argument starts with -T, then enter "test mode"
        if( *argc > 1 && strncmp((*argv)[1], "-T", 2)==0 ) {
            int nMPI=1, nOMP=1;
            std::string MPIxOMP = (*argv)[1] + 2;
            (*argv)++;
            (*argc)--;
            // Parse the optional argument MPIxOMP
            if( MPIxOMP.size() > 0 ) {
                for( unsigned int i=0; i<MPIxOMP.size(); i++ ) {
                    if( !std::isdigit(MPIxOMP[i]) && MPIxOMP[i]!='x' ) {
                        ERROR("The option -T requires the format -TMPIxOMP where MPI and OMP are integers");
                    }
                }
                unsigned int xpos = MPIxOMP.find("x");
                if( xpos == std::string::npos || xpos != MPIxOMP.rfind("x")) {
                    ERROR("The option -T requires the format -TMPIxOMP where MPI and OMP are integers");
                }
                nMPI = std::stoi( MPIxOMP.substr( 0, xpos ).c_str() );
                nOMP = std::stoi( MPIxOMP.substr( xpos+1, MPIxOMP.size()-xpos-1 ).c_str() );
                if( nMPI<1 || nOMP<1 ) {
                    ERROR("The option -T requires the format -TMPIxOMP where MPI and OMP are >=1");
                }
            }
            // Return a test-mode SmileiMPI
            return new SmileiMPI_test( nMPI, nOMP, argc, argv );
        
        
        // Otherwise, return the usual SmileiMPI
        } else {
            return new SmileiMPI( argc, argv );
        }
    
    };

};

#endif