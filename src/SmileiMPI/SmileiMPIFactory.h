#ifndef SMILEIMPIFACTORY_H
#define SMILEIMPIFACTORY_H

#include <iostream>
#include <string>

#include "SmileiMPI.h"
#include "SmileiMPI_test.h"

class SmileiMPIFactory {
public:

    static SmileiMPI* create(int* argc, char*** argv) {
    
        MESSAGE("                   _            _");
        MESSAGE(" ___           _  | |        _  \\ \\   Version : " << __VERSION);
        MESSAGE("/ __|  _ __   (_) | |  ___  (_)  | |   ");
        MESSAGE("\\__ \\ | '  \\   _  | | / -_)  _   | |");
        MESSAGE("|___/ |_|_|_| |_| |_| \\___| |_|  | |  ");
        MESSAGE("                                /_/    ");
        MESSAGE("");
        
        // If first argument is -T, then enter "test mode"
        if( *argc > 1 && strcmp((*argv)[1], "-T") ) {
            int nMPI=1, nOMP=1;
            (*argv)++;
            (*argc)--;
            // Parse the optional argument MPIxOMP
            if( *argc > 1 ) {
                std::string MPIxOMP = (*argv)[1];
                (*argv)++;
                (*argc)--;
                for( unsigned int i=0; i<MPIxOMP.size(); i++ ) {
                    if( !std::isdigit(MPIxOMP[i]) && MPIxOMP[i]!='x' ) {
                        ERROR("The option provided after the -T argument requires the format MPIxOMP where MPI and OMP are integers");
                    }
                }
                unsigned int xpos = MPIxOMP.find("x");
                if( xpos == std::string::npos || xpos != MPIxOMP.rfind("x")) {
                    ERROR("The option provided after the -T argument requires the format MPIxOMP where MPI and OMP are integers");
                }
                nMPI = std::stoi( MPIxOMP.substr( 0, xpos ).c_str() );
                nOMP = std::stoi( MPIxOMP.substr( xpos+1, MPIxOMP.size()-xpos-1 ).c_str() );
                if( nMPI<1 || nOMP<1 ) {
                    ERROR("The option provided after the -T argument requires the format MPIxOMP where MPI and OMP are >=1");
                }
            }
            MESSAGE("    ----- TEST MODE -----");
            // Return a test-mode SmileiMPI
            return new SmileiMPI_test( nMPI, nOMP );
        
        
        // Otherwise, return the usual SmileiMPI
        } else {
            return new SmileiMPI( argc, argv );
        }
    
    };

};

#endif