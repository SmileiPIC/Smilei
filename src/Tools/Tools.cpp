
#include "Tools.h"
#include <cstring>
#include <stdio.h>

#include <fcntl.h>
#include <iomanip>
#include <unistd.h>
#include <sstream>

void Tools::printMemFootPrint( std::string tag )
{
    char filename[80];
    char sbuf[1024];
    //long lmem;
    pid_t numpro;
    
    numpro = getpid();
    
    sprintf( filename, "/proc/%ld/status", ( long )numpro );
    int fd = open( filename, O_RDONLY, 0 );
    int num_read=read( fd, sbuf, ( sizeof sbuf )-1 );
    close( fd );
    
    if( !num_read ) {
        return;
    }
    
    // Peak resident set size
    int VmRSS = atoi( strstr( sbuf, "VmRSS:" )+6 );
    // Peak virtual memory usage
    // int VmSize = atoi( strstr( sbuf, "VmSize:" )+7 );
    
    std::cout << "=== Mem usage === " <<  std::setw( 20 ) << tag << "\t=== " << std::setw( 6 )
              << "\t VmRSS  << " << ( int )( ( double )VmRSS/1024. ) << " Mb" << std::endl;
              
}

double Tools::getMemFootPrint(int type_of_memory)
{
    char filename[80];
    char sbuf[1024];
    pid_t numpro;
    
    numpro = getpid();
    
    sprintf( filename, "/proc/%ld/status", ( long )numpro );
    
    if( ! fileExists( filename ) ) {
        return 0;
    }
    
    int fd = open( filename, O_RDONLY, 0 );
    int num_read=read( fd, sbuf, ( sizeof sbuf )-1 );
    close( fd );
    
    if( !num_read ) {
        return 0;
    }
    
    // Peak resident set size
    int Vm = 0;
    if (type_of_memory == 0){
        Vm = atoi( strstr( sbuf, "VmRSS:" )+6 );
    } else if (type_of_memory == 1){
        Vm = atoi( strstr( sbuf, "VmHWM:" )+6 );
    } 
    
    // Return RSS in GB
    return ( double )Vm/1024./1024.;
}


std::string Tools::printBytes( uint64_t nbytes )
{
    std::ostringstream t( "" );
    if( nbytes < 1024 ) {
        t << nbytes << " bytes";
    } else {
        double nb = nbytes;
        if( nbytes < 1048576 ) {
            t << std::fixed << std::setprecision( 2 ) << ( nb/1024. ) << " K";
        } else if( nbytes < 1073741824 ) {
            t << std::fixed << std::setprecision( 2 ) << ( nb/1048576. ) << " M";
        } else if( nbytes < 1099511627776 ) {
            t << std::fixed << std::setprecision( 2 ) << ( nb/1073741824. ) << " G";
        } else {
            t << std::fixed << std::setprecision( 2 ) << ( nb/1099511627776. ) << " T";
        }
    }
    return t.str();
}

//! Wrapper to get the thread number
int Tools::getOMPThreadNum()
{
#ifdef _OPENMP
    return ::omp_get_thread_num();
#else
    return 0;
#endif
}

// ---------------------------------------------------------------------------------------------------------------------
//! This function returns true/flase whether the file exists or not
//! \param file file name to test
// ---------------------------------------------------------------------------------------------------------------------
bool Tools::fileExists( const std::string &filename )
{
    std::ifstream file( filename.c_str() );
    return !file.fail();
}




std::string Tools::xyz = "XYZ";
