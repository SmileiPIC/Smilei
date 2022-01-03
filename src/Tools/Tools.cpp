
#include "Tools.h"
#include <cstring>
#include <stdio.h>

#include <fcntl.h>
#include <iomanip>
#include <unistd.h>
#include <sstream>

void Tools::printMemFootPrint( std::string tag )
{

    int val[4];
    
    char filename[80];
    char sbuf[1024];
    char *S;
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
    S=strstr( sbuf, "VmRSS:" )+6;
    val[1] = ( int )atoi( S );
    // Peak virtual memory usage
    S=strstr( sbuf, "VmSize:" )+6;
    val[2] = atoi( S );
    
    std::cout << "=== Mem usage === " <<  std::setw( 20 ) << tag << "\t=== " << std::setw( 6 )
              << "\t VmRSS  << " << ( int )( ( double )val[1]/1024. ) << " Mb" << std::endl;
              
}

double Tools::getMemFootPrint(int type_of_memory)
{

    int val[4];
    
    char filename[80];
    char sbuf[1024];
    char *S;
    //long lmem;
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
    if (type_of_memory == 0){
        S=strstr( sbuf, "VmRSS:" )+6;
    } else if (type_of_memory == 1){
        S=strstr( sbuf, "VmHWM:" )+6;
    } 

    val[1] = ( int )atoi( S );
    
    // Return RSS in GB
    return ( double )val[1]/1024./1024.;
    
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
