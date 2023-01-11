/*!
 \page macros Code macros

 The c++ macros used in the code should be placed in the \ref Tools.h file.

 \section caveats Warning, Error and Debug
 All these macros will print to the standard error a tag , the name, line of the source file that caused the call
 as well as the name of the function calling the macro. The arguments contained in the parenthesis will then be
 appended. Arguments can be chained together in c++ stream style (using <tt><<</tt> operator)

 The macro <tt>WARNING("text")</tt> is the most basic and is itended for warnings that should always be present in the code.

 The macro <tt>ERROR("text")</tt> is used to print an error and close the program.

 The macro <tt>DEBUG("text")</tt> can be used in two ways: using just an argument, it will display a debug message
 (similar to <tt>WARNING("text")</tt> ) but it can be used in the form <tt>DEBUG(N,"text")</tt> in this case N is a number and
 represents the debug level starting at which the dubug must be displayed.
 The debug level can be changed int the namelist vie the key <tt>debug</tt>.

 */

#ifndef TOOLS_H
#define TOOLS_H

#include <csignal>
#include <cstdlib>

#include <iostream>
#include <sstream>

#include <fstream>

#include <mpi.h>

#define ERROR_STYLE "\033[1;31m"
#define CAREFUL_STYLE "\033[1;36m"
#define FOOTER_STYLE "\033[0m"

// Headers when OMP activated
#ifdef _OMP
#include <omp.h>

#define __header(__msg,__txt) {std::cout << "\t[" << __msg << "](" << omp_get_thread_num() << ") " __FILE__ << ":" << __LINE__ << " (" \
<< __FUNCTION__ << ") " << __txt << std::endl;}

#define __header_custom_text_on_unix(__msg,__txt,__tc) { std::cout << "\033[;"<< __tc << "m" << "\n[" << __msg << "](" << omp_get_thread_num() \
<< ") " __FILE__ << ":" << __LINE__ << " (" << __FUNCTION__ << ") " << __txt << "\033[0m" << std::endl;}

#define __header_error(__msg,__txt) {std::string line = " "; for (int __ic=0; __ic < 80 ; __ic++) line += "-"; \
std::cerr << ERROR_STYLE << line << "\n [" << __msg << "](" << omp_get_thread_num() \
<< ") " __FILE__ << ":" << __LINE__ << " (" << __FUNCTION__ << ") " << __txt << "\n" << line << FOOTER_STYLE << std::endl;}

// Headers when OMP not activated
#else
#define __header(__msg,__txt) std::cout << "\t[" << __msg << "] " << __FILE__ << ":" << __LINE__ << " (" \
<< __FUNCTION__ << ") " << __txt << std::endl

#define __header_custom_text_on_unix(__msg,__txt,__tc) { std::cout << "\033[;"<< __tc <<"m" << "[" << __msg << "] " \
<< __FILE__ << ":" << __LINE__ << " (" << __FUNCTION__ << ") " << __txt << "\033[0m" << std::endl; }

#define __header_error(__msg,__txt) {std::string line = " "; for (int __ic =0; __ic < 80 ; __ic++) line += "-"; \
std::cerr << ERROR_STYLE << line << "\n [" << __msg << "] " << __FILE__ << ":" << __LINE__ << " (" \
<< __FUNCTION__ << ") " << __txt << "\n" << line << FOOTER_STYLE << std::endl;}

#endif

// Header for careful messages
#define __header_careful(__txt) {  \
std::cout << CAREFUL_STYLE << "CAREFUL: " << __txt << "\n" << FOOTER_STYLE << std::endl;}

#define MESSAGE1(__txt)  {int __rk; MPI_Comm_rank( MPI_COMM_WORLD, &__rk ); if (__rk==0) { std::cout << " ";  std::cout << __txt << std::endl;};}
#define MESSAGE2(__val,__txt) {int __rk; MPI_Comm_rank( MPI_COMM_WORLD, &__rk ); if (__rk==0) {for (int __i=0;__i<__val;__i++) std::cout << "\t";}; MESSAGE1(__txt);}
//#define MESSAGE1(__txt)  {;}
//#define MESSAGE2(__val,__txt) {MESSAGE1(__txt);}

#define MESSAGE3(arg1,arg2,arg3,...) arg3
#define MESSAGE4(...) MESSAGE3(__VA_ARGS__,MESSAGE2,MESSAGE1,)
#define MESSAGE(...) MESSAGE4(__VA_ARGS__)(__VA_ARGS__)

#define __PRINTLINE(__num) {MESSAGE(std::string(__num,'-'))}

#define TITLE(...) {MESSAGE(std::endl); MESSAGE(__VA_ARGS__); __PRINTLINE(80);}

// ATTENTION: this costs a lot! use with care!
#define MESSAGEALL1(__txt)  {int __rk; MPI_Comm_rank( MPI_COMM_WORLD, &__rk ); \
int __sz; MPI_Comm_size( MPI_COMM_WORLD, &__sz ); \
for (int __i=0;__i<__sz;__i++ ) {\
if (__i==__rk) {std::cout << "Proc [" << __i << "] " <<__txt << std::endl;} MPI_Barrier( MPI_COMM_WORLD );}}


#define MESSAGEALL2(__val,__txt) {for (int __i=0;__i<__val;__i++) std::cout << "\t"; MESSAGEALL1(__txt);}
#define MESSAGEALL3(arg1,arg2,arg3,...) arg3
#define MESSAGEALL4(...) MESSAGEALL3(__VA_ARGS__,MESSAGEALL2,MESSAGEALL1,)
#define MESSAGEALL(...) MESSAGEALL4(__VA_ARGS__)(__VA_ARGS__)

// ATTENTION: this costs a lot! use with care!
#define WARNINGALL(__txt) {int __rk; MPI_Comm_rank( MPI_COMM_WORLD, &__rk ); int __sz; MPI_Comm_size( MPI_COMM_WORLD, &__sz ); for (int __i=0;__i<__sz;__i++ ) {if (__i==__rk) {__header("WARNING proc "<<__i, __txt);} MPI_Barrier( MPI_COMM_WORLD );}}


#define PMESSAGE1(rank, __txt)  {std::cout << __FILE__ << ":" << __LINE__ << "[Process " << rank << "] : " << __txt << std::endl;}
#define PMESSAGE2(__val,rank,__txt) {for (int __i=0;__i<__val;__i++) std::cout << "\t"; std::cout << "[Process " << rank << "], " << "[" << __val << "] " << __txt << std::endl;}
#define PMESSAGE3(arg1,arg2,arg3,arg4,...) arg4
#define PMESSAGE4(...) PMESSAGE3(__VA_ARGS__,PMESSAGE2,PMESSAGE1,)
#define PMESSAGE(...) PMESSAGE4(__VA_ARGS__)(__VA_ARGS__)

#ifdef  __DEBUG

#define WARNING(__txt) {int __rk; MPI_Comm_rank( MPI_COMM_WORLD, &__rk ); if (__rk==0) {__header_custom_text_on_unix("WARNING", __txt, 33);}}

#define DEBUG(__txt) {__header("DEBUG", __txt);}

#define ERRORWITHCUSTOMSIGNAL(__txt, __sig) {int __rk; MPI_Comm_rank( MPI_COMM_WORLD, &__rk ); __header("ERROR proc "<<__rk, __txt); raise(__sig);}

#define DEBUGEXEC(...) __VA_ARGS__
#define RELEASEEXEC(...)

#define HEREIAM(__txt) {const int __num_minus=10; int __rk; MPI_Comm_rank( MPI_COMM_WORLD, &__rk ); for(int __i=0;__i<__num_minus;__i++) {std::cout << "-";}; std::cout << "> " << __FILE__ << ":" << __LINE__ << " (" << __FUNCTION__ << ")[" << __rk << "](" << omp_get_thread_num() << ") " << __txt << " <" ; for(int __i=0;__i<__num_minus;__i++) {std::cout << "-";}; std::cout << std::endl; }

#else // not DEBUG

#define WARNING(__txt) {int __rk; MPI_Comm_rank( MPI_COMM_WORLD, &__rk ); if (__rk==0) {__header_custom_text_on_unix("WARNING", __txt, 33);}}

#define DEBUG(...)
#define DEBUGEXEC(...)
#define RELEASEEXEC(...) __VA_ARGS__

#define ERRORWITHCUSTOMSIGNAL(__txt, __sig) {__header_error("ERROR", __txt); raise(__sig);}

#define HEREIAM(...)

#endif // __DEBUG

#define ERROR(__txt) {ERRORWITHCUSTOMSIGNAL(__txt, SIGABRT);}

#define ERROR_NAMELIST(__txt,__link) {      \
    int __rk;                               \
    MPI_Comm_rank( MPI_COMM_WORLD, &__rk ); \
    if (__rk==0) {                          \
        std::string __link_message = "";    \
        if (std::string(__link) != "") {    \
            __link_message = "\n\n Find out more: " + std::string(__link); \
        };                                  \
        ERRORWITHCUSTOMSIGNAL("\n A probem was found in the namelist:\n > " << __txt << __link_message, SIGABRT); \
    };                                      \
}

#define CAREFUL(__val, __txt) {      \
    int __rk;                               \
    MPI_Comm_rank( MPI_COMM_WORLD, &__rk ); \
    if (__rk==0) {                          \
        for (int __i=0;__i<__val;__i++) std::cout << "\t"; \
        __header_careful(__txt);            \
    };                                      \
}


class Tools
{
public:
    static void printMemFootPrint( std::string tag );
    static double getMemFootPrint(int type_of_memory);

    //! Converts a number of Bytes in a readable string in KiB, MiB, GiB or TiB
    static std::string printBytes( uint64_t nbytes );

    //! This function returns true/flase whether the file exists or not
    //! \param file file name to test
    static bool fileExists( const std::string &filename ) ;

    static std::string xyz;

    //! Concatenate several strings
    template<class T1, class T2, class T3=std::string, class T4=std::string>
    static std::string merge( T1 s1, T2 s2, T3 s3="", T4 s4="" )
    {
        std::ostringstream tot( "" );
        tot << s1 << s2 << s3 << s4;
        return tot.str();
    }

    //! Wrapper to get the thread number
    static int getOMPThreadNum()
    {
#ifdef _OPENMP
        return omp_get_thread_num();
#else
        return 0;
#endif
    }

};

#define LINK_NAMELIST "https://smileipic.github.io/Smilei/namelist.html"

#if defined(WIN32) || defined(_WIN32)
#define PATH_SEPARATOR "\\"
#else
#define PATH_SEPARATOR "/"
#endif

//////////////////////////////////////
/// @def SMILEI_UNUSED
///
/// Prevent a compiler warning for an unused variable.
/// It is usefull when a lot of conditional compilation is used (#if XXX then compile code for ARM #else for x86 #endif).
/// This has no performance overhead.
///
/// Example usage:
///
///     int a;
///     SMILEI_UNUSED(a); // Without this macro, it would trigger a "variable unused" warning
///     return;
///
//////////////////////////////////////
#define SMILEI_UNUSED( an_arg ) static_cast<void>( an_arg )

#endif
