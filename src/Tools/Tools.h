/*!
 \page macros Code macros
 
 The c++ macros used in the code should be paced in the \ref Tools.h file.
 
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

#include <mpi.h>

#define __header(__msg,__txt) std::cerr << "\t[" << __msg << "] " << __FILE__ << ":" << __LINE__ << " (" \
<< __FUNCTION__ << ") " << __txt << std::endl

#define MESSAGE1(__txt)  {int __rk; MPI_Comm_rank( MPI_COMM_WORLD, &__rk ); if (__rk==0) { std::cout << " ";  std::cout << __txt << std::endl;};}
#define MESSAGE2(__val,__txt) {for (int __i=0;__i<__val;__i++) std::cout << "\t"; MESSAGE1(__txt);}

#define MESSAGE3(arg1,arg2,arg3,...) arg3
#define MESSAGE4(...) MESSAGE3(__VA_ARGS__,MESSAGE2,MESSAGE1,)
#define MESSAGE(...) MESSAGE4(__VA_ARGS__)(__VA_ARGS__)

// ATTENTION: this costs a lot! use with care!
#define MESSAGEALL1(__txt)  {int __rk; MPI_Comm_rank( MPI_COMM_WORLD, &__rk ); \
int __sz; MPI_Comm_size( MPI_COMM_WORLD, &__sz ); \
for (int __i=0;__i<__sz;__i++ ) {\
if (__i==__rk) {std::cout << "Proc [" << __i << "] " <<__txt << std::endl;} MPI_Barrier( MPI_COMM_WORLD );}}


#define MESSAGEALL2(__val,__txt) {for (int __i=0;__i<__val;__i++) std::cout << "\t"; MESSAGEALL1(__txt);}
#define MESSAGEALL3(arg1,arg2,arg3,...) arg3
#define MESSAGEALL4(...) MESSAGEALL3(__VA_ARGS__,MESSAGEALL2,MESSAGEALL1,)
#define MESSAGEALL(...) MESSAGEALL4(__VA_ARGS__)(__VA_ARGS__)

#define WARNING(__txt) {int __rk; MPI_Comm_rank( MPI_COMM_WORLD, &__rk ); if (__rk==0) {__header("WARNING", __txt);}}

// ATTENTION: this costs a lot! use with care!
#define WARNINGALL(__txt) {int __rk; MPI_Comm_rank( MPI_COMM_WORLD, &__rk ); int __sz; MPI_Comm_size( MPI_COMM_WORLD, &__sz ); for (int __i=0;__i<__sz;__i++ ) {if (__i==__rk) {__header("WARNING proc "<<__i, __txt);} MPI_Barrier( MPI_COMM_WORLD );}}


#define PMESSAGE1(rank, __txt)  {std::cout << __FILE__ << ":" << __LINE__ << "[Process " << rank << "] : " << __txt << std::endl;}
#define PMESSAGE2(__val,rank,__txt) {for (int __i=0;__i<__val;__i++) std::cout << "\t"; std::cout << "[Process " << rank << "], " << "[" << __val << "] " << __txt << std::endl;}
#define PMESSAGE3(arg1,arg2,arg3,arg4,...) arg4
#define PMESSAGE4(...) PMESSAGE3(__VA_ARGS__,PMESSAGE2,PMESSAGE1,)
#define PMESSAGE(...) PMESSAGE4(__VA_ARGS__)(__VA_ARGS__)


#ifdef  __DEBUG

//#warning "DEBUG MODE "
extern int debug_level;

#define DEBUG1(__txt) {if(debug_level>=0) __header("DEBUG", __txt);}
#define DEBUG2(__val,__txt) if(((debug_level<0) && __val==-debug_level) || ((debug_level>=0) && __val<=debug_level)) __header("DEBUG "<<__val, __txt)
#define DEBUG3(arg1,arg2,arg3,...) arg3
#define DEBUG4(...) DEBUG3(__VA_ARGS__,DEBUG2,DEBUG1,)
#define DEBUG(...) DEBUG4(__VA_ARGS__)(__VA_ARGS__)

#define ERROR(__txt) {int __rk; MPI_Comm_rank( MPI_COMM_WORLD, &__rk ); __header("ERROR proc "<<__rk, __txt); MPI_Finalize(); raise(SIGSEGV);}

#define DEBUGEXEC(...) __VA_ARGS__
#define RELEASEEXEC(...)

#define HEREIAM(__txt) {const int __num_minus=40; int __rk; MPI_Comm_rank( MPI_COMM_WORLD, &__rk ); for(int __i=0;__i<__num_minus;__i++) {std::cerr << "-";}; std::cerr << "> " << __rk << " " << __FILE__ << ":" << __LINE__ << " (" << __FUNCTION__ << ") " << __txt << " <" ; for(int __i=0;__i<__num_minus;__i++) {std::cerr << "-";}; std::cerr << std::endl; }

#else
#define DEBUG(...)
#define DEBUGEXEC(...)
#define RELEASEEXEC(...) __VA_ARGS__

#define ERROR(__txt) {__header("ERROR", __txt); MPI_Finalize(); exit(EXIT_FAILURE);}

#define HEREIAM(...)

#endif

class Tools {
 public:
  static void printMemFootPrint(std::string tag);
};

#endif
