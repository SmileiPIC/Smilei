/*! \mainpage SMILEI overview
 
 \section intro Introduction
 
 This is the introduction.
 
 \section install_sec Installation
 
 \subsection Pre-requisites
 
 \subsection compiling
 
 The compilation wrks fine using <tt>g++</tt> and intel <tt>icpc</tt> compiler. 
 
 Two basic compilation mode: "release" and debug (accessible via the make commands <tt>make release</tt> and 
 <tt>make debug</tt>). 
 
 "debug" mode is without optimization and prints debug messages from <tt>DEBUG()</tt> macro (see \ref macros page).
 
 "release" mode is with <tt>-O3</tt> option and suppress the debug messages from <tt>DEBUG()</tt> macro.
 
 By default the code is compiled in "release" mode. 
 
 */

/*! @file Pic.h 
 
 @brief Pic.h 
 
 @author tommaso vinci
 @date 2013-02-15
 */

//! verbosity level use verbose keyword in data parameter
#ifdef  __DEBUG
unsigned int debug_level = 10;
#endif

//! main function
int main (int argc, char* argv[]);
