/*! \mainpage SMILEI overview
 
 \section intro Introduction
  
 This is the introduction.
 
 You can add interesting cite commands like this \cite Grech2011 (you need to put the BIBTEX record in doc/smilei.bib)
 
 \section install_sec Installation
 
 \subsection Pre-requisites
 
 \subsection compiling
 
 This is (o should be) a massive parallel code which heavily uses MPI routines. 
 The compilation works fine using <tt>g++</tt> and intel <tt>icpc</tt> compiler and openmpi. 
 
 Two basic compilation mode: "release" and debug (accessible via the make commands <tt>make release</tt> and 
 <tt>make debug</tt>). 
 
 "debug" mode is without optimization and prints debug messages from <tt>DEBUG()</tt> macro (see \ref macros page).
 
 "release" mode is with <tt>-O3</tt> option and suppress the debug messages from <tt>DEBUG()</tt> macro.
 
 By default the code is compiled in "release" mode. 
 
 \subsection Dependencies
 The code needs hdf5 libraries installed with parallel support (on mac with macport you need to install the package with:

 <tt>sudo port install hdf5-18 +openmpi</tt>
 
 */

/*! @file Pic.h 
 
 @brief Pic.h 
 
 @author tommaso vinci
 @date 2013-02-15
 */

#include "Tools.h"
#include <string>

//! verbosity level use verbose keyword in data parameter
#ifdef  __DEBUG
unsigned int debug_level = 10;
#endif

//! main function
int main (int argc, char* argv[]);

//! header of the standard message
void startingMessage(std::string inputfile);
