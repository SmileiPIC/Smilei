// You can add interesting cite commands like this \cite Grech2011 (you need to put the BIBTEX record in doc/smilei.bib)

/*! \mainpage Overview Home

 \image html smileiLogo.png

 
 \section intro Introduction

 
 The Particle-In-Cell (PIC) code SMILEI is an open-source project developed by the PIC community at the 
 Plateau de Saclay to support the development of the Apollon laser within the CILEX framework. 
 
 SMILEI stands for Simulation of Matter Irradiated by Light at Extreme Intensities, and is developed through a 
 collaboration between various teams at Ecole Polytechnique, at the CEA/Saclay and with strong support from the 
 Maison de la Simulation and IDRIS on the numerical side.
 
 \section download_sec Download

 Currently available upon registration on [llrgit](https://llrgit.in2p3.fr/smilei/) 
 to the community of the Plateau de Saclay, SMILEI is intended as an open-source code.
 Alternatively the open [github](https://github.com/SmileiPIC/Smilei) website. 
 
 */

/*! @file Smilei.h

 @brief Smilei.h

 @date 2013-02-15
 */

#include <string>
#include <vector>

#include "Tools.h"

class Params;
class SmileiMPI;
class VectorPatch;
class Timer;

//! main function
int main (int argc, char* argv[]);

std::vector<Timer> initialize_timers(SmileiMPI* smpi);
