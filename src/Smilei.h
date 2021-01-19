#include <string>
#include <vector>

#include "Tools.h"
#include "RadiationTables.h"
#include "MultiphotonBreitWheelerTables.h"

class Params;
class SmileiMPI;
class VectorPatch;
class Region;
class Checkpoint;
class SimWindow;
class OpenPMDparams;
class Timer;

//! main function
int main( int argc, char *argv[] );
int executeTestMode( VectorPatch &vecPatches,
		     Region &region,
                     SmileiMPI *smpi,
                     SimWindow *simWin, Params &params, Checkpoint &checkpoint, OpenPMDparams &openPMD, RadiationTables * radiation_tables_ );


std::vector<Timer> initialize_timers( SmileiMPI *smpi );
