/*
-----------------------------------------------------------------------
TEST PARTICLE DIAGNOSTICS
-----------------------------------------------------------------------
*/



#ifndef DiagnosticTestParticles_H
#define DiagnosticTestParticles_H

#include <cmath>

#include "Params.h"
#include "Species.h"
#include "Particles.h"
#include "SmileiMPI.h"
#include "H5.h"


// Class for the test-particles diagnostics
class DiagnosticTestParticles {

public:

    DiagnosticTestParticles(unsigned int, int, Params&, std::vector<Species*>&);
    
    ~DiagnosticTestParticles();
    
    //! ID of this diagnostic
    int diagnostic_id;
    
    //! number of the test species used
    int species_number;
    
    //! Initializes the diag (creates HDF5 file), after all species have been created
    void init( std::vector<Species*>, SmileiMPI* );
    
    //! Runs the diag (writes to file) at each timestep
    void run( int, SmileiMPI* );

private:
    
    //! HDF5 file transfer protocol
    hid_t transfer;
    //! HDF5 memory space (dimensions of the array in memory)
    hid_t mem_space;
    //! HDF5 file space (dimensions of the array in file)
    hsize_t dims[2];
    
    //! Pointer to the test species used
    Species* species;
    //! Pointer to the test particles used
    Particles * particles;
    
    //! Number of spatial dimensions
    int nDim_particle;
    
    //! Number of  timesteps between each output
    int every;
    
    //! Adds one row in a HDF5 file, within a given dataspace
    template <class T> void append(hid_t, std::string, std::vector<T>, int, hid_t, SmileiMPI*, hsize_t *);
    
};




#endif
