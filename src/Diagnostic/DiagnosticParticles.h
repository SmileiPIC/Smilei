/*
-----------------------------------------------------------------------
PARTICLE DIAGNOSTICS                -    F. Perez - 03/2015
-----------------------------------------------------------------------
  see the doc for instructions and examples
*/

#ifndef DiagnosticParticles_H
#define DiagnosticParticles_H

#include <cmath>

#include "Species.h"
#include "Particles.h"
#include "H5.h"
#include "TimeSelection.h"

class Patch;

// Class for each axis of the particle diags
struct DiagnosticParticlesAxis {

    //! quantity of the axis (e.g. 'x', 'px', ...)
    std::string type;
    
    //! starting point for the axis binning
    double min;
    //! ending point for the axis binning
    double max;
    //! number of bins for the axis binning
    int nbins;
    
    //! determines whether linear scale or log scale
    bool logscale;
    
    //! determines whether particles beyond min and max are counted in the first and last bin
    bool edge_inclusive;
    
};


// Class for the particles diagnostics
class DiagnosticParticles {
    friend class SmileiMPI;
    friend class DiagsVectorPatch;

public:

    DiagnosticParticles(unsigned int, Params& params, Patch* patch, std::vector<Species*>& vecSpecies);
    
    ~DiagnosticParticles();
        
    void run(int, std::vector<Species*>&);

    void createFile( unsigned int n_diag_particles );
    void write(int timestep);
    void clean();

protected:
    

private:
    
     //! this is the hdf5 file name
    std::string filename;
    
    //! quantity to be summed into the output array
    std::string output;
    
    //! period (in timesteps) for outputs
    unsigned int every;
//    //! Time selection
//    TimeSelection * timeSelection;
    
    //! number of timesteps during which outputs are averaged
    unsigned int time_average;
    
    //! list of the species that will be accounted for
    std::vector<unsigned int> species;
    
    //! vector of axes
    std::vector<DiagnosticParticlesAxis> axes;
    
    //! vector for saving the output array for time-averaging
    std::vector<double> data_sum;
    
    int output_size;
    
};




#endif
