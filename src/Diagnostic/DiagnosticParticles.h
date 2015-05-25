/*
-----------------------------------------------------------------------
PARTICLE DIAGNOSTICS                -    F. Perez - 03/2015
-----------------------------------------------------------------------
  During the simulation, each particle diagnostic collects the data from particles
  into a N-dimensional histogram.
  Each histogram axis can be: x, y, z, px, py, pz, p, gamma, ekin, vx, vy, vz, v or charge.
  In each bin of the histogram, several things may be summed: the weights (density), 
    the weight*charge (charge density) or weight*charge*velocity (current density).
  Examples:
      +----------+------------+---------------------+
      |   Rank   |   type     |         Axes        |
      +----------+------------+---------------------+
      |   1-D    |  density   |        'ekin'       | => energy distribution.
      |   2-D    |  density   |     'x' and 'y'     | => density map.
      |   2-D    |  x-current |     'x' and 'y'     | => x-current map.
      |   2-D    |  density   |     'x' and 'px'    | => phase space.
      |   3-D    |  density   | 'x', 'y' and 'ekin' | => density map for several energy ranges.
      +----------+------------+---------------------+

In the input (namelist) file, each diagnostics are provided as follows:


# DIAGNOSTICS ON PARTICLES - project the particles on a N-D arbitrary grid
# ---------------------------------------------------------------------------------
# output = density, charge_density, or current_density_[xyz]
#              => parameter that describes what quantity is obtained 
# every        => an integer : number of time-steps between each output
# time_average => an integer greater than 0 : number of time-steps to average
# species      => a list of one or several species whose data will be used
# axis   = type min max nsteps [logscale] [edge_inclusive]
#              => `type` can be x, y, z, px, py, pz, p, gamma, ekin, vx, vy, vz, v or charge
#              => the data is binned for `type` between `min` and `max`, in `nsteps` bins
#              => "logscale" sets the binning scale to logarithmic
#              => "edge_inclusive" forces the particles outside (`min`,`max`) to be counted in the extrema bins
#   example : axis = x 0 1 30
#   example : axis = px -1 1 100 
# >>>> MANY AXES CAN BE ADDED IN A SINGLE DIAGNOSTIC <<<<

# EXAMPLE
diag_particles
	output = density
	every = 5
	time_average = 1
	species = electron1
	axis = x    0    1    30
	axis = y    0    1    30
end

# EXAMPLE
diag_particles
	output = density
	every = 5
	time_average = 1
	species = electron1
	axis = ekin  0.0001  0.1 100 logscale
end

*/




#ifndef DiagnosticParticles_H
#define DiagnosticParticles_H

#include <cmath>

#include "Species.h"
#include "Particles.h"
#include "SmileiMPI.h"
#include "H5.h"

// Class for each axis of the particle diags
class DiagnosticParticlesAxis {

public:
    
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

public:

    DiagnosticParticles(unsigned int, std::string, unsigned int, unsigned int, std::vector<unsigned int>, std::vector<DiagnosticParticlesAxis*>);
    
    ~DiagnosticParticles();
    
    void close();
    
    void run(int, std::vector<Species*>&, SmileiMPI*);
    
    int diagnostic_id;
    
private:
    
     //! this is the hdf5 file id (we need to keep it to close at the right time)
    hid_t fileId;
    
    //! quantity to be summed into the output array
    std::string output;
    
    //! period (in timesteps) for outputs
    unsigned int every;
    
    //! number of timesteps during which outputs are averaged
    unsigned int time_average;
    
    //! list of the species that will be accounted for
    std::vector<unsigned int> species;
    
    //! vector of axes
    std::vector<DiagnosticParticlesAxis*> axes;
    
    //! vector for saving the output array for time-averaging
    std::vector<double> data_sum;
    
    int output_size;
    
};




#endif
