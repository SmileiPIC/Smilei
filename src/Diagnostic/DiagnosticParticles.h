/*
-----------------------------------------------------------------------
PARTICLE DIAGNOSTICS                -    F. Perez - 03/2015
-----------------------------------------------------------------------
  During the simulation, each particle diagnostic collects the data from particles
  into a N-dimensional histogram.
  Each histogram axis can be: x, y, z, px, py, pz, p, gamma, ekin, vx, vy, vz, v or charge.
  In each bin of the histogram, several things may be summed: the weights (density), 
    weight*charge (charge density), weight*charge*velocity (current density),
    or weight*momentum (momentum density)
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
# ------------------------------------------------------------------------
# output       = string: "density", "charge_density", etc.
#                parameter that describes what quantity is obtained 
# every        = integer > 0: number of time-steps between each output
# time_average = integer > 0: number of time-steps to average
# species      = list of strings, one or several species whose data will be used
# axes         = list of axes
# Each axis is a list: (_type_ _min_ _max_ _nsteps_ ["logscale"] ["edge_inclusive"])
#   _type_ is a string, one of the following options:
#      x, y, z, px, py, pz, p, gamma, ekin, vx, vy, vz, v or charge
#   The data is discretized for _type_ between _min_ and _max_, in _nsteps_ bins
#   The optional "logscale" sets the scale to logarithmic
#   The optional "edge_inclusive" forces the particles that are outside (_min_,_max_)
#     to be counted in the extrema bins
#   Example : axes = ("x", 0, 1, 30)
#   Example : axes = ("px", -1, 1, 100, "edge_inclusive")

# EXAMPLE
DiagParticles(
	output = "density",
	every = 5,
	time_average = 1,
	species = ["electron1"],
	axes = [
		 ["x", 0,  1, 30],
		 ["y", 0,  1, 30]
	]
)

# EXAMPLE
DiagParticles(
	output = "density",
	every = 5,
	time_average = 1,
	species = ["electron1"],
	axes = [
		 ["ekin", 0.0001, 0.1, 100, "logscale"]
	]
)

*/




#ifndef DiagnosticParticles_H
#define DiagnosticParticles_H

#include <cmath>

#include "Species.h"
#include "Particles.h"
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
    
    void run(int, std::vector<Species*>&);
    
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
