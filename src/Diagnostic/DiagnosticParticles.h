#ifndef DIAGNOSTICPARTICLES_H
#define DIAGNOSTICPARTICLES_H

#include "Diagnostic.h"

#include "Params.h"
#include "Species.h"
#include "Patch.h"
#include "SmileiMPI.h"


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
    
    //! List of coefficients (a,b,c) for a "composite" type of the form "ax+by+cz"
    std::vector<double> coefficients;
    
    //! Pointer to a function that goes through the particles and find where they 
    void (*axis_binning)(Species *, std::vector<double>&, unsigned int, DiagnosticParticlesAxis&);
};


class DiagnosticParticles : public Diagnostic {
    friend class SmileiMPI;

public :
    
    //! Default constructor
    DiagnosticParticles( Params &params, SmileiMPI* smpi, Patch* patch, int diagId );
    //! Cloning constructor
    DiagnosticParticles( DiagnosticParticles* );
    //! Default destructor
    ~DiagnosticParticles() override;
    
    void openFile( Params& params, SmileiMPI* smpi, bool newfile ) override;
    
    void closeFile() override;
    
    bool prepare( int timestep ) override;
    
    void run( Patch* patch, int timestep ) override;
    
    void write(int timestep, SmileiMPI* smpi) override;
    
    //! Clear the array
    void clear();
    
     //! Get memory footprint of current diagnostic
    int getMemFootPrint() override {
        int size = output_size*sizeof(double);
        // + data_array + index_array +  axis_array
        // + nparts_max * (sizeof(double)+sizeof(int)+sizeof(double)) 
        return size;
    }
    
private :

    //! number of timesteps during which outputs are averaged
    int time_average;
    
    //! list of the species that will be accounted for
    std::vector<unsigned int> species;
    
    //! vector of axes
    std::vector<DiagnosticParticlesAxis> axes;
    
    //! quantity to be summed into the output array
    std::string output;
    
    //! vector for saving the output array for time-averaging
    std::vector<double> data_sum;
    
    int output_size;
    
    //! Pointer to a function that goes through the particles and add their contribution to the data_array
    void (*data_filling)(Species *, std::vector<double>&, unsigned int);
};

#endif

