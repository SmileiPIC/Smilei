#ifndef DIAGNOSTICPARTICLES_H
#define DIAGNOSTICPARTICLES_H

#include "Diagnostic.h"

#include "Params.h"
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
    
};


class DiagnosticParticles : public Diagnostic {
    friend class SmileiMPI;

public :
    
    //! Default constructor
    DiagnosticParticles( Params &params, SmileiMPI* smpi, Patch* patch, int diagId );
    //! Cloning constructor
    DiagnosticParticles( DiagnosticParticles* );
    //! Default destructor
    ~DiagnosticParticles() ;
    
    void openFile( Params& params, SmileiMPI* smpi, bool newfile ) ;
    
    void closeFile() ;
    
    bool prepare( int timestep ) ;
    
    void run( Patch* patch, int timestep ) ;
    
    void write(int timestep, SmileiMPI* smpi) ;
    
    //! Clear the array
    void clear();
     
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

};

#endif

