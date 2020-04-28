#ifndef DIAGNOSTICPARTICLEBINNING_H
#define DIAGNOSTICPARTICLEBINNING_H

#include "Diagnostic.h"

#include "Histogram.h"

class DiagnosticParticleBinning : public Diagnostic
{
    friend class SmileiMPI;
    
public :

    //! Default constructor
    DiagnosticParticleBinning(
        Params &params,
        SmileiMPI *smpi,
        Patch *patch,
        int diagId,
        std::string diagName = "ParticleBinning",
        bool time_accumulate = false,
        PyObject *deposited_quantity = nullptr
    );
    //! Default destructor
    ~DiagnosticParticleBinning() override;
    
    virtual void openFile( Params &params, SmileiMPI *smpi, bool newfile ) override;
    
    void closeFile() override;
    
    bool prepare( int timestep ) override;
    
    virtual void run( Patch *patch, int timestep, SimWindow *simWindow ) override;
    
    virtual bool writeNow( int timestep );
    
    void write( int timestep, SmileiMPI *smpi ) override;
    
    //! Clear the array
    virtual void clear();
    
    virtual std::vector<std::string> excludedAxes() {
        std::vector<std::string> excluded_axes( 0 );
        excluded_axes.push_back( "a" );
        excluded_axes.push_back( "b" );
        excluded_axes.push_back( "theta" );
        excluded_axes.push_back( "phi" );
        return excluded_axes;
    }
    
    //! Get memory footprint of current diagnostic
    int getMemFootPrint() override
    {
        int size = output_size*sizeof( double );
        // + data_array + index_array +  axis_array
        // + nparts_max * (sizeof(double)+sizeof(int)+sizeof(double))
        return size;
    };
    
    //! Get disk footprint of current diagnostic
    uint64_t getDiskFootPrint( int istart, int istop, Patch *patch ) override;
    
protected:
    
    //! True for Screen only
    bool time_accumulate;
    
    //! number of timesteps during which outputs are averaged
    int time_average;
    
    //! list of the species that will be accounted for
    std::vector<unsigned int> species;
    
    //! vector for saving the output array for time-averaging
    std::vector<double> data_sum;
    
    //! Histogram object
    Histogram *histogram;
    
    unsigned int output_size;
    
    int total_axes;
    std::vector<hsize_t> dims;
    
//    //! Minimum and maximum spatial coordinates that are useful for this diag
//    std::vector<double> spatial_min, spatial_max;
};

#endif

