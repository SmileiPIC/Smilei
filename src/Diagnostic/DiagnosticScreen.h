#ifndef DIAGNOSTICSCREEN_H
#define DIAGNOSTICSCREEN_H

#include "Diagnostic.h"

#include "Params.h"
#include "Species.h"
#include "Patch.h"
#include "SmileiMPI.h"
#include "Histogram.h"

class DiagnosticScreen : public Diagnostic
{
    friend class SmileiMPI;
    
public :

    //! Default constructor
    DiagnosticScreen( Params &params, SmileiMPI *smpi, Patch *patch, int diagId );
    //! Cloning constructor
    DiagnosticScreen( DiagnosticScreen * );
    //! Default destructor
    ~DiagnosticScreen() override;
    
    void openFile( Params &params, SmileiMPI *smpi, bool newfile ) override;
    
    void closeFile() override;
    
    bool prepare( int timestep ) override;
    
    void run( Patch *patch, int timestep, SimWindow *simWindow ) override;
    
    void write( int timestep, SmileiMPI *smpi ) override;
    
    //! Clear the array
    void clear();
    
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
    
    //! vector for saving the output array
    std::vector<double> data_sum;
    
    //! Id of this diag
    int screen_id;
    
private :

    //! list of the species that will be accounted for
    std::vector<unsigned int> species;
    
    //! Histogram object
    Histogram *histogram;
    
    unsigned int output_size;
    
    std::string screen_shape;
    //! Relates to the shape of the screen (plane=0, sphere=1)
    int screen_type;
    
    //! Screen reference point (plane point or sphere center)
    std::vector<double> screen_point;
    
    //! Screen reference vector (plane normal or sphere radius)
    std::vector<double> screen_vector;
    std::vector<double> screen_unitvector;
    //! norm of the vector
    double screen_vectornorm;
    //! Vectors forming an orthogonal base with the screen vector
    std::vector<double> screen_vector_a, screen_vector_b;
    
    //! How to account for the particle direction: "both", "canceling", "forward" or "backward"
    std::string direction;
    int direction_type;
    
    //! Copy of the timestep
    double dt;
};

#endif

