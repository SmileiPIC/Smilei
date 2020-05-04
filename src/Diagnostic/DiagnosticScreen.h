#ifndef DIAGNOSTICSCREEN_H
#define DIAGNOSTICSCREEN_H

#include "DiagnosticParticleBinningBase.h"

class DiagnosticScreen : public DiagnosticParticleBinningBase
{
    friend class SmileiMPI;
    friend class Checkpoints;
    
public :

    //! Default constructor
    DiagnosticScreen( Params &params, SmileiMPI *smpi, Patch *patch, int diagId );
    //! Default destructor
    ~DiagnosticScreen();
    
    bool prepare( int timestep ) override;
    
    void run( Patch *patch, int timestep, SimWindow *simWindow ) override;
    
    bool writeNow( int timestep ) override;
    
    //! Clear the array
    void clear() override;
    
    static std::vector<std::string> excludedAxes( int idiag ) {
        std::string shape = "";
        PyTools::extract( "shape", shape, "DiagScreen", idiag );
        std::vector<std::string> excluded_axes( 0 );
        if( shape == "plane" ) {
            excluded_axes.push_back( "theta_yx" );
            excluded_axes.push_back( "theta_zx" );
        } else {
            excluded_axes.push_back( "a" );
            excluded_axes.push_back( "b" );
        }
        return excluded_axes;
    }
    
    std::vector<double> * getData() {
        return &data_sum;
    }
    
private :

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

