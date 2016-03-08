#ifndef DIAGPARTICLES_H
#define DIAGPARTICLES_H

#include "Diag.h"

#include "Params.h"
#include "Patch.h"
#include "SmileiMPI.h"


class DiagParticles : public Diag {

public :

   DiagParticles( Params &params, SmileiMPI* smpi, Patch* patch, int diagId );
   DiagParticles() {};
   ~DiagParticles();

   virtual void openFile( bool newfile );
   virtual void closeFile();

   virtual void prepare( Patch* patch, int timestep );

   virtual void run( Patch* patch, int timestep );

   virtual void write(int timestep);

private :
    void clean();

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
    
    hid_t fileId_;

};

#endif

