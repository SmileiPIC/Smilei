#ifndef DIAG_H
#define DIAG_H

#include "Params.h"
#include "Patch.h"
#include "SmileiMPI.h"


class Diag {

public :

   Diag( Params &params, SmileiMPI* smpi, Patch* patch, int diagId ) {};
   Diag() {};
   ~Diag() {};

   virtual void openFile( Params& params, SmileiMPI* smpi, VectorPatch& vecPatches, bool newfile ) = 0;
   virtual void closeFile() = 0;

   virtual void prepare( Patch* patch, int timestep ) = 0;

   virtual void run( Patch* patch, int timestep ) = 0;

   virtual void write(int timestep) = 0;

    //! Time selection
    TimeSelection * timeSelection;

    //! this is the file name
    std::string filename;
    std::string type_;
protected :
    int probeId_;


};

#endif

