#ifndef DIAGTRACK_H
#define DIAGTRACK_H

#include "Diag.h"

#include "Params.h"
#include "Patch.h"
#include "SmileiMPI.h"


class DiagTrack : public Diag {

public :

   DiagTrack( Params &params, SmileiMPI* smpi, Patch* patch, int diagId );
   DiagTrack() {};
   ~DiagTrack();

   virtual void openFile( Params& params, SmileiMPI* smpi, VectorPatch& vecPatches, bool newfile );
   void setFileSize( Params& params, SmileiMPI* smpi, VectorPatch& vecPatches );

   virtual void closeFile();

   virtual void prepare( Patch* patch, int timestep );

   virtual void run( Patch* patch, int timestep );

   virtual void write(int timestep);

   void setFile( hid_t fid );

   void setGlobalNbrParticles(int totNbrParts) {
       nbrParticles_ = totNbrParts;
   }

   hid_t getFileId() {
     return fileId_;
   }

private :
   //! Pointer to the species used
   Species* species;
   int speciesId_;

   //! Size of the diag (number of particles)
   int nbrParticles_;
    
   //! HDF5 file transfer protocol
   hid_t transfer;
   //! HDF5 file space (dimensions of the array in file)
   hsize_t dims[2];
    
   //! Number of spatial dimensions
   int nDim_particle;
 
   // iterator for dataset extension
   int iter;

    //! hdf5 file ID
    hid_t fileId_;

    template <class T> void append( hid_t fid, std::string name, T & property,  hid_t  mem_space, int nParticles, hid_t type, std::vector<hsize_t> &locator);
};

#endif

