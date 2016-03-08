#ifndef DIAGPROBES_H
#define DIAGPROBES_H

#include <hdf5.h>

#include "Diag.h"

#include "Params.h"
#include "Patch.h"
#include "SmileiMPI.h"


class DiagProbes : public Diag {

public :

    DiagProbes( Params &params, SmileiMPI* smpi, Patch* patch, int diagId );
    DiagProbes() {};
    ~DiagProbes();

    virtual void openFile( bool newfile );
    virtual void closeFile();

    virtual void prepare( Patch* patch, int timestep );

    virtual void run( Patch* patch, int timestep );

    virtual void write(int timestep);


    void setFile( hid_t masterFileId, Patch* patch, Params& params, VectorPatch& vecPatches );
    void setFile( hid_t masterFileId );
    void writePositionIn( Params &params );
    void writePositions( int ndim_Particles, int probeDim, hid_t group_id );
    void compute(unsigned int timestep, ElectroMagn* EMfields, Interpolator* interp);

private :
    int probeId_;

    //! fake particles acting as probes
    Particles probeParticles;

    //! each probe will write in a buffer
    Field2D* probesArray;

    int probesStart;

    //! number of fake particles for each probe diagnostic
    unsigned int nPart_total;

    //! Number of fields to save
    int nFields;
    
    //! List of fields to save
    std::vector<std::string> fieldname;
    
    //! Indices in the output array where each field goes
    std::vector<unsigned int> fieldlocation;
    
    //! hdf5 file ID
    hid_t fileId_;

    //! return name of the probe based on its number
    std::string probeName();

    //! E local fields for the projector
    LocalFields Eloc_fields;
    //! B local fields for the projector
    LocalFields Bloc_fields;
    //! J local fields for the projector
    LocalFields Jloc_fields;
    //! Rho local field for the projector
    double Rloc_fields;
    
};

#endif

