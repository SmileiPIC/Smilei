#ifndef DIAGNOSTICPROBES_H
#define DIAGNOSTICPROBES_H

#include <hdf5.h>

#include "Diagnostic.h"

#include "Params.h"
#include "Patch.h"
#include "SmileiMPI.h"

#include "Field2D.h"


class DiagnosticProbes : public Diagnostic {
    friend class SmileiMPI;

public :
    
    //! Default constructor
    DiagnosticProbes( Params &params, SmileiMPI* smpi, int n_probe );
    //! Default destructor
    ~DiagnosticProbes() ;
    
    void openFile( Params& params, SmileiMPI* smpi, bool newfile ) ;
    
    void closeFile() ;
    
    bool prepare( int timestep ) ;
    
    void run( SmileiMPI* smpi, VectorPatch& vecPatches, int timestep ) ;
    
    void init(Params& params, SmileiMPI* smpi, VectorPatch& vecPatches) ;    
    
    virtual bool needsRhoJs(int timestep) ;
    
    //! Creates the probe's particles (or "points")
    void createPoints(SmileiMPI* smpi, VectorPatch& vecPatches, bool createFile);
    
    //! Probes position storage at initialization and for moving window
    Field2D* posArray;
    
    //! True if load balancing or moving window have taken place 
    bool patchesHaveMoved;
    
    //! If the window has moved, then x_moved contains the movement
    double x_moved;
    
private :
    
    //! Index of the probe diagnostic
    int probe_n;
    
    //! Dimension of the probe grid
    unsigned int dimProbe;
    
    //! Dimension of particle coordinates
    unsigned int nDim_particle;
    
    //! Number of points in each dimension
    std::vector<unsigned int> vecNumber; 
    
    //! List of the coordinates of the probe vertices
    std::vector< std::vector<double> > allPos;
    
    //! Matrix containing the probe's coordinate system
    std::vector<double> axes;
    
    //! Inverse matrix from the probe's coordinate system
    std::vector<double> axesInverse;
    
    //! number of points for this probe
    unsigned int nPart_total;
    
    //! Actual number of points (without those outside the box)
    unsigned int nPart_total_actual;
    
    //! number of point for this probe, in the current MPI process
    unsigned int nPart_MPI;
    
    //! Number of fields to save
    int nFields;
    
    //! List of fields to save
    std::vector<std::string> fieldname;
    
    //! Indices in the output array where each field goes
    std::vector<unsigned int> fieldlocation;
    
    //! Variable to store the status of a dataset (whether it exists or not)
    htri_t status;
    
    //! Temporary buffer to write probes
    Field2D* probesArray;
    
    //! Array to locate the current patch in the local buffer
    std::vector<unsigned int> offset_in_MPI;
    
    //! Array to locate the current patch in the file
    std::vector<unsigned int> offset_in_file;
    
    //! True if this diagnostic requires the pre-calculation of the particle J & Rho
    bool hasRhoJs;
    
};



class ProbeParticles {
public :
    ProbeParticles() {};
    ProbeParticles( ProbeParticles* probe ) { offset_in_file=probe->offset_in_file; }
    ~ProbeParticles() {};
    
    Particles particles;
    int offset_in_file;
};





#endif

