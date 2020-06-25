#ifndef DIAGNOSTICPROBES_H
#define DIAGNOSTICPROBES_H

#include "Diagnostic.h"

#include "Field2D.h"


class DiagnosticProbes : public Diagnostic
{
    friend class SmileiMPI;
    
public :

    //! Default constructor
    DiagnosticProbes( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, int n_probe );
    //! Default destructor
    ~DiagnosticProbes() override;
    
    void openFile( Params &params, SmileiMPI *smpi ) override;
    
    void closeFile() override;
    
    bool prepare( int timestep ) override;
    
    void run( SmileiMPI *smpi, VectorPatch &vecPatches, int timestep, SimWindow *simWindow, Timers &timers ) override;
    
    void init( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches ) override;
    
    virtual bool needsRhoJs( int timestep ) override;
    
    //! Creates the probe's particles (or "points")
    void createPoints( SmileiMPI *smpi, VectorPatch &vecPatches, bool createFile, double x_moved );
    
    //! Get memory footprint of current diagnostic
    int getMemFootPrint() override
    {
        return nPart_MPI * (
                   // Size of the simili particles structure
                   ( nDim_particle+3+1 )*sizeof( double ) + sizeof( short )
                   // eval probesArray (even if temporary)
                   + 10*sizeof( double )
               );
    }
    
    //! Get disk footprint of current diagnostic
    uint64_t getDiskFootPrint( int istart, int istop, Patch *patch ) override;
    
private :

    //! Index of the probe diagnostic
    int probe_n;
    
    //! Dimension of the probe grid
    unsigned int dimProbe;
    
    //! Dimension of particle coordinates
    unsigned int nDim_particle;
    //! Dimension of field coordinates
    unsigned int nDim_field;
    
    //! Geometry
    std::string geometry;
    
    //! Number of points in each dimension
    std::vector<unsigned int> vecNumber;
    
    //! Coordinates of the probe grid origin
    std::vector<double> origin;
    
    //! List of the coordinates of the probe corners
    std::vector< std::vector<double> > corners;
    
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
    
    //! Special case of species-related fields: indices of fields in ElectroMagn
    std::vector<std::vector<unsigned int> > species_field_index;
    //! Special case of species-related fields: location of fields in output
    std::vector<std::vector<unsigned int> > species_field_location;
    
    //! Variable to store the status of a dataset (whether it exists or not)
    bool has_dataset;
    std::string dataset_name;
    
    //! Temporary buffer to write probes
    Field2D *probesArray;
    
    //! Array to locate the current patch in the local buffer
    std::vector<unsigned int> offset_in_MPI;
    
    //! Array to locate the current patch in the file
    std::vector<unsigned int> offset_in_file;
    
    //! True if this diagnostic requires the pre-calculation of the particle J & Rho
    bool hasRhoJs;
    
    //! Last iteration when points were re-calculated
    unsigned int last_iteration_points_calculated;
    
    //! Whether positions have been written in the file
    bool positions_written;
    
    //! patch size
    std::vector<double> patch_size;
};



class ProbeParticles
{
public :
    ProbeParticles() {};
    ProbeParticles( ProbeParticles *probe )
    {
        offset_in_file=probe->offset_in_file;
    }
    ~ProbeParticles() {};
    
    Particles particles;
    int offset_in_file;
};





#endif
