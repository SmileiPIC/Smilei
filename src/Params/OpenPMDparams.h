#ifndef OPENPMDPARAMS_H
#define OPENPMDPARAMS_H

#include "Params.h"
#include "H5.h"

#define SMILEI_NUNITS 9
#define SMILEI_UNIT_NONE     0
#define SMILEI_UNIT_EFIELD   1
#define SMILEI_UNIT_BFIELD   2
#define SMILEI_UNIT_CURRENT  3
#define SMILEI_UNIT_DENSITY  4
#define SMILEI_UNIT_POSITION 5
#define SMILEI_UNIT_MOMENTUM 6
#define SMILEI_UNIT_CHARGE   7
#define SMILEI_UNIT_TIME     8

class OpenPMDparams
{

public:
    OpenPMDparams( Params & );
    virtual ~OpenPMDparams() {};
    
    //! current version of the OpenPMD standard
    std::string version;
    
    //! current extension of the OpenPMD standard
    uint32_t extension;
    
    //! field axes labels
    DividedString axisLabels;
    //! Spacing of the fields grid
    std::vector<double> gridSpacing;
    //! Offset of the grid for all fields
    std::vector<double> gridGlobalOffset;
    //! Offset of the grid for each field
    std::vector<double> gridOffset;
    //! Patch size
    std::vector<unsigned int> patchSize;
    //! Units of each quantity (contains length, mass, time, current, temperature, amount, intensity)
    std::vector<std::vector<double> > unitDimension;
    //! Conversion factor to SI for each quantity
    std::vector<double> unitSI;
    //! field solver description
    std::string fieldSolver, fieldSolverParameters;
    //! boundary conditions names
    DividedString fieldBoundary, fieldBoundaryParameters, particleBoundary, particleBoundaryParameters;
    //! current smoothing description
    std::string currentSmoothing, currentSmoothingParameters;
    //! position of the fields within a cell
    std::vector<double> position;
    
    //! Returns a time string in the openPMD format
    std::string getLocalTime();
    
    //! Write the attributes for the root of the HDF5 file
    void writeRootAttributes( hid_t, std::string, std::string );
    
    //! Write the attributes for the basePath
    void writeBasePathAttributes( hid_t, unsigned int );
    
    //! Write the attributes for the meshesPath
    void writeMeshesAttributes( hid_t );
    
    //! Write the attributes for the particlesPath
    void writeParticlesAttributes( hid_t );
    
    //! Write the attributes for a field in the meshesPath
    void writeFieldAttributes( hid_t, std::vector<unsigned int> subgrid_start= {}, std::vector<unsigned int> subgrid_step= {} );
    
    //! Write the attributes for the particlesPath
    void writeSpeciesAttributes( hid_t );
    
    //! Write the attributes for a record
    void writeRecordAttributes( hid_t, unsigned int );
    
    //! Write the attributes for a field record
    void writeFieldRecordAttributes( hid_t );
    
    //! Write the attributes for a component
    void writeComponentAttributes( hid_t, unsigned int );
    
    
private:
    Params *params;
    
};

#endif

