#ifndef OPENPMDPARAMS_H
#define OPENPMDPARAMS_H

#include "Params.h"
#include "H5.h"

class OpenPMDparams {

public:
    OpenPMDparams(Params&);
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
    //! Units of each field
    std::vector<std::vector<double> > unitDimension;
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
    void writeMeshesPathAttributes( hid_t );
    
    //! Write the attributes for a field in the meshesPath
    void writeFieldAttributes( hid_t, unsigned int );
    
    //! Write the attributes for a field record (x, y or z)
    void writeFieldRecordAttributes( hid_t);
    
private:
    Params * params;
    
};

#endif

